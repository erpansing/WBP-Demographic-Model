##-----------------------------------------------------------##
##                                                           ##
##                 WBP Demographic Model                     ##
##                    Libby Pansing                          ##
##-----------------------------------------------------------##

# Note... having issues with masking of dplyr functions when loading
# plyr later in the script. IF you want to re-source the other
# scripts, you'll need to unload dplyr and plyr and reload plyr.
# Or restart R.
library(dplyr)
library(popbio)
library(tidyr)

## Import relevant parameter estimates calculated in other scripts.

## Required files:
## 1) GRIN parameter estimates.R
## 2) WBP Survival Estimates RMark Nest.R
## 3) 2017 YNP Data.Rda

source("/Users/elizabethpansing/Documents/Test/Field-et-al.-Model/WBP Demographic Model/GRIN parameter estimates.R") 


## Begin by defining the distributions from which survival and
## transition rates will be drawn to create stochastic matrixes


##-----------------------------------------------------------##
##              Empirically derived survival                 ##
##              & transition distributions                   ##
##-----------------------------------------------------------##


##-----------------------------------------------------------##
##                          SEED1                            ##
##-----------------------------------------------------------##

s_SEED1 <- 0                         # Assume seeds either transition to SEED2, 
                                     # transition to CS (first year seedling)
                                     # or die


t1_SEED1    <- function(size){        # survival probability of seeds
  rbeta(n = size,                     # Drawn from a beta to give a 
        shape1 = SEED1_survive_alpha, # probability of seeds transitioning
        shape2 = SEED1_survive_beta)  # to SEED2 stage
} 



t2_SEED1 <- function(size){           # Germination probability of seeds
  rbeta(n = size,                     # Drawn from a beta to give a 
        shape1 = SEED1_germ_alpha,    # probability of seeds transitioning
        shape2 = SEED1_germ_beta)     # to SEED2 stage
}

##-----------------------------------------------------------##
##                          SEED2                            ##
##-----------------------------------------------------------##

s_SEED2 <- 0                         # Assume seeds either transition to 
                                     # CS (first year seedling) or die

t3_SEED2 <- function(size){
  # Germination probability of seeds
  rbeta(n = size,                     # Drawn from a beta to give a 
        shape1 = SEED2_germ_alpha,    # probability of seeds transitioning
        shape2 = SEED2_germ_beta)     # to from SEED2 to CS stage (i.e., germination)
}                                     

##-----------------------------------------------------------##
##                            CS                             ##
##                     First Year Seedling                   ##
##-----------------------------------------------------------##

s_CS <- 0                            # Assume seeds transition to 
                                     # SD or die

t_CS     <- function(size){          # Survival probability of first 
  rbeta(n = size,                    # year seedlings (cotyledon seedlings)
        shape1 = CS_survive_alpha,   # Drawn from a beta to give prob
        shape2 = CS_survive_beta)    # of transitioning to SD stage
}


##-----------------------------------------------------------##
##                          SD                               ##
##-----------------------------------------------------------##

s_SD     <- function(size){         # survival probability of seedlings
  rbeta(n = size,                   # Drawn from a beta to give prob 
        shape1 = SD_survive_alpha,  # surviving any given year
        shape2 = SD_survive_beta)
}

##-----------------------------------------------------------##
##                           SAP                             ##
##-----------------------------------------------------------##

s_SAP     <- function(){         # Survival probability of saplings
  return(0.8)                    # Still looking for a distribution to 
}                                # use, so for now assuming constant

##-----------------------------------------------------------##
##                            MA                             ##
##-----------------------------------------------------------##

s_MA     <- function(){         # Survival rate of reproductively mature adults
  0.99                          # Assume limited death from senescence because 
}                               # of long lived nature of wbp (lifespan up to 1200 yrs)
                                # For now, assuming constant.

##-----------------------------------------------------------##
##                          DEFINE                           ##
##                         SURVIVAL                          ##
##                          VECTOR                           ##
##-----------------------------------------------------------##

survival_vector <- function(size){   #survival vector
  c(
    s_SD(size = size),
    s_SAP(),
    s_MA())
}

survival_vector(size = 1)

##-----------------------------------------------------------##
##                          DEFINE                           ##
##                       RESIDENCE TIME                      ##
##                          VECTOR                           ##
##-----------------------------------------------------------##

residence_SEED1 <- 1    # # years as seedling (SEED1)

residence_SEED2 <- 1    # # years as seedling (SEED2)

residence_CS    <- 1    # # years as seedling (CS)

residence_SD    <- 29   # # years as seedling (SD)

residence_SAP   <- 20   # # years as sapling (SAP)

residence_MA   <- Inf  # # years as reproductively mature
 
##-----------------------------------------------------------##
##                          DEFINE                           ##
##                       RESIDENCE TIME                      ##
##                          VECTOR                           ##
##-----------------------------------------------------------##


residence_vector <- 
  c(residence_SD,
    residence_SAP,
    residence_MA)


residence_vector


##-----------------------------------------------------------##
##                         FERTILITY                         ##
##-----------------------------------------------------------##

No_caches <- function(t, size){    #Seed production in wbp is periodic
  result <- NULL                   # with masting years every ~ 4 years
                                   # so define cone production as a function
  for(i in 1:size){                # of time described by cos with normally distributed error
    value <- ((12.5*cos(1.5 * t) + 14 + rnorm(1, sd = 3.5)) * 45)/ 3
    if( value >= 0){               # Max values and expected values from
      result[i] <- value           # IGBST cone monitoring since 1980
    } else if(value < 0){          # 
      result[i] <- value - value   # # caches assumes 45 seeds/cone
    }                              # and 3 seeds/cache. All available seeds cachek
  }                                # Assumes 45% of caches created are left for regeneration.
  return(result) 
} 


##-----------------------------------------------------------##
##                  GET MATRIX ELEMENTS                      ##
##-----------------------------------------------------------##
si <- function(size){                  # Gives probability of surviving and staying in the
  (1 - (1/residence_vector)) *         # same life stage for those life stages
    survival_vector(size = size)       # that have residence time > 1 (i.e., persist in the
}                                      # same life stage for > 1 year)

ti <- function(size) {                 # Gives probability of surviving and transitioning
  (1/residence_vector) *               # to the next life stage for those life stages
    survival_vector(size = size)       # that have residence time > 1 (i.e., persist in the
}                                      # same life stage for > 1 year)


S <- function(t){
  matrix(c(               0,              0,             0,        0,        0, No_caches(t,1) *si(1)[3],
                          t1_SEED1(1),    0,             0,        0,        0,                        0,
                          t2_SEED1(1), t3_SEED2(1),      0,        0,        0,                        0,
                          0,              0,       t_CS(1), si(1)[1],        0,                        0,
                          0,              0,             0,  ti(1)[1], si(1)[2],                       0,
                          0,              0,             0,        0,  ti(1)[2],               si(1)[3]),
         byrow = T, nrow = 6, 
         dimnames = list(c("SEED1", "SEED2", "CS", "SD", "SAP", "MA"),
                         c("SEED1", "SEED2", "CS", "SD", "SAP", "MA"))) 
}

t <- 1
S(t) 


##-----------------------------------------------------------##
##              Function that projects pop                   ##
##              sizes and incorporates fire                  ##
##-----------------------------------------------------------##

library(plyr)

n <- c(62, 580 + 38, 79, 65, 91,  353) # Arbitrary starting pop size vectors


project <- function(projection_time, n0, reps = 100){       # stochastic projection function 
                                                            # that tracks stage based pop sizes 
                                                            # over time for reps number of iterations
  
  results <- stochastic_matrixes <-                          #Create null matrix that will hold stage
    array(0, dim = c(projection_time, length(n0) + 1, reps)) # based population sizes
  
  for(j in 1:reps){       # Iterate through i years (projection_time) of population growth j times (iterations)
    # Incorporate FIRE
    if(j == 1){           
      intervals <- NULL   # Create null vectors that will hold fire intervals for each iteration
      iteration <- NULL   # Create null vectors that will show the iteration during which each fire occured
    } else if(j != 1){
      intervals <- intervals
      iteration <- iteration
    }
    
    interval <- rgamma(4, shape = fire_alpha, rate = fire_beta)  %>%     
      round(., 0) %>%     # Select fire years from a gamma distribution (assumes mean fire return interval = 233 yrs)
      cumsum(.)           # representing the waiting times btwn fires
                          # cumsum gives the time = t years during which fires should occur in projection
    interval <- interval[-which(interval > projection_time)]  #trim to only include those within the projection time
    
    ## Creates dataframe result that tracks the fire years for each iteration.
    intervals <- append(intervals, interval, after = length(intervals))
    iteration <- append(iteration, rep(j, length(interval)), 
                        after = length(iteration))
    
    
    pops <- matrix(0, nrow = projection_time , ncol = length(n)) # Creates empty matrix to hold population sizes
    
    for(i in 1:projection_time){           # get population
      t <-  i                              # time counter
      tSinceFire <- ifelse(i == 1, 1, tSinceFire)
      fire <- ifelse(t %in% interval, TRUE, FALSE)  # Determine whether this iteration experiences a stand replacing fire
      
      
      # Now, multiple possibilities
      # 1) There's a fire. This kills the population. Assumes stand replacing burn that impacts entire population
      #    And that no regeneration occurs the year of the fire
      # 2) It's a year after a fire, in which case we assume input from an outside population (i.e., system
      #    isn't entirely closed) and that the outside population is on a similar masting schedule as our 
      #    population
      # 3) There's no fire and it's >1 year after fire. In this case, there are no modifications and the 
      #    system proceeds as normal.
      
      if(fire == T){                  
        
        tSinceFire <- 1 
        pops[i,] <- c(0, 0, 0, 0, 0, 0)    # Assuming stand replacing burn with no survival and no regeneration. 
        n <- pops[i,]                      # Most fires go out with first snow. e.g., Romme 1982
        
      } else if(fire == F & tSinceFire == 2 & t != 2){
        tSinceFire <- tSinceFire + 1
        pops[i,] <- c(No_caches(size = 1, t = t), 0, 0, 0, 0, 0)
        n <- pops[i,]
        
      } else if(fire == F & tSinceFire == 2 & t == 2){
        tSinceFire <- tSinceFire + 1
        
        mat <- S(t = t)
        pops[i,] <- t(mat %*% as.matrix(n, nrow = length(n), ncol = 1))
        n <- as.matrix(pops[i,], nrow = length(pops[i,]), ncol = 1)
        
      } else if(fire == F & tSinceFire != 2){
        tSinceFire <- tSinceFire + 1
        
        mat <- S(t = t)
        pops[i,] <- t(mat %*% as.matrix(n, nrow = length(n), ncol = 1))
        n <- as.matrix(pops[i,], nrow = length(pops[i,]), ncol = 1)
      }
      
    }
    
    pops <- cbind(pops, rep(1:projection_time))  # Appends matrix to keep track of time during iteration
    
    
    results[, ,j] <- pops                        # Combines iterations into a j dimensional array
    
  }
  
  pop_sizes <- plyr::adply(results, 3)           # Changes array to dataframe so easier to manipulate later
  colnames(pop_sizes) <- c("Iteration", "SEED1", "SEED2", "CS", "SD", "SAP", "MA", "t")
  
  fire_intervals <- cbind(iteration, intervals)  # Dataframe keeping track of fire return intervals
  
  results <- list(pop_sizes = pop_sizes, fire_intervals = fire_intervals)
  
  return(results)
}


##-----------------------------------------------------------##
##                 Population projection                     ##
##-----------------------------------------------------------##

projection <- project(projection_time = 500, n0 = n, reps = 50) 

projection <- gather(projection$pop_sizes, Stage, Count, -Iteration, -t) %>%  
  filter(., !Stage == "SEED1") %>%          # Pop sizes in dataframe format
  filter(., !Stage == "SEED2") %>%          # and excluding seed numbers (most don't think of seeds)
  group_by(., Iteration, t) %>%             # as a part of the population, so presenting numbers as 
  summarise_at(., vars(Count), funs(sum)) %>% # number of living trees (i.e., post germination) is more
  ungroup(.)                                # intuitive


ggplot(data = projection, aes(x = t, y = Count, col = Iteration)) +  # plot pop sizes for all iterations.
  geom_line(lwd = 1) +
  theme(legend.position="none") +
  theme(axis.title.x=element_text( size=18, vjust=0)) +
  theme(axis.text.x=element_text(size=18))  +
  theme(axis.title.y=element_text( size=18, vjust=2.75, face = "bold")) +
  theme(axis.text.y=element_text(size = 18))


hist(projection$Count[projection$t == 500], breaks = 20,
     main = "", xlab = "Population size at time = 500 years")

##-----------------------------------------------------------##
##              Get stochastic lambda and                    ##
##                    elasticities                           ##
##-----------------------------------------------------------##

## Create a list of 10,000 matrixes to use in stochastic lambda 
## and elasticity analyses

reps <- 10000
stochastic_matrixes <- array(0, dim = c(6,6,reps))

for(i in 1:reps){
  stochastic_matrixes[, , i] <- S(i)
}

A <- list(NULL)

for(i in 1:reps){
  mat <- stochastic_matrixes[,,i]
  A[[i]] <- mat
}

rm(stochastic_matrixes)

## Estimate elasticities
stoch.elast<- stoch.sens(A, tlimit = 500)
stoch.elast

## Estimate stochastic lambda
sgr <- stoch.growth.rate(A)
sgr_real <- exp(sgr$approx)



