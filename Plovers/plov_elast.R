## plov_elast.m
## 4/8/2009 mbw
## code to calculate elasticity of lambda to low and high level parms
## x is montana, y is e colo, z is south park

rm(list = ls())
# define constants
#ratenames = {'adult survival' 'fledgling survival' 'chick survival' 'nest survival' 'prob egg hatch' '# eggs hatched'};
m.x <- c (0.87, 0.6303, 0.2060, 0.46, 0.91, 3) # mean values for montana
  
  # Is 3 the max number of offspring? Bc max value for r._[,6] is 3....
  # Model assumes only survival and transition and no survival in same lifestage, yeah?
  # So I can sample 10,000 matrixes rather than 10,000 vectors, yeah?
  
m.y <- c(0.87, 0.6302, 0.6719, 0.36, 0.92, 3) # mean values for e colorado
m.z <- c(0.87, 0.6295, 0.2051, 0.63, 0.90, 3) # mean values for south park
v.x <- c(0.05^2, 0.0025^2, 0.06^2, 0.06^2, 0.11^2) # variance for montana estimates
v.y <- c(0.05^2, 0.0025^2, 0.01^2, 0.01^2, 0.08^2) # variance for e colo estimates
v.z <- c(0.05^2, 0.0025^2, 0.06^2, 0.16^2, 0.12^2) # variance for SP estimates

reps<-10000; # define iteration number
# load vectors of random parameter estimates using method of moments for
# beta parameters from empirical estimates given above
r.x <- matrix(0, nrow = reps, ncol  = length(m.x))
r.y <- matrix(0, nrow = reps, ncol  = length(m.y))
r.z <- matrix(0, nrow = reps, ncol  = length(m.z))

for(i in 1:5){
r.x[,i] <- rbeta(shape1 = -(m.x[i]*(v.x[i]+m.x[i]^2-m.x[i]))/v.x[i],
                 shape2 = (v.x[i]+m.x[i]^2-m.x[i])*(m.x[i]-1)/v.x[i],
                 n      = reps)
r.y[,i] = rbeta(shape1 = -(m.x[i]*(v.x[i]+m.x[i]^2-m.x[i]))/v.x[i],
                shape2 = (v.x[i]+m.x[i]^2-m.x[i])*(m.x[i]-1)/v.x[i],
                n      = reps)
r.z[,i] = rbeta(shape1 = -(m.x[i]*(v.x[i]+m.x[i]^2-m.x[i]))/v.x[i],
                shape2 = (v.x[i]+m.x[i]^2-m.x[i])*(m.x[i]-1)/v.x[i],
                n      = reps)
}

r.x[,6] = rbinom(size = 3, p = r.x[,5], n = reps)  
r.y[,6] = rbinom(size = 3, p = r.y[,5], n = reps)
r.z[,6] = rbinom(size = 3, p = r.y[,5], n = reps)

## iterate through the random vectors, build submatrices and track reps

elem.x <- matrix(0, nrow = reps, ncol = 4)
elem.y <- matrix(0, nrow = reps, ncol = 4)
elem.z <- matrix(0, nrow = reps, ncol = 4)

# for(i in 1:reps){
# elem.x[i,] = c(r.x[i,3]*r.x[i,4]*r.x[i,6], # chicksurvival * fledgesurvival * # eggs hatched
#                r.x[i,1]*r.x[i,4]*r.x[i,6], # adultsurvival * fledgesurvival * # eggs hatched
#                r.x[i,3]*r.x[i,2],          # chicksurvival * fledgesurvival
#                r.x[i,1])                   # adult survival
# elem.y[i,] = c(r.y[i,3]*r.y[i,4]*r.x[i,6],
#                r.y[i,1]*r.y[i,4]*r.y[i,6],
#                r.y[i,3]*r.y[i,2],
#                r.y[i,1])
# elem.z[i,] <- c(r.z[i,3]*r.z[i,4]*r.x[i,6],
#                 r.z[i,1]*r.z[i,4]*r.z[i,6], 
#                 r.z[i,3]*r.z[i,2],
#                 r.z[i,1])
# }

# eigenanalysis, sensitivity, elasticity
(lams,lambda.x(i,1)] = eigenall(x);   

## This is function from eigenall.m file


#sensmx = v*w'/(v'*w); # matrix element sensitivities
#elastmx = (sensmx.*x)/lambda.x(i); # matrix element elasticities
[lams,lambda.y(i,1)] = eigenall(y);
#sensmx = v*w'/(v'*w);
#elastmx = (sensmx.*y)/lambda.y(i);
[lams,lambda.z(i,1)] = eigenall(z);
#sensmx = v*w'/(v'*w);
#elastmx = (sensmx.*z)/lambda.z(i);
end
## do summary computations
# compute mean of matrices
mnmx.x = reshape(mean(elem.x),2,2);
mnmx.y = reshape(mean(elem.y),2,2);
mnmx.z = reshape(mean(elem.z),2,2);
# compute covariance among matrix elements
cvmx.x = cov(elem.x);
cvmx.y = cov(elem.y);
cvmx.z = cov(elem.z);
# compute sensitivity and elasticity for mean of matrices
[lams,lamb.mean(1),W,w,V,v] = eigenall(mnmx.x);
sens.x.mn = v*w'/(v'*w);
elas.x.mn = (sens.x.mn.*mnmx.x)/lamb.mean(1);
[lams,lamb.mean(2),W,w,V,v] = eigenall(mnmx.y);
sens.y.mn = v*w'/(v'*w);
elas.y.mn = (sens.y.mn.*mnmx.y)/lamb.mean(2);
[lams,lamb.mean(3),W,w,V,v] = eigenall(mnmx.z);
sens.z.mn = v*w'/(v'*w);
elas.z.mn = (sens.z.mn.*mnmx.z)/lamb.mean(3);

## compute lambdas at max vals
# max adult surv (1)
a11=mean(r.x(:,3))*mean(r.x(:,4))*mean(r.x(:,6));
a12=max(r.x(:,1))*mean(r.x(:,4))*mean(r.x(:,6));
a21=mean(r.x(:,3))*mean(r.x(:,2)); a22=max(r.x(:,1));
a=[a11 a12; a21 a22]; [lams,lamb.max(1,1)] = eigenall(a);
a11=mean(r.y(:,3))*mean(r.y(:,4))*mean(r.y(:,6));
a12=max(r.y(:,1))*mean(r.y(:,4))*mean(r.y(:,6));
a21=mean(r.y(:,3))*mean(r.y(:,2)); a22=max(r.y(:,1));
a=[a11 a12; a21 a22]; [lams,lamb.max(2,1)] = eigenall(a);
a11=mean(r.z(:,3))*mean(r.z(:,4))*mean(r.z(:,6));
a12=max(r.z(:,1))*mean(r.z(:,4))*mean(r.z(:,6));
a21=mean(r.z(:,3))*mean(r.z(:,2)); a22=max(r.z(:,1));
a=[a11 a12; a21 a22]; [lams,lamb.max(3,1)] = eigenall(a);

# max fledge surv (2)
a11=mean(r.x(:,3))*mean(r.x(:,4))*mean(r.x(:,6));
a12=mean(r.x(:,1))*mean(r.x(:,4))*mean(r.x(:,6));
a21=mean(r.x(:,3))*max(r.x(:,2)); a22=mean(r.x(:,1));
a=[a11 a12; a21 a22]; [lams,lamb.max(1,2)] = eigenall(a);
a11=mean(r.y(:,3))*mean(r.y(:,4))*mean(r.y(:,6));
a12=mean(r.y(:,1))*mean(r.y(:,4))*mean(r.y(:,6));
a21=mean(r.y(:,3))*max(r.y(:,2)); a22=mean(r.y(:,1));
a=[a11 a12; a21 a22]; [lams,lamb.max(2,2)] = eigenall(a);
a11=mean(r.z(:,3))*mean(r.z(:,4))*mean(r.z(:,6));
a12=mean(r.z(:,1))*mean(r.z(:,4))*mean(r.z(:,6));
a21=mean(r.z(:,3))*max(r.z(:,2)); a22=mean(r.z(:,1));
a=[a11 a12; a21 a22]; [lams,lamb.max(3,2)] = eigenall(a);

# max chick surv (3)
a11=max(r.x(:,3))*mean(r.x(:,4))*mean(r.x(:,6));
a12=mean(r.x(:,1))*mean(r.x(:,4))*mean(r.x(:,6));
a21=max(r.x(:,3))*mean(r.x(:,2)); a22=mean(r.x(:,1));
a=[a11 a12; a21 a22]; [lams,lamb.max(1,3)] = eigenall(a);
a11=max(r.y(:,3))*mean(r.y(:,4))*mean(r.y(:,6));
a12=mean(r.y(:,1))*mean(r.y(:,4))*mean(r.y(:,6));
a21=max(r.y(:,3))*mean(r.y(:,2)); a22=mean(r.y(:,1));
a=[a11 a12; a21 a22]; [lams,lamb.max(2,3)] = eigenall(a);
a11=max(r.z(:,3))*mean(r.z(:,4))*mean(r.z(:,6));
a12=mean(r.z(:,1))*mean(r.z(:,4))*mean(r.z(:,6));
a21=max(r.z(:,3))*mean(r.z(:,2)); a22=mean(r.z(:,1));
a=[a11 a12; a21 a22]; [lams,lamb.max(3,3)] = eigenall(a);

# max nest surv (4)
a11=mean(r.x(:,3))*max(r.x(:,4))*mean(r.x(:,6));
a12=mean(r.x(:,1))*max(r.x(:,4))*mean(r.x(:,6));
a21=mean(r.x(:,3))*mean(r.x(:,2)); a22=mean(r.x(:,1));
a=[a11 a12; a21 a22]; [lams,lamb.max(1,4)] = eigenall(a);
a11=mean(r.y(:,3))*max(r.y(:,4))*mean(r.y(:,6));
a12=mean(r.y(:,1))*max(r.y(:,4))*mean(r.y(:,6));
a21=mean(r.y(:,3))*mean(r.y(:,2)); a22=mean(r.y(:,1));
a=[a11 a12; a21 a22]; [lams,lamb.max(2,4)] = eigenall(a);
a11=mean(r.z(:,3))*max(r.z(:,4))*mean(r.z(:,6));
a12=mean(r.z(:,1))*max(r.z(:,4))*mean(r.z(:,6));
a21=mean(r.z(:,3))*mean(r.z(:,2)); a22=mean(r.z(:,1));
a=[a11 a12; a21 a22]; [lams,lamb.max(3,4)] = eigenall(a);

# max number of eggs (6)
a11=mean(r.x(:,3))*mean(r.x(:,4))*max(r.x(:,6));
a12=mean(r.x(:,1))*mean(r.x(:,4))*max(r.x(:,6));
a21=mean(r.x(:,3))*mean(r.x(:,2)); a22=mean(r.x(:,1));
a=[a11 a12; a21 a22]; [lams,lamb.max(1,5)] = eigenall(a);
a11=mean(r.y(:,3))*mean(r.y(:,4))*max(r.y(:,6));
a12=mean(r.y(:,1))*mean(r.y(:,4))*max(r.y(:,6));
a21=mean(r.y(:,3))*mean(r.y(:,2)); a22=mean(r.y(:,1));
a=[a11 a12; a21 a22]; [lams,lamb.max(2,5)] = eigenall(a);
a11=mean(r.z(:,3))*mean(r.z(:,4))*max(r.z(:,6));
a12=mean(r.z(:,1))*mean(r.z(:,4))*max(r.z(:,6));
a21=mean(r.z(:,3))*mean(r.z(:,2)); a22=mean(r.z(:,1));
a=[a11 a12; a21 a22]; [lams,lamb.max(3,5)] = eigenall(a);

## compute stochastic lambda, sensitivity/elasticity of st.lambda to matrix elements
[sens.x.st,elas.x.st,stlam(1)] = stochsens(mnmx.x,cvmx.x);
[sens.y.st,elas.y.st,stlam(2)] = stochsens(mnmx.y,cvmx.y);
[sens.z.st,elas.z.st,stlam(3)] = stochsens(mnmx.z,cvmx.z);
elast.stoch(:,:,1)=elas.x.st.mn;
elast.stoch(:,:,2)=elas.y.st.mn;
elast.stoch(:,:,3)=elas.z.st.mn;
elast.mean(:,:,1)=elas.x.mn;
elast.mean(:,:,2)=elas.y.mn;
elast.mean(:,:,3)=elas.z.mn;
for i = 1:5;
lamdiff(1,i)=lamb.max(1,i)-lamb.mean(1);
lamdiff(2,i)=lamb.max(2,i)-lamb.mean(2);
lamdiff(3,i)=lamb.max(3,i)-lamb.mean(3);
end 

save results elem lambda elast lamb mnmx cvmx r m v stlam lamdiff
clear all

## compute regression stats of lambda on vital rates
load results lambda r
for i = 1:6;
r2x=regstats(lambda.x,r.x(:,i),'linear','rsquare');
r2y=regstats(lambda.y,r.y(:,i),'linear','rsquare');
r2z=regstats(lambda.z,r.z(:,i),'linear','rsquare');
r2(1,i)=r2x.rsquare;
r2(2,i)=r2y.rsquare;
r2(3,i)=r2z.rsquare;
end

save reg_stats r2
clear all