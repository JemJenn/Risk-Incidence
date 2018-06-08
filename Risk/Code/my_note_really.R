## Bivariate normal matrix generator
## This script generates matrices whose elements contain the probabilities that somebody exhibiting test responsiveness
## corresponding to a particular delay between infection and testing ("lying on a particular curve") will become infectious
## according to the curve with the same rank (eg the second curve); will become infectious in the same order relative
## to the rest of the population..
## This probability, known as a correlation coefficient, will be calculated using a  bivariate normal distribution. the inputs to the
## bivariate normal distribution will be the distances from their respective means of the position variables [ie the difference between
## the individual test sensitivity curve's test delay and the population mean test delay, and the difference between the individual infectiousness curve's delay
## and the population-level mean infectiousness delay.] The index of a particular sensitivity/infectiousness curve (which has an associated delay)
## gives the column/row number will give 
##

# couldn't I just have a matrix that reflects the probability of a certain "deviance" ie how "far"
# in time-delay is the second test from the curve that correlates "best" to the first?
library(calibrate)
library(scatterplot3d)
library(rgl)
library(mvtnorm)

generate_positions_cumulative_normal = function(n,mean_center_position,sd_size){
  myseq <- seq(1/(2*n),1-1/(2*n),1/n)
  positions<-qnorm( myseq,mean=mean_center_position,sd=sd_size )
  return(positions)
}

likelihood_by_DDI_mixed = function(set_of_positive_curves,set_of_negative_curves,matrix_of_correlations,times){ 
  # accepts sets of curves, meant to represent probabilities of achieving certain test results (names "positive/negative" aren't strict)
  # outputs the likelihood for each time point that the HDI (hypothetical date of infection) was on that time point,
  # given a certain structure of correlation between the population-distribution of delay profiles for each test.
  # the notation here ("forward" or "backward")) assumes a negative test followed by a positive test
  
  if(ncol(set_of_negative_curves) != ncol(set_of_positive_curves)) {stop ("Huh?? Why are there different numbers of curves?")}
  if(nrow(set_of_negative_curves) != nrow(set_of_positive_curves)) {stop ("Huh?? Why are there different numbers of time steps?")}
  if(ncol(set_of_negative_curves) != nrow(matrix_of_correlations)) {stop ("The correlation matrix has a different number of curves to to the data. Re-compute with same inputs")}
  
  likelihoods_per_time_forward = rep(0,length(times))
  likelihoods_per_time_backward = rep(0,length(times))
  
  for (tpos in seq(1:length(times))){
    cumu_likely_forward <- 0
    
    for (curvenumber in seq(1:ncol(set_of_negative_curves))){ #chronological order makes a trivial difference - there is 
      #just a product per person
      for (other_curve in seq(1:ncol(set_of_positive_curves))){
        likelihood_forward <- set_of_positive_curves[tpos,curvenumber]*set_of_negative_curves[tpos,other_curve]*the_matrix[curvenumber,other_curve]
        likelihood_backward <- "The same thing unless I don't actually understand this task"
        cumu_likely_forward <- cumu_likely_forward + likelihood_forward 
      }
    }
    
    cumu_likely_backward <- 0
    for (curvenumber in seq(1:ncol(set_of_positive_curves))){ #chronological order makes a trivial difference - there is 
      #just a product per person
      for (other_curve in seq(1:ncol(set_of_negative_curves))){
        likelihood_backward <- set_of_negative_curves[tpos,curvenumber]*set_of_positive_curves[tpos,other_curve]*the_matrix[curvenumber,other_curve]
        cumu_likely_backward <- cumu_likely_backward + likelihood_backward
      }
    }
    likelihoods_per_time_forward[tpos] <- cumu_likely_forward   #ncol(set_of_negative_curves)
    likelihoods_per_time_backward[tpos] <- cumu_likely_backward
  }
  
  # check how much the two directions differ  
  print (max(abs(likelihoods_per_time_backward - likelihoods_per_time_forward)))
  

  # check whether the tiny difference is likely to effect result
  # refactored_likelihood <- (likelihoods_per_time_backward + likelihoods_per_time_backward)*1/2*max(abs(likelihoods_per_time_backward - likelihoods_per_time_forward))
  # plot(timeaxis,abs(likelihoods_per_time_backward - likelihoods_per_time_forward))
  # lines(timeaxis,refactored_likelihood)
  # 
  
  return(likelihoods_per_time_forward)
}


# Parameters

n = 10 + 50

mean1 <- 10
mean2 <- 15

sd1 <- 5
sd2 <- 5

sd_correlation <- .05

# test positions
test1 <- generate_positions_cumulative_normal(n = n,mean_center_position = mean1,sd_size = sd1)
test2 <- generate_positions_cumulative_normal(n = n,mean_center_position = mean2,sd_size = sd2)


# create the correlation matrix

the_matrix <- matrix(nrow=length(test1),ncol=length(test2))
# other_mat <-  matrix(nrow=length(test1),ncol=length(test2))
normalise <- 0
for (i in 1:nrow(the_matrix)){
  for (j in 1:ncol(the_matrix)){
    the_matrix[i,j] <- dnorm(abs((test1[i]-mean1)/sd1 - (test2[j]-mean2)/sd2),sd=sd_correlation,mean=0)
    normalise <- normalise + the_matrix[i,j]
    # other_mat[i,j] <- i +j
  }
}
the_matrix <- the_matrix / normalise


#plot the correlation matrix in two dimensions

plot(test1,test2,type='n')
for (i in 1:nrow(the_matrix)){
  for (j in 1:ncol(the_matrix)){
    points(test1[i],test2[j],cex=the_matrix[i,j]*(5/max(the_matrix))*2/sqrt(n))
    # print(the_matrix[i,j])
    textxy(test1[i],test2[j],paste("(",i,",",j,")"))
  }
}

# plot the correlation matrix in three dimensions
persp(test1,test2,the_matrix, theta = 110, phi = 10, ticktype = "detailed", expand = 0.5, shade = .2)

# check the relative deviance from symetry in the correlation matrix (due to numerical error)

max(abs(the_matrix - t(the_matrix))/(the_matrix/2+t(the_matrix)/2))

# 
hold_on

likelihoods_mixed <- likelihood_by_DDI_mixed(set_of_positive_curves = plotdata_positive_background,set_of_negative_curves = plotdata_negative_background,matrix_of_correlations = the_matrix,times = timeaxis)

lines(timeaxis,likelihoods_mixed,lwd = 2,col=3)
# scatter3D(x=test1,y=test2,z=the_matrix)
# plot3d(test1,test2,the_matrix)

# the next step in the risk estimation is to calculate the likelihood of a particular time position by adding the normalised likelihoods
# of all potential curve combinations.
# for (focus on each curve for one test). check each other curve, and calculate the [implied to be the average] likelihood of a person being on the 
# focus curve for the focus test while being on the other curve for the other test. Add these together, then add all the likelihoods obtained from the
# focus curves (and divide by the number of focus curves) to get the total likelihood of a point.



