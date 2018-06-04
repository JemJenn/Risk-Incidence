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

n = 15

mean1<-10
mean2<-15

sd1 <- .2
sd2 <- 5



test1 <- generate_positions_cumulative_normal(n = n,mean_center_position = mean1,sd_size = sd1)
test2 <- generate_positions_cumulative_normal(n = n,mean_center_position = mean2,sd_size = sd2)

the_matrix <- matrix(nrow=length(test1),ncol=length(test2))
other_mat <-  matrix(nrow=length(test1),ncol=length(test2))

for (i in 1:nrow(the_matrix)){
  for (j in 1:ncol(the_matrix)){
    the_matrix[i,j] <- dnorm(abs((test1[i]-mean1)/sd1 - (test2[j]-mean2)/sd2),sd=5,mean=0)
    other_mat[i,j] <- i +j
  }
}


# the_matrix
plot(the_matrix[,1],the_matrix[1,])
plot(other_mat)
other_mat
plot(test1,test2,type='n')
# textxy(.003975,.003975,"yay")



for (i in 1:nrow(the_matrix)){
  for (j in 1:ncol(the_matrix)){
    points(test1[i],test2[j],cex=the_matrix[i,j]*(5/max(the_matrix))*2/sqrt(n))
    # print(the_matrix[i,j])
    textxy(the_matrix[i,j],paste("(",i,",",j,")"))
  }
}

persp(test1,test2,the_matrix, theta = 110, phi = 10, ticktype = "detailed", expand = 0.5, shade = .2)

scatter3D(x=test1,y=test2,z=the_matrix)
plot3d(test1,test2,the_matrix)
plot3d(iris)

