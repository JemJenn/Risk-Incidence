library(calibrate)
library(scatterplot3d)
library(rgl)
library(mvtnorm)

generate_positions_cumulative_normal = function(n,mean_center_position,sd_size){
  myseq <- seq(1/(2*n),1-1/(2*n),1/n)
  positions<-qnorm( myseq,mean=mean_center_position,sd=sd_size )
  return(positions)
}


timeaxis <- seq(0,200,1)

n <- 10
mean_center_position <- 100
donation_time <- 150
sd_size <- 2

infectiousness_shape <- 2
infectiousness_scale <- 1

detectability_shape <- 5
detectability_scale <- 5

#constant_shift

infectious_positions <- generate_positions_cumulative_normal(n,mean_center_position,sd_size)
detectable_positions <- infectious_positions - 14



# then generate the curves for each set of positions, with much slower delays for the detectability
# curves.

set_of_infectious_curves <- matrix(nrow=length(timeaxis),ncol=n)
for(i in seq(1,n)){
  set_of_infectious_curves[,i] <- individual_negative_likelihood_weibul(times=timeaxis,scale=infectiousness_scale,shape=infectiousness_shape,delay = infectious_positions[i],test_time=donation_time)
}


set_of_detectable_curves <- matrix(nrow=length(timeaxis),ncol=n)
for(i in seq(1,n)){
  set_of_detectable_curves[,i] <- individual_positive_likelihood_weibul(times=timeaxis,scale=detectability_scale,shape=detectability_shape,delay = detectable_positions[i],test_time=donation_time)
}

##copying some plotting code from likelihood_calculator.R

plot(timeaxis,set_of_detectable_curves[,1],type='l',xlim=c(timeaxis[1],timeaxis[length(timeaxis)]),ylim=c(-.002,1.002),xaxt='n',yaxt='n',xaxs='i',yaxs='i',bty='l',xlab='',ylab='',col='green') #clarify label in comment

yaxis_pos <- c(0,0.5,1)
yaxis_names <- c('0',"",'1')

xaxis_pos <- c(donation_time)
xaxis_names <- c(expression('donation time'))

zero_pos <- c(0)
zero_name <- c(expression('0'['']))
#shift axes
axis(side=2, at=yaxis_pos, labels= yaxis_names,tck=-0.037, padj=.437)
axis(side=1, at=xaxis_pos, labels= xaxis_names,padj=-.35,hadj=-.137)
#axis(side=1, at=zero_pos, labels=zero_name,padj=-0.45,hadj=0.37)

# goto_do

#points(plotdata[,1],plotdata[,3])
for (i in seq(1:n)){
  lines(timeaxis,set_of_detectable_curves[,i],col=col_negative,lwd=2)
}
#points(plotdata[,1],plotdata[,3])
for (i in seq(1:n)){
  lines(timeaxis,set_of_infectious_curves[,i],col=3,lwd=2)
}

# positive_mean_background <- rowMeans(plotdata_positive_background)
# negative_mean_background <- rowMeans(plotdata_negative_background)
# 
# lines(timeaxis,negative_mean_background, lwd=lwd_means, col=col_negative)
# lines(timeaxis,positive_mean_background, lwd=lwd_means, col=col_positive)