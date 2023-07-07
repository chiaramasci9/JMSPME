JMSPME_discrete_mass_plot<-function(mod){
  case=mod$case
  knots=mod$Knots
  weights=mod$weights
  if(case=='int'){
    n_ran = 1
    random_param = 1
  }
  
  if(case=='slope'){
    n_ran = 1
    random_param = 3
  }
  
  if(case=='inteslope'){
    n_ran = 2
    random_param = 13
  }
  
  
  if(n_ran==2){
    par(mfrow=c(2,1))
    print(plot(knots[[1]][[1]],knots[[1]][[2]],lwd=100*weights[[1]],xlab='Intercept',ylab='Slope',pch=19))
    print(plot(knots[[2]][[1]],knots[[2]][[2]],lwd=100*weights[[2]],xlab='Intercept',ylab='Slope',pch=19))
  }
  else{
    if(random_param==3){
      par(mfrow=c(2,1))
      print(plot(knots[[1]],lwd=100*weights[[1]],xlab='Index',ylab='Slope',pch=19))
      print(plot(knots[[2]],lwd=100*weights[[2]],xlab='Index',ylab='Slope',pch=19))
    }
    else{
      par(mfrow=c(2,1))
      print(plot(knots[[1]],lwd=100*weights[[1]],xlab='Index',ylab='Intercept',pch=19))
      print(plot(knots[[2]],lwd=100*weights[[2]],xlab='Index',ylab='Intercept',pch=19))
    }
  }
}
