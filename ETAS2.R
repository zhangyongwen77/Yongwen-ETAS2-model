
ETAS_T<-function(data, params, sequence.length, seed, bvalue, M0) 
{
  set.seed(seed)
  
  if(is.null(data))
  {
    data.time<-Inf 
    data.mag<-0
    
  }
  
  data.mag<-data$mag-M0
  data.time<-data$time
  
  ti<-max(data.time)
  Rmax <- conditional.intensity(data.mag=data.mag,
                                data.time=data.time, eval.time=ti, params=params)
  
  
  repeat
  {
    
    if(Rmax > 0) tau <- rexp(1, rate = Rmax)
    else tau <- Inf
    ti <- ti + tau
    if (ti > sequence.length) break
    #if(length(aftershock.times)>)
    rate <- conditional.intensity(data.mag=data.mag, data.time=data.time, eval.time=ti, params=params)
    if (runif(1, 0, 1) <= rate/Rmax)
    {
      # select a magnitude
      
      
      new.mag <- rexp(1, bvalue * log(10))
      
      data.time <- c(data.time, ti)
      data.mag <- c(data.mag, new.mag)
      
      ti<-max(data.time)
      Rmax <- conditional.intensity(data.mag=data.mag,
                                    data.time=data.time, eval.time=ti, params=params)
    }else{
      Rmax <- rate
    }
  }
  
  
  final.events <- list(times=data.time,
                       magnitudes=data.mag+M0)
  return(final.events)
  
}


conditional.intensity <- function(data.mag, data.time, eval.time, params)
{
  mu <- params[1]
  A <- params[2]
  alpha1 <- params[3]
  alpha2 <- params[6]
  CC <- params[4]
  P <- params[5] 
  coss <- params[10] 
  
  
  if (length(data.time) > 0)
  {
    if (length(data.time) > coss)
    {
      triggering.ci1 <- A*sum(exp(alpha1 * data.mag[(length(data.mag)-coss+1):length(data.mag)]) *
                                 (1 + (eval.time - data.time[(length(data.mag)-coss+1):length(data.mag)])/CC)^(-P))
      
      triggering.ci2 <- A*sum(exp(alpha2 * data.mag[1:(length(data.mag)-coss)]) *
                                 (1 + (eval.time - data.time[1:(length(data.mag)-coss)])/CC)^(-P))
      
      triggering.ci<-triggering.ci1+triggering.ci2
    }else{
      triggering.ci <- A*sum(exp(alpha1 * data.mag) *
                                (1 + (eval.time - data.time)/CC)^(-P))
    }
  }else triggering.ci <- 0
  ci <- mu + triggering.ci
  return(ci)
}


aftershock.xy<-function(x0,y0,data.mag,data.time,lmax,xx,yy,bg,params,seed)
{ 
  
  set.seed(seed)
  data.x<-x0
  data.y<-y0
  
  mu <- params[1]
  A <- params[2]
  alpha1 <- params[3] 
  alpha2 <- params[6]
  CC <- params[4]
  P <- params[5] 
  Q<-params[7]
  D<-params[8]
  gamma<-params[9]
  coss<-params[10]
  
  
  for(t0 in (lmax+1):(length(data.mag)))
  {
    
    triggering.ci1 <- A*exp(alpha1 * data.mag[(t0-coss+1):(t0-1)]) *
      (1 + (data.time[t0] - data.time[(t0-coss+1):(t0-1)])/CC)^(-P)
    
    triggering.ci2 <- A*exp(alpha2 * data.mag[(t0-lmax):(t0-coss)]) *
      (1 + (data.time[t0] - data.time[(t0-lmax):(t0-coss)])/CC)^(-P)
    
    triggering.ci<-c(triggering.ci2,triggering.ci1,mu)
    
    index<-sample(1:length(triggering.ci),1,replace=T,prob = triggering.ci)
    
    if(index<(lmax+1))
    {
      
      sigma<-D*exp(gamma*data.mag[t0])
      index2<-t0-lmax+index-1
      R<-sigma*sqrt(runif(1)^(1/(1-Q))-1)
      theta<-runif(1,0,2*pi)
      rx<-R*cos(theta)+data.x[index2]
      ry<-R*sin(theta)+data.y[index2]
      
    }else{
      index3<-sample(1:length(bg),1,replace = T,prob = bg)
      rx<-xx[index3]
      ry<-yy[index3]
    }
    
    data.x<-c(data.x,rx)
    data.y<-c(data.y,ry)
  }
  
  return(list(x=data.x,y=data.y))
  
}