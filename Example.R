


source(file = "~/ETAS2.R")               # change your Path
load(file = "~/Initial_data.RData")


#===============ETAS1 model===================

lmax<-4000
history_data<-italy_catlog[1:lmax,]        #input historical data from the real catalog 

st<-max(history_data$time)+10000           # maximum time length of simulation  
mu=0.2
A=6.26
c=0.007
alpha1<-1.5
alpha2<-1.5
p<-1.13
Q<-2
D<-0.03
gamma<-0.48
coss<-200     # Cossover k
M0<-3         # Magnitude threshold
params<- c(mu,A,alpha1,c,p,alpha2,Q,D,gamma,coss)

i<-1
sim_catlog<-ETAS_T(data=history_data, params=params, sequence.length=st, seed=i, bvalue=1, M0=M0)
data.xy<-aftershock.xy(history_data$x,history_data$y,sim_catlog$magnitudes-M0,sim_catlog$times,lmax,italy_bg_density$lon,italy_bg_density$lat,italy_bg_density$prb, params, seed=i)
sim_catlog[["lon"]]<-data.xy$x
sim_catlog[["lat"]]<-data.xy$y


sim_catlog$times<-sim_catlog$times[(lmax+1):length(sim_catlog$times)]
sim_catlog$magnitudes<-sim_catlog$magnitudes[(lmax+1):length(sim_catlog$magnitudes)]
sim_catlog$lon<-sim_catlog$lon[(lmax+1):length(sim_catlog$lon)]
sim_catlog$lat<-sim_catlog$lat[(lmax+1):length(sim_catlog$lat)]

#===========================ETAS2 model=============================================


lmax<-4000
history_data<-italy_catlog[1:lmax,]       #input history data from real catalog 

st<-max(history_data$time)+10000
mu=0.2
A=3.26
c=0.007
alpha1<-2.0
alpha2<-1.4
p<-1.13
Q<-2
D<-0.03
gamma<-0.48
coss<-200         #Cossover k
M0<-3            # Magnitude threshold 
params<- c(mu,A,alpha1,c,p,alpha2,Q,D,gamma,coss)

i<-1
sim_catlog<-ETAS_T(data=history_data, params=params, sequence.length=st, seed=i, bvalue=1, M0=M0)
data.xy<-aftershock.xy(history_data$x,history_data$y,sim_catlog$magnitudes-M0,sim_catlog$times,lmax,italy_bg_density$lon,italy_bg_density$lat,italy_bg_density$prb, params, seed=i)
sim_catlog[["lon"]]<-data.xy$x
sim_catlog[["lat"]]<-data.xy$y


sim_catlog$times<-sim_catlog$times[(lmax+1):length(sim_catlog$times)]
sim_catlog$magnitudes<-sim_catlog$magnitudes[(lmax+1):length(sim_catlog$magnitudes)]
sim_catlog$lon<-sim_catlog$lon[(lmax+1):length(sim_catlog$lon)]
sim_catlog$lat<-sim_catlog$lat[(lmax+1):length(sim_catlog$lat)]



  
