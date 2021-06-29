#'Predation under IPCC RCP2.6 or RCP8.5 scenarios
#'
#'@description This function allows simulating the effect of the IPCC RCP2.6 or RCP8.5 scenarios
#'             (2014) on the abundances of two species interacting through predation. The prey is an
#'             ectotherm population, and the predator is not affected by temperature.
#'
#'
#'
#'@param y_ini Initial population values (must be written with its name: N).
#'@param temp_ini Initial temperature.
#'@param temp_cmin Minimum critical temperature (prey).
#'@param temp_cmax Maximum critical temperature (prey).
#'@param ro Population growth rate at optimum temperature (prey).
#'@param lambda Marginal loss by non-thermodependent intraspecific competition (prey).
#'@param e Efficiency with which food is converted into population growth (predator).
#'@param mp Mortality rate (predator).
#'@param q Maximum per capita consumption rate (predator).
#'@param a Mean saturation rate (predator).
#'@param RCP Representative concentration trajectories (RCP2.6 and RCP8.5 scenarios).
#'@param time_start Start of time sequence.
#'@param time_end End of time sequence.
#'@param leap Time sequence step.
#'
#'@details Three scenarios can be evaluated for a predation interaction where the prey is an ectotherm
#'         population. The temperature trends correspond to IPCC projections under the RCP2.6 or RCP8.5
#'         scenarios. In each input vector, the parameters for the three simulations must be specified.
#'
#'
#'@return (1) A data.frame with columns having the simulated trends.
#'@return (2) A four-panel figure where (a), (b), and (c) show the population abundance curves for each simulation.
#'        The brown curve corresponds to the abundance of prey and the green curve to predators. Panel (d)
#'        shows the temperature trend curves used for each simulation, green, blue, and black, respectively.
#'
#'@references IPCC. (2014): Climate Change 2014: Synthesis Report. Contribution of Working Groups I,
#'            II and III to the Fifth Assessment Report of the Intergovernmental Panel on Climate
#'            Change [Core Writing Team, R.K. Pachauri and L.A. Meyer (eds.)]. IPCC, Geneva,
#'            Switzerland, 151 pp.
#'
#'@export
#'@examples
#'
#'#######################################################################
#'   #Example 1: Different thermal tolerance ranges (scenario RCP2.6).
#'#######################################################################
#'
#'temp_cmin <- 18
#'
#'# Temperature that occurs before the minimum simulation time.
#'temp_i <- 22
#'
#'time_end <- 2100
#'
#'# Temperature that occurs in the maximum time of the simulation.
#'temp_max <- get_RCP2.6(time_end)+temp_i
#'
#'# Simulation thermal range.
#'RS <- temp_max-temp_cmin
#'
#'temp_cmax1 <- 4/3*RS+temp_cmin
#'temp_cmax2 <- 2/3*RS+temp_cmin
#'temp_cmax3 <- 1/3*RS+temp_cmin
#'temp_ini <- (temp_cmin+temp_cmax3)/2
#'
#'predation(y_ini = c(V = 800, V = 800, V = 800,
#'                    P = 600, P = 600, P = 600),
#'          temp_ini = rep(temp_ini,3),
#'          temp_cmin = rep(temp_cmin,3),
#'          temp_cmax = c(temp_cmax1,temp_cmax2,temp_cmax3),
#'          ro = rep(0.7,3),
#'          lambda = rep(0.0004,3),
#'          e = rep(0.9,3),
#'          mp = rep(0.1,3),
#'          q = rep(0.7,3),
#'          a = rep(1000,3),
#'          RCP = 2.6,
#'          time_start = 2005,
#'          time_end = time_end,
#'          leap = 1/50)
#'\donttest{
#'#######################################################################
#'   #Example 2: Different thermal tolerance ranges (scenario RCP8.5).
#'#######################################################################
#'
#'temp_cmin <- 18
#'
#'# Temperature that occurs before the minimum simulation time.
#'temp_i <- 22
#'
#'time_end <- 2100
#'
#'# Temperature that occurs in the maximum time of the simulation.
#'temp_max <- get_RCP8.5(time_end)+temp_i
#'
#'# Simulation thermal range.
#'RS <- temp_max-temp_cmin
#'
#'temp_cmax1 <- 4/3*RS+temp_cmin
#'temp_cmax2 <- 2/3*RS+temp_cmin
#'temp_cmax3 <- 1/3*RS+temp_cmin
#'temp_ini <- (temp_cmin+temp_cmax3)/2
#'
#'predation(y_ini = c(V = 800, V = 800, V = 800,
#'                    P = 600, P = 600, P = 600),
#'          temp_ini = rep(temp_ini,3),
#'          temp_cmin = rep(temp_cmin,3),
#'          temp_cmax = c(temp_cmax1,temp_cmax2,temp_cmax3),
#'          ro = rep(0.7,3),
#'          lambda = rep(0.0004,3),
#'          e = rep(0.9,3),
#'          mp = rep(0.1,3),
#'          q = rep(0.7,3),
#'          a = rep(1000,3),
#'          RCP = 8.5,
#'          time_start = 2005,
#'          time_end = time_end,
#'          leap = 1/50)
#'
#'#######################################################################
#'   #Example 3: Different conversion efficiencies (scenario RCP2.6).
#'#######################################################################
#'
#'e1 <- 0.2
#'e2 <- 2*e1
#'e3 <- 2*e2
#'
#'predation(y_ini = c(V = 800, V = 800, V = 800,
#'                    P = 400, P = 400, P = 400),
#'          temp_ini = rep(22,3),
#'          temp_cmin = rep(20,3),
#'          temp_cmax = rep(35,3),
#'          ro = rep(0.9,3),
#'          lambda = rep(0.0006,3),
#'          e = c(e1,e2,e3),
#'          mp = rep(0.1,3),
#'          q = rep(0.7,3),
#'          a = rep(800,3),
#'          RCP = 2.6,
#'          time_start = 2005,
#'          time_end = 2100,
#'          leap = 1/50)
#'
#'#######################################################################
#'   #Example 4: Different conversion efficiencies (scenario RCP8.5).
#'#######################################################################
#'
#'e1 <- 0.2
#'e2 <- 2*e1
#'e3 <- 2*e2
#'
#'predation(y_ini = c(V = 800, V = 800, V = 800,
#'                    P = 400, P = 400, P = 400),
#'          temp_ini = rep(22,3),
#'          temp_cmin = rep(20,3),
#'          temp_cmax = rep(35,3),
#'          ro = rep(0.9,3),
#'          lambda = rep(0.0006,3),
#'          e = c(e1,e2,e3),
#'          mp = rep(0.1,3),
#'          q = rep(0.7,3),
#'          a = rep(800,3),
#'          RCP = 8.5,
#'          time_start = 2005,
#'          time_end = 2100,
#'          leap = 1/50)
#'}



predation <- function(y_ini = c(V = 400, V = 400, V = 400,
                                P = 200, P = 200, P = 200),
                      temp_ini = rep(25,3),
                      temp_cmin = rep(18,3),
                      temp_cmax = c(25,28,32),
                      ro = rep(0.7,3),
                      lambda = rep(0.00005,3),
                      e = rep(0.3,3),
                      mp = rep(0.08,3),
                      q = rep(0.7,3),
                      a = rep(800,3),
                      RCP = 2.6,
                      time_start = 2005,
                      time_end = 2100,
                      leap = 1/50){

times<- seq(time_start, time_end, leap)

if(time_end<=2100){
  if(time_start<=time_end){

if(temp_cmin[1]<temp_cmax[1] && temp_cmin[2]<temp_cmax[2] && temp_cmin[3]<temp_cmax[3] ){


if(temp_cmin[1]<=temp_ini[1] && temp_ini[1]<=temp_cmax[1] && temp_cmin[2]<=temp_ini[2] &&
   temp_ini[2]<=temp_cmax[2] && temp_cmin[3]<=temp_ini[3] && temp_ini[3]<=temp_cmax[3]){




##########################################################
# Respuesta funcional
H2<- function(V,q,a){

  (q*V/(V+a))
}
##########################################################

##########################################################
# Optimum growing temperature
##########################################################

temp_op1<- (temp_cmax[1]+temp_cmin[1])/3+sqrt(((temp_cmax[1]+
            temp_cmin[1])/3)^2-(temp_cmax[1]*temp_cmin[1])/3)

temp_op2<- (temp_cmax[2]+temp_cmin[2])/3+sqrt(((temp_cmax[2]+
            temp_cmin[2])/3)^2-(temp_cmax[2]*temp_cmin[2])/3)

temp_op3<- (temp_cmax[3]+temp_cmin[3])/3+sqrt(((temp_cmax[3]+
            temp_cmin[3])/3)^2-(temp_cmax[3]*temp_cmin[3])/3)

##########################################################
# Parameters
##########################################################

parms1<-c(temp_cmin[1],temp_ini[1],temp_cmax[1],temp_op1,ro[1], lambda[1])
parms2<-c(temp_cmin[2],temp_ini[2],temp_cmax[2],temp_op2,ro[2], lambda[2])
parms3<-c(temp_cmin[3],temp_ini[3],temp_cmax[3],temp_op3,ro[3], lambda[3])

##############################################

if(RCP==2.6) {

##########################################################
# Model for each trend
##########################################################
model1 <- function (times, y,parms1) {
  with(as.list(c(y)), {
    T1 <- get_RCP2.6(times)+temp_ini[1]  # IPCC1
    r1<- rate_TPC(T1,ro[1],temp_cmin[1],temp_cmax[1],temp_op1)
    dV <-   r1 * V * (1 - lambda[1]*(V / r1))-H2(V,q[1],a[1])*P
    dP<- e[1]*H2(V,q[1],a[1])*P-mp[1]*P

    return(list(c(dV,dP)))
  })
}
###############################################################

model2 <- function (times, y,parms2) {
  with(as.list(c(y)), {
    T2 <- get_RCP2.6(times)+temp_ini[2]  # IPCC1
    r2<- rate_TPC(T2,ro[2],temp_cmin[2],temp_cmax[2],temp_op2)
    dV <-   r2 * V * (1 - lambda[2]*(V / r2))-H2(V,q[2],a[2])*P
    dP<- e[2]*H2(V,q[2],a[2])*P-mp[2]*P

    return(list(c(dV,dP)))
  })
}
###############################################################

model3 <- function (times, y,parms3) {
  with(as.list(c(y)), {
    T3 <- get_RCP2.6(times)+temp_ini[3]  # IPCC1
    r3<- rate_TPC(T3,ro[3],temp_cmin[3],temp_cmax[3],temp_op3)
    dV <-   r3 * V * (1 - lambda[3]*(V / r3))-H2(V,q[3],a[3])*P
    dP<- e[3]*H2(V,q[3],a[3])*P-mp[3]*P

    return(list(c(dV,dP)))
  })
}
###############################################################
y_ini1<-c(y_ini[1],y_ini[4])
y_ini2<-c(y_ini[2],y_ini[5])
y_ini3<-c(y_ini[3],y_ini[6])

###############################################################
# Solution
##############################################################

out1 <- ode(y=y_ini1, times, model1, parms1,method = "ode45")
out2 <- ode(y=y_ini2, times, model2, parms2,method = "ode45")
out3 <- ode(y=y_ini3, times, model3, parms3,method = "ode45")
#############################################################


###############################################################
# Abundance
##############################################################

data1<-data.frame('x'=times,'y'=out1[,2] )
data2<-data.frame('x'=times,'y'=out1[,3] )

dat1<-data.frame('x'=times,'y'=out2[,2] )
dat2<-data.frame('x'=times,'y'=out2[,3] )

da1<-data.frame('x'=times,'y'=out3[,2] )
da2<-data.frame('x'=times,'y'=out3[,3] )

T1 <- get_RCP2.6(times)+temp_ini[1]
T2 <- get_RCP2.6(times)+temp_ini[2]
T3 <- get_RCP2.6(times)+temp_ini[3]

d1<-data.frame('x'=times,'y'=T1)
d2<-data.frame('x'=times,'y'=T2)
d3<-data.frame('x'=times,'y'=T3)

r1<- rate_TPC(T1,ro[1],temp_cmin[1],temp_cmax[1],temp_op1)
r2<- rate_TPC(T2,ro[2],temp_cmin[2],temp_cmax[2],temp_op2)
r3<- rate_TPC(T3,ro[3],temp_cmin[3],temp_cmax[3],temp_op3)

K1=r1/lambda[1]
K2=r2/lambda[2]
K3=r3/lambda[3]

cap1<-data.frame('x'=times,'y'=K1 )
cap2<-data.frame('x'=times,'y'=K2 )
cap3<-data.frame('x'=times,'y'=K3 )

###############################################################
# Data
###############################################################

Data<- data.frame(times,out1[,2],out1[,3],K1,out2[,2],out2[,3],K2,
                  out3[,2],out3[,3],K3)
names(Data)<- c("Time","Abundance species-1","Abundance species-2",
                "Carrying capacity scenario 1","Abundance species-1",
                "Abundance species-2","Carrying capacity scenario 2",
                "Abundance species-2","Abundance species-2","Carrying
                capacity scenario 3")
u<- formattable(Data, align = c("l", rep("r", NCOL(Data))))
print(u)
################################################################

times_new1<-vector(mode = "numeric", length = 0)
times_new2<-vector(mode = "numeric", length = 0)
times_new3<-vector(mode = "numeric", length = 0)

for (i in 2: length(times)){

  if(out1[i-1,2]>=0 && ( out1[i,2])<0){
    times_new1[i-1]<- times[i-1]
  }else{
    times_new1[i-1]<- 0
  }
}


for (i in 2: length(times)){

  if(out2[i-1,2]>=0 && ( out2[i,2])<0){
    times_new2[i-1]<- times[i-1]
  }else{
    times_new2[i-1]<- 0
  }
}


for (i in 2: length(times)){

  if(out3[i-1,2]>=0 && ( out3[i,2])<0){
    times_new3[i-1]<- times[i-1]
  }else{
    times_new3[i-1]<- 0
  }
}

index1<- which(times_new1!=0)[1]
index2<- which(times_new2!=0)[1]
index3<- which(times_new3!=0)[1]

index1<- as.integer(index1)
index2<- as.integer(index2)
index3<- as.integer(index3)

if(!is.na(as.integer(index1))== FALSE){
  times_sup11<- times[length(times)]
}else{
  times_sup11<- times[index1]
}
if(!is.na(as.integer(index2))== FALSE){
  times_sup21<- times[length(times)]
}else{
  times_sup21<- times[index2]
}

if(!is.na(as.integer(index3))== FALSE){
  times_sup31<- times[length(times)]
}else{
  times_sup31<- times[index3]
}

times_new7<-vector(mode = "numeric", length = 0)
times_new8<-vector(mode = "numeric", length = 0)
times_new9<-vector(mode = "numeric", length = 0)

for (i in 2: length(times)){

  if(( temp_cmax[1]-T1[i-1])>=0 && ( temp_cmax[1]-T1[i])<0){
    times_new7[i-1]<- times[i-1]

  }else if(( temp_cmax[1]-T1[i-1])<=0 && ( temp_cmax[1]-T1[i])>0){

    times_new7[i-1]<- times[i-1]
  }else{
    times_new7[i-1]<- 0
  }
}


for (i in 2: length(times)){

  if(( temp_cmax[2]-T2[i-1])>=0 && ( temp_cmax[2]-T2[i])<0){
    times_new8[i-1]<- times[i-1]

  }else if(( temp_cmax[2]-T2[i-1])<=0 && ( temp_cmax[2]-T2[i])>0){

    times_new8[i-1]<- times[i-1]
  }else{
    times_new8[i-1]<- 0
  }
}


for (i in 2: length(times)){

  if(( temp_cmax[3]-T3[i-1])>=0 && ( temp_cmax[3]-T3[i])<0){
    times_new9[i-1]<- times[i-1]

  }else if(( temp_cmax[3]-T3[i-1])<=0 && ( temp_cmax[3]-T3[i])>0){

    times_new9[i-1]<- times[i-1]
  }else{
    times_new9[i-1]<- 0
  }
}

index7<- which(times_new7!=0)[1]
index8<- which(times_new8!=0)[1]
index9<- which(times_new9!=0)[1]

index7<- as.integer(index7)
index8<- as.integer(index8)
index9<- as.integer(index9)

if(!is.na(as.integer(index7))== FALSE){
  times_sup12<- times[length(times)]
}else{
  times_sup12<- times[index7]
}
if(!is.na(as.integer(index8))== FALSE){
  times_sup22<- times[length(times)]
}else{
  times_sup22<- times[index8]
}

if(!is.na(as.integer(index9))== FALSE){
  times_sup32<- times[length(times)]
}else{
  times_sup32<- times[index9]
}


if(times_sup11<= times_sup12){
  times_sup1<-times_sup11
}else{
  times_sup1<-times_sup12
}

if(times_sup21<= times_sup22){
  times_sup2<-times_sup21
}else{
  times_sup2<-times_sup22
}

if(times_sup31<= times_sup32){
  times_sup3<-times_sup31
}else{
  times_sup3<-times_sup32
}
###############################################################
# Carrying capacity
##############################################################

times_new4<-vector(mode = "numeric", length = 0)
times_new5<-vector(mode = "numeric", length = 0)
times_new6<-vector(mode = "numeric", length = 0)

for (i in 2: length(times)){

  if(K1[i-1]>=0 && ( K1[i])<0){
    times_new4[i-1]<- times[i-1]
  }else{
    times_new4[i-1]<- 0
  }
}


for (i in 2: length(times)){

  if(K2[i-1]>=0 && ( K2[i])<0){
    times_new5[i-1]<- times[i-1]
  }else{
    times_new5[i-1]<- 0
  }
}


for (i in 2: length(times)){

  if(K3[i-1]>=0 && ( K3[i])<0){
    times_new6[i-1]<- times[i-1]
  }else{
    times_new6[i-1]<- 0
  }
}

index4<- which(times_new4!=0)[1]
index5<- which(times_new5!=0)[1]
index6<- which(times_new6!=0)[1]

index4<- as.integer(index4)
index5<- as.integer(index5)
index6<- as.integer(index6)

if(!is.na(as.integer(index4))== FALSE){
  times_sup4<- times[length(times)]
}else{
  times_sup4<- times[index4]
}
if(!is.na(as.integer(index5))== FALSE){
  times_sup5<- times[length(times)]
}else{
  times_sup5<- times[index5]
}

if(!is.na(as.integer(index6))== FALSE){
  times_sup6<- times[length(times)]
}else{
  times_sup6<- times[index6]
}


###############################################################
# Plots
##############################################################

data<-rbind(data1, data2, cap1)

p1<- ggplot(data, aes(x=.data$x, y=.data$y)) +
theme_bw()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_ribbon(data=subset(cap1,times>times[1] & times<times_sup4),aes(x=.data$x,
                                ymax=.data$y),ymin=0,alpha=0.3, fill="brown") +
geom_vline(xintercept = times_sup1, size=.5, color="brown",linetype="dashed")+
geom_line(data =subset(data1,times>times[1] & times<times_sup1), color = "brown")+
geom_line(data =subset(data2,times>times[1] & times<times_sup1), color = "green4")+
labs(x = "Time",y="Abundance")+
theme(plot.title = element_text(size=40))+
theme(plot.title = element_text(hjust = 0.5))+
theme(axis.title.y = element_text(size = rel(1), angle = 90))+
theme(axis.title.x = element_text(size = rel(1), angle = 00))+
labs(tag = "(a)")



dat<-rbind(dat1, dat2, cap2)

p2<- ggplot(dat, aes(x=.data$x, y=.data$y)) +
theme_bw()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_ribbon(data=subset(cap2,times>times[1] & times<times_sup5),aes(x=.data$x,
                                ymax=.data$y),ymin=0,alpha=0.3, fill="brown") +
geom_vline(xintercept = times_sup2, size=.5, color="brown",linetype="dashed")+
geom_line(data =subset(dat1,times>times[1] & times<times_sup2), color = "brown")+
geom_line(data =subset(dat2,times>times[1] & times<times_sup2), color = "green4")+
labs(x = "Time",y="Abundance")+
theme(plot.title = element_text(size=40))+
theme(plot.title = element_text(hjust = 0.5))+
theme(axis.title.y = element_text(size = rel(1), angle = 90))+
theme(axis.title.x = element_text(size = rel(1), angle = 00))+
labs(tag = "(b)")


da<-rbind(da1, da2, cap3)

p3<- ggplot(da, aes(x=.data$x, y=.data$y)) +
theme_bw()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_ribbon(data=subset(cap3,times>times[1] & times<times_sup6),aes(x=.data$x,
                                ymax=.data$y),ymin=0,alpha=0.3, fill="brown") +
geom_vline(xintercept = times_sup3, size=.5, color="brown",linetype="dashed")+
geom_line(data =subset(da1,times>times[1] & times<times_sup3), color = "brown")+
geom_line(data =subset(da2,times>times[1] & times<times_sup3), color = "green4")+
labs(x = "Time",y="Abundance")+
theme(plot.title = element_text(size=40))+
theme(plot.title = element_text(hjust = 0.5))+
theme(axis.title.y = element_text(size = rel(1), angle = 90))+
theme(axis.title.x = element_text(size = rel(1), angle = 00))+
labs(tag = "(c)")


d<-rbind(d1, d2, d3)

p4<- ggplot(d, aes(x=.data$x, y=.data$y)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(xintercept = times_sup1, size=.5, color="green",linetype="dashed")+
  geom_vline(xintercept = times_sup2, size=.5, color="blue",linetype="dashed")+
  geom_vline(xintercept = times_sup3, size=.5, color="black",linetype="dashed")+
  geom_line(data =subset(d1,times>times[1] & times<times_sup1), color = "green")+
  geom_line(data =subset(d2,times>times[1] & times<times_sup2), color = "blue")+
  geom_line(data =subset(d3,times>times[1] & times<times_sup3), color = "black")+
  labs(x = "Time",y="Temperature")+
  theme(plot.title = element_text(size=40))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.y = element_text(size = rel(1), angle = 90))+
  theme(axis.title.x = element_text(size = rel(1), angle = 00))+
  labs(tag = "(d)")

plot_grid(p1, p2,p3,p4)

} else if(RCP==8.5) {

  RCP8.5 <- function(date,a,b) {a * exp(b * date)}
  values <- c(0.61, 2, 3.7)
  x<- c(2005,2065,2100)
  y<- values
  df <- data.frame(x, y)

  m<- nls(y ~ exp(loga + b * x), df, start = list( loga = log(2), b = 0.005),control = list (maxiter = 500))
  y_est<-predict(m,df$x)


##########################################################
  # Model for each trend
##########################################################
  model1 <- function (times, y,parms1) {
  with(as.list(c(y)), {
  T1<- RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[1]    #IPCC2
  r1<- rate_TPC(T1,ro[1],temp_cmin[1],temp_cmax[1],temp_op1)
  dV <-   r1 * V * (1 - lambda[1]*(V / r1))-H2(V,q[1],a[1])*P
  dP<- e[1]*H2(V,q[1],a[1])*P-mp[1]*P

      return(list(c(dV,dP)))
    })
  }
###############################################################

  model2 <- function (times, y,parms2) {
  with(as.list(c(y)), {
  T2<-  RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[2]    #IPCC2
  r2<- rate_TPC(T2,ro[2],temp_cmin[2],temp_cmax[2],temp_op2)
  dV <-   r2 * V * (1 - lambda[2]*(V / r2))-H2(V,q[2],a[2])*P
  dP<- e[2]*H2(V,q[2],a[2])*P-mp[2]*P

      return(list(c(dV,dP)))
    })
  }
###############################################################

  model3 <- function (times, y,parms3) {
  with(as.list(c(y)), {
  T3<- RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[3]    #IPCC2
  r3<- rate_TPC(T3,ro[3],temp_cmin[3],temp_cmax[3],temp_op3)
  dV <-   r3 * V * (1 - lambda[3]*(V / r3))-H2(V,q[3],a[3])*P
  dP<- e[3]*H2(V,q[3],a[3])*P-mp[3]*P

      return(list(c(dV,dP)))
    })
  }
###############################################################
  y_ini1<-c(y_ini[1],y_ini[4])
  y_ini2<-c(y_ini[2],y_ini[5])
  y_ini3<-c(y_ini[3],y_ini[6])

###############################################################
  # Solution
##############################################################

  out1 <- ode(y=y_ini1, times, model1, parms1,method = "ode45")
  out2 <- ode(y=y_ini2, times, model2, parms2,method = "ode45")
  out3 <- ode(y=y_ini3, times, model3, parms3,method = "ode45")
#############################################################

###############################################################
  # Abundance
##############################################################

  data1<-data.frame('x'=times,'y'=out1[,2] )
  data2<-data.frame('x'=times,'y'=out1[,3] )

  dat1<-data.frame('x'=times,'y'=out2[,2] )
  dat2<-data.frame('x'=times,'y'=out2[,3] )

  da1<-data.frame('x'=times,'y'=out3[,2] )
  da2<-data.frame('x'=times,'y'=out3[,3] )

  T1<-  RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[1]
  T2<-  RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[2]
  T3<-  RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[3]

  d1<-data.frame('x'=times,'y'=T1)
  d2<-data.frame('x'=times,'y'=T2)
  d3<-data.frame('x'=times,'y'=T3)

  r1<- rate_TPC(T1,ro[1],temp_cmin[1],temp_cmax[1],temp_op1)
  r2<- rate_TPC(T2,ro[2],temp_cmin[2],temp_cmax[2],temp_op2)
  r3<- rate_TPC(T3,ro[3],temp_cmin[3],temp_cmax[3],temp_op3)

  K1=r1/lambda[1]
  K2=r2/lambda[2]
  K3=r3/lambda[3]


  cap1<-data.frame('x'=times,'y'=K1 )
  cap2<-data.frame('x'=times,'y'=K2 )
  cap3<-data.frame('x'=times,'y'=K3 )

###############################################################
  # Data
###############################################################

Data<- data.frame(times,out1[,2],out1[,3],K1,out2[,2],out2[,3],
                  K2,out3[,2],out3[,3],K3)
names(Data)<- c("Time","Abundance species-1","Abundance species-2",
                "Carrying capacity scenario 1","Abundance species-2",
                "Abundance species-2","Carrying capacity scenario 2",
                "Abundance species-1","Abundance species-2","Carrying
                capacity scenario 3")
u<- formattable(Data, align = c("l", rep("r", NCOL(Data))))
print(u)

###############################################################
  times_new1<-vector(mode = "numeric", length = 0)
  times_new2<-vector(mode = "numeric", length = 0)
  times_new3<-vector(mode = "numeric", length = 0)

  for (i in 2: length(times)){

    if(out1[i-1,2]>=0 && ( out1[i,2])<0){
      times_new1[i-1]<- times[i-1]
    }else{
      times_new1[i-1]<- 0
    }
  }


  for (i in 2: length(times)){

    if(out2[i-1,2]>=0 && ( out2[i,2])<0){
      times_new2[i-1]<- times[i-1]
    }else{
      times_new2[i-1]<- 0
    }
  }


  for (i in 2: length(times)){

    if(out3[i-1,2]>=0 && ( out3[i,2])<0){
      times_new3[i-1]<- times[i-1]
    }else{
      times_new3[i-1]<- 0
    }
  }

  index1<- which(times_new1!=0)[1]
  index2<- which(times_new2!=0)[1]
  index3<- which(times_new3!=0)[1]

  index1<- as.integer(index1)
  index2<- as.integer(index2)
  index3<- as.integer(index3)

  if(!is.na(as.integer(index1))== FALSE){
    times_sup11<- times[length(times)]
  }else{
    times_sup11<- times[index1]
  }
  if(!is.na(as.integer(index2))== FALSE){
    times_sup21<- times[length(times)]
  }else{
    times_sup21<- times[index2]
  }

  if(!is.na(as.integer(index3))== FALSE){
    times_sup31<- times[length(times)]
  }else{
    times_sup31<- times[index3]
  }


  times_new7<-vector(mode = "numeric", length = 0)
  times_new8<-vector(mode = "numeric", length = 0)
  times_new9<-vector(mode = "numeric", length = 0)

  for (i in 2: length(times)){

    if(( temp_cmax[1]-T1[i-1])>=0 && ( temp_cmax[1]-T1[i])<0){
      times_new7[i-1]<- times[i-1]

    }else if(( temp_cmax[1]-T1[i-1])<=0 && ( temp_cmax[1]-T1[i])>0){

      times_new7[i-1]<- times[i-1]
    }else{
      times_new7[i-1]<- 0
    }
  }


  for (i in 2: length(times)){

    if(( temp_cmax[2]-T2[i-1])>=0 && ( temp_cmax[2]-T2[i])<0){
      times_new8[i-1]<- times[i-1]

    }else if(( temp_cmax[2]-T2[i-1])<=0 && ( temp_cmax[2]-T2[i])>0){

      times_new8[i-1]<- times[i-1]
    }else{
      times_new8[i-1]<- 0
    }
  }


  for (i in 2: length(times)){

    if(( temp_cmax[3]-T3[i-1])>=0 && ( temp_cmax[3]-T3[i])<0){
      times_new9[i-1]<- times[i-1]

    }else if(( temp_cmax[3]-T3[i-1])<=0 && ( temp_cmax[3]-T3[i])>0){

      times_new9[i-1]<- times[i-1]
    }else{
      times_new9[i-1]<- 0
    }
  }

  index7<- which(times_new7!=0)[1]
  index8<- which(times_new8!=0)[1]
  index9<- which(times_new9!=0)[1]

  index7<- as.integer(index7)
  index8<- as.integer(index8)
  index9<- as.integer(index9)

  if(!is.na(as.integer(index7))== FALSE){
    times_sup12<- times[length(times)]
  }else{
    times_sup12<- times[index7]
  }
  if(!is.na(as.integer(index8))== FALSE){
    times_sup22<- times[length(times)]
  }else{
    times_sup22<- times[index8]
  }

  if(!is.na(as.integer(index9))== FALSE){
    times_sup32<- times[length(times)]
  }else{
    times_sup32<- times[index9]
  }


  if(times_sup11<= times_sup12){
    times_sup1<-times_sup11
  }else{
    times_sup1<-times_sup12
  }

  if(times_sup21<= times_sup22){
    times_sup2<-times_sup21
  }else{
    times_sup2<-times_sup22
  }

  if(times_sup31<= times_sup32){
    times_sup3<-times_sup31
  }else{
    times_sup3<-times_sup32
  }


  ###############################################################
  # Carrying capacity
  ##############################################################

  times_new4<-vector(mode = "numeric", length = 0)
  times_new5<-vector(mode = "numeric", length = 0)
  times_new6<-vector(mode = "numeric", length = 0)

  for (i in 2: length(times)){

    if(K1[i-1]>=0 && ( K1[i])<0){
      times_new4[i-1]<- times[i-1]
    }else{
      times_new4[i-1]<- 0
    }
  }


  for (i in 2: length(times)){

    if(K2[i-1]>=0 && ( K2[i])<0){
      times_new5[i-1]<- times[i-1]
    }else{
      times_new5[i-1]<- 0
    }
  }


  for (i in 2: length(times)){

    if(K3[i-1]>=0 && ( K3[i])<0){
      times_new6[i-1]<- times[i-1]
    }else{
      times_new6[i-1]<- 0
    }
  }

  index4<- which(times_new4!=0)[1]
  index5<- which(times_new5!=0)[1]
  index6<- which(times_new6!=0)[1]

  index4<- as.integer(index4)
  index5<- as.integer(index5)
  index6<- as.integer(index6)

  if(!is.na(as.integer(index4))== FALSE){
    times_sup4<- times[length(times)]
  }else{
    times_sup4<- times[index4]
  }
  if(!is.na(as.integer(index5))== FALSE){
    times_sup5<- times[length(times)]
  }else{
    times_sup5<- times[index5]
  }

  if(!is.na(as.integer(index6))== FALSE){
    times_sup6<- times[length(times)]
  }else{
    times_sup6<- times[index6]
  }


###############################################################
  # Plots
##############################################################

data<-rbind(data1, data2, cap1)

p1<- ggplot(data, aes(x=.data$x, y=.data$y)) +
theme_bw()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_ribbon(data=subset(cap1,times>times[1] & times<times_sup4),aes(x=.data$x,
                                ymax=.data$y),ymin=0,alpha=0.3, fill="brown") +
geom_vline(xintercept = times_sup1, size=.5, color="brown",linetype="dashed")+
geom_line(data =subset(data1,times>times[1] & times<times_sup1), color = "brown")+
geom_line(data =subset(data2,times>times[1] & times<times_sup1), color = "green4")+
labs(x = "Time",y="Abundance")+
theme(plot.title = element_text(size=40))+
theme(plot.title = element_text(hjust = 0.5))+
theme(axis.title.y = element_text(size = rel(1), angle = 90))+
theme(axis.title.x = element_text(size = rel(1), angle = 00))+
labs(tag = "(a)")


dat<-rbind(dat1, dat2, cap2)

p2<- ggplot(dat, aes(x=.data$x, y=.data$y)) +
theme_bw()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_ribbon(data=subset(cap2,times>times[1] & times<times_sup5),aes(x=.data$x,
                               ymax=.data$y),ymin=0,alpha=0.3, fill="brown") +
geom_vline(xintercept = times_sup2, size=.5, color="brown",linetype="dashed")+
geom_line(data =subset(dat1,times>times[1] & times<times_sup2), color = "brown")+
geom_line(data =subset(dat2,times>times[1] & times<times_sup2), color = "green4")+
labs(x = "Time",y="Abundance")+
theme(plot.title = element_text(size=40))+
theme(plot.title = element_text(hjust = 0.5))+
theme(axis.title.y = element_text(size = rel(1), angle = 90))+
theme(axis.title.x = element_text(size = rel(1), angle = 00))+
labs(tag = "(b)")


da<-rbind(da1, da2, cap3)

p3<- ggplot(da, aes(x=.data$x, y=.data$y)) +
theme_bw()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_ribbon(data=subset(cap3,times>times[1] & times<times_sup6),aes(x=.data$x,ymax=.data$y),ymin=0,alpha=0.3, fill="brown") +
geom_vline(xintercept = times_sup3, size=.5, color="brown",linetype="dashed")+
geom_line(data =subset(da1,times>times[1] & times<times_sup3), color = "brown")+
geom_line(data =subset(da2,times>times[1] & times<times_sup3), color = "green4")+
labs(x = "Time",y="Abundance")+
theme(plot.title = element_text(size=40))+
theme(plot.title = element_text(hjust = 0.5))+
theme(axis.title.y = element_text(size = rel(1), angle = 90))+
theme(axis.title.x = element_text(size = rel(1), angle = 00))+
labs(tag = "(c)")


  d<-rbind(d1, d2, d3)

  p4<- ggplot(d, aes(x=.data$x, y=.data$y)) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    geom_vline(xintercept = times_sup1, size=.5, color="green",linetype="dashed")+
    geom_vline(xintercept = times_sup2, size=.5, color="blue",linetype="dashed")+
    geom_vline(xintercept = times_sup3, size=.5, color="black",linetype="dashed")+
    geom_line(data =subset(d1,times>times[1] & times<times_sup1), color = "green")+
    geom_line(data =subset(d2,times>times[1] & times<times_sup2), color = "blue")+
    geom_line(data =subset(d3,times>times[1] & times<times_sup3), color = "black")+
    labs(x = "Time",y="Temperature")+
    theme(plot.title = element_text(size=40))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.title.y = element_text(size = rel(1), angle = 90))+
    theme(axis.title.x = element_text(size = rel(1), angle = 00))+
    labs(tag = "(d)")

  plot_grid(p1, p2,p3,p4)

} else {

  stop("The available trends are associated with the RCP2.6 and RCP8.5 scenarios.")


}

}else{
  stop("The initial study temperature must be within the thermal tolerance range")
}

}else{

  stop("The minimum critical temperature must be less than the maximum critical temperature")
}

  }else{

    stop("time_start must be less than time_end ")
  }

}else{

  stop("The maximum simulation time is the year 2100 ")
}

}




