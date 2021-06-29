#'Age structure under IPCC RCP2.6 or RCP8.5 scenarios
#'
#'@description This function allows to simulate the effect of IPCC (2014) RCP2.6
#'             or RCP8.5 scenarios on the abundances of an ectotherm population
#'             considering three stages as age structure.
#'
#'
#'
#'@param y_ini Initial population values (must be written with its name: N).
#'@param temp_ini Initial temperature (K).
#'@param temp_cmin Minimum critical temperature (K).
#'@param temp_cmax Maximum critical temperature (K).
#'@param ro Population growth rate at optimum temperature.
#'@param lambda1,lambda2,lambda3 Marginal loss by non-thermodependent intraspecific competition.
#'@param alpha1,alpha2 Stage 1 to 2 and 2 to 3 transition coefficient respectively.
#'@param d2,d3 Mortality rate at a reference temperature.
#'@param Ad2,Ad3 Arrhenius constant which quantifies the temperature sensitivity of mortality.
#'@param Tr2,Tr3 Reference temperature (K).
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
#'@return (2) Figure of four panels in which (a), (b) and (c) show the population abundance curves for
#'            each age stage. Panel (d) shows the temperature trend curves used for each simulation,
#'            colored brown, green and blue, respectively.
#'
#'@references IPCC. (2014): Climate Change 2014: Synthesis Report. Contribution of Working Groups I,
#'            II and III to the Fifth Assessment Report of the Intergovernmental Panel on Climate
#'            Change [Core Writing Team, R.K. Pachauri and L.A. Meyer (eds.)]. IPCC, Geneva,
#'            Switzerland, 151 pp.
#'@export
#'@examples
#'\donttest{
#'#######################################################################
#'   #Example 1: Different thermal tolerance ranges (scenario RCP2.6).
#'#######################################################################
#'
#'temp_cmin <- 291
#'
#'# Temperature that occurs before the minimum simulation time.
#'temp_i <- 295
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
#'age_structure(y_ini = c(N1 = 800, N1 = 800, N1 = 800,
#'                        N2 = 600, N2 = 600, N2 = 600,
#'                        N3 = 400, N3 = 400, N3 = 400),
#'              temp_ini = rep(temp_ini,3),
#'              temp_cmin = rep(temp_cmin,3),
#'              temp_cmax = c(temp_cmax1,temp_cmax2,temp_cmax3),
#'              ro = rep(0.7,3),
#'              lambda1 = c(0.00002,0,0),
#'              lambda2 = c(0,0.00004,0.00003),
#'              lambda3 = c(0,0.00003,0.00004),
#'              alpha1 = rep(0.3,3),
#'              alpha2 = rep(0.4,3),
#'              d2 = rep(0.004,3),
#'              d3 = rep(0.005,3),
#'              Ad2 = rep(0.5,3),
#'              Ad3 = rep(0.75,3),
#'              Tr2 = rep(298,3),
#'              Tr3 = rep(298,3),
#'              RCP = 2.6,
#'              time_start = 2005,
#'              time_end = time_end,
#'              leap = 1/50)
#'
#'#######################################################################
#'   #Example 2: Different thermal tolerance ranges (scenario RCP8.5).
#'#######################################################################
#'
#'temp_cmin <- 291
#'
#'# Temperature that occurs before the minimum simulation time.
#'temp_i <- 295
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
#'age_structure(y_ini = c(N1 = 800, N1 = 800, N1 = 800,
#'                        N2 = 600, N2 = 600, N2 = 600,
#'                        N3 = 400, N3 = 400, N3 = 400),
#'              temp_ini = rep(temp_ini,3),
#'              temp_cmin = rep(temp_cmin,3),
#'              temp_cmax = c(temp_cmax1,temp_cmax2,temp_cmax3),
#'              ro = rep(0.7,3),
#'              lambda1 = c(0.00002,0,0),
#'              lambda2 = c(0,0.00004,0.00003),
#'              lambda3 = c(0,0.00003,0.00004),
#'              alpha1 = rep(0.3,3),
#'              alpha2 = rep(0.4,3),
#'              d2 = rep(0.004,3),
#'              d3 = rep(0.003,3),
#'              Ad2 = rep(0.5,3),
#'              Ad3 = rep(0.6,3),
#'              Tr2 = rep(298,3),
#'              Tr3 = rep(298,3),
#'              RCP = 8.5,
#'              time_start = 2005,
#'              time_end = time_end,
#'              leap = 1/50)
#'}




age_structure<- function(y_ini = c(N1 = 800, N1 = 800, N1 = 800,
                        N2 = 600, N2 = 600, N2 = 600,
                        N3 = 400, N3 = 400, N3 = 400),
              temp_ini = rep(25+273.15,3),
              temp_cmin = rep(18+273.15,3),
              temp_cmax = c(25+273.15,28+273.15,35+273.15),
              ro = rep(0.7,3),
              lambda1 = rep(0.0004,3),
              lambda2 = rep(0.0004,3),
              lambda3 = rep(0.0004,3),
              alpha1 = rep(0.1,3),
              alpha2 = rep(0.7,3),
              d2 = rep(0.005,3),
              d3 = rep(0.5,3),
              Ad2 = rep(0.5,3),
              Ad3 = rep(0.75,3),
              Tr2 = rep(298,3),
              Tr3 = rep(298,3),
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
          # Optimum growing temperature
##########################################################

temp_op1<- (temp_cmax[1]+temp_cmin[1])/3+sqrt(((temp_cmax[1]+
            temp_cmin[1])/3)^2-(temp_cmax[1]*temp_cmin[1])/3)
temp_op2<- (temp_cmax[2]+temp_cmin[2])/3+sqrt(((temp_cmax[2]+
            temp_cmin[2])/3)^2-(temp_cmax[2]*temp_cmin[2])/3)
temp_op3<- (temp_cmax[3]+temp_cmin[3])/3+sqrt(((temp_cmax[3]+
            temp_cmin[3])/3)^2-(temp_cmax[3]*temp_cmin[3])/3)


##########################################################
if(RCP==2.6) {

  ##########################################################
  # Parameters
  ##########################################################
  parms1<-c(temp_cmin[1],temp_ini[1],temp_cmax[1],temp_op1,ro[1],
            lambda1[1],lambda1[2],lambda1[3],lambda2[1],lambda2[2],
            lambda2[3],lambda3[1],lambda3[2],lambda3[3],alpha1[1],
            alpha2[1],d2[1],d3[1])
  parms2<-c(temp_cmin[2],temp_ini[2],temp_cmax[2],temp_op2,ro[2],
            lambda1[1],lambda1[2],lambda1[3],lambda2[1],lambda2[2],
            lambda2[3],lambda3[1],lambda3[2],lambda3[3],alpha1[2],
            alpha2[2],d2[2],d3[2])
  parms3<-c(temp_cmin[3],temp_ini[3],temp_cmax[3],temp_op3,ro[3],
            lambda1[1],lambda1[2],lambda1[3],lambda2[1],lambda2[2],
            lambda2[3],lambda3[1],lambda3[2],lambda3[3],alpha1[3],
            alpha2[3],d2[3],d3[3])
##########################################################
            # Model for each trend
##########################################################
model1 <- function (times, y,parms1) {
    with(as.list(c(y)), {
  T1<- get_RCP2.6(times)+temp_ini[1]  # IPCC1
  r1<- rate_TPC(T1,ro[1],temp_cmin[1],temp_cmax[1],temp_op1)
  D2<- d2[1]*exp(Ad2[1]*(1/Tr2[1]-1/T1))
  D3<- d3[1]*exp(Ad3[1]*(1/Tr3[1]-1/T1))
  dN1 <-  r1*N3-lambda1[1]*N1*N1-lambda1[2]*N1*N2-lambda1[3]*N1*N3-alpha1[1]*N1
  dN2 <-  alpha1[1]*N1-D2*N2-lambda2[1]*N2*N1-lambda2[2]*N2*N2-lambda2[3]*N2*N3-alpha2[1]*N2
  dN3 <-  alpha2[1]*N2-D3*N3-lambda3[1]*N3*N1-lambda3[2]*N3*N2-lambda3[3]*N3*N3

return(list(c(dN1,dN2,dN3)))
  })
}
###############################################################
model2 <- function (times, y,parms2) {
   with(as.list(c(y)), {
  T2 <- get_RCP2.6(times)+temp_ini[2]  # IPCC1
  r2<- rate_TPC(T2,ro[2],temp_cmin[2],temp_cmax[2],temp_op2)
  D2<- d2[2]*exp(Ad2[2]*(1/Tr2[2]-1/T2))
  D3<- d3[2]*exp(Ad3[2]*(1/Tr3[2]-1/T2))
  dN1 <-  r2*N3-lambda1[1]*N1*N1-lambda1[2]*N1*N2-lambda1[3]*N1*N3-alpha1[2]*N1
  dN2 <-  alpha1[2]*N1-D2*N2-lambda2[1]*N2*N1-lambda2[2]*N2*N2-lambda2[3]*N2*N3-alpha2[2]*N2
  dN3 <-  alpha2[2]*N2-D3*N3-lambda3[1]*N3*N1-lambda3[2]*N3*N2-lambda3[3]*N3*N3

return(list(c(dN1,dN2,dN3)))
    })
}
###############################################################
model3 <- function (times, y,parms3) {
   with(as.list(c(y)), {
 T3 <- get_RCP2.6(times)+temp_ini[3]  # IPCC1
 r3<- rate_TPC(T3,ro[3],temp_cmin[3],temp_cmax[3],temp_op3)
 D2<- d2[3]*exp(Ad2[3]*(1/Tr2[3]-1/T3))
 D3<- d3[3]*exp(Ad3[3]*(1/Tr3[3]-1/T3))
 dN1 <-  r3*N3-lambda1[1]*N1*N1-lambda1[2]*N1*N2-lambda1[3]*N1*N3-alpha1[3]*N1
 dN2 <-  alpha1[3]*N1-D2*N2-lambda2[1]*N2*N1-lambda2[2]*N2*N2-lambda2[3]*N2*N3-alpha2[3]*N2
 dN3 <-  alpha2[3]*N2-D3*N3-lambda3[1]*N3*N1-lambda3[2]*N3*N2-lambda3[3]*N3*N3

return(list(c(dN1,dN2,dN3)))
  })
}
###############################################################
y_ini1<-c(y_ini[1],y_ini[4],y_ini[7])
y_ini2<-c(y_ini[2],y_ini[5],y_ini[8])
y_ini3<-c(y_ini[3],y_ini[6],y_ini[9])
###############################################################
            # Solution
##############################################################
out1 <- ode(y=y_ini1, times, model1, parms1, method = "ode45")
out2 <- ode(y=y_ini2, times, model2, parms2, method = "ode45")
out3 <- ode(y=y_ini3, times, model3, parms3, method = "ode45")
#############################################################
###############################################################
            # Abundance
##############################################################
data1<-data.frame('x'=times,'y'=out1[,2])
data2<-data.frame('x'=times,'y'=out1[,3])
data3<-data.frame('x'=times,'y'=out1[,4])

dat1<-data.frame('x'=times,'y'=out2[,2])
dat2<-data.frame('x'=times,'y'=out2[,3])
dat3<-data.frame('x'=times,'y'=out2[,4])

da1<-data.frame('x'=times,'y'=out3[,2])
da2<-data.frame('x'=times,'y'=out3[,3])
da3<-data.frame('x'=times,'y'=out3[,4])

T1 <- get_RCP2.6(times)+temp_ini[1]
T2 <- get_RCP2.6(times)+temp_ini[2]
T3 <- get_RCP2.6(times)+temp_ini[3]

d1<-data.frame('x'=times,'y'=T1)
d2<-data.frame('x'=times,'y'=T2)
d3<-data.frame('x'=times,'y'=T3)

r1<- rate_TPC(T1,ro[1],temp_cmin[1],temp_cmax[1],temp_op1)
r2<- rate_TPC(T2,ro[2],temp_cmin[2],temp_cmax[2],temp_op2)
r3<- rate_TPC(T3,ro[3],temp_cmin[3],temp_cmax[3],temp_op3)

###############################################################
            # Data
###############################################################
Data<- data.frame(times,out1[,2],out1[,3], out1[,4],out2[,2],out2[,3],
                  out2[,4],out3[,2],out3[,3],out3[,4])
names(Data)<- c("Time","Abundance Stage-1,Scenario-1","Abundance Stage-2,Scenario-1",
                "Abundance Stage-3,Scenario-1","Abundance Stage-1,Scenario-2",
                "Abundance Stage-2,Scenario-2","Abundance Stage-3,Scenario-2",
                "Abundance Stage-1,Scenario-3","Abundance Stage-2,Scenario-3",
                "Abundance Stage-3,Scenario-3")
        u<- formattable(Data, align = c("l", rep("r", NCOL(Data))))
print(u)

times_new11<-vector(mode = "numeric", length = 0)
times_new12<-vector(mode = "numeric", length = 0)
times_new13<-vector(mode = "numeric", length = 0)
times_new21<-vector(mode = "numeric", length = 0)
times_new22<-vector(mode = "numeric", length = 0)
times_new23<-vector(mode = "numeric", length = 0)
times_new31<-vector(mode = "numeric", length = 0)
times_new32<-vector(mode = "numeric", length = 0)
times_new33<-vector(mode = "numeric", length = 0)

for (i in 2: length(times)){

  if(is.na(out1[i,2])) {times_new11[i-1]<- 0
  } else {
    if(out1[i-1,2]>=0 && ( out1[i,2]<0 || is.na(out1[i,2]))){
      times_new11[i-1]<- times[i-1]
    }else{
      times_new11[i-1]<- 0
    }
  }
}


for (i in 2: length(times)){
  if(is.na(out1[i,3])) {times_new12[i-1]<- 0
  } else {
    if(out1[i-1,3]>=0 && ( out1[i,3]<0 || is.na(out1[i,3]))){
      times_new12[i-1]<- times[i-1]
    }else{
      times_new12[i-1]<- 0
    }
  }
}


for (i in 2: length(times)){
  if(is.na(out1[i,4])) {times_new13[i-1]<- 0
  } else {
    if(out1[i-1,4]>=0 && ( out1[i,4]<0 || is.na(out1[i,4]))){
      times_new13[i-1]<- times[i-1]
    }else{
      times_new13[i-1]<- 0
    }
  }

}
index11<- which(times_new11!=0)[1]
index12<- which(times_new12!=0)[1]
index13<- which(times_new13!=0)[1]

index11<- as.integer(index11)
index11<- as.integer(index12)
index13<- as.integer(index13)

if(!is.na(as.integer(index11))== FALSE){
  times_super11<- times[length(times)]
}else{
  times_super11<- times[index11]
}
if(!is.na(as.integer(index12))== FALSE){
  times_super12<- times[length(times)]
}else{
  times_super12<- times[index12]
}

if(!is.na(as.integer(index13))== FALSE){
  times_super13<- times[length(times)]
}else{
  times_super13<- times[index13]
}
###################################################################

for (i in 2: length(times)){

  if(is.na(out2[i,2])) {times_new21[i-1]<- 0
  } else {
    if(out2[i-1,2]>=0 && ( out2[i,2]<0 || is.na(out2[i,2]))){
      times_new21[i-1]<- times[i-1]
    }else{
      times_new21[i-1]<- 0
    }
  }
}


for (i in 2: length(times)){
  if(is.na(out2[i,3])) {times_new22[i-1]<- 0
  } else {
    if(out2[i-1,3]>=0 && ( out2[i,3]<0 || is.na(out2[i,3]))){
      times_new22[i-1]<- times[i-1]
    }else{
      times_new22[i-1]<- 0
    }
  }
}


for (i in 2: length(times)){
  if(is.na(out2[i,4])) {times_new23[i-1]<- 0
  } else {
    if(out2[i-1,4]>=0 && ( out2[i,4]<0 || is.na(out2[i,4]))){
      times_new23[i-1]<- times[i-1]
    }else{
      times_new23[i-1]<- 0
    }
  }

}
index21<- which(times_new21!=0)[1]
index22<- which(times_new22!=0)[1]
index23<- which(times_new23!=0)[1]

index21<- as.integer(index21)
index22<- as.integer(index22)
index23<- as.integer(index23)

if(!is.na(as.integer(index21))== FALSE){
  times_super21<- times[length(times)]
}else{
  times_super21<- times[index21]
}
if(!is.na(as.integer(index22))== FALSE){
  times_super22<- times[length(times)]
}else{
  times_super22<- times[index22]
}

if(!is.na(as.integer(index23))== FALSE){
  times_super23<- times[length(times)]
}else{
  times_super23<- times[index23]
}
##################################################################


for (i in 2: length(times)){

  if(is.na(out3[i,2])) {times_new31[i-1]<- 0
  } else {
    if(out3[i-1,2]>=0 && ( out3[i,2]<0 || is.na(out3[i,2]))){
      times_new31[i-1]<- times[i-1]
    }else{
      times_new31[i-1]<- 0
    }
  }
}


for (i in 2: length(times)){
  if(is.na(out3[i,3])) {times_new32[i-1]<- 0
  } else {
    if(out3[i-1,3]>=0 && ( out3[i,3]<0 || is.na(out3[i,3]))){
      times_new32[i-1]<- times[i-1]
    }else{
      times_new32[i-1]<- 0
    }
  }
}


for (i in 2: length(times)){
  if(is.na(out3[i,4])) {times_new33[i-1]<- 0
  } else {
    if(out3[i-1,4]>=0 && ( out3[i,4]<0 || is.na(out3[i,4]))){
      times_new33[i-1]<- times[i-1]
    }else{
      times_new33[i-1]<- 0
    }
  }

}

index31<- which(times_new31!=0)[1]
index32<- which(times_new32!=0)[1]
index33<- which(times_new33!=0)[1]

index31<- as.integer(index31)
index32<- as.integer(index32)
index33<- as.integer(index33)

if(!is.na(as.integer(index31))== FALSE){
  times_super31<- times[length(times)]
}else{
  times_super31<- times[index31]
}
if(!is.na(as.integer(index32))== FALSE){
  times_super32<- times[length(times)]
}else{
  times_super32<- times[index32]
}

if(!is.na(as.integer(index33))== FALSE){
  times_super33<- times[length(times)]
}else{
  times_super33<- times[index33]
}
##################################################################

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
  times_sup1<- times[length(times)]
}else{
  times_sup2<- times[index7]
}
if(!is.na(as.integer(index8))== FALSE){
  times_sup2<- times[length(times)]
}else{
  times_sup2<- times[index8]
}

if(!is.na(as.integer(index9))== FALSE){
  times_sup3<- times[length(times)]
}else{
  times_sup3<- times[index9]
}




###############################################################
            # Plots
##############################################################

data<-rbind(data1, data2, data3)

p1<- ggplot(data, aes(x=.data$x, y=.data$y)) +
theme_bw()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_vline(xintercept = times_super11, size=.5, color="brown",linetype="dashed")+
geom_vline(xintercept = times_super12, size=.5, color="green4",linetype="dashed")+
geom_vline(xintercept = times_super13, size=.5, color="blue",linetype="dashed")+
geom_line(data =subset(data1,times>times[1] & times<times_super11), color = "brown")+
geom_line(data =subset(data2,times>times[1] & times<times_super12), color = "green4")+
geom_line(data =subset(data3,times>times[1] & times<times_super13), color = "blue")+
labs(x = "Time",y="Abundance")+
theme(plot.title = element_text(size=40))+
theme(plot.title = element_text(hjust = 0.5))+
theme(axis.title.y = element_text(size = rel(1), angle = 90))+
theme(axis.title.x = element_text(size = rel(1), angle = 00))+
labs(tag = "(a)")


dat<-rbind(dat1, dat2, dat3)

p2<- ggplot(dat, aes(x=.data$x, y=.data$y)) +
theme_bw()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_vline(xintercept = times_super21, size=.5, color="brown",linetype="dashed")+
geom_vline(xintercept = times_super22, size=.5, color="green4",linetype="dashed")+
geom_vline(xintercept = times_super23, size=.5, color="blue",linetype="dashed")+
geom_line(data =subset(dat1,times>times[1] & times<times_super21), color = "brown")+
geom_line(data =subset(dat2,times>times[1] & times<times_super22), color = "green4")+
geom_line(data =subset(dat3,times>times[1] & times<times_super23), color = "blue")+
labs(x = "Time",y="Abundance")+
theme(plot.title = element_text(size=40))+
theme(plot.title = element_text(hjust = 0.5))+
theme(axis.title.y = element_text(size = rel(1), angle = 90))+
theme(axis.title.x = element_text(size = rel(1), angle = 00))+
labs(tag = "(b)")


da<-rbind(da1, da2, da3)

p3<- ggplot(da, aes(x=.data$x, y=.data$y)) +
theme_bw()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_vline(xintercept = times_super31, size=.5, color="brown",linetype="dashed")+
geom_vline(xintercept = times_super32, size=.5, color="green4",linetype="dashed")+
geom_vline(xintercept = times_super33, size=.5, color="blue",linetype="dashed")+
geom_line(data =subset(da1,times>times[1] & times<times_super31), color = "brown")+
geom_line(data =subset(da2,times>times[1] & times<times_super32), color = "green4")+
geom_line(data =subset(da3,times>times[1] & times<times_super33), color = "blue")+
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
geom_vline(xintercept = times_sup1, size=.5, color="brown",linetype="dashed")+
geom_vline(xintercept = times_sup2, size=.5, color="green4",linetype="dashed")+
geom_vline(xintercept = times_sup3, size=.5, color="blue",linetype="dashed")+
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
  # Parameters
  ##########################################################
  parms1<-c(temp_cmin[1],temp_ini[1],temp_cmax[1],temp_op1,ro[1],
            lambda1[1],lambda1[2],lambda1[3],lambda2[1],lambda2[2],
            lambda2[3],lambda3[1],lambda3[2],lambda3[3],alpha1[1],
            alpha2[1],d2[1],d3[1],exp(coef(m)[1]),coef(m)[2])
  parms2<-c(temp_cmin[2],temp_ini[2],temp_cmax[2],temp_op2,ro[2],
            lambda1[1],lambda1[2],lambda1[3],lambda2[1],lambda2[2],
            lambda2[3],lambda3[1],lambda3[2],lambda3[3],alpha1[2],
            alpha2[2],d2[2],d3[2],exp(coef(m)[1]),coef(m)[2])
  parms3<-c(temp_cmin[3],temp_ini[3],temp_cmax[3],temp_op3,ro[3],
            lambda1[1],lambda1[2],lambda1[3],lambda2[1],lambda2[2],
            lambda2[3],lambda3[1],lambda3[2],lambda3[3],alpha1[3],
            alpha2[3],d2[3],d3[3],exp(coef(m)[1]),coef(m)[2])

##########################################################
            # Model for each trend
##########################################################
model1 <- function (times, y,parms1) {
    with(as.list(c(y)), {
 T1<- RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[1]    #IPCC2
 r1<- rate_TPC(T1,ro[1],temp_cmin[1],temp_cmax[1],temp_op1)
 D2<- d2[1]*exp(Ad2[1]*(1/Tr2[1]-1/T1))
 D3<- d3[1]*exp(Ad3[1]*(1/Tr3[1]-1/T1))
 dN1 <-  r1*N3-lambda1[1]*N1*N1-lambda1[2]*N1*N2-lambda1[3]*N1*N3-alpha1[1]*N1
 dN2 <-  alpha1[1]*N1-D2*N2-lambda2[1]*N2*N1-lambda2[2]*N2*N2-lambda2[3]*N2*N3-alpha2[1]*N2
 dN3 <-  alpha2[1]*N2-D3*N3-lambda3[1]*N3*N1-lambda3[2]*N3*N2-lambda3[3]*N3*N3
return(list(c(dN1,dN2,dN3)))
      })
  }
###############################################################
model2 <- function (times, y,parms2) {
    with(as.list(c(y)), {
 T2<-  RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[2]    #IPCC2
 r2<- rate_TPC(T2,ro[2],temp_cmin[2],temp_cmax[2],temp_op2)
 D2<- d2[2]*exp(Ad2[2]*(1/Tr2[2]-1/T2))
 D3<- d3[2]*exp(Ad3[2]*(1/Tr3[2]-1/T2))
 dN1 <-  r2*N3-lambda1[1]*N1*N1-lambda1[2]*N1*N2-lambda1[3]*N1*N3-alpha1[2]*N1
 dN2 <-  alpha1[2]*N1-D2*N2-lambda2[1]*N2*N1-lambda2[2]*N2*N2-lambda2[3]*N2*N3-alpha2[2]*N2
 dN3 <-  alpha2[2]*N2-D3*N3-lambda3[1]*N3*N1-lambda3[2]*N3*N2-lambda3[3]*N3*N3
return(list(c(dN1,dN2,dN3)))
    })
}
###############################################################
model3 <- function (times, y,parms3) {
    with(as.list(c(y)), {
 T3<-  RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[3]    #IPCC2
 r3<- rate_TPC(T3,ro[3],temp_cmin[3],temp_cmax[3],temp_op3)
 D2<- d2[3]*exp(Ad2[3]*(1/Tr2[3]-1/T3))
 D3<- d3[3]*exp(Ad3[3]*(1/Tr3[3]-1/T3))
 dN1 <-  r3*N3-lambda1[1]*N1*N1-lambda1[2]*N1*N2-lambda1[3]*N1*N3-alpha1[3]*N1
 dN2 <-  alpha1[3]*N1-D2*N2-lambda2[1]*N2*N1-lambda2[2]*N2*N2-lambda2[3]*N2*N3-alpha2[3]*N2
 dN3 <-  alpha2[3]*N2-D3*N3-lambda3[1]*N3*N1-lambda3[2]*N3*N2-lambda3[3]*N3*N3

 return(list(c(dN1,dN2,dN3)))
    })
}
###############################################################
y_ini1<-c(y_ini[1],y_ini[4],y_ini[7])
y_ini2<-c(y_ini[2],y_ini[5],y_ini[8])
y_ini3<-c(y_ini[3],y_ini[6],y_ini[9])
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
data3<-data.frame('x'=times,'y'=out1[,4] )

dat1<-data.frame('x'=times,'y'=out2[,2] )
dat2<-data.frame('x'=times,'y'=out2[,3] )
dat3<-data.frame('x'=times,'y'=out2[,4] )

da1<-data.frame('x'=times,'y'=out3[,2] )
da2<-data.frame('x'=times,'y'=out3[,3] )
da3<-data.frame('x'=times,'y'=out3[,4] )

T1<-  RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[1]
T2<-  RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[2]
T3<-  RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[3]

d1<-data.frame('x'=times,'y'=T1)
d2<-data.frame('x'=times,'y'=T2)
d3<-data.frame('x'=times,'y'=T3)

r1<- rate_TPC(T1,ro[1],temp_cmin[1],temp_cmax[1],temp_op1)
r2<- rate_TPC(T2,ro[2],temp_cmin[2],temp_cmax[2],temp_op2)
r3<- rate_TPC(T3,ro[3],temp_cmin[3],temp_cmax[3],temp_op3)

###############################################################
            # Data
###############################################################
Data<- data.frame(times,out1[,2],out1[,3],out1[,4],out2[,2],out2[,3],
                  out2[,4],out3[,2],out3[,3],out3[,4])
names(Data)<- c("Time","Abundance Stage-1,Scenario-1","Abundance Stage-2,Scenario-1",
                "Abundance Stage-3,Scenario-1","Abundance Stage-1,Scenario-2",
                "Abundance Stage-2,Scenario-2","Abundance Stage-3,Scenario-2",
                "Abundance Stage-1,Scenario-3","Abundance Stage-2,Scenario-3",
                "Abundance Stage-3,Scenario-3")
    u<- formattable(Data, align = c("l", rep("r", NCOL(Data))))
print(u)
###############################################################
            times_new11<-vector(mode = "numeric", length = 0)
            times_new12<-vector(mode = "numeric", length = 0)
            times_new13<-vector(mode = "numeric", length = 0)
            times_new21<-vector(mode = "numeric", length = 0)
            times_new22<-vector(mode = "numeric", length = 0)
            times_new23<-vector(mode = "numeric", length = 0)
            times_new31<-vector(mode = "numeric", length = 0)
            times_new32<-vector(mode = "numeric", length = 0)
            times_new33<-vector(mode = "numeric", length = 0)

            for (i in 2: length(times)){

              if(is.na(out1[i,2])) {times_new11[i-1]<- 0
              } else {
                if(out1[i-1,2]>=0 && ( out1[i,2]<0 || is.na(out1[i,2]))){
                  times_new11[i-1]<- times[i-1]
                }else{
                  times_new11[i-1]<- 0
                }
              }
            }


            for (i in 2: length(times)){
              if(is.na(out1[i,3])) {times_new12[i-1]<- 0
              } else {
                if(out1[i-1,3]>=0 && ( out1[i,3]<0 || is.na(out1[i,3]))){
                  times_new12[i-1]<- times[i-1]
                }else{
                  times_new12[i-1]<- 0
                }
              }
            }


            for (i in 2: length(times)){
              if(is.na(out1[i,4])) {times_new13[i-1]<- 0
              } else {
                if(out1[i-1,4]>=0 && ( out1[i,4]<0 || is.na(out1[i,4]))){
                  times_new13[i-1]<- times[i-1]
                }else{
                  times_new13[i-1]<- 0
                }
              }

            }
            index11<- which(times_new11!=0)[1]
            index12<- which(times_new12!=0)[1]
            index13<- which(times_new13!=0)[1]

            index11<- as.integer(index11)
            index11<- as.integer(index12)
            index13<- as.integer(index13)

            if(!is.na(as.integer(index11))== FALSE){
              times_super11<- times[length(times)]
            }else{
              times_super11<- times[index11]
            }
            if(!is.na(as.integer(index12))== FALSE){
              times_super12<- times[length(times)]
            }else{
              times_super12<- times[index12]
            }

            if(!is.na(as.integer(index13))== FALSE){
              times_super13<- times[length(times)]
            }else{
              times_super13<- times[index13]
            }
###################################################################

            for (i in 2: length(times)){

              if(is.na(out2[i,2])) {times_new21[i-1]<- 0
              } else {
                if(out2[i-1,2]>=0 && ( out2[i,2]<0 || is.na(out2[i,2]))){
                  times_new21[i-1]<- times[i-1]
                }else{
                  times_new21[i-1]<- 0
                }
              }
            }


            for (i in 2: length(times)){
              if(is.na(out2[i,3])) {times_new22[i-1]<- 0
              } else {
                if(out2[i-1,3]>=0 && ( out2[i,3]<0 || is.na(out2[i,3]))){
                  times_new22[i-1]<- times[i-1]
                }else{
                  times_new22[i-1]<- 0
                }
              }
            }


            for (i in 2: length(times)){
              if(is.na(out2[i,4])) {times_new23[i-1]<- 0
              } else {
                if(out2[i-1,4]>=0 && ( out2[i,4]<0 || is.na(out2[i,4]))){
                  times_new23[i-1]<- times[i-1]
                }else{
                  times_new23[i-1]<- 0
                }
              }

            }
            index21<- which(times_new21!=0)[1]
            index22<- which(times_new22!=0)[1]
            index23<- which(times_new23!=0)[1]

            index21<- as.integer(index21)
            index22<- as.integer(index22)
            index23<- as.integer(index23)

            if(!is.na(as.integer(index21))== FALSE){
              times_super21<- times[length(times)]
            }else{
              times_super21<- times[index21]
            }
            if(!is.na(as.integer(index22))== FALSE){
              times_super22<- times[length(times)]
            }else{
              times_super22<- times[index22]
            }

            if(!is.na(as.integer(index23))== FALSE){
              times_super23<- times[length(times)]
            }else{
              times_super23<- times[index23]
            }
##################################################################


            for (i in 2: length(times)){

              if(is.na(out3[i,2])) {times_new31[i-1]<- 0
              } else {
                if(out3[i-1,2]>=0 && ( out3[i,2]<0 || is.na(out3[i,2]))){
                  times_new31[i-1]<- times[i-1]
                }else{
                  times_new31[i-1]<- 0
                }
              }
            }


            for (i in 2: length(times)){
              if(is.na(out3[i,3])) {times_new32[i-1]<- 0
              } else {
                if(out3[i-1,3]>=0 && ( out3[i,3]<0 || is.na(out3[i,3]))){
                  times_new32[i-1]<- times[i-1]
                }else{
                  times_new32[i-1]<- 0
                }
              }
            }


            for (i in 2: length(times)){
              if(is.na(out3[i,4])) {times_new33[i-1]<- 0
              } else {
                if(out3[i-1,4]>=0 && ( out3[i,4]<0 || is.na(out3[i,4]))){
                  times_new33[i-1]<- times[i-1]
                }else{
                  times_new33[i-1]<- 0
                }
              }

            }

            index31<- which(times_new31!=0)[1]
            index32<- which(times_new32!=0)[1]
            index33<- which(times_new33!=0)[1]

            index31<- as.integer(index31)
            index32<- as.integer(index32)
            index33<- as.integer(index33)

            if(!is.na(as.integer(index31))== FALSE){
              times_super31<- times[length(times)]
            }else{
              times_super31<- times[index31]
            }
            if(!is.na(as.integer(index32))== FALSE){
            times_super32<- times[length(times)]
            }else{
              times_super32<- times[index32]
            }

            if(!is.na(as.integer(index33))== FALSE){
              times_super33<- times[length(times)]
            }else{
              times_super33<- times[index33]
            }
##################################################################

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
              times_sup1<- times[length(times)]
            }else{
              times_sup2<- times[index7]
            }
            if(!is.na(as.integer(index8))== FALSE){
              times_sup2<- times[length(times)]
            }else{
              times_sup2<- times[index8]
            }

            if(!is.na(as.integer(index9))== FALSE){
              times_sup3<- times[length(times)]
            }else{
              times_sup3<- times[index9]
            }



###############################################################
            # Plots
##############################################################

data<-rbind(data1, data2, data3)

p1<- ggplot(data, aes(x=.data$x, y=.data$y)) +
theme_bw()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_vline(xintercept = times_super11, size=.5, color="brown",linetype="dashed")+
geom_vline(xintercept = times_super12, size=.5, color="green4",linetype="dashed")+
geom_vline(xintercept = times_super13, size=.5, color="blue",linetype="dashed")+
geom_line(data =subset(data1,times>times[1] & times<times_super11), color = "brown")+
geom_line(data =subset(data2,times>times[1] & times<times_super12), color = "green4")+
geom_line(data =subset(data3,times>times[1] & times<times_super13), color = "blue")+
labs(x = "Time",y="Abundance")+
theme(plot.title = element_text(size=40))+
theme(plot.title = element_text(hjust = 0.5))+
theme(axis.title.y = element_text(size = rel(1), angle = 90))+
theme(axis.title.x = element_text(size = rel(1), angle = 00))+
labs(tag = "(a)")



dat<-rbind(dat1, dat2, dat3)

p2<- ggplot(dat, aes(x=.data$x, y=.data$y)) +
theme_bw()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_vline(xintercept = times_super21, size=.5, color="brown",linetype="dashed")+
geom_vline(xintercept = times_super22, size=.5, color="green4",linetype="dashed")+
geom_vline(xintercept = times_super23, size=.5, color="blue",linetype="dashed")+
geom_line(data =subset(dat1,times>times[1] & times<times_super21), color = "brown")+
geom_line(data =subset(dat2,times>times[1] & times<times_super22), color = "green4")+
geom_line(data =subset(dat3,times>times[1] & times<times_super23), color = "blue")+
labs(x = "Time",y="Abundance")+
theme(plot.title = element_text(size=40))+
theme(plot.title = element_text(hjust = 0.5))+
theme(axis.title.y = element_text(size = rel(1), angle = 90))+
theme(axis.title.x = element_text(size = rel(1), angle = 00))+
labs(tag = "(b)")

da<-rbind(da1, da2, da3)

p3<- ggplot(da, aes(x=.data$x, y=.data$y)) +
theme_bw()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_vline(xintercept = times_super31, size=.5, color="brown",linetype="dashed")+
geom_vline(xintercept = times_super32, size=.5, color="green4",linetype="dashed")+
geom_vline(xintercept = times_super33, size=.5, color="blue",linetype="dashed")+
geom_line(data =subset(da1,times>times[1] & times<times_super31), color = "brown")+
geom_line(data =subset(da2,times>times[1] & times<times_super32), color = "green4")+
geom_line(data =subset(da3,times>times[1] & times<times_super33), color = "blue")+
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




