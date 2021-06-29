#'Cooling pulse-2
#'
#' @description This function allows simulating the effect of an environmental
#'              cooling pulse on the abundance of ectotherm populations. After
#'              the pulse, the temperature stabilizes at a temperature q units
#'              lower than the starting temperature (temp_ini).
#'
#'
#'@param y_ini Initial population values (must be written with its name: N).
#'@param temp_ini Initial temperature.
#'@param temp_cmin Minimum critical temperature.
#'@param temp_cmax Maximum critical temperature.
#'@param ro Population growth rate at optimum temperature.
#'@param lambda Marginal loss by non-thermodependent intraspecific competition.
#'@param temp_peak Peak pulse temperature.
#'@param time_peak Time when temp_peak is reached.
#'@param q Difference between temp_ini and the stabilized temperature.
#'@param time_start Start of time sequence.
#'@param time_end End of time sequence.
#'@param leap Time sequence step.
#'
#'
#'@details Three populations and/or scenarios can be simulated simultaneously.
#'         The temperature trend is determined by a Gaussian function, whose
#'         main parameters are the mean and standard deviation. In each input
#'         vector, the parameters for the three simulations must be specified
#'         (finite numbers for initial population abundance). The simulations
#'         are obtained by a model that incorporates the effects of temperature
#'         over time, which leads to a non-autonomous ODE approach. This is
#'         function uses the ODE solver implemented in the package deSolve
#'         (Soetaert et al., 2010).
#'
#'@return (1) A data.frame with columns having the simulated trends.
#'@return (2) A two-panel figure in which (a) shows the population abundance curves
#'            represented by solid lines and the corresponding carrying capacities
#'            are represented by shaded areas. In (b) the temperature trend is shown.
#'            The three simultaneous simulations are depicted by different colors, i.e.
#'            1st brown, 2nd green and 3rd blue.
#'
#'@references Soetaert, K., Petzoldt, T., & Setzer, R. (2010). Solving Differential
#'            Equations in R: Package deSolve. Journal of Statistical Software, 33(9),
#'            1 - 25. doi:http://dx.doi.org/10.18637/jss.v033.i09
#'
#'@export
#'@examples
#'
#'#######################################################################
#'   #Example 1: Different initial population abundances.
#'#######################################################################
#'
#' cooling_pulse2(y_ini = c(N = 100, N = 200, N = 400),
#'                temp_ini = rep(28,3),
#'                temp_cmin = rep(18,3),
#'                temp_cmax = rep(30,3),
#'                ro = rep(0.7,3),
#'                lambda = rep(0.00005,3),
#'                temp_peak = rep(19,3),
#'                time_peak = rep(2060,3),
#'                q = rep(1,3),
#'                time_start = 2005,
#'                time_end = 2100,
#'                leap = 1/12)
#'
#'#######################################################################
#'   #Example 2: Different thermal tolerance ranges.
#'#######################################################################
#'
#'temp_cmin3 <- 18
#'temp_cmin2 <- 10/9*temp_cmin3
#'temp_cmin1 <- 10/9*temp_cmin2
#'
#'temp_cmax1 <- 32.4
#'temp_cmax2 <- 10/9*temp_cmax1
#'temp_cmax3 <- 10/9*temp_cmax2
#'
#'cooling_pulse2(y_ini = c(N = 100, N = 100, N = 100),
#'               temp_ini = rep(30,3),
#'               temp_cmin = c(temp_cmin1,temp_cmin2,temp_cmin3),
#'               temp_cmax = c(temp_cmax1,temp_cmax2,temp_cmax3),
#'               ro = rep(0.7,3),
#'               lambda = rep(0.00005,3),
#'               temp_peak = rep(19,3),
#'               time_peak = rep(2060,3),
#'               q = rep(1,3),
#'               time_start = 2005,
#'               time_end = 2100,
#'               leap = 1/12)
#'\donttest{
#'#######################################################################
#'   #Example 3: Different relationships between initial environmental
#'   #           temperature and optimum temperature
#'#######################################################################
#'
#'temp_cmin <- 18
#'temp_cmax <- 30
#'
#'# Temperature at which performance is at its maximum value.
#'temp_op <- (temp_cmax+temp_cmin)/3+sqrt(((temp_cmax+temp_cmin)/3)^2-
#'            (temp_cmax*temp_cmin)/3)
#'
#'temp_ini1<- (temp_cmin+temp_op)/2
#'temp_ini2 <- temp_op
#'temp_ini3 <-(temp_op+temp_cmax)/2
#'
#'cooling_pulse2(y_ini = c(N = 100, N = 100, N = 100),
#'               temp_ini = c(temp_ini1,temp_ini2,temp_ini3),
#'               temp_cmin = rep(temp_cmin,3),
#'               temp_cmax = rep(temp_cmax,3) ,
#'               ro = rep(0.7,3),
#'               lambda = rep(0.00005,3),
#'               temp_peak = rep(19,3),
#'               time_peak = rep(2060,3),
#'               q = rep(1,3),
#'               time_start = 2005,
#'               time_end = 2100,
#'               leap = 1/12)
#'
#'
#'#######################################################################
#'   #Example 4: Different peaks of temperature.
#'#######################################################################
#'
#'temp_peak1 <- 16
#'temp_peak2 <- 5/4*temp_peak1
#'temp_peak3 <- 5/4*temp_peak2
#'
#'cooling_pulse2(y_ini = c(N = 100, N = 100, N = 100),
#'               temp_ini = rep(28,3),
#'               temp_cmin = rep(18,3),
#'               temp_cmax = rep(30,3),
#'               ro = rep(0.7,3),
#'               lambda = rep(0.00005,3),
#'               temp_peak = c(temp_peak1,temp_peak2,temp_peak3),
#'               time_peak = rep(2060,3),
#'               q = rep(1,3),
#'               time_start = 2005,
#'               time_end = 2100,
#'               leap = 1/12)
#'
#'#######################################################################
#'   #Example 5: Different marginal losses by a non-thermodependent
#'   #           component of intraspecific competition.
#'#######################################################################
#'
#'lambda3 <- 0.01
#'lambda2 <- 1/2*lambda3
#'lambda1 <- 1/2*lambda2
#'
#'cooling_pulse2(y_ini = c(N = 100, N = 100, N = 100),
#'               temp_ini = rep(28,3),
#'               temp_cmin = rep(18,3),
#'               temp_cmax = rep(30,3),
#'               ro = rep(0.7,3),
#'               lambda = c(lambda1,lambda2,lambda3),
#'               temp_peak = rep(25,3),
#'               time_peak = rep(2060,3),
#'               q = rep(1,3),
#'               time_start = 2005,
#'               time_end = 2100,
#'               leap = 1/12)
#'}
###################################################

cooling_pulse2<- function(y_ini = c(N = 400, N = 400, N = 400),
                          temp_ini = rep(35,3),
                          temp_cmin = c(15,18,20),
                          temp_cmax = rep(40,3),
                          ro = rep(0.7,3),
                          lambda = rep(0.00005,3),
                          temp_peak = rep(25,3),
                          time_peak = rep(2060,3),
                          q = rep(5,3),
                          time_start = 2005,
                          time_end = 2100,
                          leap = 1/12){

  times<- seq(time_start, time_end, leap)

if(temp_cmin[1]<temp_cmax[1] && temp_cmin[2]<temp_cmax[2] && temp_cmin[3]<temp_cmax[3] ){

if(temp_cmin[1]<=temp_ini[1] && temp_ini[1]<=temp_cmax[1] && temp_cmin[2]<=temp_ini[2] &&
   temp_ini[2]<=temp_cmax[2] && temp_cmin[3]<=temp_ini[3] && temp_ini[3]<=temp_cmax[3]){

if(temp_peak[1]<=(temp_ini[1]-q[1]) && temp_peak[2]<=(temp_ini[2]-q[2]) && temp_peak[3]<=(temp_ini[3]-q[3])){

if(time_start<=time_peak[1] && time_peak[1]<=time_end && time_start<=time_peak[2] &&
   time_peak[2]<=time_end  && time_start<=time_peak[3] && time_peak[3]<=time_end){


  P2E <- function (times,temp_ini,a,b,q) {
    T <-  temp_ini-(q*(times-time_start)^{2})/((times-time_start)^{2}-b*(times-time_start)+a)

  }



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
  # Temperature parameters
  ##########################################################
  b1<- 2*(time_peak[1]-time_start)*(1+q[1]/(temp_peak[1]-temp_ini[1]))
  a1<- (b1*(time_peak[1]-time_start))/2

  b2<- 2*(time_peak[2]-time_start)*(1+q[2]/(temp_peak[2]-temp_ini[2]))
  a2<- (b2*(time_peak[2]-time_start))/2

  b3<- 2*(time_peak[3]-time_start)*(1+q[3]/(temp_peak[3]-temp_ini[3]))
  a3<- (b3*(time_peak[3]-time_start))/2


  ##########################################################
  # Time
  ##########################################################

time_op11=suppressWarnings((b1*(temp_ini[1]-temp_op1)+sqrt(b1^2*(temp_ini[1]-temp_op1)^2-
4*a1*((temp_ini[1]-temp_op1)-q[1])*(temp_ini[1]-temp_op1)))/(2*((temp_ini[1]-temp_op1)-q[1]))+time_start)
time_op12=suppressWarnings((b1*(temp_ini[1]-temp_op1)-sqrt(b1^2*(temp_ini[1]-temp_op1)^2-
4*a1*((temp_ini[1]-temp_op1)-q[1])*(temp_ini[1]-temp_op1)))/(2*((temp_ini[1]-temp_op1)-q[1]))+time_start)
time_cmin11=suppressWarnings((b1*(temp_ini[1]-temp_cmin[1])+sqrt(b1^2*(temp_ini[1]-temp_cmin[1])^2-
4*a1*((temp_ini[1]-temp_cmin[1])-q[1])*(temp_ini[1]-temp_cmin[1])))/(2*((temp_ini[1]-temp_cmin[1])-q[1]))+time_start)
time_cmin12=suppressWarnings((b1*(temp_ini[1]-temp_cmin[1])-sqrt(b1^2*(temp_ini[1]-temp_cmin[1])^2-
4*a1*((temp_ini[1]-temp_cmin[1])-q[1])*(temp_ini[1]-temp_cmin[1])))/(2*((temp_ini[1]-temp_cmin[1])-q[1]))+time_start)
suppressWarnings(sqrt(b1^2*(temp_ini[1]-temp_cmin[1])^2-4*a1*((temp_ini[1]-temp_cmin[1])-q[1])*(temp_ini[1]-temp_cmin[1])))
  ##########################################################

time_op21=suppressWarnings((b2*(temp_ini[2]-temp_op2)+sqrt(b2^2*(temp_ini[2]-temp_op2)^2-
4*a2*((temp_ini[2]-temp_op2)-q[2])*(temp_ini[2]-temp_op2)))/(2*((temp_ini[2]-temp_op2)-q[2]))+time_start)
time_op22=suppressWarnings((b2*(temp_ini[2]-temp_op2)-sqrt(b2^2*(temp_ini[2]-temp_op2)^2-
4*a2*((temp_ini[2]-temp_op2)-q[2])*(temp_ini[2]-temp_op2)))/(2*((temp_ini[2]-temp_op2)-q[2]))+time_start)
time_cmin21=suppressWarnings((b2*(temp_ini[2]-temp_cmin[2])+sqrt(b2^2*(temp_ini[2]-temp_cmin[2])^2-
4*a2*((temp_ini[2]-temp_cmin[2])-q[2])*(temp_ini[2]-temp_cmin[2])))/(2*((temp_ini[2]-temp_cmin[2])-q[2]))+time_start)
time_cmin22=suppressWarnings((b2*(temp_ini[2]-temp_cmin[2])-sqrt(b2^2*(temp_ini[2]-temp_cmin[2])^2-
4*a2*((temp_ini[2]-temp_cmin[2])-q[2])*(temp_ini[2]-temp_cmin[2])))/(2*((temp_ini[2]-temp_cmin[2])-q[2]))+time_start)
suppressWarnings(sqrt(b2^2*(temp_ini[2]-temp_cmin[2])^2-4*a2*((temp_ini[2]-temp_cmin[2])-q[2])*(temp_ini[2]-temp_cmin[2])))
  ##########################################################
time_op31=suppressWarnings((b3*(temp_ini[3]-temp_op3)+sqrt(b3^2*(temp_ini[3]-temp_op3)^2-
4*a3*((temp_ini[3]-temp_op3)-q[3])*(temp_ini[3]-temp_op3)))/(2*((temp_ini[3]-temp_op3)-q[3]))+time_start)
time_op32=suppressWarnings((b3*(temp_ini[3]-temp_op3)-sqrt(b3^2*(temp_ini[3]-temp_op3)^2-
4*a3*((temp_ini[3]-temp_op3)-q[3])*(temp_ini[3]-temp_op3)))/(2*((temp_ini[3]-temp_op3)-q[3]))+time_start)
time_cmin31=suppressWarnings((b3*(temp_ini[3]-temp_cmin[3])+sqrt(b3^2*(temp_ini[3]-temp_cmin[3])^2-
4*a3*((temp_ini[3]-temp_cmin[3])-q[3])*(temp_ini[3]-temp_cmin[3])))/(2*((temp_ini[3]-temp_cmin[3])-q[3]))+time_start)
time_cmin32=suppressWarnings((b3*(temp_ini[3]-temp_cmin[3])-sqrt(b3^2*(temp_ini[3]-temp_cmin[3])^2-
4*a3*((temp_ini[3]-temp_cmin[3])-q[3])*(temp_ini[3]-temp_cmin[3])))/(2*((temp_ini[3]-temp_cmin[3])-q[3]))+time_start)
  suppressWarnings(sqrt(b3^2*(temp_ini[3]-temp_cmin[3])^2-4*a3*((temp_ini[3]-temp_cmin[3])-q[3])*(temp_ini[3]-temp_cmin[3])))
  #########################################################

  ##########################################################
  # Time limits
  ##########################################################
  tm<-c(time_cmin11[1],time_cmin12[1],time_cmin21[1],time_cmin22[1],time_cmin31[1],time_cmin32[1])
  tm_new <- tm
  tm_new[is.nan(tm_new)] <- times[length(times)]

  if(times[length(times)]<tm_new[1]){
    tm_new[1]=times[length(times)]
  }
  if(times[length(times)]<tm_new[2]){
    tm_new[2]=times[length(times)]
  }
  if(times[length(times)]<tm_new[3]){
    tm_new[3]=times[length(times)]
  }

  if(times[length(times)]<tm_new[4]){
    tm_new[4]=times[length(times)]
  }

  if(times[length(times)]<tm_new[5]){
    tm_new[5]=times[length(times)]
  }

  if(times[length(times)]<tm_new[6]){
    tm_new[6]=times[length(times)]
  }

  if(tm_new[1]<=tm_new[2]){

    tm_new1<-tm_new[1]
  }else{
    tm_new1<-tm_new[2]
  }


  if(tm_new[3]<=tm_new[4]){

    tm_new2<-tm_new[3]
  }else{
    tm_new2<-tm_new[4]
  }

  if(tm_new[5]<=tm_new[6]){

    tm_new3<-tm_new[5]
  }else{
    tm_new3<-tm_new[6]
  }






  ##########################################################
  # Parameters
  ##########################################################
parms1<-c(temp_cmin[1],temp_ini[1],temp_cmax[1],temp_op1, ro[1], q[1], a1, b1, lambda[1])
parms2<-c(temp_cmin[2],temp_ini[2],temp_cmax[2],temp_op2, ro[2], q[2], a2, b2, lambda[2])
parms3<-c(temp_cmin[3],temp_ini[3],temp_cmax[3],temp_op3, ro[3], q[3], a3, b3, lambda[3])

  ##########################################################
  # Model for each trend
  ##########################################################
model1 <- function (times, y,parms1) {
      with(as.list(c(y)), {
  T <- P2E(times,temp_ini[1],a1,b1,q[1])
  r1<- rate_TPC(T,ro[1],temp_cmin[1],temp_cmax[1],temp_op1)
  dN <-   r1 * N * (1 - lambda[1]*(N / r1))
  list(dN,T,r1) })
    }
###############################################################
model2 <- function (times, y,parms2) {
      with(as.list(c(y)), {
 T <- P2E(times,temp_ini[2],a2,b2,q[2])
 r2<- rate_TPC(T,ro[2],temp_cmin[2],temp_cmax[2],temp_op2)
 dN <-   r2 * N * (1 - lambda[2]*(N / r2))
 list(dN,T,r2) })
    }
###############################################################
model3 <- function (times, y,parms3) {
      with(as.list(c(y)), {
  T <- P2E(times,temp_ini[3],a3,b3,q[3])
  r3<- rate_TPC(T,ro[3],temp_cmin[3],temp_cmax[3],temp_op3)
  dN <-   r3 * N * (1 - lambda[3]*(N / r3))
  list(dN,T,r3)})
    }
###############################################################
  # Solution
##############################################################
out1 <- ode(y=y_ini[1], times, model1, parms1, method = "ode45")
out2 <- ode(y=y_ini[2], times, model2, parms2, method = "ode45")
out3 <- ode(y=y_ini[3], times, model3, parms3, method = "ode45")
#############################################################
###############################################################
    # Temperature trend
##############################################################

    da1<-data.frame('x'=times,'y'=out1[,3] )
    da2<-data.frame('x'=times,'y'=out2[,3] )
    da3<-data.frame('x'=times,'y'=out3[,3] )

###############################################################
    # Abundance
##############################################################

    data1<-data.frame('x'=times,'y'=out1[,2] )
    data2<-data.frame('x'=times,'y'=out2[,2] )
    data3<-data.frame('x'=times,'y'=out3[,2] )

###############################################################
    # Carrying capacity
  ##############################################################
    K1=out1[,4]/lambda[1]
    K2=out2[,4]/lambda[2]
    K3=out3[,4]/lambda[3]

    dat1<-data.frame('x'=times,'y'=K1 )
    dat2<-data.frame('x'=times,'y'=K2 )
    dat3<-data.frame('x'=times,'y'=K3 )

###############################################################
# Data
###############################################################
Data<- data.frame(times,out1[,3],out1[,2],K1,out2[,3],out2[,2],K2,
                  out3[,3],out3[,2],K3)
names(Data)<- c("Time","Temperature Scenario 1","Abundance scenario 1",
                "Carrying capacity scenario 1","Temperature scenario 2",
                "Abundance scenario 2","Carrying capacity scenario 2",
                "Temperature scenario 3","Abundance scenario 3","Carrying
                capacity scenario 3")
    u<- formattable(Data, align = c("l", rep("r", NCOL(Data))))
    print(u)

###############################################################
    # Plots
##############################################################
  data<-rbind(data1,data2,data3,dat1,dat2,dat3,da1,da2,da3)

    p1<- ggplot(data, aes(x=.data$x, y=.data$y)) +
            theme_bw()+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
            geom_ribbon(data=subset(dat1,times>times[1] & times<tm_new1),aes(x=.data$x,ymax=.data$y),ymin=0,alpha=0.3, fill="brown") +
            geom_ribbon(data=subset(dat2,times>times[1] & times<tm_new2),aes(x=.data$x,ymax=.data$y),ymin=0,alpha=0.3, fill="green4") +
            geom_ribbon(data=subset(dat3,times>times[1] & times<tm_new3),aes(x=.data$x,ymax=.data$y),ymin=0,alpha=0.3,fill="blue") +
            geom_vline(xintercept = tm_new1, size=.5, color="brown",linetype="dashed")+
            geom_vline(xintercept = tm_new2, size=.5, color="green4",linetype="dashed")+
            geom_vline(xintercept = tm_new3, size=.5, color="blue",linetype="dashed")+
            geom_line(data =subset(data1,times>times[1] & times<tm_new1), color = "brown")+
            geom_line(data =subset(data2,times>times[1] & times<tm_new2), color = "green4")+
            geom_line(data =subset(data3,times>times[1] & times<tm_new3), color = "blue")+
            labs(x = "Time",y="Abundance")+
            theme(plot.title = element_text(size=40))+
            theme(plot.title = element_text(hjust = 0.5))+
            theme(axis.title.y = element_text(size = rel(1), angle = 90))+
            theme(axis.title.x = element_text(size = rel(1), angle = 00))+
            labs(tag = "(a)")


    p2 <- ggplot(data, aes(x=.data$x, y=.data$y)) +
            theme_bw()+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
            geom_vline(xintercept = tm_new1, size=.5, color="brown",linetype="dashed")+
            geom_vline(xintercept = tm_new2, size=.5, color="green4",linetype="dashed")+
            geom_vline(xintercept = tm_new3, size=.5, color="blue",linetype="dashed")+
            geom_line(data =subset(da1,times>times[1] & times<tm_new1), color = "brown")+
            geom_line(data =subset(da2,times>times[1] & times<tm_new2), color = "green4")+
            geom_line(data =subset(da3,times>times[1] & times<tm_new3), color = "blue")+
            labs(x = "Time",y="Temperature")+
            theme(axis.title.y = element_text(size = rel(1), angle = 90))+
            theme(axis.title.x = element_text(size = rel(1), angle = 00))+
            labs(tag = "(b)")

    plot_grid(p1, p2)

        }else{
          stop("It is recommended that the time in which the peak temperature is reached is within the time sequence")
        }

      }else{
        stop("The peak temperature must be less than or equal to the temp_ini-q")
      }


    }else{
      stop("The initial study temperature must be within the thermal tolerance range")
    }

  }else{

    stop("The minimum critical temperature must be less than the maximum critical temperature")
  }

 }


