#' IPCC RCP8.5 scenario
#'
#' @description This function allows simulating the effect of an increase in environmental
#'              temperature according to the IPCC RCP8.5 scenario (2014) on the abundance of
#'              ectotherm populations.
#'
#'
#'@param y_ini Initial population values (must be written with its name: N).
#'@param temp_ini Initial temperature.
#'@param temp_cmin Minimum critical temperature.
#'@param temp_cmax Maximum critical temperature.
#'@param ro Population growth rate at optimum temperature.
#'@param lambda Marginal loss by non-thermodependent intraspecific competition.
#'@param time_start Start of time sequence.
#'@param time_end End of time sequence.
#'@param leap Time sequence step.
#'
#'@details Three populations can be simulated simultaneously. The temperature trend is determined
#'         by a projection of the change in global mean surface temperature according to the IPCC
#'         RCP8.5 scenario. In each input vector, the parameters for the three simulations must be
#'         specified (finite numbers for initial population abundance). The simulations are obtained
#'         by a model that incorporates the effects of temperature over time, which leads to a
#'         non-autonomous ODE approach. This is function uses the ODE solver implemented in the package
#'         deSolve (Soetaert et al., 2010).
#'
#'@references IPCC. (2014): Climate Change 2014: Synthesis Report. Contribution of Working Groups I,
#'            II and III to the Fifth Assessment Report of the Intergovernmental Panel on Climate
#'            Change [Core Writing Team, R.K. Pachauri and L.A. Meyer (eds.)]. IPCC, Geneva,
#'            Switzerland, 151 pp.
#'
#'@references Soetaert, K., Petzoldt, T., & Setzer, R. (2010). Solving Differential Equations in R: Package
#'            deSolve. Journal of Statistical Software, 33(9), 1 - 25.
#'            doi:http://dx.doi.org/10.18637/jss.v033.i09
#'
#'@return (1) A data.frame with columns having the simulated trends.
#'@return (2) A two-panel figure in which (a) shows the population abundance curves represented by solid
#'            lines and the corresponding carrying capacities are represented by shaded areas. In (b) the
#'            temperature trend is shown. The three simultaneous simulations are depicted by different colors,
#'            i.e. 1st brown, 2nd green and 3rd blue.
#'
#'@export
#'@examples
#'
#'#######################################################################
#'   #Example 1: Different initial population abundances.
#'#######################################################################
#'
#'IPCC_RCP8_5(y_ini = c(N = 100, N = 200, N = 400),
#'            temp_ini = rep(26,3),
#'            temp_cmin = rep(18,3),
#'            temp_cmax = rep(30,3),
#'            ro = rep(0.7,3),
#'            lambda = rep(0.00005,3),
#'            time_start = 2005,
#'            time_end = 2100,
#'            leap = 1/12)
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
#'IPCC_RCP8_5(y_ini = c(N = 100, N = 100, N = 100),
#'            temp_ini = rep(30,3),
#'            temp_cmin = c(temp_cmin1,temp_cmin2,temp_cmin3),
#'            temp_cmax = c(temp_cmax1,temp_cmax2,temp_cmax3),
#'            ro = rep(0.7,3),
#'            lambda = rep(0.00005,3),
#'            time_start = 2005,
#'            time_end = 2100,
#'            leap = 1/12)
#'\donttest{
#'#######################################################################
#'   #Example 3: Different relationships between initial environmental
#'   #           temperature and optimum temperature.
#'#######################################################################
#'
#'temp_cmin <- 18
#'temp_cmax <- 30
#'
#'# Temperature at which performance is at its maximum value.
#'temp_op <- (temp_cmax+temp_cmin)/3+sqrt(((temp_cmax+temp_cmin)/3)^2-
#'            (temp_cmax*temp_cmin)/3)
#'
#'temp_ini1 <- (temp_cmin+temp_op)/2
#'temp_ini2 <- temp_op
#'temp_ini3 <- (temp_op+temp_cmax)/2
#'
#'IPCC_RCP8_5(y_ini = c(N = 100, N = 100, N = 100),
#'            temp_ini = c(temp_ini1,temp_ini2,temp_ini3),
#'            temp_cmin = rep(temp_cmin,3),
#'            temp_cmax = rep(temp_cmax,3),
#'            ro = rep(0.7,3),
#'            lambda = rep(0.00005,3),
#'            time_start = 2005,
#'            time_end = 2100,
#'            leap = 1/12)
#'
#'#######################################################################
#'   #Example 4:  Different marginal losses by a non-thermodependent
#'   #            component of intraspecific competition.
#'#######################################################################
#'
#'lambda3 <- 0.01
#'lambda2 <- 1/2*lambda3
#'lambda1 <- 1/2*lambda2
#'
#'IPCC_RCP8_5(y_ini = c(N = 100, N = 100,N = 100),
#'            temp_cmin = rep(18,3),
#'            temp_ini = rep(25,3),
#'            temp_cmax = rep(30,3),
#'            ro = rep(0.7,3),
#'            lambda = c(lambda1,lambda2,lambda3),
#'            time_start = 2005,
#'            time_end = 2100,
#'            leap = 1/12)
#'}
###################################################


IPCC_RCP8_5<- function(y_ini = c(N = 400, N = 400, N = 400),
                       temp_ini = rep(20,3),
                       temp_cmin = rep(18,3),
                       temp_cmax = c(25,28,30),
                       ro = rep(0.7,3),
                       lambda = rep(0.00005,3),
                       time_start = 2005,
                       time_end = 2100,
                       leap = 1/12){


  times<- seq(time_start, time_end, leap)

if(time_end<=2100){
    if(time_start<=time_end){

if(temp_cmin[1]<temp_cmax[1] && temp_cmin[2]<temp_cmax[2] && temp_cmin[3]<temp_cmax[3] ){


if(temp_cmin[1]<=temp_ini[1] && temp_ini[1]<=temp_cmax[1] && temp_cmin[2]<=temp_ini[2] &&
   temp_ini[2]<=temp_cmax[2] && temp_cmin[3]<=temp_ini[3] && temp_ini[3]<=temp_cmax[3]){



  RCP8.5 <- function(date,a,b) {a * exp(b * date)}
  values <- c(0.61, 2, 3.7)
  x<- c(2005,2065,2100)
  y<- values
  df <- data.frame(x, y)

  m<- nls(y ~ exp(loga + b * x), df, start = list( loga = log(2), b = 0.005),control = list (maxiter = 500))
  y_est<-predict(m,df$x)

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
  # Time
##########################################################

time_cmax1= suppressWarnings(1/coef(m)[2]*log((temp_cmax[1]-temp_ini[1])/exp(coef(m)[1])))
time_cmax2= suppressWarnings(1/coef(m)[2]*log((temp_cmax[2]-temp_ini[2])/exp(coef(m)[1])))
time_cmax3= suppressWarnings(1/coef(m)[2]*log((temp_cmax[3]-temp_ini[3])/exp(coef(m)[1])))

##########################################################

    if(time_cmax1<=times[length(times)]){
      time_ext1<- time_cmax1
    }else{
      time_ext1<- times[length(times)]
    }

    if(time_cmax2<=times[length(times)]){
      time_ext2<- time_cmax2
    }else{
      time_ext2<- times[length(times)]
    }


    if(time_cmax3<=times[length(times)]){
      time_ext3<- time_cmax3
    }else{
      time_ext3<- times[length(times)]
    }

##########################################################
    # Parameters
##########################################################

parms1<-c(temp_cmin[1],temp_ini[1],temp_cmax[1],temp_op1,ro[1], lambda[1])
parms2<-c(temp_cmin[2],temp_ini[2],temp_cmax[2],temp_op2,ro[2], lambda[2])
parms3<-c(temp_cmin[3],temp_ini[3],temp_cmax[3],temp_op3,ro[3], lambda[3])

##########################################################
##########################################################
      # Model for each trend
##########################################################
model1 <- function (times, y,parms1) {
with(as.list(c(y)), {
T <- RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[1]
r1<- rate_TPC(T,ro[1],temp_cmin[1],temp_cmax[1],temp_op1)
dN <-   r1 * N * (1 - lambda[1]*(N / r1))
list(dN,T,r1) })
}

###############################################################
###############################################################
model2 <- function (times, y,parms2) {
with(as.list(c(y)), {
T <- RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[2]
r2<- rate_TPC(T,ro[2],temp_cmin[2],temp_cmax[2],temp_op2)
dN <-   r2 * N * (1 - lambda[2]*(N / r2))
list(dN,T,r2)})
      }
###############################################################
###############################################################
model3 <- function (times, y,parms3) {
with(as.list(c(y)), {
T <- RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[3]
r3<- rate_TPC(T,ro[3],temp_cmin[3],temp_cmax[3],temp_ini[3])
dN <-   r3 * N * (1 - lambda[3]*(N / r3))
list(dN,T,r3)})
      }
###############################################################

###############################################################
      # Solution
###############################################################

out1 <- ode(y=y_ini[1], times, model1, parms1,method = "ode23")
out2 <- ode(y=y_ini[2], times, model2, parms2,method = "ode23")
out3 <- ode(y=y_ini[3], times, model3, parms3,method = "ode23")

##############################################################
##############################################################
      # Temperature trend
##############################################################

      da1<-data.frame('x'=times,'y'=out1[,3] )
      da2<-data.frame('x'=times,'y'=out2[,3] )
      da3<-data.frame('x'=times,'y'=out3[,3] )

##############################################################
      # Abundance
##############################################################

      data1<-data.frame('x'=times,'y'=out1[,2] )
      data2<-data.frame('x'=times,'y'=out2[,2] )
      data3<-data.frame('x'=times,'y'=out3[,2] )

##############################################################
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

Data<- data.frame(times,out1[,3],out1[,2],K1,out2[,3],out2[,2],
                  K2,out3[,3],out3[,2],K3)
names(Data)<- c("Time","Temperature Scenario 1","Abundance scenario 1",
                "Carrying capacity scenario 1","Temperature scenario 2",
                "Abundance scenario 2","Carrying capacity scenario 2",
                "Temperature scenario 3","Abundance scenario 3","Carrying
                capacity scenario 3")
u<- formattable(Data, align = c("l", rep("r", NCOL(Data))))
print(u)

###############################################################

###############################################################
      # Plots
##############################################################

      data<-rbind(data1, data2, data3)
      dat<-rbind(dat1,dat2,dat3)
      da<-rbind(da1, da2, da3)

p1 <- ggplot(data, aes(x=.data$x, y=.data$y)) +
theme_bw()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_ribbon(data=subset(dat1,times>times[1] & times<time_ext1),aes(x=.data$x,
                              ymax=.data$y),ymin=0,alpha=0.3, fill="brown") +
geom_ribbon(data=subset(dat2,times>times[1] & times<time_ext2),aes(x=.data$x,
                             ymax=.data$y),ymin=0,alpha=0.3, fill="green4") +
geom_ribbon(data=subset(dat3,times>times[1] & times<time_ext3),aes(x=.data$x,
                               ymax=.data$y),ymin=0,alpha=0.3, fill="blue") +
geom_vline(xintercept = time_ext1, size=.5, color="brown",linetype="dashed")+
geom_vline(xintercept = time_ext2, size=.5, color="green4",linetype="dashed")+
geom_vline(xintercept = time_ext3, size=.5, color="blue",linetype="dashed")+
geom_line(data =subset(data1,times>times[1] & times<time_ext1), color = "brown")+
geom_line(data =subset(data2,times>times[1] & times<time_ext2), color = "green4")+
geom_line(data =subset(data3,times>times[1] & times<time_ext3), color = "blue")+
labs(x = "Time",y="Abundance")+
theme(plot.title = element_text(size=40))+
theme(plot.title = element_text(hjust = 0.5))+
theme(axis.title.y = element_text(size = rel(1), angle = 90))+
theme(axis.title.x = element_text(size = rel(1), angle = 00))+
labs(tag = "(a)")


p2 <- ggplot(da, aes(x=.data$x, y=.data$y)) +
theme_bw()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_vline(xintercept = time_ext1, size=.5, color="brown",linetype="dashed")+
geom_vline(xintercept = time_ext2, size=.5, color="green4",linetype="dashed")+
geom_vline(xintercept = time_ext3, size=.5, color="blue",linetype="dashed")+
geom_line(data =subset(da1,times>times[1] & times<time_cmax1), color = "brown")+
geom_line(data =subset(da2,times>times[1] & times<time_cmax2), color = "green4")+
geom_line(data =subset(da3,times>times[1] & times<time_cmax3), color = "blue")+
labs(x = "Time",y="Temperature")+
theme(axis.title.y = element_text(size = rel(1), angle = 90))+
theme(axis.title.x = element_text(size = rel(1), angle = 00))+
labs(tag = "(b)")


plot_grid(p1, p2)


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



