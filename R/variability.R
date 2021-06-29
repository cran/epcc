#' Increasing variable temperature trend
#'
#' @description This function allows simulating the effect of increasing variable temperature
#'              trend on the abundance of ectotherm populations.
#'
#'
#'
#'@param y_ini Initial population values (must be written with its name: N).
#'@param temp_ini Initial temperature.
#'@param temp_cmin Minimum critical temperature.
#'@param temp_cmax Maximum critical temperature.
#'@param ro Population growth rate at optimum temperature.
#'@param lambda Marginal loss by non-thermodependent intraspecific competition.
#'@param mean Average noise.
#'@param sd Standard deviation of noise.
#'@param time_start Start of time sequence.
#'@param time_end End of time sequence.
#'@param leap Time sequence step.
#'
#'@details Three populations and/or scenarios can be simulated simultaneously.
#'         The temperature trend considers the IPCC RCP8.5 but adding random noise
#'         following a normal distribution, which is characterized by the mean and
#'         standard deviation. In each input vector, the parameters for the three
#'         simulations must be specified (finite numbers for the initial population
#'         abundance). The simulations are obtained by a model that incorporates the
#'         effects of temperature over time, which leads to a non-autonomous ODE approach.
#'         This is function uses the ODE solver implemented in the package deSolve
#'         (Soetaert et al., 2010).
#'
#'
#'
#'@return (1) A data.frame with columns having the simulated trends.
#'@return (2) A two-panel figure in which (a) shows the population abundance curves represented
#'            by solid lines and the corresponding carrying capacities are represented by shaded areas.
#'            In (b) the temperature trend is shown. The three simultaneous simulations are depicted by
#'            different colors, i.e. 1st brown, 2nd green and 3rd blue.
#'
#'
#'@references Soetaert, K., Petzoldt, T., & Setzer, R. (2010). Solving Differential Equations in R: Package
#'            deSolve. Journal of Statistical Software, 33(9), 1 - 25.
#'            doi:http://dx.doi.org/10.18637/jss.v033.i09
#'@export
#'@examples
#'
#'#######################################################################
#'   #Example 1: Different initial population abundances.
#'#######################################################################
#'
#'variability(y_ini = c(N = 100, N = 200, N = 400),
#'            temp_ini = rep(22,3),
#'            temp_cmin = rep(18,3),
#'            temp_cmax = rep(35,3),
#'            ro = rep(0.7,3),
#'            lambda = rep(0.00005,3),
#'            mean = rep(2,3),
#'            sd = rep(3,3),
#'            time_start = 2005,
#'            time_end = 2100,
#'            leap = 1/12)
#'\donttest{
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
#'variability(y_ini = c(N = 100, N = 100, N = 100),
#'            temp_ini = rep(26,3),
#'            temp_cmin = c(temp_cmin1,temp_cmin2,temp_cmin3),
#'            temp_cmax = c(temp_cmax1,temp_cmax2,temp_cmax3),
#'            ro = rep(0.7,3),
#'            lambda = rep(0.00005,3),
#'            mean = rep(2,3),
#'            sd = rep(2,3),
#'            time_start = 2005,
#'            time_end = 2100,
#'            leap = 1/12)
#'
#'#######################################################################
#'   #Example 3: Different relationships between initial environmental
#'   #           temperature and optimum temperature.
#'#######################################################################
#'
#'temp_cmin <- 18
#'temp_cmax <- 40
#'
#'# Temperature at which performance is at its maximum value.
#'temp_op <- (temp_cmax+temp_cmin)/3+sqrt(((temp_cmax+temp_cmin)/3)^2-
#'            (temp_cmax*temp_cmin)/3)
#'
#'temp_ini1 <- (temp_cmin+temp_op)/2
#'temp_ini2 <- temp_op
#'temp_ini3 <- (temp_op+temp_cmax)/2
#'
#'variability(y_ini = c(N = 100, N = 100, N = 100),
#'            temp_ini = c(temp_ini1,temp_ini2,temp_ini3),
#'            temp_cmin = rep(temp_cmin,3),
#'            temp_cmax = rep(temp_cmax,3),
#'            ro = rep(0.7,3),
#'            lambda = rep(0.00005,3),
#'            mean = rep(1,3),
#'            sd = rep(0.5,3),
#'            time_start = 2005,
#'            time_end = 2100,
#'            leap = 1/12)
#'
#'#######################################################################
#'   #Example 4: Different marginal losses by a non-thermodependent
#'   #           component of intraspecific competition.
#'#######################################################################
#'
#'lambda3 <- 0.01
#'lambda2 <- 1/2*lambda3
#'lambda1 <- 1/2*lambda2
#'
#'variability(y_ini = c(N = 100, N = 100, N = 100),
#'            temp_ini = rep(26,3),
#'            temp_cmin = rep(18,3),
#'            temp_cmax = rep(35,3),
#'            ro = rep(0.7,3),
#'            lambda = c(lambda1,lambda2,lambda3),
#'            mean = rep(1,3),
#'            sd = rep(2,3),
#'            time_start = 2005,
#'            time_end = 2100,
#'            leap = 1/12)
#'
#'#######################################################################
#'   #Example 5: Different noise means.
#'#######################################################################
#'
#'mean3 <- 2
#'mean2 <- 1/2*mean3
#'mean1 <- 1/2*mean2
#'
#'variability(y_ini = c(N = 100, N = 100, N = 100),
#'            temp_ini = rep(22,3),
#'            temp_cmin = rep(18,3),
#'            temp_cmax = rep(35,3),
#'            ro = rep(0.7,3),
#'            lambda = rep(0.00005,3),
#'            mean = c(mean1,mean2,mean3),
#'            sd = rep(2,3),
#'            time_start = 2005,
#'            time_end = 2100,
#'            leap = 1/12)
#'
#'#######################################################################
#'   #Example 6: Different noise standard deviations.
#'#######################################################################
#'
#'sd3 <- 4
#'sd2 <- 1/2*sd3
#'sd1 <- 1/2*sd2
#'
#'variability(y_ini = c(N = 100, N = 100, N = 100),
#'            temp_ini = rep(22,3),
#'            temp_cmin = rep(18,3),
#'            temp_cmax = rep(35,3),
#'            ro = rep(0.7,3),
#'            lambda = rep(0.00005,3),
#'            mean = rep(3,3),
#'            sd = c(sd1,sd2,sd3),
#'            time_start = 2005,
#'            time_end = 2100,
#'            leap = 1/12)
#'}
###################################################


variability<- function(y_ini = c(N =400, N = 400, N = 400),
                       temp_ini = rep(20,3),
                       temp_cmin = rep(18,3),
                       temp_cmax = c(25,28,32),
                       ro = rep(0.7,3),
                       lambda = rep(0.00005,3),
                       mean = rep(5,3),
                       sd = rep(2,3),
                       time_start = 2005,
                       time_end = 2100,
                       leap = 1/12){

  times<- seq(time_start, time_end, leap)


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

temp_op1<- ( temp_cmax[1]+temp_cmin[1])/3+sqrt((( temp_cmax[1]+temp_cmin[1])/3)^2-
                                                   ( temp_cmax[1]*temp_cmin[1])/3)

temp_op2<- ( temp_cmax[2]+temp_cmin[2])/3+sqrt((( temp_cmax[2]+temp_cmin[2])/3)^2-
                                                   ( temp_cmax[2]*temp_cmin[2])/3)

temp_op3<- ( temp_cmax[3]+temp_cmin[3])/3+sqrt((( temp_cmax[3]+temp_cmin[3])/3)^2-
                                                   ( temp_cmax[3]*temp_cmin[3])/3)


##########################################################
    # Parameters
##########################################################
parms1<-c(temp_cmin[1],temp_ini[1], temp_cmax[1],temp_op1,ro[1], lambda[1])
parms2<-c(temp_cmin[2],temp_ini[2], temp_cmax[2],temp_op2,ro[2], lambda[2])
parms3<-c(temp_cmin[3],temp_ini[3], temp_cmax[3],temp_op3,ro[3], lambda[3])
##########################################################


##########################################################
    # Model for each trend
##########################################################
model1 <- function (times, y,parms1 ) {
    with(as.list(c(y)), {
    set.seed(times/length(times))
    cd1<- rnorm(length(times),mean= mean[1],sd=sd[1])
    T <- RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[1]+cd1
    r1<- rate_TPC(T,ro[1],temp_cmin[1], temp_cmax[1],temp_op1)
    dN <-   r1 * N * (1 - lambda[1]*(N / r1))
    list(dN,T,r1,cd1) })
    }
#########################################################
model2 <- function (times, y,parms2) {
    with(as.list(c(y)), {
    set.seed(times/length(times))
    cd2<- rnorm(length(times),mean= mean[2],sd=sd[2])
    T <- RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[2]+cd2
    r2<- rate_TPC(T,ro[2],temp_cmin[2], temp_cmax[2],temp_op2)
    dN <-   r2 * N * (1 - lambda[2]*(N / r2))
    list(dN,T,r2,cd2)})
    }
#########################################################
model3 <- function (times, y,parms3) {
    with(as.list(c(y)), {
    set.seed(times/length(times))
    cd3<-rnorm(length(times), mean=mean[3],sd=sd[3])
    T <- RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[3]+cd3
    r3<- rate_TPC(T,ro[3],temp_cmin[3], temp_cmax[3],temp_op3)
    dN <-   r3 * N * (1 - lambda[3]*(N / r3))
    list(dN,T,r3,cd3)})
    }
##########################################################

##########################################################
    # Solution
##########################################################

out1 <- ode(y=y_ini[1], times, model1, parms1 ,method = "ode45")
out2 <- ode(y=y_ini[2], times, model2, parms2, method = "ode45")
out3 <- ode(y=y_ini[3], times, model3, parms3, method = "ode45")

#############################################################

##############################################################
    # Time limits
##############################################################

times_new1<-vector(mode = "numeric", length = 0)
times_new2<-vector(mode = "numeric", length = 0)
times_new3<-vector(mode = "numeric", length = 0)

for (i in 2: length(times)){

if((out1[i-1,3]-temp_cmin[1])>=0 && (out1[i,3]-temp_cmin[1])<0){
        times_new1[i-1]<- times[i-1]

}else if((out1[i-1,3]-temp_cmin[1])<=0 && (out1[i,3]-temp_cmin[1])>0){

 times_new1[i-1]<- times[i-1]
 }else{
  times_new1[i-1]<- 0
      }
    }


for (i in 2: length(times)){

if((out2[i-1,3]-temp_cmin[2])>=0 && (out2[i,3]-temp_cmin[2])<0){
        times_new2[i-1]<- times[i-1]

}else if((out2[i-1,3]-temp_cmin[2])<=0 && (out2[i,3]-temp_cmin[2])>0){

  times_new2[i-1]<- times[i-1]
   }else{
     times_new2[i-1]<- 0
      }
    }

for (i in 2: length(times)){

if((out3[i-1,3]-temp_cmin[3])>=0 && (out3[i,3]-temp_cmin[3])<0){
        times_new3[i-1]<- times[i-1]

}else if((out3[i-1,3]-temp_cmin[3])<=0 && (out3[i,3]-temp_cmin[3])>0){

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

    if(!is.na(as.integer(index1))==FALSE){
      times_inf1<- times[length(times)]
    }else{
      times_inf1<- times[index1]
    }
    if(!is.na(as.integer(index2))== FALSE){
      times_inf2<- times[length(times)]
    }else{
      times_inf2<- times[index2]
    }

    if(!is.na(as.integer(index3))== FALSE){
      times_inf3<- times[length(times)]
    }else{
      times_inf3<- times[index3]
    }


#############################################################
#############################################################

    times_new4<-vector(mode = "numeric", length = 0)
    times_new5<-vector(mode = "numeric", length = 0)
    times_new6<-vector(mode = "numeric", length = 0)

for (i in 2: length(times)){

if(( temp_cmax[1]-out1[i-1,3])>=0 && ( temp_cmax[1]-out1[i,3])<0){
    times_new4[i-1]<- times[i-1]

}else if(( temp_cmax[1]-out1[i-1,3])<=0 && ( temp_cmax[1]-out1[i,3])>0){

      times_new4[i-1]<- times[i-1]
      }else{
        times_new4[i-1]<- 0
      }
    }


for (i in 2: length(times)){

if(( temp_cmax[2]-out2[i-1,3])>=0 && ( temp_cmax[2]-out2[i,3])<0){
        times_new5[i-1]<- times[i-1]

}else if(( temp_cmax[2]-out2[i-1,3])<=0 && ( temp_cmax[2]-out2[i,3])>0){

        times_new5[i-1]<- times[i-1]
      }else{
        times_new5[i-1]<- 0
      }
    }


for (i in 2: length(times)){

if(( temp_cmax[3]-out3[i-1,3])>=0 && ( temp_cmax[3]-out3[i,3])<0){
      times_new6[i-1]<- times[i-1]

}else if(( temp_cmax[3]-out3[i-1,3])<=0 && ( temp_cmax[3]-out3[i,3])>0){

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
      times_sup1<- times[length(times)]
    }else{
      times_sup1<- times[index4]
    }
    if(!is.na(as.integer(index5))== FALSE){
      times_sup2<- times[length(times)]
    }else{
      times_sup2<- times[index5]
    }

    if(!is.na(as.integer(index6))== FALSE){
      times_sup3<- times[length(times)]
    }else{
      times_sup3<- times[index6]
    }

#############################################################
#############################################################

    if(times_inf1<= times_sup1){
      times_ext1<-times_inf1
    }else{
      times_ext1<-times_sup1
    }

    if(times_inf2<=times_sup2){
      times_ext2<-times_inf2
    }else{
      times_ext2<-times_sup2
    }

    if(times_inf3<=times_sup3){
      times_ext3<-times_inf3
    }else{
      times_ext3<-times_sup3
    }


###############################################################
    # Temperature trend
##############################################################

    da1<-data.frame('x'=times,'y'=out1[,3] )
    da2<-data.frame('x'=times,'y'=out2[,3] )
    da3<-data.frame('x'=times,'y'=out3[,3] )

###############################################################
    # Abundance
###############################################################

    data1<-data.frame('x'=times,'y'=out1[,2] )
    data2<-data.frame('x'=times,'y'=out2[,2] )
    data3<-data.frame('x'=times,'y'=out3[,2] )

###############################################################
    # Carrying capacity
###############################################################

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
    # Plots
###############################################################

  data<-rbind(data1,data2,data3)
  dat<-rbind(dat1,dat2,dat3)
  da<-rbind(da1,da2,da3)

p1 <- ggplot(data, aes(x=.data$x, y=.data$y)) +
theme_bw()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_ribbon(data=subset(dat1,times>times[1] & times<times_ext1),aes(x=.data$x,
                                ymax=.data$y),ymin=0,alpha=0.3, fill="brown") +
geom_ribbon(data=subset(dat2,times>times[1] & times<times_ext2),aes(x=.data$x,
                                ymax=.data$y),ymin=0,alpha=0.3, fill="green4") +
geom_ribbon(data=subset(dat3,times>times[1] & times<times_ext3),aes(x=.data$x,
                                ymax=.data$y),ymin=0,alpha=0.3, fill="blue") +
geom_vline(xintercept = times_ext1, size=.5, color="brown",linetype="dashed")+
geom_vline(xintercept = times_ext2, size=.5, color="green4",linetype="dashed")+
geom_vline(xintercept = times_ext3, size=.5, color="blue",linetype="dashed")+
geom_line(data =subset(data1,times>times[1] & times<times_ext1), color = "brown")+
geom_line(data =subset(data2,times>times[1] & times<times_ext2), color = "green4")+
geom_line(data =subset(data3,times>times[1] & times<times_ext3), color = "blue")+
labs(x = "Time",y="Abundance")+
theme(plot.title = element_text(size=40))+
theme(plot.title = element_text(hjust = 0.5))+
theme(axis.title.y = element_text(size = rel(1), angle = 90))+
theme(axis.title.x = element_text(size = rel(1), angle = 00))+
labs(tag = "(a)")


p2 <- ggplot(da, aes(x=.data$x, y=.data$y)) +
theme_bw()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_vline(xintercept = times_ext1, size=.5, color="brown",linetype="dashed")+
geom_vline(xintercept = times_ext2, size=.5, color="green4",linetype="dashed")+
geom_vline(xintercept = times_ext3, size=.5, color="blue",linetype="dashed")+
geom_line(data =subset(da1,times>times[1] & times<times_ext1), color = "brown")+
geom_line(data =subset(da2,times>times[1] & times<times_ext2), color = "green4")+
geom_line(data =subset(da3,times>times[1] & times<times_ext3), color = "blue")+
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


}








