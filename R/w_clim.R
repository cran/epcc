#' Temperature trend obtained from WorldClim.
#'
#' @description This function allows simulating the effect of temperature trends on the
#'              abundance of ectotherm populations in different geographic locations.
#'              Temperature data is obtained from WorldClim.
#'
#'
#'
#'@param y_ini Initial population values (must be written with its name: N).
#'@param temp_cmin Minimum critical temperature.
#'@param temp_cmax Maximum critical temperature.
#'@param ro Population growth rate at optimum temperature.
#'@param lambda Marginal loss by non-thermodependent intraspecific competition.
#'@param lat,lon Geographical coordinates of some place of interest of study, in decimal degrees.
#'@param s Bioclimatic variable.
#'@param res Spatial resolution.
#'@param time_start Start of time sequence.
#'@param time_end End of time sequence.
#'@param leap Sequence increase.
#'
#'@details Three populations and/or scenarios can be simulated simultaneously. The temperature trends
#'         are obtained by data extracted from WorldClim for the years 2000, 2050 and 2070 at a specific
#'         location. The function internally calls the function getData of the raster package (Hijmans, 2020) to
#'         obtain the bioclimatic variable of interest a given spatial resolution. An exponential
#'         expression is fitted using the nls function. In each input vector, the parameters for
#'         the three simulations must be specified (finite numbers for the initial population abundance).
#'         The simulations are obtained by a model that incorporates the effects of temperature over time,
#'         which leads to a non-autonomous ODE approach. This is function uses the ODE solver implemented
#'         in the package deSolve (Soetaert et al., 2010). In the first three examples, three geographic locations are
#'         considered for Macrolophus pygmaeus as reported in Sánchez et al. (2012).
#'
#'
#'
#'@return (1) A data.frame with columns having the simulated trends.
#'@return (2) A two-panel figure in which (a) shows the population abundance curves represented by solid lines
#'            and the corresponding carrying capacities are represented by shaded areas. In (b) the temperature
#'            trend is shown. The three simultaneous simulations are depicted by different colors, i.e. 1st brown,
#'            2nd green and 3rd blue.
#'
#'@references Hijmans, R.J. (2020). Package `raster' (Version 3.3-13). pp. 1-249.
#'
#'@references Rezende, E. L., & Bozinovic, F. (2019). Thermal performance across levels of biological
#'            organization. Philosophical Transactions of the Royal Society B: Biological Sciences,
#'            374(1778), 20180549.doi:10.1098/rstb.2018.0549
#'
#'@references Sanchez, J. A., Spina, M. L., & Perera, O. P. (2012). Analysis of the population
#'            structure of Macrolophus pygmaeus (Rambur) (Hemiptera: Miridae) in the Palaearctic
#'            region using microsatellite markers. Ecology and Evolution, 2(12), 3145-3159.
#'            doi:10.1002/ece3.420
#'
#'@references Soetaert, K., Petzoldt, T., & Setzer, R. (2010). Solving Differential Equations in R:
#'            Package deSolve. Journal of Statistical Software, 33(9), 1 - 25.
#'            doi:http://dx.doi.org/10.18637/jss.v033.i09
#'
#'@export
#'@examples
#'\dontrun{
#'#######################################################################
#'   #Example 1: Different initial population abundances.
#'#######################################################################
#'
#'w_clim(y_ini = c(N = 100, N = 200, N = 400),
#'       temp_cmin = rep(18,3),
#'       temp_cmax = rep(30,3),
#'       ro = rep(0.7,3),
#'       lambda = rep(0.00005,3),
#'       lat = rep(-33,3),
#'       lon = rep(-71,3),
#'       s = 5,
#'       res = 5,
#'       time_start = 2000,
#'       time_end = 2070,
#'       leap = 1/12)
#'
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
#'w_clim(y_ini = c(N = 100, N = 100, N = 100),
#'       temp_cmin = c(temp_cmin1,temp_cmin2,temp_cmin3),
#'       temp_cmax = c(temp_cmax1,temp_cmax2,temp_cmax3),
#'       ro = rep(0.7,3),
#'       lambda = rep(0.00005,3),
#'       lat = rep(-33,3),
#'       lon = rep(-71,3),
#'       s = 5,
#'       res = 5,
#'       time_start = 2000,
#'       time_end = 2070,
#'       leap = 1/12)
#'
#'#######################################################################
#'   #Example 3: Different latitudes.
#'#######################################################################
#'
#'lat1 <- -10
#'lat2 <- -33
#'lat3 <- -42
#'
#'w_clim(y_ini = c(N = 100, N = 100, N = 100),
#'       temp_cmin = rep(18,3),
#'       temp_cmax = rep(40,3),
#'       ro = rep(0.7,3),
#'       lambda = rep(0.00005,3),
#'       lat = c(lat1,lat2,lat3),
#'       lon = rep(-71,3),
#'       s = 5,
#'       res = 5,
#'       time_start = 2000,
#'       time_end = 2070,
#'       leap = 1/12)
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
#' w_clim(y_ini = c(N = 100, N = 100, N = 100),
#'        temp_cmin = rep(18,3),
#'        temp_cmax = rep(30,3),
#'        ro = rep(0.7,3),
#'        lambda = c(lambda1,lambda2,lambda3),
#'        lat = rep(-33,3),
#'        lon = rep(-71,3),
#'        s = 5,
#'        res = 5,
#'        time_start = 2000,
#'        time_end = 2070,
#'        leap = 1/12)
#'
#########################################################################
#'
#'#######################################################################
#'   #Application example I: Bioclimatic variable
#'   #                       (Annual Mean Temperature).
#'#######################################################################
#'
#'#We consider a population of Macrolophus pygmaeus in three different
#'#locations, and its intrinsic growth rate is adjusted to data obtained
#'#from Rezende and Bozinovic (2019).
#'
#'
#'github_link <- "https://github.com/Victor-Saldana/epcc/raw/main/M_pygmaeus.xlsx"
#'library(httr)
#'temp_file <- tempfile(fileext = ".xlsx")
#'req <- GET(github_link,
#'           authenticate(Sys.getenv("GITHUB_PAT"), ""),
#'           write_disk(path = temp_file))
#'M_pygmaeus <- readxl::read_excel(temp_file)
#'
#'TPC <- rate_adjustment(data = M_pygmaeus)
#'
#'#locality 1
#'lat1 <- 38.1827778
#'lon1 <- -1.7380555
#'
#'#locality 2
#'lat2 <- 41.01384
#'lon2 <- 28.94966
#'
#'#locality 3
#'lat3 <- 39.7213889
#'lon3 <- 21.63416638888889
#'
#' w_clim(y_ini = c(N = 100, N = 100, N = 100),
#'        temp_cmin = rep(TPC$temp_cmin,3),
#'        temp_cmax = rep(TPC$temp_cmax,3),
#'        ro = rep(TPC$ro,3),
#'        lambda = rep(0.00005,3),
#'        lat = c(lat1,lat2,lat3),
#'        lon = c(lon1,lon2,lon3),
#'        s = 1,
#'        res = 5,
#'        time_start = 2000,
#'        time_end = 2070,
#'        leap = 1/12)
#'
#'#######################################################################
#'   #Application example II: Bioclimatic variable
#'   #                        (Max Temperature of Warmest Month).
#'#######################################################################
#'
#'#We consider a population of Macrolophus pygmaeus in three different
#'#locations, and its intrinsic growth rate is adjusted to data obtained
#'#from Rezende and Bozinovic (2019).
#'
#'
#'github_link <- "https://github.com/Victor-Saldana/epcc/raw/main/M_pygmaeus.xlsx"
#'library(httr)
#'temp_file <- tempfile(fileext = ".xlsx")
#'req <- GET(github_link,
#'           authenticate(Sys.getenv("GITHUB_PAT"), ""),
#'           write_disk(path = temp_file))
#'M_pygmaeus <- readxl::read_excel(temp_file)
#'
#'TPC <- rate_adjustment(data = M_pygmaeus)
#'
#'#locality 1
#'lat1 <- 38.1827778
#'lon1 <- -1.7380555
#'
#'#locality 2
#'lat2 <- 41.01384
#'lon2 <- 28.94966
#'
#'#locality 3
#'lat3 <- 39.7213889
#'lon3 <- 21.63416638888889
#'
#' w_clim(y_ini = c(N = 100, N = 100, N = 100),
#'        temp_cmin = rep(TPC$temp_cmin,3),
#'        temp_cmax = rep(TPC$temp_cmax,3),
#'        ro = rep(TPC$ro,3),
#'        lambda = rep(0.00005,3),
#'        lat = c(lat1,lat2,lat3),
#'        lon = c(lon1,lon2,lon3),
#'        s = 5,
#'        res = 5,
#'        time_start = 2000,
#'        time_end = 2070,
#'        leap = 1/12)
#'
#'#######################################################################
#'   #Application example III: Bioclimatic variable
#'   #                         (Mean Temperature of Warmest Quarter).
#'#######################################################################
#'
#'#We consider a population of Macrolophus pygmaeus in three different
#'#locations, and its intrinsic growth rate is adjusted to data obtained
#'#from Rezende and Bozinovic (2019).
#'
#'
#'github_link <- "https://github.com/Victor-Saldana/epcc/raw/main/M_pygmaeus.xlsx"
#'library(httr)
#'temp_file <- tempfile(fileext = ".xlsx")
#'req <- GET(github_link,
#'           authenticate(Sys.getenv("GITHUB_PAT"), ""),
#'           write_disk(path = temp_file))
#'M_pygmaeus <- readxl::read_excel(temp_file)
#'
#'TPC <- rate_adjustment(data = M_pygmaeus)
#'
#'#locality 1
#'lat1 <- 38.1827778
#'lon1 <- -1.7380555
#'
#'#locality 2
#'lat2 <- 41.01384
#'lon2 <- 28.94966
#'
#'#locality 3
#'lat3 <- 39.7213889
#'lon3 <- 21.63416638888889
#'
#' w_clim(y_ini = c(N = 100, N = 100, N = 100),
#'        temp_cmin = rep(TPC$temp_cmin,3),
#'        temp_cmax = rep(TPC$temp_cmax,3),
#'        ro = rep(TPC$ro,3),
#'        lambda = rep(0.00005,3),
#'        lat = c(lat1,lat2,lat3),
#'        lon = c(lon1,lon2,lon3),
#'        s = 10,
#'        res = 5,
#'        time_start = 2000,
#'        time_end = 2070,
#'        leap = 1/12)
#'}


w_clim<- function(y_ini = c(N = 400, N = 400, N = 400),
                       temp_cmin = rep(18,3),
                       temp_cmax = c(25,28,32),
                       ro = rep(0.7,3),
                       lambda = rep(0.00005,3),
                       lat = rep(-33,3),
                       lon = rep(-71,3),
                       s = 10,
                       res = 5,
                       time_start = 2000,
                       time_end = 2070,
                       leap = 1/12){

times<- seq(time_start, time_end, leap)

if(time_end<=2100){
if(time_start<=time_end){

if(temp_cmin[1]<temp_cmax[1] && temp_cmin[2]<temp_cmax[2] && temp_cmin[3]<temp_cmax[3] ){



wclimCurr <- getData("worldclim", var="bio", res=res)
wclimCurr <- wclimCurr[[s]]
plot(wclimCurr,cex.axis=1.2, tcl=-0.7,las=1)
axis(side=1,at=c(-200,200),col="black",lwd=3)
axis(side=2,at=c(-200,200),col="black",lwd=3)
axis(side=3,at=c(-200,200),col="black",lwd=3)
axis(side=4,at=c(-200,200),col="black",lwd=3)

wclim50 <- getData("CMIP5", var="bio", rcp=85, res=res, model="AC", year="50")
wclim50 <- wclim50[[s]]
wclim70 <- getData("CMIP5", var="bio", rcp=85, res=res, model="AC", year="70")
wclim70 <- wclim70[[s]]


coord1 <- data.frame(x=lon[1], y=lat[1])
coord2 <- data.frame(x=lon[2], y=lat[2])
coord3 <- data.frame(x=lon[3], y=lat[3])

point1 <- SpatialPoints(coord1, proj4string=wclimCurr@crs)
point2 <- SpatialPoints(coord2, proj4string=wclimCurr@crs)
point3 <- SpatialPoints(coord3, proj4string=wclimCurr@crs)
plot(point1, add=TRUE,col="brown")
plot(point2, add=TRUE,col="red")
plot(point3, add=TRUE,col="black")

valueCurr1 <- extract(wclimCurr, point1)
valueCurr2 <- extract(wclimCurr, point2)
valueCurr3 <- extract(wclimCurr, point3)

value50_1 <- extract(wclim50, point1)
value70_1 <- extract(wclim70, point1)

value50_2 <- extract(wclim50, point2)
value70_2 <- extract(wclim70, point2)

value50_3 <- extract(wclim50, point3)
value70_3 <- extract(wclim70, point3)

values1 <- c(valueCurr1, value50_1, value70_1)
values1 <- values1/10

if (is.na(values1[1])== FALSE){

values2 <- c(valueCurr2, value50_2, value70_2)
values2 <- values2/10

if (is.na(values2[1])== FALSE){

values3 <- c(valueCurr3, value50_3, value70_3)
values3 <- values3/10

if (is.na(values3[1])== FALSE){

x<- c(2000,2050,2070)
y1<- values1
y2<- values2
y3<- values3
df1 <- data.frame(x, y1)
df2 <- data.frame(x, y2)
df3 <- data.frame(x, y3)
f <- function(times,a,b) {a * exp(b * times)}

m1<- nls(y1 ~ exp(loga + b * x), df1, start = list( loga = log(2),
                        b = 0.005),control = list (maxiter = 500))
y_est1<-predict(m1,df1$x)

m2<- nls(y2 ~ exp(logw + z * x), df2, start = list(logw = log(2),
                        z = 0.005), control = list (maxiter = 500))
y_est2<-predict(m2,df2$x)

m3<- nls(y3 ~ exp(loge + g * x), df3, start = list(loge = log(2),
                        g = 0.005),control = list (maxiter = 500))
y_est3<-predict(m3,df3$x)

a1=exp(coef(m1)[1])
b1=coef(m1)[2]

a2=exp(coef(m2)[1])
b2=coef(m2)[2]

a3=exp(coef(m3)[1])
b3=coef(m3)[2]
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
  # Time
##########################################################
time_op1= suppressWarnings(1/b1*log(temp_op1/a1))
time_cmax1= suppressWarnings(1/b1*log(temp_cmax[1]/a1))

#########################################################
time_op2= suppressWarnings(1/b2*log(temp_op2/a2))
time_cmax2= suppressWarnings(1/b2*log(temp_cmax[2]/a2))

#########################################################
time_op3= suppressWarnings(1/b3*log(temp_op3/a3))
time_cmax3= suppressWarnings(1/b3*log(temp_cmax[3]/a3))

#########################################################
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

##############################################################
    # Parameters
##############################################################

parms1<-c(temp_cmin[1],temp_cmax[1],temp_op1,ro[1], lambda[1])
parms2<-c(temp_cmin[2],temp_cmax[2],temp_op2,ro[2], lambda[2])
parms3<-c(temp_cmin[3],temp_cmax[3],temp_op3,ro[3], lambda[3])

###############################################################

###############################################################
      # Model for each trend
###############################################################

model1 <- function (times, y,parms1) {
  with(as.list(c(y)), {
  T1<-f(times,a=exp(coef(m1)[1]), b=coef(m1)[2])
  r1<- rate_TPC(T1,ro[1],temp_cmin[1],temp_cmax[1],temp_op1)
  dN <-   r1 * N * (1 - lambda[1]*(N / r1))
  list(dN,T1,r1) })
}

###############################################################

model2 <- function (times, y,parms2) {
  with(as.list(c(y)), {
  T2<-f(times,a=exp(coef(m2)[1]), b=coef(m2)[2])
  r2<- rate_TPC(T2,ro[2],temp_cmin[2],temp_cmax[2],temp_op2)
  dN <-   r2 * N * (1 - lambda[2]*(N / r2))
  list(dN,T2,r2)})
}

###############################################################

model3 <- function (times, y,parms3) {
with(as.list(c(y)), {
T3<-f(times,a=exp(coef(m3)[1]), b=coef(m3)[2])
r3<- rate_TPC(T3,ro[3],temp_cmin[3],temp_cmax[3],temp_op3)
dN <-   r3 * N * (1 - lambda[3]*(N / r3))
list(dN,T3,r3)})
}

###############################################################
###############################################################
  # Solution
##############################################################

out1 <- ode(y=y_ini[1], times, model1, parms1,method = "ode45")
out2 <- ode(y=y_ini[2], times, model2, parms2,method = "ode45")
out3 <- ode(y=y_ini[3], times, model3, parms3,method = "ode45")

##############################################################
##############################################################
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

data<-rbind(data1,data2,data3,dat1,dat2,dat3,da1,da2,da3)

p1<- ggplot(data, aes(x=.data$x, y=.data$y)) +
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


p2<-ggplot(data, aes(x=.data$x, y=.data$y)) +
theme_bw()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_vline(xintercept = time_ext1, size=.5, color="brown",linetype="dashed")+
geom_vline(xintercept = time_ext2, size=.5, color="green4",linetype="dashed")+
geom_vline(xintercept = time_ext3, size=.5, color="blue",linetype="dashed")+
geom_line(data =subset(da1,times>times[1] & times<time_ext1), color = "brown")+
geom_line(data =subset(da2,times>times[1] & times<time_ext2), color = "green4")+
geom_line(data =subset(da3,times>times[1] & times<time_ext3), color = "blue")+
labs(x = "Time",y="Temperature")+
theme(axis.title.y = element_text(size = rel(1), angle = 90))+
theme(axis.title.x = element_text(size = rel(1), angle = 00))+
labs(tag = "(b)")

plot_grid(p1, p2)


}else{

stop("No information is recorded for the entered coordinates (number 3), enter new coordinates.")
   }


 }else{

stop("No information is recorded for the entered coordinates (number 2), enter new coordinates.")
      }

}else{

stop("No information is recorded for the entered coordinates (number 1), enter new coordinates.")
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





