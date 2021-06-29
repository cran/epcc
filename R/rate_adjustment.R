#'Intrinsic growth rate adjustment
#'
#' @description This function allows you to adjust the intrinsic growth rate of an ectothermic
#'             population using temperature and growth rate data obtained empirically, using a
#'             cubic TPC (Saldaña et al., 2019).
#'
#'@param data database where the first column shows the ambient
#'            temperature (TA) and the second column contains the
#'            intrinsic growth rate (R) associated with them.
#'
#'
#'@details This function allows you to adjust the intrinsic growth rate of an ectotherm
#'         population by providing as input the values of environmental or body temperature
#'         together with the growth rate, using the nls2 function you can find the parameters
#'         temp_cmax, temp_cmax and ro, which are necessary to adjust the curve through a cubic polynomial
#'         which fulfills the essential conditions of a TPC.
#'
#'@return A figure showing the fitting curve corresponding to the intrinsic growth rate of an ectothermic
#'        population, the empirically obtained temperature and growth rate data points, in addition,
#'        called for side effects.
#'
#'@export
#'@examples
#'
#'######################################################################
#'   #Example 1: We consider a population of Macrolophus pygmaeus whose
#'   #intrinsic growth rate is adjusted to the data obtained from Rezende
#'   #and Bozinovic (2019).
#'#######################################################################
#'
#'github_link <- "https://github.com/Victor-Saldana/epcc/raw/main/M_pygmaeus.xlsx"
#'library(httr)
#'temp_file <- tempfile(fileext = ".xlsx")
#'req <- GET(github_link,
#'           authenticate(Sys.getenv("GITHUB_PAT"), ""),
#'           write_disk(path = temp_file))
#'M_pygmaeus <- readxl::read_excel(temp_file)
#'TPC <- rate_adjustment(data = M_pygmaeus)
#'
#'######################################################################
#'   #Example 2: We consider a population of Eretmocerus furuhashii whose
#'   #intrinsic growth rate is adjusted to the data obtained from Rezende
#'   #and Bozinovic (2019).
#'#######################################################################
#'
#'github_link <- "https://github.com/Victor-Saldana/epcc/raw/main/E_furuhashii.xlsx"
#'library(httr)
#'temp_file <- tempfile(fileext = ".xlsx")
#'req <- GET(github_link,
#'           authenticate(Sys.getenv("GITHUB_PAT"), ""),
#'           write_disk(path = temp_file))
#'E_furuhashii <- readxl::read_excel(temp_file)
#'TPC <- rate_adjustment(data = E_furuhashii)
#'
#'######################################################################
#'   #Example 3: We consider a population of Trichogramma pretoisum whose
#'   #intrinsic growth rate is adjusted to the data obtained from Rezende
#'   #and Bozinovic (2019).
#'#######################################################################
#'
#'github_link <- "https://github.com/Victor-Saldana/epcc/raw/main/T_pretoisum.xlsx"
#'library(httr)
#'temp_file <- tempfile(fileext = ".xlsx")
#'req <- GET(github_link,
#'           authenticate(Sys.getenv("GITHUB_PAT"), ""),
#'           write_disk(path = temp_file))
#'T_pretoisum <- readxl::read_excel(temp_file)
#'TPC <- rate_adjustment(data = T_pretoisum)
#'
#'
#'@references Rezende, E. L., & Bozinovic, F. (2019). Thermal performance across levels of biological
#'            organization. Philosophical Transactions of the Royal Society B: Biological Sciences,
#'            374(1778), 20180549.doi:10.1098/rstb.2018.0549
#'
#'@references Saldaña-Núñez, V.N., Córdova-Lepe, F.D. & Moreno-Gómez, F.N. (2019). Population dynamics in the face of
#'            climate change: Analysis of a cubic thermal performance curve in ectotherms. J. Phys.: Conf.
#'            Ser. 1329 012007.   doi:10.1088/1742-6596/1329/1/012007
#'





rate_adjustment<-function(data = data){


  rate <- function(T,temp_cmin,temp_cmax,ro){ro*T*(T-temp_cmin)*(temp_cmax-T)/(((temp_cmin+temp_cmax+sqrt((temp_cmin+temp_cmax)^2-
          3*temp_cmin*temp_cmax))/3)*(((temp_cmin+temp_cmax+sqrt((temp_cmin+temp_cmax)^2-3*temp_cmin*temp_cmax))/3)-temp_cmin)*
          (temp_cmax-((temp_cmin+temp_cmax+sqrt((temp_cmin+temp_cmax)^2-3*temp_cmin*temp_cmax))/3)))}

  m <- nls2(data$R~rate(data$TA,temp_cmin,temp_cmax,ro),data=data,start=list(temp_cmin=10,temp_cmax=30,ro=0.2))

  temp_cmin<- coefficients(m)[1]
  temp_cmax<- coefficients(m)[2]
  ro<-coefficients(m)[3]
  Top<- (temp_cmin+temp_cmax+sqrt((temp_cmin+temp_cmax)^2-3*temp_cmin*temp_cmax))/3

  s<- seq(temp_cmin,temp_cmax)

  r<- ro*s*(s-temp_cmin)*(temp_cmax-s)/(((temp_cmin+temp_cmax+sqrt((temp_cmin+temp_cmax)^2-
      3*temp_cmin*temp_cmax))/3)*(((temp_cmin+temp_cmax+sqrt((temp_cmin+temp_cmax)^2-
      3*temp_cmin*temp_cmax))/3)-temp_cmin)*(temp_cmax-((temp_cmin+temp_cmax+
      sqrt((temp_cmin+temp_cmax)^2-3*temp_cmin*temp_cmax))/3)))



  plot(s,r,xlab="Temperature", ylab="r(T)",
  main="Intrinsic growth rate", xlim = c(temp_cmin, temp_cmax),
  ylim = c(0, ro), type="l",lty=c(1),cex.axis=1.5,
  tcl=-0.7,las=0, cex.main=2,bty="n",cex=1.5,lwd=2,cex.lab=1.5,lwd.ticks=2)
  par(new=TRUE)
  plot(data$TA,data$R,xlab = "",ylab = "", xlim = c(temp_cmin, temp_cmax), ylim = c(0, ro), axes=FALSE,lwd=2)
  axis(side=1,at=c(-20,100),col="black",lwd=3)
  axis(side=2,at=c(-20,100),col="black",lwd=3)
  as.list(c(temp_cmin,temp_cmax,ro))
}



