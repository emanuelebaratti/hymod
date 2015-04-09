\name{hymodf}
\alias{hymodf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An R-package for the hydrological model HyMOD (FORTRAN implementation)
}
\description{
FORTRAN implementation of the hydrological model HyMOD (see \code{hymodr}). The function combines the R user friendly interface with a fast processing unit.
}
\usage{
hymodf(param, area, tdelta, e, p,
       w_initial=0, wslow_initial = 0.0, wquick_initial = 0.0)
}
\arguments{
  \item{param}{
vector of parameters: 
\enumerate{
\item \code{C_max} {maximum water storage capacity of the catchment [mm]}
\item \code{Beta} {degree of spatial variability in the water storage capacities \code{C_max} [/]}
\item \code{Alpha} {partition parameter distributing the flow between slow and quick component [/]}
\item \code{k_slow} {storage coefficient for slow response (i.e. the linear slow reservoir) [1/time_step]}
\item \code{k_quick} {storage coefficient for quick response (i.e. the three linear quick reservoirs) [1/time_step]}
}
}
  \item{area}{watershed area [km^2]}
  \item{tdelta}{conversion factor (used to express the quantities in mm). 
  It has to be set \code{tdelta}=3600 for hourly data, \code{tdelta}=86400 for daily data}
  \item{e}{vector of potential evapotranspiration [mm/time_step]}
  \item{p}{vector of precipitation for each time step [mm/time_step]}
  \item{w_initial}{initial water (storage) in the catchment per unit area [mm]. Defaults is 0 mm} 
  \item{wslow_initial}{initial value for slow response storage [mm]. Defaults is 0 mm} 
  \item{wquick_initial}{initial value for quick response storage [mm]. Defaults is 0 mm}
}
\value{
The outfut of the function hymodf returns a list containing the following elements:
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{w_tot}{simulated water (storage) in the catchment [mm/time_step]}
 \item{ep_tot}{simulated actual evapotranspiration [mm/time_step]} 
 \item{q_slow}{simulated slow flow response (i.e. component of the total runoff derives from the single linear reservoir) [m^3 s^-1]}
 \item{q_quick}{simulated quick flow response (i.e. component of the total runoff derives from the three linear reservoirs) [m^3 s^-1]}
 \item{q_tot}{total outflow of the catchment [m^3 s^-1]}
}
\examples{
## Load data
data(sievefornacina)

## Catchment area [km^2]
area<-830  

## Load data
data(sievefornacina)

## Precipitation [mm/hour]
precipit<-sievefornacina[,2]

## Potential evapotranspiration [mm/hour]
evapo<-sievefornacina[,3]

## Observed discharge [m^3/s]
qoss<-sievefornacina[,4]

## Trial parameter set
param<-c(355,0.2860,0.1591,0.0017,0.2006) 

## Time for run the model 10 times
comput_time_fortran<-vector()
for ( i in 1:10){
## Simulated output
comput_time_fortran[i]<-system.time({
out_sim<-hymodf(param=param,
                area=area,
                tdelta=3600,
                e=evapo,
                p=precipit)
})[3]

}

comput_time_R<-vector()
for ( i in 1:10){
  ## Simulated output
  comput_time_R[i]<-system.time({
    out_sim<-hymodr(param=param,
                    area=area,
                    tdelta=3600,
                    e=evapo,
                    p=precipit)
  })[3]
  
}

## Plot
ylimiti<-range(comput_time_fortran,comput_time_R)
plot(comput_time_R,col="blue",ylab="Elapsed Time [s]",
xlab="Iteration number",ylim=ylimiti)
points(comput_time_fortran,col="red")
title(main="Sieve@Fornacina - Hourly data \n Elapsed Time for 10 run of the model",
cex.main=1.2,col.main="darkblue")
legend("left",legend=c("R-implementation","FORTRAN implementation"),
pch=c(1,1),col=c("blue","red"),ncol=1,bty="n")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
