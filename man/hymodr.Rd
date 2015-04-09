\name{hymodr}
\alias{hymodr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An R-package for the hydrological model HyMOD (R implementation)
}
\description{
R implementation of the hydrological model HyMOD.
}
\usage{
hymodr(param, area, tdelta, e, p,
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

\details{
HyMOD is a five-parameter conceptual rainfall-runoff model that was introduced by Boyle 2000 [1]. The model is based on the probability-distributed soil storage capacity principle introduced by R. J. Moore 1985 [2].
In HyMOD, the rainfall-runoff process is represented through a nonlinear tank connected with three identical linear tanks in parallel representing the surface flow and a slow-flow tank representing groundwater flow. The model requires the optimization of five parameters: \code{C_max} (the maximum storage capacity within the watershed), \code{Beta} (the degree of spatial variability of the soil moisture capacity within the watershed), \code{Alpha} (a factor for partitioning the flow between two series of tanks) and the two residence time parameters of quick-flow and slow-flow tanks, \code{k_quick} and \code{k_slow} respectively. 
Given its relatively simplicity but robustness, the model is widely used in the literature. The input data consist of precipitation and potential evapotranspiration at the given time scale.

[1] Boyle D.P. (2000), Multicriteria calibration of hydrological models, Ph.D. dissertation, Dep. of Hydrol. and Water Resour., Univ of Arizona, Tucson.

[2] Moore, R.J., (1985), The probability-distributed principle and runoff production at point and basin scale, Hydrol. Sci. J., 30(2), 273-297.
}



\value{
The outfut of the function hymodr returns a list containing the following elements:
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{w_tot}{simulated water (storage) in the catchment [mm/time_step]}
 \item{ep_tot}{simulated actual evapotranspiration [mm/time_step]} 
 \item{q_slow}{simulated slow flow response (i.e. component of the total runoff derives from the single linear reservoir) [m^3 s^-1]}
 \item{q_quick}{simulated quick flow response (i.e. component of the total runoff derives from the three linear reservoirs) [m^3 s^-1]}
 \item{q_tot}{total outflow of the catchment [m^3 s^-1]}
}
\examples{

    #################################
    ## Hourly data
    #################################

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

    ## Simulated output
    out_sim<-hymodf(param=param,
                    area=area,
                    tdelta=3600,
                    e=evapo,
                    p=precipit)
    
    ## Plot
    time<-sievefornacina[,1]
    plot(time,qoss,ty="l",lwd=0.7,col="blue",ylab="Discharge [m^3 s-1]")
    title(main="Sieve@Fornacina - Hourly data",cex.main=1.2,col.main="darkblue")
    lines(time,out_sim$q_tot,col="green4")
    legend("topleft",legend=c("Observed","Simulated"),lty=c(1,1),
    col=c("blue","green4"),bty="n",ncol=1)



    #################################
    ## Daily data
    #################################
    
    ## Load data
    data(arnosubbiano)
    
    ## Catchment area [km^2]
    area<-751  
    
    ## Precipitation [mm/day]
    precipit<-arnosubbiano[,2]
    
    ## Potential arnosubbiano [mm/day]
    evapo<-arnosubbiano[,3]
    
    ## Observed discharge [m^3/s]
    qoss<-arnosubbiano[,4]
    
    ## Trial parameter set
    param<-c(601,0.1872,0.3778,0.03956, 0.74985)

    ## Simulated output
    out_sim<-hymodr(param=param,
                     area=area,
                     tdelta=86400,
                     e=evapo,
                     p=precipit)
    
    ## Plot
    time<-arnosubbiano[,1]
    plot(time,qoss,ty="l",lwd=0.7,col="blue",ylab="Discharge [m^3 s-1]")
    title(main="Arno@Subbiano - Daily data",cex.main=1.2,col.main="darkblue")
    lines(time,out_sim$q_tot,col="green4")
    legend("topleft",legend=c("Observed","Simulated"),lty=c(1,1),
    col=c("blue","green4"),bty="n",ncol=1)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
