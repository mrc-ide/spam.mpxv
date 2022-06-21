# `spam.mpxv`: Stochastic pair approximation model for Monkeypox virus

Population structure is modelled as a dynamic network based on
behavioural surveillance of MSM. Each individual has a mixture of
long-duration and short-duration(casual) contacts as well as one-off
contacts.

Degree distributions and partnership durations are based on: [Anderson
et al., Epidemiology
2021](http://dx.doi.org/10.1097/EDE.0000000000001390)

Relative contact rates for different contact types are based on [Jenness
et al., JID 2016](https://doi.org/10.1093%2Finfdis%2Fjiw223)

The model is designed within the ‘edge-based compartmental’(EBCM)
framework laid out in

-   Miller JC, Volz EM. Model hierarchies in edge-based compartmental
    modeling for infectious disease spread. J Math Biol. 2013;67:
    869–899. <doi:10.1007/s00285-012-0572-3>
-   Volz E, Meyers LA. Susceptible–infected–recovered epidemics in
    dynamic contact networks. Proceedings of the Royal Society B. 2007.
    Available:
    <https://royalsocietypublishing.org/doi/abs/10.1098/rspb.2007.1159>

A system of stochastic differential equations were designed around the
EBCM by introducing Guassian noise on 1) the force of infection and 2)
the average degree among new infected individuals.

The model is designed to work with the `pomp` package for simulation and
model fitting.

## Example simulation

    library( spam.mpxv )
    set.seed( 1111 ) 
    st0 <- Sys.time() 
    # Model m1.0 implements the first basic model with three contact types 
    # Default parameters are not calibrated and only for demonstration
    s1 <- pomp::simulate( m3.0 )
    st1 <- Sys.time() 
    print( st1 - st0 )

    ## Time difference of 0.01085997 secs

Reproduction number

    calcR0.3.0(beta = 2.25, gamma1 = 1/4) 

    ## [1] 1.564376

    plot( s1 ) 

![](README_files/figure-markdown_strict/plot-1.png)![](README_files/figure-markdown_strict/plot-2.png)![](README_files/figure-markdown_strict/plot-3.png)

Example simulated data:

    s1d <- as.data.frame( s1 ) 
    print( tail( s1d ))

    ##    time Ytravel Yendog Yunk    thetaf         MSEf          MEf      MSSf         MSIf          MIf    thetag         MSEg          MEg
    ## 85   85       0      5    0 0.9998131 9.170162e-05 0.0001514393 0.9993615 1.511577e-05 6.324420e-05 0.9995270 0.0002672865 0.0002543847
    ## 86   86       0      5    0 0.9998024 9.173478e-05 0.0001544219 0.9993407 1.606797e-05 6.636306e-05 0.9995097 0.0002803447 0.0002666551
    ## 87   87       1      5    0 0.9997917 9.982899e-05 0.0001650752 0.9993119 1.635344e-05 6.907504e-05 0.9994924 0.0002969041 0.0002824125
    ## 88   88       0      5    0 0.9997757 8.973335e-05 0.0001581575 0.9992994 1.746761e-05 7.244068e-05 0.9994786 0.0002866366 0.0002712892
    ## 89   89       0      2    0 0.9997650 9.114120e-05 0.0001611069 0.9992778 1.653199e-05 7.410020e-05 0.9994579 0.0002924269 0.0002764542
    ## 90   90       0      2    0 0.9997650 9.570996e-05 0.0001671799 0.9992526 1.647858e-05 7.571351e-05 0.9994337 0.0002984528 0.0002818882
    ##         MSSg         MSIg          MIg    thetah          MEh          MIh        E        I     newI cutf cutg cuth cuts  seedrate
    ## 85 0.9987887 5.548255e-05 0.0001009370 0.9830379 0.0003860357 0.0001774621 72.12746 31.61915 9.015933   35  137   69   75 0.8608633
    ## 86 0.9987367 5.864850e-05 0.0001075009 0.9826876 0.0003658050 0.0001813510 73.11153 32.85330 9.138941   37  142   72   75 0.7536780
    ## 87 0.9986794 6.173658e-05 0.0001139575 0.9819873 0.0003957859 0.0001817389 76.97259 34.26155 9.621573   39  147   78   75 0.7162264
    ## 88 0.9986469 6.522938e-05 0.0001207697 0.9816373 0.0003871005 0.0001857774 77.35101 35.36504 9.668877   42  151   81   75 0.7354071
    ## 89 0.9985993 6.554290e-05 0.0001244884 0.9812875 0.0004145907 0.0001877206 78.68214 36.35905 9.835267   44  157   84   75 0.8342455
    ## 90 0.9985503 6.647204e-05 0.0001279231 0.9808213 0.0003547351 0.0001926143 79.84687 37.25014 9.980859   44  164   85   78 0.7713144

How long does it take to run a particle filter on the simulated data?

    print(
        system.time(pfilter(s1,Np=100))
    )

    ##    user  system elapsed 
    ##   0.884   0.000   0.879

## Version 0.0.3 2022-06-20

-   Revise natural history parameters, SEIR
