#~ 4.0
#~ - including vaccination prob prop to degree

#~ 3.0
#~ - adding latent period 

#~ 2.1
#~ - simplifying measurement model 

#~ m2: 
#~ -seeding rate gaussian process 
#~ -accumulate variables for fitting 

#~  m1.1: 
#~ Continuous seeding 
#~ Measurement model for seeded cases (travel associated) and endogenous



f <- function(x) (2809 + 2000*x + 95*x^2) / 4904 
fp <- function(x) (2000 + 2*95*x) / 4904
fpp <- function(x) (2*95) / 4904
fppp <- function(x) 0

g <- function(x) (2943 + 1009*x + 477*x^2 + 475*x^3)/4904
gp <- function(x) (1009 + 2*477*x^1 + 3*475*x^2)/4904
gpp <- function(x) (2*477 + 2*3*475*x^1)/4904
gppp <- function(x) (2*3*475)/4904

hshape <- 0.26
hrate <- 1.85 * 7

h <- function(x) (1-log(x)/hrate)^(-hshape)
hp <- function(x) hshape *(1-log(x)/hrate)^(-hshape-1) / (  (hrate*x) )
hpp <- function(x) hshape*((hrate-log(x))/hrate)^(-hshape) * (-hrate+hshape+log(x)+1) / (x^2 * (hrate-log(x))^2)
hppp <- function(x) hshape*((hrate-log(x))/hrate)^(-hshape)  * (hshape^2+3*hshape+2*(hrate-log(x))^2 -3*(hrate-log(x))*(hshape+1)+2 ) / (x^3*(hrate-log(x))^3) 

#' @export
calcR0.4.0 <- function(beta = 0.75, gamma1 = 1/4)
{
	rf <- beta * 1.5/7# factor 1.5/7 is act rate per day; factor 1.5 b/c more contacts per week in long partnerships 
	rg <- beta * 1/7
	pf = rf / (rf + gamma1)
	pg = rg / (rg + gamma1)
	#NOTE no eta terms b/c of separate timescales
	(1+hpp(1)/hp(1)) * beta / gamma1 + pg*gpp(1)/gp(1) + pf*fpp(1)/fp(1) 
	Rff = pf*fpp(1)/fp(1) 
	Rgg = pg*gpp(1)/gp(1)
	Rhh = (1+hpp(1)/hp(1)) * beta / gamma1
	Rgf = Rhf = pf * fp(1) 
	Rfg = Rhg = pg * gp(1) 
	Rfh = Rgh = (beta / gamma1)*hp(1) 
	K = matrix( c( 
		  Rff, Rfg, Rfh
		, Rgf, Rgg, Rgh
		, Rhf, Rhg, Rhh 
	) , nrow = 3, byrow=TRUE )
	eigen(K)$values[1] 
}



.update_theta_vacc4.0 <- function(  theta_vacc, vacc_amt  )
{ #TODO should modify h( thetah x ), not h(x), depending on if we allow those infected to be vaccinated 
	p0 <- h( theta_vacc ) 
	p1 <- max(0.01, p0 - vacc_amt )
	of <- function( lntheta ){
		(h(exp(lntheta)) - p1)^2
	}
	o = optimise( lower = -1e4, upper = 0 , f = of )
	exp( o$minimum  )
}

poispgf <- function(x, z ) exp( z*(x-1))
FGH365 <- function(x,etaf,etag) f(poispgf(x,etaf*365)) * g(poispgf(x,etag*365)) * h(poispgf(x,365))
.update_theta_vacc4.0.365 <- function( theta_vacc, vacc_amt, etaf = 1/200, etag=1/100  )
{
	p0 <- FGH365( theta_vacc, etaf, etag ) 
	p1 <- max(0.01, p0 - vacc_amt )
	of <- function( lntheta ){
		(  FGH365(exp(lntheta),etaf,etag) - p1)^2
	}
	o = optimise( lower = -1e4, upper = 0 , f = of )
	exp( o$minimum  )
}



#' @export
rstep4.0 <- function( thetaf
		, MSEf, MEf
		, MSSf, MSIf, MIf
	, thetag
		, MSEg, MEg
		, MSSg, MSIg, MIg
	, thetah
		, MEh
		, MIh
	, E, I, newI, cutf, cutg, cuth, cuts 
	, seedrate
	, theta_vacc ##
	
	, beta,  gamma0, gamma1, etaf, etag,  N, i0, delta, seedrate0, seedrate_sd
	, vacc_freq, vacc_amt, vacc_start_day, vacc_fin_day ##
	
	, time
	, ...) 
{
	# vacc model 
	## every vacc_freq days inoculate proportion vacc_amt with prob prop to degree
	if (  (time>=vacc_start_day)  &  (time<=vacc_fin_day)  &  (( (time-vacc_start_day) %% vacc_freq) == 0) ){
		theta_vacc <- unname( .update_theta_vacc4.0(  theta_vacc, vacc_amt  ) )
		MSEf <- MSEf * (1-vacc_amt)
		MSSf <- MSSf * (1-vacc_amt)^2
		MSIf <- MSIf * (1-vacc_amt)
		MSEg <- MSEg * (1-vacc_amt)
		MSSg <- MSSg * (1-vacc_amt)^2
		MSIg <- MSIg * (1-vacc_amt)
	}
	.thetah <- thetah * theta_vacc 
	
	seedrate <- max(0, seedrate + rnorm (1, 0, sd = seedrate_sd  ) )
	
	MSf  <- thetaf * fp( thetaf ) / fp(1) 
	MSg  <- thetag * gp( thetag ) / gp(1) 
	MSh  <- .thetah * hp( .thetah ) / hp(1) 
	
	# beta is transm probability per contact * contact rate inflation factor. Transm rates:
	#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5142082/
	rf <- beta * 1.5/7# factor 1.5/7 is act rate per day; factor 1.5 b/c more contacts per week in long partnerships 
	rg <- beta * 1/7
	
	transmf <- rpois(1,MSIf*N*fp(1)*rf )#
	dthetaf <- -thetaf * transmf / (MSf*N*fp(1)) 
	dSf <- fp(thetaf)*dthetaf #note prop to transm
	meanfield_delta_si_f = u1f <- (thetaf * fpp(thetaf) / fp(thetaf) )
	u2f <- (thetaf * fpp(thetaf) + thetaf^2 * fppp(thetaf ) ) / fp(thetaf)
	vf = u2f - u1f^2
	if ( transmf > 0 )
		delta_si_f <- mean( rnorm( transmf, meanfield_delta_si_f , sd = sqrt(vf)) )
	else
		delta_si_f <- 0 
	
	transmg <- rpois(1, MSIg*N*gp(1)*rg)
	dthetag <- -thetag * transmg / (MSg*N*gp(1)) 
	dSg <- gp(thetag)*dthetag #note prop to transm
	meanfield_delta_si_g = u1g <- (thetag * gpp(thetag) / gp(thetag) )
	u2g <- (thetag * gpp(thetag) + thetag^2 * gppp(thetag ) ) / gp(thetag)
	vg = u2g - u1g^2
	if ( transmg > 0 )
		delta_si_g <- mean( rnorm( transmg, meanfield_delta_si_g , sd = sqrt(vg)) )
	else
		delta_si_g <- 0 
	
	
	transmseed <- rpois( 1, seedrate ) 
	transmh <- rpois(1, beta*MIh*N*hp(1) ) 
	dthetah <- -.thetah * (transmh + transmseed ) / (N*hp(1) ) # seeding happens here 
	dSh <- hp(.thetah)*dthetah #note prop to transm
	meanfield_delta_si_h = u1h <- (1 + .thetah * hpp(.thetah) / hp(.thetah) ) #note + 1 for mfsh
	u2h <- hppp(.thetah)*.thetah^2/hp(.thetah) + 2*.thetah*hpp(.thetah)/hp(.thetah) + u1h 
	vh = u2h - u1h^2
	if ( transmh > 0 )
		delta_si_h <- mean( rnorm( transmh, meanfield_delta_si_h , sd = sqrt(vh)) )
	else
		delta_si_h <- 0 
	
	
	dMSEf <- -gamma1*MSEf + 
		+2 * etaf * MSf * MEf + 
		-etaf * MSEf + 
		+ (-dSf) * (delta_si_f/fp(1)) * (MSSf / MSf ) +
		+ (-dSg) * (thetaf*fp(thetaf)/f(thetaf)/fp(1)) * (MSSf / MSf) +
		+ (-dSh) * (thetaf*fp(thetaf)/f(thetaf)/fp(1)) * (MSSf / MSf) 	
	dMSIf <- -rf * MSIf +
		-gamma1*MSIf + 
		+gamma0*MSEf + 
		+2 * etaf * MSf * MIf + 
		-etaf * MSIf + 
		+ (-dSf) * (delta_si_f/fp(1)) * ( -MSIf/MSf) +
		+ (-dSg) * (thetaf*fp(thetaf)/f(thetaf)/fp(1)) * ( - MSIf/MSf) +
		+ (-dSh) * (thetaf*fp(thetaf)/f(thetaf)/fp(1)) * ( - MSIf/MSf) 	
	
	dMSEg <- -gamma0*MSEg + 
		+2 * etag * MSg * MEg + 
		-etag * MSEg + 
		+ (-dSg) * (delta_si_g/gp(1)) * (MSSg / MSg ) +
		+ (-dSf) * (thetag*gp(thetag)/g(thetag)/gp(1)) * (MSSg / MSg) +
		+ (-dSh) * (thetag*gp(thetag)/g(thetag)/gp(1)) * (MSSg / MSg) 
	dMSIg <- -rg * MSIg +
		-gamma1*MSIg + 
		+gamma0*MSEg + 
		+2 * etag * MSg * MIg + 
		-etag * MSIg + 
		+ (-dSg) * (delta_si_g/gp(1)) * ( - MSIg/MSg) +
		+ (-dSf) * (thetag*gp(thetag)/g(thetag)/gp(1)) * ( - MSIg/MSg) +
		+ (-dSh) * (thetag*gp(thetag)/g(thetag)/gp(1)) * ( - MSIg/MSg) 
	
	
	
	
	
	dMSSf <- +1 * etaf * MSf^2 + #2 or 1? 
		- etaf * MSSf + 
		- (-dSf) * (delta_si_f/fp(1)) * MSSf / MSf +  #2 or 1 ? 
		- ((-dSg)+(-dSh)) * (thetaf*fp(thetaf)/f(thetaf)/fp(1)) * MSSf / MSf 
	dMSSg <- +1 * etag * MSg^2 + #2 or 1? 
		- etag * MSSg + 
		- (-dSg) * (delta_si_g/gp(1)) * MSSg / MSg  + #2 or 1 ? 
		- ((-dSf)+(-dSh)) * (thetag*gp(thetag)/g(thetag)/gp(1)) * MSSg / MSg 
	
	
	
	
	dMEf <- -gamma0 * MEf + 
		+ (-dSf) * (delta_si_f/fp(1)) +
		+ ((-dSg)+(-dSh)) * (thetaf*fp(thetaf)/f(thetaf)/fp(1))
	dMEg <- -gamma0 * MEg + 
		+ (-dSg) * (delta_si_g/gp(1)) +
		+ ((-dSf)+(-dSh)) * (thetag*gp(thetag)/g(thetag)/gp(1))
	dMEh <- -gamma0 * MEh + 
		+ (-dSh) * (delta_si_h/hp(1)) +
		+ ((-dSf)+(-dSg)) * (.thetah*hp(.thetah)/h(.thetah)/hp(1))
		
	dMIf <- -gamma1 * MIf + 
		+ gamma0 * MEf 
	dMIg <- -gamma1 * MIg + 
		+ gamma0 * MEg 
	dMIh <- -gamma1 * MIh + 
		+ gamma0 * MEh
	
	# infected, infectious and not detected: 
	newE <- transmf + transmg + transmh + transmseed 
	# exposed, infectious, diagnosed or undiagnosed
	E <- max(0, E + newE - gamma0*E)
	I <- max(0, I + gamma0*E - gamma1*I)
	
	# new case detections
	newI <- newI + gamma0*E # will accumulate between obvs
	# cumulative transm
	cutf <- cutf + transmf
	cutg <- cutg + transmg
	cuth <- cuth + transmh
	cuts <- cuts + transmseed 
	
	# time step one day , implicitly included 
	c( 
		thetaf = max(1e-9,min(1, thetaf + dthetaf)  )
		, MSEf = max(0, MSEf + dMSEf )
		, MEf = max(0, MEf + dMEf )
		, MSSf = max(0, MSSf + dMSSf)
		, MSIf = max(0, MSIf + dMSIf)
		, MIf  = max(0, MIf + dMIf)
		
		, thetag = max(1e-9, min(1, thetag + dthetag)   )
		, MSEg = max(0, MSEg + dMSEg )
		, MEg = max(0, MEg + dMEg )
		, MSSg = max(0, MSSg + dMSSg)
		, MSIg = max(0, MSIg + dMSIg)
		, MIg  = max(0, MIg + dMIg)
		
		, thetah = max(1e-9, min(1, thetah + dthetah)   )
		, MEh = max(0, MEh + dMEh )
		, MIh  = max(0, MIh + dMIh)
		
		, E = E
		, I = I 
		, newI = newI
		, cutf = cutf 
		, cutg = cutg 
		, cuth = cuth 
		, cuts = cuts 
		, seedrate = seedrate 
		, theta_vacc = theta_vacc 
	) 

}

#' @export
rmeas4.0 <- function(newI,MSIf,MSIg,MIh,seedrate,N,beta,delta, ...){
	rf <- beta * 1.5/7# factor 1.5/7 is act rate per day; factor 1.5 b/c more contacts per week in long partnerships 
	rg <- beta * 1/7
	mftransm <- MSIf*N*fp(1)*rf + MSIg*N*gp(1)*rg + beta*MIh*N*hp(1)
	Y <- rbinom( length( newI), size =ceiling(newI), prob = delta )
	Ytravel <- rbinom( length(Y), size = Y, prob = seedrate / (seedrate + mftransm)  )  
	c( Ytravel = Ytravel,  Yendog = Y - Ytravel, Yunk = 0  )
}

#' @export 
dmeas4.0 <- function(time, Ytravel, Yendog, Yunk, I, newI, MSIf,  MSIg, MIh, seedrate, N, beta, delta, gamma,  ..., log)
{
	Y <- Ytravel + Yendog + Yunk 
	if ( any( is.na( Y )) )
		return( ifelse(log, 0, 1) )
	t1 <- ifelse( log, -Inf, 0 )
	if ( newI >= Y )
		t1 <- suppressWarnings(  dbinom(Y , size = ceiling(newI) , prob = delta, log = log ) )
	if ( is.na( t1 ))
		t1 <- ifelse( log, -Inf, 0 )
	
	rf <- beta * 1.5/7# factor 1.5/7 is act rate per day; factor 1.5 b/c more contacts per week in long partnerships 
	rg <- beta * 1/7
	mftransm <- MSIf*N*fp(1)*rf + MSIg*N*gp(1)*rg + beta*MIh*N*hp(1)
		
	t2 <- ifelse(log,  0, 1)  
	if ( (Ytravel + Yendog) > 0 ) 
		t2 <- dbinom( Ytravel , size = Ytravel + Yendog, prob = seedrate / (seedrate + mftransm) , log = log )
	if ( is.na( t2 ))
		t2 <- ifelse( log, -Inf, 0 )
	
	rv = ifelse( log, t1 + t2, t1 * t2 ) 
	
	rv
}

#' @export 
rinit4.0 <- function(i0,N,seedrate0, ... ){
	xinit <- i0 / N 
	c(
		thetaf = 1 
			, MSEf = 0
			, MEf = 0
			, MSSf = 1
			, MSIf = 0
			, MIf = 0
		, thetag = 1 
			, MSEg = 0 
			, MEg = 0 
			, MSSg = 1
			, MSIg = 0
			, MIg = 0
		, thetah = 1 - xinit 
			, MEh = xinit / 2
			, MIh = xinit / 2 #* (1+hpp(1)/hp(1)) / hp(1)
		, E = xinit*N / 2
		, I = xinit*N / 2
		, newI = 0 
		, cutf = 0 
		, cutg = 0 
		, cuth = 0 
		, cuts = 0 
		, seedrate = seedrate0
		, theta_vacc = 1 
	)
}

#' @export 
m4.0 <- pomp::pomp(
	t0 = 1
	, data = NULL
	, times = seq( 1, 90 ) 
	, rprocess = pomp::discrete_time(step.fun = rstep4.0)
	, rmeasure = rmeas4.0
	, dmeasure = dmeas4.0
	, rinit = rinit4.0
	, statenames = c('thetaf'
			, 'MSEf', 'MEf'
			, 'MSSf', 'MSIf', 'MIf'
		, 'thetag'
			, 'MSEg', 'MEg' 
			, 'MSSg', 'MSIg', 'MIg'
		, 'thetah'
			, 'MEh', 'MIh'
		, 'E', 'I', 'newI', 'cutf', 'cutg', 'cuth', 'cuts'
		, 'seedrate'
		, 'theta_vacc'
		)
	, obsnames = c('Ytravel', 'Yendog', 'Yunk')
	, paramnames = c( 'beta', 'gamma0', 'gamma1', 'etaf', 'etag',  'N', 'i0', 'delta', 'seedrate0', 'seedrate_sd' 
		,  'vacc_freq', 'vacc_amt', 'vacc_start_day', 'vacc_fin_day' )
	, accumvars = c('newI')
	, params = c( beta = 2.25
		, gamma0 = 1/8
		, gamma1 = 1/4
		, etaf = 1/200 ## Anderson Epidemiology 2021
		, etag = 1/100 #Anderson Epidemiology 2021
		, N =  .49 * 56.3 * (1-.24 ) * 1e6 * 0.02 # male*england(millions)*adult*msm
		, i0 = 0
		, delta = .4 #40pc detection rate
		, seedrate0 = 0.75
		, seedrate_sd = 0.05 
		, vacc_freq = 7
		, vacc_amt = 0.04
		, vacc_start_day = 91
		, vacc_fin_day = 91+5*7 ##
		)
		
)
