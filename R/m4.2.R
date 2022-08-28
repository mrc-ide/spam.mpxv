#~ 4.2
#~ - duplicating state variables for seeds so that they do not go directly to I compartment 
#~ - time variable case ascertainment 

#~ 4.1 
#~ - beta(t) 
#~ - seedrate AR(1) 

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



.update_theta_vacc4.2 <- function(  theta_vacc, vacc_amt  )
{ # note allows previous infections to receive vacc 
	p0 <- f(theta_vacc) * g(theta_vacc) * h( theta_vacc ) 
	p1 <- max(0.01, p0 - vacc_amt )
	of <- function( lntheta ){
		(f(exp(lntheta))*g(exp(lntheta))*h(exp(lntheta)) - p1)^2
	}
	o = optimise( lower = -1e4, upper = 0 , f = of )
	exp( o$minimum  )
}



poispgf <- function(x, z ) exp( z*(x-1))




#' @export
rstep4.2 <- function( thetaf
		, MSEf, MEf
		, MSSf, MSIf, MIf
	, thetag
		, MSEg, MEg
		, MSSg, MSIg, MIg
	, thetah
		, MEh
		, MIh
	, E, I, newI
	, Eseed, newIseed 
	, cutf, cutg, cuth, cuts 
	, seedrate
	, dseedrate 
	, theta_vacc # targetted vacc surv function
	, S_vacc # random vacc surv function
	, beta # now a state variable 
	
	, beta0, beta_freq, beta_sd,  gamma0, gamma1, etaf, etag,  N, i0, delta0, delta1, delta_slope, seedrate0, seedrate_sd
	, vacc_freq, vacc_amt, vacc_start_day, vacc_fin_day, vacc_targetted ##
	
	, time
	, ...) 
{
	N <- max( N, 1e4 ) # impose minimum pop size 
	# vacc model 
	## every vacc_freq days inoculate proportion vacc_amt with prob prop to degree
	if (  (time>=vacc_start_day)  &  (time<=vacc_fin_day)  &  (( (time-vacc_start_day) %% vacc_freq) == 0) ){
		amt_targetted <- vacc_targetted * vacc_amt / N  
		amt_random <- (1-vacc_targetted) * vacc_amt / N
		
		S_vacc <- S_vacc * (1-amt_random)
		MSEf <- MSEf * (1-amt_random) 
		MSSf <- MSSf * (1-amt_random)^2
		MSIf <- MSIf * (1-amt_random)
		MSEg <- MSEg * (1-amt_random)
		MSSg <- MSSg * (1-amt_random)^2
		MSIg <- MSIg * (1-amt_random)
		
		theta_vacc0 <- theta_vacc 
		theta_vacc  <- unname( .update_theta_vacc4.2( theta_vacc, amt_targetted  ) )
		red_f <- (theta_vacc*fp(theta_vacc )) / (theta_vacc0*fp(theta_vacc0) )
		red_g <- (theta_vacc*gp(theta_vacc )) / (theta_vacc0*gp(theta_vacc0) )
		MSEf <- MSEf * red_f ##? TODO should be greater for targetted vacc on g and f 
		MSSf <- MSSf * red_f^2
		MSIf <- MSIf * red_f
		MSEg <- MSEg * red_g
		MSSg <- MSSg * red_g^2
		MSIg <- MSIg * red_g
	}
	.thetaf <- thetaf * theta_vacc 
	.thetag <- thetag * theta_vacc 
	.thetah <- thetah * theta_vacc 
	
	if ((time %%  beta_freq) == 0 )  {
		dseedrate <- rnorm (1, dseedrate, sd = seedrate_sd  ) 
		#seedrate <- max(0,  rnorm (1, seedrate, sd = seedrate_sd*beta_freq  ) )
	}
	seedrate <- max(0, seedrate + dseedrate  )
	
	if ((time %%  beta_freq) == 0 )  
		beta <- max(0, rnorm(1,  beta, sd = beta_sd ) )
	
	MSf  <- .thetaf * S_vacc * fp( .thetaf ) / fp(1) 
	MSg  <- .thetag * S_vacc * gp( .thetag ) / gp(1) 
	MSh  <- .thetah * S_vacc * hp( .thetah ) / hp(1) 
	
	# beta is transm probability per contact * contact rate inflation factor. Transm rates:
	#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5142082/
	rf <- beta * 1.5/7# factor 1.5/7 is act rate per day; factor 1.5 b/c more contacts per week in long partnerships 
	rg <- beta * 1/7
	
	tratef <- max(0, MSIf*N*fp(1)*rf )
	if ( is.na( tratef ) | is.infinite( tratef ) ){
		print( 'NA tratef' )
		for (x in ls()){
			print(x)
			print( get(x) )
		}
		print( c( MSIf , N , fp(1), rf ))
		
		tratef <- 0 
	}
	transmf <- rpois(1, tratef )#
	dthetaf <- -.thetaf * transmf / (MSf*N*fp(1)) 
	dSf <- fp(.thetaf)*dthetaf #note prop to transm
	meanfield_delta_si_f = u1f <- (.thetaf * fpp(.thetaf) / fp(.thetaf) )
	u2f <- (.thetaf * fpp(.thetaf) + .thetaf^2 * fppp(.thetaf ) ) / fp(.thetaf)
	vf = u2f - u1f^2
	if ( transmf > 0 )
		delta_si_f <- mean( rnorm( transmf, meanfield_delta_si_f , sd = sqrt(vf)) )
	else
		delta_si_f <- 0 
	
	trateg <- max(0, MSIg*N*gp(1)*rg )
	if ( is.na( trateg ) | is.infinite(trateg)  ) trateg <- 0 
	transmg <- rpois(1, trateg)
	dthetag <- -.thetag * transmg / (MSg*N*gp(1)) 
	dSg <- gp(.thetag)*dthetag #note prop to transm
	meanfield_delta_si_g = u1g <- (.thetag * gpp(.thetag) / gp(.thetag) )
	u2g <- (.thetag * gpp(.thetag) + .thetag^2 * gppp(.thetag ) ) / gp(.thetag)
	vg = u2g - u1g^2
	if ( transmg > 0 )
		delta_si_g <- mean( rnorm( transmg, meanfield_delta_si_g , sd = sqrt(vg)) )
	else
		delta_si_g <- 0 
	
	
	transmseed <- rpois( 1, seedrate ) 
	
	trateh <- max(0,  beta*MIh*N*hp(1)*S_vacc )
	if ( is.na( trateh ) | is.infinite( trateh )) trateh <- 0 
	transmh <- rpois(1, trateh ) # note S_vacc here, because there is no MSI in MFSH model 
	dthetah <- -.thetah * (transmh + transmseed ) / (N*hp(1) ) # seeding happens here 
	dSh <- hp(.thetah)*dthetah #note prop to transm
	# may also need this separated into seed and non-seed components: 
		#dSh0 <- hp(.thetah)*(-.thetah * transmh   / (N*hp(1) )) #
		#dSh_seed <- hp(.thetah)*(-.thetah * transmseed  / (N*hp(1) )) #
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
		+ (-dSg) * (.thetaf*fp(.thetaf)/f(.thetaf)/fp(1)) * (MSSf / MSf) +
		+ (-dSh) * (.thetaf*fp(.thetaf)/f(.thetaf)/fp(1)) * (MSSf / MSf) 	
	dMSIf <- -rf * MSIf +
		-gamma1*MSIf + 
		+gamma0*MSEf + 
		+2 * etaf * MSf * MIf + 
		-etaf * MSIf + 
		+ (-dSf) * (delta_si_f/fp(1)) * ( -MSIf/MSf) +
		+ (-dSg) * (.thetaf*fp(.thetaf)/f(.thetaf)/fp(1)) * ( - MSIf/MSf) +
		+ (-dSh) * (.thetaf*fp(.thetaf)/f(.thetaf)/fp(1)) * ( - MSIf/MSf) 	
	
	
	dMSEg <- -gamma0*MSEg + 
		+2 * etag * MSg * MEg + 
		-etag * MSEg + 
		+ (-dSg) * (delta_si_g/gp(1)) * (MSSg / MSg ) +
		+ (-dSf) * (.thetag*gp(.thetag)/g(.thetag)/gp(1)) * (MSSg / MSg) +
		+ (-dSh) * (.thetag*gp(.thetag)/g(.thetag)/gp(1)) * (MSSg / MSg) 
	dMSIg <- -rg * MSIg +
		-gamma1*MSIg + 
		+gamma0*MSEg + 
		+2 * etag * MSg * MIg + 
		-etag * MSIg + 
		+ (-dSg) * (delta_si_g/gp(1)) * ( - MSIg/MSg) +
		+ (-dSf) * (.thetag*gp(.thetag)/g(.thetag)/gp(1)) * ( - MSIg/MSg) +
		+ (-dSh) * (.thetag*gp(.thetag)/g(.thetag)/gp(1)) * ( - MSIg/MSg) 
	
	
	
	
	
	dMSSf <- +1 * etaf * MSf^2 + #2 or 1? 
		- etaf * MSSf + 
		- (-dSf) * (delta_si_f/fp(1)) * MSSf / MSf +  #2 or 1 ? 
		- ((-dSg)+(-dSh)) * (.thetaf*fp(.thetaf)/f(.thetaf)/fp(1)) * MSSf / MSf 
	dMSSg <- +1 * etag * MSg^2 + #2 or 1? 
		- etag * MSSg + 
		- (-dSg) * (delta_si_g/gp(1)) * MSSg / MSg  + #2 or 1 ? 
		- ((-dSf)+(-dSh)) * (.thetag*gp(.thetag)/g(.thetag)/gp(1)) * MSSg / MSg 
	
	
	
	
	dMEf <- -gamma0 * MEf + 
		+ (-dSf) * (delta_si_f/fp(1)) +
		+ ((-dSg)+(-dSh)) * (.thetaf*fp(.thetaf)/f(.thetaf)/fp(1))
	dMEg <- -gamma0 * MEg + 
		+ (-dSg) * (delta_si_g/gp(1)) +
		+ ((-dSf)+(-dSh)) * (.thetag*gp(.thetag)/g(.thetag)/gp(1))
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
	I <- max(0, I + gamma0*E - gamma1*I)
	# new case detections. some subset of these will be detected each week
	newI <- newI + gamma0*E # will accumulate between obvs
	E <- max(0, E + newE - gamma0*E)
	
	# seed state variables 
	# exposed, infectious, diagnosed or undiagnosed
	# new case detections. some subset of these will be detected each week
	newIseed <- newIseed + gamma0*Eseed # will accumulate between obvs
	Eseed <- max(0, Eseed + transmseed - gamma0*Eseed)
	
	
	# cumulative transm
	cutf <- cutf + transmf
	cutg <- cutg + transmg
	cuth <- cuth + transmh
	cuts <- cuts + transmseed 
	
	# protect outputs 
	if ( is.na( dthetaf )) dthetaf <- 0 
	if ( is.na( dthetag )) dthetag <- 0 
	if ( is.na( dthetah )) dthetah <- 0 
	if ( is.na( dMSEf )) dMSEf <- 0 
	if ( is.na( dMSEg )) dMSEg <- 0 
	if ( is.na( dMSSf )) dMSSf <- 0 
	if ( is.na( dMSSg )) dMSSg <- 0 
	if ( is.na( dMSIf )) dMSIf <- 0 
	if ( is.na( dMSIg )) dMSIg <- 0 
	if ( is.na( dMEf )) dMEf <- 0 
	if ( is.na( dMEg )) dMEg <- 0 
	if ( is.na( dMEh )) dMEh <- 0 
	if ( is.na( dMIf )) dMIf <- 0 
	if ( is.na( dMIg )) dMIg <- 0 
	if ( is.na( dMIh )) dMIh <- 0 
	
	# time step one day , implicitly included 
	rv <- c( 
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
		
		, Eseed = Eseed 
		, newIseed = newIseed 
		
		, cutf = cutf 
		, cutg = cutg 
		, cuth = cuth 
		, cuts = cuts 
		, seedrate = seedrate 
		, dseedrate = dseedrate 
		, theta_vacc = theta_vacc 
		, S_vacc = S_vacc
		, beta = beta 
	) 
#~ browser()
	rv[ is.na(rv) ] <- 0
	rv
}

#' @export
rmeas4.2 <- function(newI,MSIf,MSIg,MIh,seedrate,N,beta, delta0, delta1, delta_slope, ...){
	delta <- max(.01, min( delta1, delta0 + delta_slope * time ) )
	Yendog <- rbinom( length( newI), size =ceiling(newI), prob = delta )
	Ytravel <- rbinom( length( newIseed ), size = ceiling( newIseed ), prob = delta ) 
	c( Ytravel = Ytravel,  Yendog = Yendog, Yunk = 0  )
}

#' time-varying case ascertainment
#' @export 
dmeas4.2 <- function(time, Ytravel, Yendog, Yunk, I, newI, newIseed, MSIf,  MSIg, MIh, seedrate, N, beta, delta0, delta1, delta_slope, gamma,  ..., log)
{
	Y <- ceiling( Ytravel + Yendog + Yunk  )
	delta <- max(.01, min( delta1, delta0 + delta_slope * time ) )
	if ( any( is.na( Y )) )
		return( ifelse(log, 0, 1) )
	t1 <- ifelse( log, -Inf, 0 )
	if ( newI >= Y )
		t1 <- suppressWarnings( dbinom(Y , size = ceiling(newI) , prob = delta, log = log ) )
	else
		return( ifelse(log, -Inf, 0) )
	if ( is.na( t1 ))
		t1 <- ifelse( log, -Inf, 0 )
	
	t2 <- ifelse(log, 0, 1)  
	if ( (Ytravel + Yendog) > 0 ) 
		t2 <- dbinom( ceiling(Ytravel) , size = ceiling(Ytravel + Yendog), prob = newIseed / (newIseed + newI) , log = log )
	if ( is.na( t2 ))
		t2 <- ifelse( log, -Inf, 0 )
	
	rv = ifelse( log, t1 + t2, t1 * t2 ) 
#~ browser()
	rv
}

#' @export 
rinit4.2 <- function(i0,N,seedrate0,beta0, ... ){
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
		, Eseed = 0
		, newIseed  = 0 
		, cutf = 0 
		, cutg = 0 
		, cuth = 0 
		, cuts = 0 
		, seedrate = seedrate0
		, dseedrate = 0 
		, theta_vacc = 1 
		, S_vacc = 1
		, beta = beta0
	)
}

#' @export 
m4.2 <- pomp::pomp(
	t0 = 1
	, data = NULL
	, times = seq( 1, 90 ) 
	, rprocess = pomp::discrete_time(step.fun = rstep4.2)
	, rmeasure = rmeas4.2
	, dmeasure = dmeas4.2
	, rinit = rinit4.2
	, statenames = c('thetaf'
			, 'MSEf', 'MEf'
			, 'MSSf', 'MSIf', 'MIf'
		, 'thetag'
			, 'MSEg', 'MEg' 
			, 'MSSg', 'MSIg', 'MIg'
		, 'thetah'
			, 'MEh', 'MIh'
		, 'E', 'I', 'newI'
		, 'Eseed', 'newIseed' 
		, 'cutf', 'cutg', 'cuth', 'cuts'
		, 'seedrate'
		, 'dseedrate'
		, 'theta_vacc'
		, 'S_vacc'
		, 'beta'
		)
	, obsnames = c('Ytravel', 'Yendog', 'Yunk')	
	, paramnames = c( 'beta0', 'beta_freq', 'beta_sd'
		, 'gamma0', 'gamma1', 'etaf', 'etag',  'N', 'i0'
		, 'delta0', 'delta1', 'delta_slope' 
		, 'seedrate0', 'seedrate_sd' 
		,  'vacc_freq', 'vacc_amt', 'vacc_start_day', 'vacc_fin_day' , 'vacc_targetted')
	, accumvars = c('newI')
	, params = c( beta0 = 2.25
		, beta_freq = 7
		, beta_sd = 0.15  # 1 = 2*sqrt( (75/7)*sigma2 ) 
		, gamma0 = 1/8
		, gamma1 = 1/4
		, etaf = 1/200 ## Anderson Epidemiology 2021
		, etag = 1/100 #Anderson Epidemiology 2021
		, N =  750e3 # > 1e4 
		, i0 = 0
		, delta0 = .20
		, delta1 = .80 
		, delta_slope = 0.0
		, seedrate0 = 0.75
		, seedrate_sd = 3 #sd of random walk of daily diff in seedrate 
		, vacc_freq = 1
		, vacc_amt = 0.85*1040 # assuming about 30k doses over 30 days, 85% vacc eff 
		, vacc_start_day = 91
		, vacc_fin_day = 91+29 ##
		, vacc_targetted = .8 # prop vacc targetted vs random 
		)
)
