#~ m2: 
#~ -seeding rate gaussian process 
#~ -accumulate variables for fitting 

#~  m1.1: 
#~ Continuous seeding 
#~ Measurement model for seeded cases (travel associated) and endogenous


#' @export
rstep2.0 <- function( thetaf, MSSf, MSIf, MIf
	, thetag, MSSg, MSIg, MIg
	, thetah, MIh
	, I, C, newC, cutf, cutg, cuth, cuts 
	, seedrate
	
	, beta,  gamma, etaf, etag,  N, i0, delta, seedrate0, seedrate_sd
	
	, time
	, ...) 
{
	seedrate <- max(0, seedrate + rnorm (1, 0, sd = seedrate_sd  ) )
	
	MSf  <- thetaf * fp( thetaf ) / fp(1) 
	MSg  <- thetag * gp( thetag ) / gp(1) 
	MSh  <- thetah * hp( thetah ) / hp(1) 
	
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
	dthetah <- -thetah * (transmh + transmseed ) / (N*hp(1) ) # seeding happens here 
	dSh <- hp(thetah)*dthetah #note prop to transm
	meanfield_delta_si_h = u1h <- (1 + thetah * hpp(thetah) / hp(thetah) ) #note + 1 for mfsh
	u2h <- hppp(thetah)*thetah^2/hp(thetah) + 2*thetah*hpp(thetah)/hp(thetah) + u1h 
	vh = u2h - u1h^2
	if ( transmh > 0 )
		delta_si_h <- mean( rnorm( transmh, meanfield_delta_si_h , sd = sqrt(vh)) )
	else
		delta_si_h <- 0 
	
	
	
	dMSIf <- -rf * MSIf +
		-gamma*MSIf + 
		+2 * etaf * MSf * MIf + 
		-etaf * MSIf + 
		+ (-dSf) * (delta_si_f/fp(1)) * (MSSf / MSf - MSIf/MSf) +
		+ (-dSg) * (thetaf*fp(thetaf)/f(thetaf)/fp(1)) * (MSSf / MSf - MSIf/MSf) +
		+ (-dSh) * (thetaf*fp(thetaf)/f(thetaf)/fp(1)) * (MSSf / MSf - MSIf/MSf) 	
	dMSIg <- -rg * MSIg +
		-gamma*MSIg + 
		+2 * etag * MSg * MIg + 
		-etag * MSIg + 
		+ (-dSg) * (delta_si_g/gp(1)) * (MSSg / MSg - MSIg/MSg) +
		+ (-dSf) * (thetag*gp(thetag)/g(thetag)/gp(1)) * (MSSg / MSg - MSIg/MSg) +
		+ (-dSh) * (thetag*gp(thetag)/g(thetag)/gp(1)) * (MSSg / MSg - MSIg/MSg) 
	
	
	
	
	
	dMSSf <- +1 * etaf * MSf^2 + #2 or 1? 
		- etaf * MSSf + 
		- (-dSf) * (delta_si_f/fp(1)) * MSSf / MSf +  #2 or 1 ? 
		- ((-dSg)+(-dSh)) * (thetaf*fp(thetaf)/f(thetaf)/fp(1)) * MSSf / MSf 
	dMSSg <- +1 * etag * MSg^2 + #2 or 1? 
		- etag * MSSg + 
		- (-dSg) * (delta_si_g/gp(1)) * MSSg / MSg  + #2 or 1 ? 
		- ((-dSf)+(-dSh)) * (thetag*gp(thetag)/g(thetag)/gp(1)) * MSSg / MSg 
	
	
	
	
	dMIf <- -gamma * MIf + 
		+ (-dSf) * (delta_si_f/fp(1)) +
		+ ((-dSg)+(-dSh)) * (thetaf*fp(thetaf)/f(thetaf)/fp(1))
	dMIg <- -gamma * MIg + 
		+ (-dSg) * (delta_si_g/gp(1)) +
		+ ((-dSf)+(-dSh)) * (thetag*gp(thetag)/g(thetag)/gp(1))
	dMIh <- -gamma * MIh + 
		+ (-dSh) * (delta_si_h/hp(1)) +
		+ ((-dSf)+(-dSg)) * (thetah*hp(thetah)/h(thetah)/hp(1))
	
	# infected, infectious and not detected: 
	newI <- transmf + transmg + transmh + transmseed 
	I <- max(0, I + newI - gamma*I - delta*I)
	# infected, infectious & detected case
	C <- max(0, C + delta*I - gamma*C ) 
	# new case detections
	newC <- newC + delta*I # will accumulate between obvs
	# cumulative transm
	cutf <- cutf + transmf
	cutg <- cutg + transmg
	cuth <- cuth + transmh
	cuts <- cuts + transmseed 
	
	
	
	# time step one day , implicitly included 
	c( 
		thetaf = max(1e-9,min(1, thetaf + dthetaf)  )
		, MSSf = max(0, MSSf + dMSSf)
		, MSIf = max(0, MSIf + dMSIf)
		, MIf  = max(0, MIf + dMIf)
		
		, thetag = max(1e-9, min(1, thetag + dthetag)   )
		, MSSg = max(0, MSSg + dMSSg)
		, MSIg = max(0, MSIg + dMSIg)
		, MIg  = max(0, MIg + dMIg)
		
		, thetah = max(1e-9, min(1, thetah + dthetah)   )
		, MIh  = max(0, MIh + dMIh)
		
		, I = I 
		, C = C
		, newC = newC 
		, cutf = cutf 
		, cutg = cutg 
		, cuth = cuth 
		, cuts = cuts 
		, seedrate = seedrate 
	) 

}

#' @export
rmeas2.0 <- function(newC,MSIf,MSIg,MIh,seedrate,N,beta, ...){
	rf <- beta * 1.5/7# factor 1.5/7 is act rate per day; factor 1.5 b/c more contacts per week in long partnerships 
	rg <- beta * 1/7
	mftransm <- MSIf*N*fp(1)*rf + MSIg*N*gp(1)*rg + beta*MIh*N*hp(1)
	Y <- rpois(length(newC), newC )
	Ytravel <- rbinom( length(Y), size = Y, prob = seedrate / (seedrate + mftransm)  )  
	
	c( Ytravel = Ytravel,  Yendog = Y - Ytravel, Yunk = 0  )
}

#' @export 
dmeas2.0 <- function(Ytravel, Yendog, Yunk, newC, MSIf,  MSIg, MIh, seedrate, N, beta,  ..., log)
{
	rf <- beta * 1.5/7# factor 1.5/7 is act rate per day; factor 1.5 b/c more contacts per week in long partnerships 
	rg <- beta * 1/7
	mftransm <- MSIf*N*fp(1)*rf + MSIg*N*gp(1)*rg + beta*MIh*N*hp(1)
	Y <- Ytravel + Yendog + Yunk
	if ( newC == 0 & Y > 0 )
		t1 <- ifelse( log, -Inf, 0 )
	else
		t1 <- dpois( Y, lambda = newC, log = log )
	
	t2 <- ifelse(log,  0, 1)  
	if ( (Ytravel + Yendog) > 0 ) 
		t2 <- dbinom( Ytravel , size = Ytravel + Yendog, prob = seedrate / (seedrate + mftransm) , log = log )
	if ( is.na( t2 ))
		t2 <- ifelse( log, -Inf, 0 )
	
	ifelse( log, t1 + t2, t1 * t2 ) 
}

#' @export 
rinit2.0 <- function(i0,N,seedrate0, ... ){
	xinit <- i0 / N 
	c(
		thetaf = 1 
			, MSSf = 1
			, MSIf = 0
			, MIf = 0
		, thetag = 1 
			, MSSg = 1
			, MSIg = 0
			, MIg = 0
		, thetah = 1 - xinit 
			, MIh = xinit #* (1+hpp(1)/hp(1)) / hp(1)
		, I = xinit*N
		, C = 0
		, newC = 0 
		, cutf = 0 
		, cutg = 0 
		, cuth = 0 
		, cuts = 0 
		, seedrate = seedrate0
		
	)
}

#' @export 
m2.0 <- pomp::pomp(
	t0 = 1
	, data = NULL
	, times = seq( 1, 365 ) 
	, rprocess = pomp::discrete_time(step.fun = rstep2.0)
	, rmeasure = rmeas2.0
	, dmeasure = dmeas2.0
	, rinit = rinit2.0
	, statenames = c('thetaf', 'MSSf', 'MSIf', 'MIf'
		, 'thetag', 'MSSg', 'MSIg', 'MIg'
		, 'thetah', 'MIh'
		, 'I', 'C', 'newC', 'cutf', 'cutg', 'cuth', 'cuts'
		, 'seedrate'
		)
	, obsnames = c('Y')
	, paramnames = c( 'beta', 'gamma', 'etaf', 'etag',  'N', 'i0', 'delta', 'seedrate0', 'seedrate_sd' )
	, accumvars = c('newC')
	, params = c( beta = .75
		, gamma = 0.1
		, etaf = 1/200 ## Anderson Epidemiology 2021
		, etag = 1/100 #Anderson Epidemiology 2021
		, N =  .49 * 56.3 * (1-.24 ) * 1e6 * 0.02 # male*england(millions)*adult*msm
		, i0 = 0
		, delta = .4*.1/(1-.4)
		, seedrate0 = 0.75
		, seedrate_sd = 0.05 
		)
		
)
