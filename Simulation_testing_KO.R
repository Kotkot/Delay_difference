#########################
#
# NOTES
#
# 1.  Fixing F_equil=F_t[1] for spatial and nonspatial models increases bias for F, S, N, etc.
#
# TEMP
# A.  par()$usr gives plot region boundaries for plotting a legend
#
#########################

#########################
# Simulation or data analysis model
#########################

File = "C:/Users/Kotaro Ono/Dropbox/Postdoc_projects/side_projects/Delay_difference/"
File = "C:/Users/Kotkot/Dropbox/Postdoc_projects/side_projects/Delay_difference/"
TmbFile = paste0(File, "/executables/")

Run_simul <- function(SD_A=0.5, SD_E=0.5, Scale=250, CV_w=0.2, Accel=0.2, SD_F=0.2, ...)
{

# Directory name
  Dir_save=paste0(File,"SD_A=", SD_A, "_SD_E=", SD_E, "_Scale=", Scale, "_CV_w=", CV_w, "_Accel=", Accel, "_SD_F=", SD_F)
  dir.create(Dir_save)

# Libraries
library(INLA)
library(TMB)
newtonOption(smartsearch=TRUE)
library(RandomFields)
library(PBSmapping)
library(matlab)
library(pscl)  # zeroinfl
library(RANN) # nn2 -- nearest neighbors assignment
library(cluster)	#pam 
    
# Source helper functions
source( file=paste(File,"Fn_helpers_2014-05-09.R",sep="") )

#### Settings
# Simulation test
  #RandomSeed = as.numeric(paste(na.omit(as.numeric(strsplit(as.character(Date),"")[[1]])),collapse=""))
  RandomSeed = 0#ceiling( runif(1, min=1, max=1e6) )
  if( !exists("ThreadNum") ) ThreadNum = 1
  RepSet = 1:50
  RepSet = RepSet + max(RepSet)*(ThreadNum-1) 
# Estimation
  ModelSet = c("Nonspatial", "Spatial", "Index", "Strata")[c(1:2,4)]
  ErrorModel_CatchRates = 1   # 0: Poisson; 1: Negative binomial for counts
  ErrorModel_MeanWeight = 1   # 0: Fixed-CV; 1: Est-CV
  Smooth_F = 0      # 0: No; 1: Yes
  Fix_Q = TRUE   # 
  SpatialSimModel = "Matern"
  MeshType = c("Minimal","Recommended")[1]
  Version = c("delay_difference_v4d", "delay_difference_v6c", "delay_difference_v8e", "delay_difference_v9")[4]
  n_s = 25
# Domain
  n_t = 30
  Range_X = c(0,1000)
  Range_Y = c(0,1000)
  DomainArea = diff(Range_X) * diff(Range_Y)
# Biological - growth
  k = 3         # Age at when individual is fully recruited
  ro = 0.8      # Brody growth coefficient
  alpha_g = 0.2   # Weight at age 1
# Biological - survival
  M = 0.3  # Natural mortality 
# Biological recruitment
  RecFn = c("Ricker", "BH", "Constant")[3]
  mu_R0_total = 1e9 # Median total R0
  mu_N0 = mu_R0_total / (1-exp(-M))
# Variability
  SD_A = SD_A  # 0.5  # Spatial variation in productivity
  SD_E = SD_E  # 0.5  # Spatiotemporal variation in recruitment
  Scale = Scale		  # Controls the range of the semivariogram
# Fishing mortality 
  F_equil = 0.05  # Initial equilibrium fishing mortality
  S_bioecon = 0.4
  Accel = Accel
  SD_F = SD_F # 0.2
# Data
  n_samp_per_year = 100 
  AreaSwept = 0.0025 # 10 / mu_N0 * DomainArea  # in km^2
  q_I = 1
  CV_w = CV_w		# sampling variability in the weight data
  CV_c = 0.05		# observation error with the catch data
# Visualization
  Ngrid_sim = 1e4
  Ngrid_proj = 1e4  
  
  # Save settings
  SettingsList = list( "RandomSeed"=RandomSeed, "ThreadNum"=ThreadNum, "RepSet"=RepSet, "ModelSet"=ModelSet, "ErrorModel_CatchRates"=ErrorModel_CatchRates, "ErrorModel_MeanWeight"=ErrorModel_MeanWeight, "Smooth_F"=Smooth_F, "Fix_Q"=Fix_Q, "SpatialSimModel"=SpatialSimModel, "MeshType"=MeshType, "Version"=Version, "n_s"=n_s, "n_t"=n_t, "Range_X"=Range_X, "Range_Y"=Range_Y, "k"=k, "ro"=ro, "alpha_g"=alpha_g, "M"=M, "RecFn"=RecFn, "mu_R0_total"=mu_R0_total, "SD_A"=SD_A, "SD_E"=SD_E, "Scale"=Scale, "F_equil"=F_equil, "S_bioecon"=S_bioecon, "Accel"=Accel, "SD_F"=SD_F, "n_samp_per_year"=n_samp_per_year, "AreaSwept"=AreaSwept, "q_I"=q_I, "CV_w"=CV_w, "CV_c"=CV_c)  
    capture.output(SettingsList, file=paste(Dir_save,"SettingsList.txt",sep=""))
    save(SettingsList, file=paste(Dir_save,"SettingsList.RData",sep=""))
    # file.copy( from=paste(TmbFile,Version,".cpp",sep=""), to=paste(Dir_save,Version,".cpp",sep=""), overwrite=TRUE)
  
  # Replication loop
  RepI=1; ModelI=3
  for(RepI in 1:length(RepSet)){

    # RepFile
    RepFile = paste(Dir_save,"Rep=",RepSet[RepI],"/",sep="")
    dir.create(RepFile)
    set.seed( RandomSeed + RepSet[RepI] )
  
	### Set-up the operating model
	# Derived values
      w_a = rep(NA,length=20)
      w_a[1] = alpha_g
      for(a in 2:length(w_a)) w_a[a] = alpha_g + ro*w_a[a-1]
      w_k = w_a[k]
      MRPSB = (1 - exp(-M) - ro*exp(-M) + ro*exp(-2*M) ) / ( w_k - (w_k-alpha_g)*exp(-M) )  # Maximum recruits per spawning biomass
      mu_R0_per_area = mu_R0_total / DomainArea
      mu_S0_per_area = mu_R0_per_area / MRPSB
    if(RecFn=="Ricker"){
      log_mu_alpha = log(mu_R0_total)   # median of maximum recruits per spawning biomass per unit area; min possible = MRPSB
      beta = log( mu_R0/exp(log_mu_alpha)/mu_S0_per_area ) / (-mu_S0_per_area)
    }
    if(RecFn=="BH"){
      log_mu_alpha = log(mu_R0_total)   # median of maximum recruits per spawning biomass per unit area; min possible = MRPSB
      beta = (exp(log_mu_alpha)*mu_S0_per_area/mu_R0_per_area - 1) / mu_S0_per_area
      tmp = seq(0,mu_S0_per_area,length=1000)
      tmp2 = exp(log_mu_alpha)*tmp / (1 + beta*tmp)
      plot( x=tmp, y=tmp2) 
    }
    if(RecFn=="Constant"){
      log_mu_alpha = log(mu_R0_total / DomainArea) # Expected recruitment per unit area
    }

    # Simulate random variability
    Eps_F = rnorm(n_t, mean=0, sd=SD_F)
    Eps_C = rnorm(n_t, mean=0, sd=CV_c)  
    Year_Range = c(1,n_t)
    F_t = rep(NA, n_t)

    # Define fields
      Nu=1
      if(SpatialSimModel=="Matern"){
        # matern: C(r=h/Scale) = Var * 2^{1-Nu} * gamma(Nu)^{-1} * (sqrt{2*Nu}*r)^Nu * besselK(sqrt{2*Nu} * r)
        model_A <- RMmatern(nu=Nu, var=SD_A^2, scale=Scale)
        model_E <- RMmatern(nu=Nu, var=SD_E^2, scale=Scale)
      }
      # Calculate range at which correlation is 13% of maximum (comparable to range defined in Lindgren and Rue 2013, immediately before Eq. 4)
      h = seq(-3*diff(Range_X),3*diff(Range_X),length=10000)
      r = abs(h)/Scale
      if(SpatialSimModel=="Matern") C = SD_A^2 * 2^(1-Nu) * gamma(Nu)^(-1) * (sqrt(2*Nu)*r)^Nu * besselK(sqrt(2*Nu) * r, nu=Nu)
      if(SpatialSimModel=="Gaussian") C = SD_A^2 * exp(-r^2)
      par(mfrow=c(1,2))
      plot(model_A, ylim=c(0,1), xlim=c(-3,3)*diff(Range_X))
      plot(x=h, y=C, ylim=c(0,1), type="l")
      RangeTrue = abs(h[which.min( (C-0.10*SD_A)^2 )])
  
    # locations of sampling stations
    if( n_samp_per_year < n_s*2 ) error("Increase n_samp_per_year for K-means algorithm")
    if(TRUE){
      x_stations = runif(n_samp_per_year, min=Range_X[1], max=Range_X[2])
      y_stations = runif(n_samp_per_year, min=Range_Y[1], max=Range_Y[2])
      loc_samples = cbind("long"=x_stations, "lat"=y_stations)
    }
    if(FALSE){
      x_stations = seq(Range_X[1],Range_X[2],length=sqrt(n_samp_per_year))
      y_stations = seq(Range_Y[1],Range_Y[2],length=sqrt(n_samp_per_year))
      loc_samples = expand.grid("long"=x_stations, "lat"=y_stations)
    }
   
    # K-means
      K = kmeans( x=loc_samples[,c('long','lat')], centers=n_s, iter.max=1e3, nstart=100) # K$tot.withinss
      K$centers = K$centers + rnorm( n_s, mean=0, sd=rep(0.001*c(diff(Range_X),diff(Range_Y)),each=n_s) )  # Necessary to ensure Kmeans aren't identical to other samples
      # K = pam( x=loc_samples[,c('long','lat')], k=n_s) # K$tot.withinss
      # K$centers = K$medoids
	  
    # Generate list of points to track
      loc_samples = cbind( loc_samples, "Station_j"=K$cluster )
      loc_stations = cbind( K$centers, "Station_j"=1:n_s )
      loc_grid = expand.grid( "long"=seq(Range_X[1],Range_X[2],length=ceiling(sqrt(Ngrid_sim))),"lat"=seq(Range_Y[1],Range_Y[2],length=ceiling(sqrt(Ngrid_sim))) )
      NN = nn2( query=loc_grid, data=loc_stations[,1:2], k=1 )
      loc_grid = cbind( loc_grid, "Station_j"=NN$nn.idx )
      loc_all = rbind( loc_stations, loc_samples, loc_grid )
                                                      
    # Get area for each location
      Voronoi_samples = calcVoronoi( cbind(X=loc_samples[,'long'], Y=loc_samples[,'lat']), xlim=range(Range_X,loc_samples[,'long']), ylim=range(Range_Y,loc_samples[,'lat']))
      Area_samples = calcArea( Voronoi_samples )[,2]
      Area_all = c( rep(0,n_s), Area_samples, rep(0,Ngrid_sim) )

    # Simulate null field
      Omega_s = RFsimulate(model=model_A, x=loc_all[,'long'], y=loc_all[,'lat'])@data[,1]
      Alpha_s = exp(Omega_s + log_mu_alpha)
    
    # Annual variation
    Epsilon_tmp_s = Epsilon_s = array(NA, dim=c(n_s+n_samp_per_year+Ngrid_sim, n_t))
    for(t in 1:n_t){
      Epsilon_s[,t] = Epsilon_tmp_s[,t] = RFsimulate(model_E, x=loc_all[,'long'], y=loc_all[,'lat'])@data[,1]
    }
      
    # Equilibrium conditions - recruitment
    if(RecFn=="Ricker"){
      mu_S_equil_per_area = log((1 - exp(-M-F_equil) - ro*exp(-M-F_equil) + ro*exp(-2*M-2*F_equil) ) / (exp(log_mu_alpha) * ( w_k - (w_k-alpha_g)*exp(-M-F_equil) ))) / (-beta)
      mu_R_equil_per_area = exp(log_mu_alpha)*mu_S_equil*exp(-beta*mu_S_equil)
    }
    if(RecFn=="BH"){
      mu_S_equil_per_area = ((exp(log_mu_alpha) * ( w_k - (w_k-alpha_g)*exp(-M-F_equil) ) / (1 - exp(-M-F_equil) - ro*exp(-M-F_equil) + ro*exp(-2*M-2*F_equil) )) - 1) / beta
      mu_R_equil_per_area = exp(log_mu_alpha)*mu_S_equil / (1+beta*mu_S_equil)
    }
    if(RecFn=="Constant"){
      mu_S_equil_per_area = (exp(log_mu_alpha) * ( w_k - (w_k-alpha_g)*exp(-M-F_equil) ) / (1 - exp(-M-F_equil) - ro*exp(-M-F_equil) + ro*exp(-2*M-2*F_equil) ))
      mu_R_equil_per_area = exp(log_mu_alpha)
      mu_S0_per_area = (exp(log_mu_alpha) * ( w_k - (w_k-alpha_g)*exp(-M-(0)) ) / (1 - exp(-M-(0)) - ro*exp(-M-(0)) + ro*exp(-2*M-2*(0)) ))
    }
    mu_N_equil_per_area = mu_R_equil_per_area / (1 - exp(-M-F_equil))
    S_equil = mu_S_equil_per_area * exp( Omega_s )
    R_equil = mu_R_equil_per_area * exp( Omega_s )
    N_equil = mu_N_equil_per_area * exp( Omega_s )
    W_equil = S_equil / N_equil
    S0 = mu_S0_per_area * exp( Omega_s )
    sum_S0 = sum(S0 * Area_all)
    if( abs(log(sum_S0 / (sum( R_equil * Area_all )/MRPSB)))>0.5  ) stop("Problem with area units")
    print( paste("Expected unfished catch rate = ", sum( mu_S_equil_per_area*AreaSwept ) ))
        
    #### Simulate B
    C_st = S_st = N_st = R_st = W_st = matrix(NA, nrow=n_s+n_samp_per_year+Ngrid_sim, ncol=n_t) # These are "per unit area"
    # Project forward
    for(t in 1:n_t){
      # Recruitment
      R_st[,t] = R_equil * exp(Epsilon_s[,t])
      # Numbers
      if(t==1) N_st[,t] = N_equil * exp(-M-F_equil) + R_st[,t]  
      if(t>=2) N_st[,t] = N_st[,t-1] * exp(-M-F_t[t-1]) + R_st[,t]  
      # Biomass
      #if(t==1) S_st[,t] = (1+ro)*exp(-M-F_equil)*S_equil - ro*exp(-2*M-F_equil-F_equil)*S_equil + w_k*R_st[,t] - (w_k-alpha_g)*exp(-F_equil-M)*R_equil
      #if(t==2) S_st[,t] = (1+ro)*exp(-M-F_t[t-1])*S_st[,t-1] - ro*exp(-2*M-F_t[t-1]-F_equil)*S_equil + w_k*R_st[,t] - (w_k-alpha_g)*exp(-F_t[t-1]-M)*R_st[,t-1]
      #if(t>=3) S_st[,t] = (1+ro)*exp(-M-F_t[t-1])*S_st[,t-1] - ro*exp(-2*M-F_t[t-1]-F_t[t-2])*S_st[,t-2] + w_k*R_st[,t] - (w_k-alpha_g)*exp(-F_t[t-1]-M)*R_st[,t-1]
      if(t==1) S_st[,t] = exp(-M-F_equil)* (alpha_g*N_equil + ro*S_equil) + w_k*R_st[,t]
      if(t>=2) S_st[,t] = exp(-M-F_t[t-1])* (alpha_g*N_st[,t-1] + ro*S_st[,t-1]) + w_k*R_st[,t]
      # Weight
      W_st[,t] = S_st[,t] / N_st[,t]
      # Fishing mortality
      if(t==1) F_t[t] = F_equil 
      if(t>=2) F_t[t] = F_t[t-1] * (sum(Area_all*S_st[,t])/sum(Area_all*S_equil)/S_bioecon)^Accel * exp( Eps_F[t] - SD_F^2/2 )
      # Catches
      C_st[,t] = F_t[t]/(M+F_t[t]) * (1-exp(-M-F_t[t])) * S_st[,t]
    }
    ## Simulate data
    # Survey catches
      Y = I_st = rpois(n=n_samp_per_year*n_t, lambda=AreaSwept*q_I*N_st[n_s+1:n_samp_per_year,] )
      Y_stations = matrix(Y, ncol=n_t, byrow=FALSE)
      NAind_stations = array( as.integer(ifelse(is.na(Y_stations),1,0)), dim=dim(Y_stations))
    # Average weight
      Wobs_st = array(rnorm(n=n_samp_per_year*n_t, mean=W_st[n_s+1:n_samp_per_year,], sd=W_st[n_s+1:n_samp_per_year,]*CV_w), dim=c(n_samp_per_year,n_t))
    # Total catch
      C_t = colSums(C_st * outer(Area_all,rep(1,ncol(C_st))))
      Cobs_t = (C_t + C_t*Eps_C)
    # Total abundance
      S_t = colSums(S_st*outer(Area_all,rep(1,ncol(C_st))))
      R_t = colSums(R_st*outer(Area_all,rep(1,ncol(C_st))))
      N_t = colSums(N_st*outer(Area_all,rep(1,ncol(C_st))))
  
    # Make INLA inputs for TMB
      if(MeshType=="Recommended"){
        mesh_stations = inla.mesh.create( loc_stations[,c('lat','long')], plot.delay=NULL, extend=list(n=8,offset=-0.15), refine=list(min.angle=26,max.edge.data=0.08,max.edge.extra=0.2) )  # loc_samp
      }
      if(MeshType=="Minimal"){
        mesh_stations = inla.mesh.create( loc_stations[,c('lat','long')], plot.delay=NULL, extend=list(n=8,offset=-0.15), refine=F )  # loc_samp
      }
      spde_stations = inla.spde2.matern(mesh_stations, alpha=2)  # Given 2D field, alpha=2 <=> Matern Nu=1 (Lindgren and Rue 2013, between Eq. 1&2)
      n_i = nrow(mesh_stations$loc)
  
    # Get area for each location
      Voronoi_s = calcVoronoi( cbind(X=loc_stations[,1], Y=loc_stations[,2]), xlim=Range_X, ylim=Range_Y)
      Area_s = calcArea( Voronoi_s )[,2]
      Area_i = c( Area_s, rep(0,n_i-n_s) )

    # Save Varanoi stuff
      png( file=paste(RepFile,"Voronoi.png",sep=""), width=6, height=3, res=200, units="in")
      par( mfrow=c(1,2), mar=c(3,3,2,0), mgp=c(2,0.5,0) )
      # Visualize samples
      plotMap( Voronoi_samples )
      points( y=loc_samples[,'lat'], x=loc_samples[,'long'], col=c("red","black")[rep.int( c(1,2), times=c(n_s,n_samp_per_year-n_s))], pch=loc_samples[,'Station_j'] )
      # Visualize stations
      plotMap( Voronoi_s )
      points( y=loc_samples[,'lat'], x=loc_samples[,'long'], col=c("red","black")[rep.int( c(1,2), times=c(n_s,n_samp_per_year-n_s))] )
      dev.off()

    # Convert data to data frame
      DF = data.frame( "I_j"=as.vector(I_st), "W_j"=as.vector(Wobs_st), "Station_j"=rep(loc_samples[,'Station_j'],n_t), "Year_j"=as.vector(col(Wobs_st)), "AreaSwept_j"=AreaSwept, "long_j"=rep(loc_samples[,'long'],n_t), "lat_j"=rep(loc_samples[,'lat'],n_t) )
      n_j = nrow(DF)
        
    # Calculate index of abundance
      PredDF = data.frame( "Year_j"=sort(unique(DF[,'Year_j'])) )
      Glm_I = zeroinfl(I_j ~ 0 + factor(Year_j) | 1, data=DF, dist=c("poisson","negbin","geometric")[2], method="Nelder-Mead", control=zeroinfl.control("EM"=TRUE))
      Glm_W = lm(W_j ~ 0 + factor(Year_j), data=DF, subset=DF[,'I_j']>0)
      #Index_hat = cbind("Mean_I"=predict(Glm_I, newdata=PredDF, type="response")/(AreaSwept*q_I), "logSD_I"=summary(Glm_I)$coef$count[1:n_t,'Std. Error'])
      #Index_hat = cbind(Index_hat, "Mean_W"=predict(Glm_W, newdata=PredDF, type="response"), "SD_W"=summary(Glm_W)$coef[1:n_t,'Std. Error'])
      Index_hat = array(-999, dim=c(n_t,4), dimnames=list(Year_Range[1]:Year_Range[2],c("Mean_I","logSD_I","Mean_W","SD_W")))
      Index_hat[unique(DF[,'Year_j']),c("Mean_I","logSD_I")] = cbind(predict(Glm_I, newdata=PredDF, type="response")/(AreaSwept*q_I)*DomainArea, summary(Glm_I)$coef$count[1:length(unique(DF[,'Year_j'])),'Std. Error'])
      Index_hat[unique(DF[,'Year_j']),c("Mean_W","SD_W")] = cbind(predict(Glm_W, newdata=PredDF, type="response"), summary(Glm_W)$coef[1:length(unique(DF[,'Year_j'])),'Std. Error'])
    
    # Plot
      png( file=paste(RepFile,"True_Timeseries.png",sep=""), width=3*3, height=2*3, res=200, units="in")
      par(mfrow=c(2,3), mar=c(3,3,2,0))
      plot( C_t, type="l", ylim=c(0,max(C_t)), main="Catch" )
      plot( F_t, type="l", ylim=c(0,max(F_t)), main="Fishing mortality" )
      plot( S_t, type="l", ylim=c(0,max(S_t)), main="Spawning biomass" )
      plot( S_t/sum_S0, type="l", ylim=c(0,1.2), main="Depletion" )
      plot( N_t, type="l", ylim=c(0,max(N_t)), main="Total abundance" )
      for(t in 1:n_t){
        points( x=t, y=Index_hat[t,'Mean_I'], col="red")
        lines(x=c(t,t), y=Index_hat[t,'Mean_I']*exp(c(-1,1)*Index_hat[t,'logSD_I']), col="red" )
      }
      plot( S_t/N_t, type="l", ylim=c(0,max(S_t/N_t)), main="Average weight" )
      for(t in 1:n_t){
        points( x=t, y=Index_hat[t,'Mean_W'], col="red")
        lines(x=c(t,t), y=Index_hat[t,'Mean_W']+c(-1,1)*Index_hat[t,'SD_W'], col="red" )
      }
      dev.off()
    
    # Save true stuff
      TrueList = list( "n_s"=n_s, "S_t"=S_t, "R_t"=R_t, "C_t"=C_t, "F_t"=F_t, "N_t"=N_t, "sum_S0"=sum_S0, "S_st"=S_st, "R_st"=R_st, "C_st"=C_st, "W_st"=W_st, "N_st"=N_st, "Index_hat"=Index_hat, "DF"=DF, "Omega_s"=Omega_s, "Alpha_s"=Alpha_s, "Epsilon_s"=Epsilon_s )
      save(TrueList, file=paste(RepFile,"TrueList.RData",sep=""))
      
  for(ModelI in 1:length(ModelSet)){
  
    # Model File
    Model = ModelSet[ModelI]
    ModelFile = paste(RepFile,"Model=",ModelSet[ModelI],"/",sep="")
    dir.create(ModelFile)
    
    # Run TMB
      if(FALSE){
        # Compile spatial
        if(TRUE){
          dyn.unload( paste(TmbFile,dynlib(Version),sep="") )
          file.remove( paste(TmbFile,Version,c(".o",".dll"),sep="") )
        }
        setwd( TmbFile )
        compile( paste(Version,".cpp",sep="") )
      }
        
      # Covariates
        X_stations = cbind( rep(1,n_s*n_t) )
        DF_input = DF
        if(ModelSet[ModelI]=="Index"){
          IndexMat = Index_hat
          DF_input[,c("I_j","W_j")] = -999
        }
        if(ModelSet[ModelI]!="Index"){
          IndexMat = array(-999, dim=dim(Index_hat) )
        }
        #DF_input[,"I_j"] = -999  # Removing I_j still results in biased C in 5th/6th year
        #DF_input[,"W_j"] = -999  # Removing W_j is much more biased C in 5th/6th year
        
      # Run spatial model
        dyn.load( paste(TmbFile,dynlib(Version),sep="") )
        if(Version=="delay_difference_v4d"){
          #dyn.load( paste(TmbFile,dynlib("delay_difference_v4d"),sep="") )
          Error_Model = ErrorModel_CatchRates
          Data = list(Error_Model=Error_Model, Smooth_F=Smooth_F, n_j=n_j, n_i=n_i, n_s=n_s, n_t=n_t, I_j=DF_input[,'I_j'], W_j=DF_input[,'W_j'], AreaSwept_j=DF_input[,'AreaSwept_j'], Station_j=DF_input[,'Station_j']-1, Year_j=DF_input[,'Year_j']-1, C_t=C_t, meshidxloc=mesh_stations$idx$loc-1, G0=spde_stations$param.inla$M0, G1=spde_stations$param.inla$M1, G2=spde_stations$param.inla$M2, Area_i=Area_i, alpha_g=alpha_g, ro=ro, w_k=w_k, M=M, k=k, CV_c=CV_c, CV_w=CV_w)
          Parameters = list(log_F_sd=log(1), log_F_equil=log(F_equil), log_F_t_input=log(F_t), log_q_I=log(q_I), beta=c(0.0), log_tau_E=0.0, log_tau_O=0.0, log_kappa=0.0,	ln_VarInfl=c(0.0,0.0), Epsilon_input=matrix(0,spde_stations$n.spde,n_t), Omega_input=rep(0,spde_stations$n.spde))
          #obj <- MakeADFun(data=Data, parameters=Parameters, random=c("Epsilon_input","Omega_input"), hessian=TRUE)
        }
        if(Version=="delay_difference_v6c"){    
          #dyn.load( paste(TmbFile,dynlib("delay_difference_v6c"),sep="") )
          Data = list(ModelType=ifelse(Model=="Spatial",1,2), ErrorModel_CatchRates=ErrorModel_CatchRates, ErrorModel_MeanWeight=ErrorModel_MeanWeight, Smooth_F=Smooth_F, n_j=n_j, n_i=n_i, n_s=n_s, n_t=n_t, I_j=DF_input[,'I_j'], W_j=DF_input[,'W_j'], AreaSwept_j=DF_input[,'AreaSwept_j'], Station_j=DF_input[,'Station_j']-1, Year_j=DF_input[,'Year_j']-1, C_t=C_t, meshidxloc=mesh_stations$idx$loc-1, G0=spde_stations$param.inla$M0, G1=spde_stations$param.inla$M1, G2=spde_stations$param.inla$M2, Area_i=Area_i, alpha_g=alpha_g, ro=ro, w_k=w_k, M=M, k=k, CV_c=CV_c, CV_w=CV_w)
          Parameters = list(log_F_sd=log(1), log_F_equil=log(F_equil), log_F_t_input=log(F_t), log_q_I=log(q_I), beta=c(0.0), log_tau_E=0.0, log_tau_O=0.0, log_kappa=0.0,	ln_VarInfl=c(0.0,0.0), log_extraCV_w=log(0.05), log_tau_N=log(1), Epsilon_input=matrix(0,spde_stations$n.spde,n_t), Omega_input=rep(0,spde_stations$n.spde), Nu_input=rep(0,n_t))
          #obj <- MakeADFun(data=Data, parameters=Parameters, random=c("Epsilon_input","Omega_input"), hessian=TRUE)
        } 
        if(Version=="delay_difference_v8e"){    
          Data = list(ModelType=ifelse(Model=="Spatial",1,2), ErrorModel_CatchRates=ErrorModel_CatchRates, ErrorModel_MeanWeight=ErrorModel_MeanWeight, Smooth_F=Smooth_F, n_j=n_j, n_i=n_i, n_s=n_s, n_t=n_t, I_j=DF_input[,'I_j'], W_j=DF_input[,'W_j'], AreaSwept_j=DF_input[,'AreaSwept_j'], Station_j=DF_input[,'Station_j']-1, Year_j=DF_input[,'Year_j']-1, C_t=C_t, IndexMat=IndexMat, meshidxloc=mesh_stations$idx$loc-1, G0=spde_stations$param.inla$M0, G1=spde_stations$param.inla$M1, G2=spde_stations$param.inla$M2, Area_i=Area_i, alpha_g=alpha_g, ro=ro, w_k=w_k, M=M, k=k, CV_c=CV_c, CV_w=CV_w)
          Parameters = list(log_F_sd=log(1), log_F_t_input=log(c(F_equil,F_t)), log_q_I=log(q_I), beta=log_mu_alpha, log_tau_E=0.0, log_tau_O=0.0, log_kappa=0.0,	ln_VarInfl=c(0.0,0.0), log_extraCV_w=log(0.05), log_tau_N=log(1), log_extraCV_Index=rep(log(0.1),2), Epsilon_input=matrix(0,spde_stations$n.spde,n_t), Omega_input=rep(0,spde_stations$n.spde), Nu_input=rep(0,n_t))
        }
		if(Version=="delay_difference_v9"){    
          Data = list(ModelType=ifelse(Model=="Spatial",1,2), ErrorModel_CatchRates=ErrorModel_CatchRates, ErrorModel_MeanWeight=ErrorModel_MeanWeight, Smooth_F=Smooth_F, n_j=n_j, n_i=n_i, n_s=n_s, n_t=n_t, I_j=DF_input[,'I_j'], W_j=DF_input[,'W_j'], AreaSwept_j=DF_input[,'AreaSwept_j'], Station_j=DF_input[,'Station_j']-1, Year_j=DF_input[,'Year_j']-1, C_t=C_t, IndexMat=IndexMat, meshidxloc=mesh_stations$idx$loc-1, G0=spde_stations$param.inla$M0, G1=spde_stations$param.inla$M1, G2=spde_stations$param.inla$M2, Area_i=Area_i, alpha_g=alpha_g, ro=ro, w_k=w_k, M=M, k=k, CV_c=CV_c, CV_w=CV_w)
          Parameters = list(log_F_sd=log(1), log_F_t_input=log(c(F_equil,F_t)), log_q_I=log(q_I), beta=log_mu_alpha, log_tau_E=0.0, log_tau_O=0.0, log_kappa=0.0,	ln_VarInfl=c(0.0,0.0), log_extraCV_w=log(0.05), Strata = rep(0, n_s), log_tau_N=log(1), log_extraCV_Index=rep(log(0.1),2), Epsilon_input=matrix(0,spde_stations$n.spde,n_t), Omega_input=rep(0,spde_stations$n.spde), Nu_input=rep(0,n_t))
        } 
      # Define random effects
        if( Model=="Spatial" ){
          if( Smooth_F==0 ) Random = c("Epsilon_input", "Omega_input")
          if( Smooth_F!=0 ) Random = c("Epsilon_input", "Omega_input", "log_F_equil")
        }
        if( Model=="Nonspatial" ){
          if( Smooth_F==0 ) Random = c("Nu_input")
          if( Smooth_F!=0 ) Random = c("Nu_input", "log_F_t_input")
        }
        if( Model=="Index" ){
          if( Smooth_F==0 ) Random = c("Nu_input")
          if( Smooth_F!=0 ) Random = c("Nu_input", "log_F_t_input")
        }
        if( Model=="Strata" ){
          if( Smooth_F==0 ) Random = c("Nu_input")
          if( Smooth_F!=0 ) Random = c("Nu_input", "log_F_t_input")
        }
      # Define fixed parameters
       Map = list()
       if( Smooth_F==0 ) Map[["log_F_sd"]] = factor(NA)
       if( ErrorModel_CatchRates==0 ) Map[["ln_VarInfl"]] = factor(c(NA,NA))
       if( ErrorModel_MeanWeight==0 ) Map[["log_extraCV_w"]] = factor(NA)
       if( Model=="Spatial" ){
         Map[["log_tau_N"]] = factor(NA)
         Map[["Nu_input"]] = factor( rep(NA,n_t) )
		 Map[["Strata"]] = factor(rep(NA,n_s))
         Map[["log_extraCV_Index"]] = factor( rep(NA,2) )
         Map[["log_F_t_input"]] = factor( c(1,1:n_t) )
       }
       if( Model=="Nonspatial" ){
         Map[["log_tau_E"]] = factor(NA)
         Map[["log_tau_O"]] = factor(NA)
 		 Map[["Strata"]] = factor(rep(NA,n_s))
         Map[["Epsilon_input"]] = factor( rep(NA,length(Parameters$Epsilon_input)) )
         Map[["Omega_input"]] = factor( rep(NA,length(Parameters$Omega_input)) )
         Map[["log_kappa"]] = factor(NA)
         Map[["log_extraCV_Index"]] = factor( rep(NA,2) )
         Map[["log_F_t_input"]] = factor( c(1,1:n_t) )
       }
       if( Model=="Index" ){
         Map[["log_tau_E"]] = factor(NA)
         Map[["log_tau_O"]] = factor(NA)
 		 Map[["Strata"]] = factor(rep(NA,n_s))
         Map[["Epsilon_input"]] = factor( rep(NA,length(Parameters$Epsilon_input)) )
         Map[["Omega_input"]] = factor( rep(NA,length(Parameters$Omega_input)) )
         Map[["log_kappa"]] = factor(NA)
         Map[["ln_VarInfl"]] = factor(c(NA,NA))
         Map[["log_extraCV_w"]] = factor(NA)
         Map[["log_F_t_input"]] = factor( c(1,1:n_t) )
       }
        if( Model=="Strata" ){
         Map[["beta"]] = factor(NA)
         Map[["log_tau_E"]] = factor(NA)
         Map[["log_tau_O"]] = factor(NA)
         Map[["Epsilon_input"]] = factor( rep(NA,length(Parameters$Epsilon_input)) )
         Map[["Omega_input"]] = factor( rep(NA,length(Parameters$Omega_input)) )
         Map[["log_kappa"]] = factor(NA)
         Map[["log_extraCV_Index"]] = factor( rep(NA,2) )
         Map[["log_F_t_input"]] = factor( c(1,1:n_t) )
       }
       if( Fix_Q==TRUE ) Map[["log_q_I"]] = factor(NA)
      # Save stuff
        MapList = list( "Map"=Map, "Random"=Random)
        capture.output(MapList, file=paste(ModelFile,"MapList.txt",sep=""))
      # Build object                                                              #    
        # no random effect for the annual variation in recruitment
			if(Model!="Spatial")
			{
				Map1 <- Map
				Map1[["log_tau_N"]] = factor(NA)
				obj <- MakeADFun(data=Data, parameters=Parameters, random=NULL, hessian=TRUE, map=Map1)
			}
		# with the annual variation in recruitment as random effect
			if(Model=="Spatial")
			{
				obj <- MakeADFun(data=Data, parameters=Parameters, random=Random, hessian=TRUE, map=Map)
			}
			
      # Separability
        if(FALSE){
          h = obj$env$spHess(random=TRUE)
          image(h)
        }
      
      # Settings
        obj$control <- list(trace=1, parscale=1, REPORT=1, reltol=1e-12, maxit=1000)
        obj$hessian <- FALSE
        obj$fn(obj$par)
      # Run optimizer
      tic()
        opt = nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, lower=-20, upper=20, control=list(trace=1, eval.max=1e4, iter.max=1e4))
        Report = obj$report()
        Sdreport = try( sdreport(obj) )
      toc()

		# Save results
		capture.output( opt, file=paste(ModelFile,"opt.txt",sep=""))
		capture.output( Report, file=paste(ModelFile,"Report.txt",sep=""))
		Save = list( "opt"=opt, "Report"=Report, "obj"=obj, "Sdreport"=Sdreport)
		save( Save, file=paste(ModelFile,"Save.RData",sep=""))

		# Re-load
		  if(FALSE){
			load(file=paste(ModelFile,"Save.RData",sep=""))
			attach(Save)
		  }
		
		# Debug
		  if(FALSE){
			# Biased estimates of strata-specific abundance
			plot(N_st[,1], Report$N_it[1:n_s,1])
			  abline(a=0,b=1)
			# 
			Which = which(Data$Year_j==0)
			Data$I_j[Which]
			Report$I_j_hat[Which]
			# Could it be a bias correction issue?
			apply(Report$Epsilon[1:n_s,], MARGIN=2, FUN=sd)
		  }
		
		# Check correlation of random fields
		  Cor_Omega = cor( Report$Omega[1:n_s], Omega_s[1:n_s] )
		  Cor_Epsilon = rep(NA, n_t)
		  for(t in 1:n_t) Cor_Epsilon[t] = cor( Report$Epsilon[1:n_s,t], Epsilon_s[1:n_s,t] )
		
		########### Simulated vs. True comparison
		# Time series
		  FUN = function(InputMat) c( InputMat[,1]-1*InputMat[,2], rev(InputMat[,1]+1*InputMat[,2]))
		  png( file=paste(ModelFile,"True_vs_Est.png",sep=""), width=9, height=6, res=200, units="in")
			par(mfrow=c(2,3), mar=c(3,3,2,0))
			# Abundance
			Mat = cbind( "Year"=unique(DF[,'Year_j']), "True"=N_t, "Est"=Report$N_t_hat)
			matplot( y=Mat[,c("True",'Est')], x=Mat[,c("Year")], type="l", col=c("black","red"), lty="solid", ylim=c(0,max(Mat[,c('True','Est')])), main="Abundance")
			if( !("condition" %in% names(attributes(Sdreport))) ) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="N_t_hat"),]), x=c(Mat[,c("Year")],rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)  
			# Biomass
			Mat = cbind( "Year"=unique(DF[,'Year_j']), "True"=S_t, "Est"=Report$S_t_hat)
			matplot( y=Mat[,c("True",'Est')], x=Mat[,c("Year")], type="l", col=c("black","red"), lty="solid", ylim=c(0,max(Mat[,c('True','Est')])), main="Spawning biomass")
			if( !("condition" %in% names(attributes(Sdreport))) ) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_t_hat"),]), x=c(Mat[,c("Year")],rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)  
			# Average weight
			Mat = cbind( "Year"=unique(DF[,'Year_j']), "True"=S_t/N_t, "Est"=Report$S_t_hat/Report$N_t_hat)
			matplot( y=Mat[,c("True",'Est')], x=Mat[,c("Year")], type="l", col=c("black","red"), lty="solid", ylim=c(0,max(Mat[,c('True','Est')])), main="Average weight")
			if( !("condition" %in% names(attributes(Sdreport))) ) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_t_hat / N_t_hat"),]), x=c(Mat[,c("Year")],rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)  
			# Recruitment
			Mat = cbind( "Year"=unique(DF[,'Year_j']), "True"=R_t, "Est"=Report$R_t_hat)
			matplot( y=Mat[,c("True",'Est')], x=Mat[,c("Year")], type="l", col=c("black","red"), lty="solid", ylim=c(0,max(Mat[,c('True','Est')])), main="Recruitment")
			if( !("condition" %in% names(attributes(Sdreport))) ) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="R_t_hat"),]), x=c(Mat[,c("Year")],rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)  
			# Fishing mortality
			Mat = cbind( "Year"=unique(DF[,'Year_j']), "True"=F_t, "Est"=Report$F_t)
			matplot( y=Mat[,c("True",'Est')], x=Mat[,c("Year")], type="l", col=c("black","red"), lty="solid", ylim=c(0,max(Mat[,c('True','Est')])), main="Fishing mortality")
			if( !("condition" %in% names(attributes(Sdreport))) ) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="F_t"),]), x=c(Mat[,c("Year")],rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)  
			# Depletion
			Mat = cbind( "Year"=unique(DF[,'Year_j']), "True"=S_t/sum_S0, "Est"=Report$S_t/Report$sum_S0)
			matplot( y=Mat[,c("True",'Est')], x=Mat[,c("Year")], type="l", col=c("black","red"), lty="solid", ylim=c(0,max(Mat[,c('True','Est')])), main="Depletion")
			if( !("condition" %in% names(attributes(Sdreport))) ) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_t_hat / sum_S0"),]), x=c(Mat[,c("Year")],rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)  
		  dev.off()

		# fields
		  Years2Show = 1:3
		  png( file=paste(ModelFile,"True_vs_Est_Recruitment.png",sep=""), width=2*2, height=(1+length(Years2Show))*2, res=200, units="in")
			par(mfrow=c(1+length(Years2Show),2), mar=c(0.2,0.2,0,0), oma=c(4,4,2,2), mgp=c(2,0.5,0), tck=-0.02)
			Zlim = range( log(c(R_st[-c(1:(n_s+n_samp_per_year)),Years2Show],Report$R_it[1:n_s,Years2Show])) ) 
			X = seq(Range_X[1],Range_X[2],length=ceiling(sqrt(Ngrid_sim)))
			Y = seq(Range_Y[1],Range_Y[2],length=ceiling(sqrt(Ngrid_sim)))
			# Expected recruitment
			image( x=X, y=Y, z=matrix(Omega_s[-c(1:(n_s+n_samp_per_year))]+log_mu_alpha, ncol=ceiling(sqrt(Ngrid_sim)), nrow=ceiling(sqrt(Ngrid_sim))), zlim=Zlim, col=Col(10), xaxt="n", yaxt="n" )
			mtext(side=3, line=0.5, outer=FALSE, text="Simulated")
			mtext(side=2, line=2, outer=FALSE, text=expression(italic(beta)+bold(Omega)))
			axis(2)
			image( x=X, y=Y, z=matrix(Report$Omega[loc_grid[,'Station_j']]+opt$par['Strata']+opt$par['beta'], ncol=ceiling(sqrt(Ngrid_sim)), nrow=ceiling(sqrt(Ngrid_sim))), zlim=Zlim, col=Col(10), xaxt="n", yaxt="n" )
			mtext(side=3, line=0.5, outer=FALSE, text="Estimated")
			# Realized recruitment
			for(t in 1:length(Years2Show)){
			  image( x=X, y=Y, z=matrix(log(R_st[-c(1:(n_s+n_samp_per_year)),Years2Show[t]]), ncol=ceiling(sqrt(Ngrid_sim)), nrow=ceiling(sqrt(Ngrid_sim))), zlim=Zlim, col=Col(10), xaxt="n", yaxt="n" )
			  if(t==length(Years2Show)) axis(1)
			  axis(2)
			  #if(t==1) mtext(side=3, line=0.5, outer=FALSE, text="True")
			  mtext(side=2, line=2, outer=FALSE, text=switch(t,expression(bold(R[1])),expression(bold(R[2])),expression(bold(R[3]))) )
			  image( x=X, y=Y, z=matrix(log(Report$R_it[cbind(loc_grid[,'Station_j'],Years2Show[t])]), ncol=ceiling(sqrt(Ngrid_sim)), nrow=ceiling(sqrt(Ngrid_sim))), zlim=Zlim, col=Col(10), xaxt="n", yaxt="n" )
			  if(t==length(Years2Show)) axis(1)
			  #if(t==1) mtext(side=3, line=0.5, outer=FALSE, text="Est")
			}
			mtext(side=1, outer=TRUE, line=2, text="Eastings (km)")
			mtext(side=4, outer=TRUE, line=0.5, text="Northings (km)")      
		  dev.off()
		
		# Check probabilities of data
		  if(FALSE){
			cbind( Report$log_pIndexMat_0, Report$log_pIndexMat_2, Report$log_pC_t )
			Report$extraCV_Index
		  }
		  
		# Check results
		  opt$par[c("beta")]
		  unlist( Report[c('Range','SigmaE','SigmaO','q_I')] )

		# Population trends
		  png( file=paste(ModelFile,"Est_S_and_N.png",sep=""), width=6, height=6, res=200, units="in")
		    par(mfrow=c(2,2), mar=c(3,3,2,0))
		    # Abundance
		    Glm = glm( I_j ~ 0 + factor(Year_j) + factor(Station_j), data=DF, family="poisson")
		    Mat = cbind( "Mean"=tapply(DF[,'I_j'],INDEX=DF[,'Year_j'],FUN=mean), "Glm"=exp(Glm$coef[1:length(unique(DF[,'Year_j']))]), "Year"=unique(DF[,'Year_j']))
		    Mat[,c('Mean','Glm')] = Mat[,c('Mean','Glm')] * outer( rep(1,nrow(Mat)), mean(Report$N_t_hat)/c(mean(Mat[,"Mean"]),mean(Mat[,"Glm"])) ) 
		    matplot( y=Report$N_t_hat, x=Year_Range[1]:Year_Range[2], type="l", col="black", lty="solid", ylim=c(0,max(Report$N_t_hat,Mat[,c('Mean','Glm')])), main="Abundance")
		    matplot( y=Mat[,c("Mean",'Glm')], x=(Year_Range[1]:Year_Range[2])[Mat[,"Year"]], col=c("red","blue"), type="p", pch=21, add=TRUE)
		    # Biomass
		    Glm = glm( I(I_j*W_j) ~ 0 + factor(Year_j) + factor(Station_j), data=DF, family="gaussian")
		    Mat = cbind( "Mean"=tapply(DF[,'W_j']*DF[,'I_j'],INDEX=DF[,'Year_j'],FUN=mean,na.rm=TRUE), "Glm"=Glm$coef[1:length(unique(DF[,'Year_j']))], "Year"=unique(DF[,'Year_j']))
		    Mat[,c('Mean','Glm')] = Mat[,c('Mean','Glm')] * outer( rep(1,nrow(Mat)), mean(Report$S_t_hat)/c(mean(Mat[,"Mean"]),mean(Mat[,"Glm"])) ) 
		    matplot( y=Report$S_t_hat, x=Year_Range[1]:Year_Range[2], type="l", col="black", lty="solid", ylim=c(0,max(Report$S_t_hat,Mat[,c('Mean','Glm')])), main="Spawning biomass")
		    matplot( y=Mat[,c("Mean",'Glm')], x=(Year_Range[1]:Year_Range[2])[Mat[,"Year"]], col=c("red","blue"), type="p", pch=21, add=TRUE)
		    # Average weight
		    Glm = glm( I(ifelse(W_j==Inf,NA,W_j)) ~ 0 + factor(Year_j) + factor(Station_j), data=DF, family="gaussian", na.action="na.omit")
		    Mat = cbind( "Mean"=tapply(DF[which(DF[,'I_j']>0),'W_j'],INDEX=DF[which(DF[,'I_j']>0),'Year_j'],FUN=mean,na.rm=TRUE), "Glm"=Glm$coef[1:length(unique(DF[,'Year_j']))], "Year"=unique(DF[,'Year_j']))
		    Mat[,c('Mean','Glm')] = Mat[,c('Mean','Glm')] * outer( rep(1,nrow(Mat)), mean(Report$S_t_hat/Report$N_t_hat)/c(mean(Mat[,"Mean"]),mean(Mat[,"Glm"])) ) 
		    matplot( y=Report$S_t_hat/Report$N_t_hat, x=Year_Range[1]:Year_Range[2], type="l", col="black", lty="solid", ylim=c(0,max(Report$S_t_hat/Report$N_t_hat,Mat[,c('Mean','Glm')])), main="Average weight")
		    matplot( y=Mat[,c("Mean",'Glm')], x=(Year_Range[1]:Year_Range[2])[Mat[,"Year"]], col=c("red","blue"), type="p", pch=21, add=TRUE)
		    # Recruitment
		    matplot( y=Report$R_t_hat, x=Year_Range[1]:Year_Range[2], type="l", col="black", lty="solid", ylim=c(0,max(Report$R_t_hat)), main="Recruits")
		  dev.off()

		# Fishing mortality
		  png( file=paste(ModelFile,"Est_Frate.png",sep=""), width=4, height=4, res=200, units="in")
		    par(mar=c(3,3,2,0), mgp=c(2,0.5,0), tck=-0.02, xaxs="i", yaxs="i")
		    matplot( y=cbind(Report$F_t,Report$Exploit_Frac), x=Year_Range[1]:Year_Range[2], type="l", col="black", lty=c("solid","dotted"), ylim=c(0,max(cbind(Report$F_t,Report$Exploit_Frac))), xlab="", ylab="")
		  dev.off()
		
		# Catch history
		  png( file=paste(ModelFile,"Pred_v_Obs_Catch.png",sep=""), width=4, height=4, res=200, units="in")
		    matplot( y=cbind(C_t,Report$C_t_hat), x=Year_Range[1]:Year_Range[2], col=c("black","red"), type=c("l","p"), pch=21, ylim=c(0,max(cbind(C_t,Report$C_t_hat))))
		  dev.off()

		# Weight -- pred vs. obs
		  png( file=paste(Dir_save,"Pred_v_Obs_MeanWeight.png",sep=""), width=4, height=4, res=200, units="in")
		    Mat = cbind( "Obs"=DF[which(DF[,'I_j']>0),'W_j'], "Pred"=Report$W_it[as.matrix(DF[which(DF[,'I_j']>0),c('Station_j','Year_j')])] )
		    # Cbind( Report$W_it[as.matrix(DF[which(DF[,'I_j']>0),c('Station_j','Year_j')])], (Report$S_it/Report$N_it)[as.matrix(DF[which(DF[,'I_j']>0),c('Station_j','Year_j')])]
		    plot(x=Mat[,"Pred"], y=Mat[,"Obs"], xlim=c(0,max(Mat,na.rm=TRUE)), ylim=c(0,max(Mat,na.rm=TRUE)), col=rgb(0,0,0,alpha=0.02), , ylab="Obs", xlab="Pred"  ) 
		    abline(a=0, b=1, col="red")                                                         # site 42 is wrong somehow
		  dev.off()
		#  (1+ro)*exp(-M-Report$F_equil)*Report$S_equil - ro*exp(-2*M-Report$F_equil-Report$F_equil)*Report$S_equil + w_k*Report$R_it[,1] - (w_k-alpha_g)*exp(-Report$F_equil-M)*Report$R_equil
		
		# Catch rates -- pred vs. obs
		  png( file=paste(Dir_save,"Pred_v_Obs_CatchRates.png",sep=""), width=4, height=4, res=200, units="in")
		    Mat = cbind( "Obs"=DF[,'I_j'], "Pred"=Report$I_j_hat )
		    plot(x=sqrt(Mat[,"Pred"]), y=sqrt(Mat[,"Obs"]), xlim=c(0,max(sqrt(Mat),na.rm=TRUE)), ylim=c(0,max(sqrt(Mat),na.rm=TRUE)), xaxt="n", yaxt="n", col=rgb(0,0,0,alpha=0.02), ylab="Obs", xlab="Pred"  ) 
		    axis(1, at=axTicks(1), labels=axTicks(1)^2)
		    axis(2, at=axTicks(2), labels=axTicks(2)^2)
		    abline(a=0, b=1, col="red")
		  dev.off()
		  			
		##### Fields
		# Productivity
		  png( file=paste(ModelFile,"Fields_Omega_LL.png",sep=""), width=3, height=6, res=200, units="in")
		    par(mfrow=c(2,1), mar=c(0,0,0,0))
		    #map("worldHires", ylim=y.lim,xlim=x.lim,col="grey90",fill=T, main="", mar=c(0,0,2.5,0))
		    plot(1, type="n", ylim=range(loc_stations[,'long']), xlim=range(loc_stations[,'lat']))
		    #points(y=loc_stations[,'lat'], x=loc_stations[,'long'], col=Col(n=10)[ceiling(f(matrix(Report$Omega[1:n_s],ncol=1))*9)+1], cex=3, pch=20)
		    Bin = Bin_Quantile( Report$Omega[1:n_s], Nregions=5 )  
		    points(y=loc_stations[,'lat'], x=loc_stations[,'long'], col=Col(Bin$Nregions)[Bin$Region], cex=3, pch=20)
		    Legend( Bin, Col )
		  dev.off()
		  # Dynamical states
		  Ncol = ceiling(sqrt(n_t+1)); Nrow = ceiling( (n_t+1)/Ncol )
		  for(FigI in 1:5){
		    if(FigI==1) Mat = (Report$Epsilon)[1:n_s,]
		    if(FigI==2) Mat = log(Report$S_it)[1:n_s,]
		    if(FigI==3) Mat = log(Report$W_it)[1:n_s,]
		    if(FigI==4) Mat = log(Report$N_it)[1:n_s,]
		    if(FigI==5) Mat = log(Report$R_it)[1:n_s,]
		    png( file=paste(ModelFile,"Fields_",c("Epsilon","S","W","N","R")[FigI],"_LL.png",sep=""), width=4*Ncol, height=2*Nrow, res=200, units="in")
		  	par(mfrow=c(Nrow,Ncol), mar=c(0,0,0,0), oma=c(2,2,0,0) )
		  	for(i in 1:n_t){
		  	  #map("worldHires", ylim=y.lim,xlim=x.lim,col="grey90",fill=T, main="", mar=c(0,0,2.5,0))
		  	  plot(1, type="n", ylim=range(loc_stations[,'long']), xlim=range(loc_stations[,'lat']))
		  	  title( (Year_Range[1]:Year_Range[2])[i] )
		  	  #points(y=loc_stations[,'lat'], x=loc_stations[,'long'], col=Col(n=10)[ceiling(f(Mat[1:n_s,])[,i]*9)+1], cex=3, pch=20)
		  	  Bin = Bin_Quantile( Mat, Nregions=5 )  
		  	  points(y=loc_stations[,'lat'], x=loc_stations[,'long'], col=Col(Bin$Nregions)[Bin$Region[,i]], cex=3, pch=20)
		  		if( (i-1)%%Ncol == 0) axis(2)
		  		if( (i-1)/Ncol >= (Nrow-1)) axis(1,las=2)
		  		box(bty="o",lwd=2)
		  	}
		  	Legend( Bin, Col )
		    dev.off()
		  }

		# Spatial residuals (red is positive bias; blue is negative bias)
		  if(ModelSet[ModelI]!="Index"){
		  Ncol = ceiling(sqrt(length(unique(DF[,'Year_j'])))); Nrow = ceiling( length(unique(DF[,'Year_j']))/Ncol )
		  for(FigI in 1:2){
			if(FigI==1){
			  Mat1 = tapply( (Report[["W_j_hat"]]-Data[["W_j"]])^2/(Report[["W_j_hat"]]*(Data$CV_w+Report[["extraCV_w"]]))^2, INDEX=list(Data[['Station_j']], Data[['Year_j']]), FUN=mean, na.rm=TRUE)
			  Mat2 = tapply( (Report[["W_j_hat"]]-Data[["W_j"]]), INDEX=list(Data[['Station_j']], Data[['Year_j']]), FUN=mean, na.rm=TRUE)
			}
			if(FigI==2){
			  Mat1 = tapply( (Report[["I_j_hat"]]-Data[["I_j"]])^2/Report[["I_j_var"]], INDEX=list(Data[['Station_j']], Data[['Year_j']]), FUN=mean)
			  Mat2 = tapply( (Report[["I_j_hat"]]-Data[["I_j"]]), INDEX=list(Data[['Station_j']], Data[['Year_j']]), FUN=mean)
			}
			# lat/lon Pearson residuals
			png( file=paste(ModelFile,"PearsonResid_",c("W","I")[FigI],"_LL_stations.png",sep=""), width=4*Ncol, height=2*Nrow, res=200, units="in")
			  par(mfrow=c(Nrow,Ncol), mar=c(0,0,0,0), oma=c(2,2,0,0) )
			  for(i in 1:length(unique(DF[,'Year_j']))){
				#map("worldHires", ylim=y.lim,xlim=x.lim,col="grey90",fill=T, main="", mar=c(0,0,2.5,0))
				plot(1, type="n", ylim=range(loc_stations[,'long']), xlim=range(loc_stations[,'lat']))
				title( unique(DF[,'Year_j'])[i] )
				points(y=loc_stations[,'lat'], x=loc_stations[,'long'], col=c("blue","red")[ifelse(Mat2[,i]>0,2,1)], cex=sqrt(Mat1[,i]), pch=21)
				points(y=loc_stations[,'lat'], x=loc_stations[,'long'], col="black", cex=sqrt(0.5), pch=20)
				if( (i-1)%%Ncol == 0) axis(2)
				if( (i-1)/Ncol >= (Nrow-1)) axis(1,las=2)
				box(bty="o",lwd=2)
			  }
			dev.off()
			# individual Pearson residuals
			RangeFn = function(Vec) range(ifelse(abs(Vec)==Inf,NA,Vec),na.rm=TRUE) 
			png( file=paste(ModelFile,"PearsonResid_",c("W","I")[FigI],".png",sep=""), width=3*Ncol, height=3*Nrow, res=200, units="in")
			  par(mfrow=c(Nrow,Ncol), mar=c(3,3,2,0) )
			  for(i in 1:length(unique(DF[,'Year_j']))){
				Which = which(Data[["Year_j"]]==unique(Data[["Year_j"]])[i])
				if(FigI==1){
				  Vec1 = ( (Report[["W_j_hat"]]-Data[["W_j"]])^2/(Report[["W_j_hat"]]*(Data$CV_w+Report[["extraCV_w"]]))^2 )
				  Vec2 = ( (Report[["W_j_hat"]]-Data[["W_j"]]) )
				}
				if(FigI==2){
				  Vec1 = ( (Report[["I_j_hat"]]-Data[["I_j"]])^2/Report[["I_j_var"]] )
				  Vec2 = ( (Report[["I_j_hat"]]-Data[["I_j"]]) )
				}
				plot(y=ifelse(Vec2[Which]>0,1,-1)*Vec1[Which], x=Data[["Station_j"]][Which], pch=21, main=unique(Data[["Year_j"]])[i], ylim=RangeFn(ifelse(Vec2>0,1,-1)*Vec1))
			  }
			dev.off()
			# individual Pearson residuals on map
			RangeFn = function(Vec) range(ifelse(abs(Vec)==Inf,NA,Vec),na.rm=TRUE) 
			png( file=paste(ModelFile,"PearsonResid_",c("W","I")[FigI],"_LL.png",sep=""), width=4*Ncol, height=2*Nrow, res=200, units="in")
			  par(mfrow=c(Nrow,Ncol), mar=c(3,3,2,0) )
			  for(i in 1:length(unique(DF[,'Year_j']))){
				Which = which(Data[["Year_j"]]==unique(Data[["Year_j"]])[i])
				if(FigI==1){
				  Vec1 = ( (Report[["W_j_hat"]]-Data[["W_j"]])^2/(Report[["W_j_hat"]]*(Data$CV_w+Report[["extraCV_w"]]))^2 )
				  Vec2 = ( (Report[["W_j_hat"]]-Data[["W_j"]]) )
				}
				if(FigI==2){
				  Vec1 = ( (Report[["I_j_hat"]]-Data[["I_j"]])^2/Report[["I_j_var"]] )
				  Vec2 = ( (Report[["I_j_hat"]]-Data[["I_j"]]) )
				}
				#map("worldHires", ylim=y.lim,xlim=x.lim,col="grey90",fill=T, main="", mar=c(0,0,2.5,0))
				plot(1, type="n", ylim=range(loc_stations[,'long']), xlim=range(loc_stations[,'lat']))
				title( unique(DF[,'Year_j'])[i] )
				points(y=DF[Which,'lat_j'], x=DF[Which,'long_j'], col=c("blue","red")[ifelse(Vec2[Which]>0,2,1)], cex=sqrt(Vec1[Which]), pch=21)
				points(y=DF[Which,'lat_j'], x=DF[Which,'long_j'],, col="black", cex=sqrt(0.5), pch=20)
				if( (i-1)%%Ncol == 0) axis(2)
				if( (i-1)/Ncol >= (Nrow-1)) axis(1,las=2)
				box(bty="o",lwd=2)
				#plot(y=ifelse(Vec2[Which]>0,1,-1)*Vec1[Which], x=Data[["Station_j"]][Which], pch=21, main=unique(Data[["Year_j"]])[i], ylim=RangeFn(ifelse(Vec2>0,1,-1)*Vec1))
			  }
			dev.off()
		    } # End FigI loop
		  } # End if(ModelI!="Index")
    }} # End ModelI, RepI loops

}




# Base case
	Run_simul(SD_A=0.5, SD_E=0.5, Scale=250, CV_w=0.2, Accel=0.2, SD_F=0.2)
# Test SD_A
	Run_simul(SD_A=0.25, SD_E=0.5, Scale=250, CV_w=0.2, Accel=0.2, SD_F=0.2)
	Run_simul(SD_A=1, SD_E=0.5, Scale=250, CV_w=0.2, Accel=0.2, SD_F=0.2)
# Test SD_E
	Run_simul(SD_A=0.5, SD_E=0.25, Scale=250, CV_w=0.2, Accel=0.2, SD_F=0.2)
	Run_simul(SD_A=0.5, SD_E=1, Scale=250, CV_w=0.2, Accel=0.2, SD_F=0.2)
	





################################
#
# Explore simulation results
#
################################

if(FALSE){

  load(file=paste(Dir_save,"SettingsList.RData",sep=""))
  attach( SettingsList )
  # n_s = 25

  # Save objects
  Timeseries = array(NA, dim=c(length(RepSet),n_t,5,length(ModelSet)+1), dimnames=list( paste("Rep",RepSet), paste("Year",1:n_t), c("N_t","S_t","R_t","F_t","D_t"), c("True",ModelSet)) )
  Omega = array(NA, dim=c(length(RepSet),n_s,length(ModelSet)+1), dimnames=list( paste("Rep",RepSet), paste("Omega",1:n_s), c("True",ModelSet)) )
  Epsilon = array(NA, dim=c(length(RepSet),n_t,n_s,length(ModelSet)+1), dimnames=list( paste("Rep",RepSet), paste("Year",1:n_t), paste("Omega",1:n_s), c("True",ModelSet)) )
  Param = array(NA, dim=c(length(RepSet),4,length(ModelSet)+1), dimnames=list( paste("Rep",RepSet), c("R2_omega","Cor_omega","R2_epsilon","Cor_epsilon"), c("True",ModelSet)) )

  RepI=1; ModelI=3
  for(RepI in 1:length(RepSet)){
  for(ModelI in 1:length(ModelSet)){
  
    # Unload previous
    if( 'TrueList' %in% search() ) detach(TrueList)
    if( 'Save' %in% search() ) detach(Save)

    # RepFile
    RepFile = paste(Dir_save,"Rep=",RepSet[RepI],"/",sep="")
    ModelFile = paste(RepFile,"Model=",ModelSet[ModelI],"/",sep="")

    # Load true values
    load( file=paste(RepFile,"TrueList.RData",sep=""))
    Match = which(names(TrueList) %in% ls() )
    if( length(Match)>=1) remove( list=names(TrueList)[Match])
    attach(TrueList)
    # Load estimated values
    load( file=paste(ModelFile,"Save.RData",sep=""))
    Match = which(names(Save) %in% ls() )
    if( length(Match)>=1) remove( list=names(Save)[Match])
    attach(Save)
                        
    #### Time series
    # Save true values
    if(ModelI==1){
      Timeseries[RepI,1:n_t,'N_t','True'] = N_t 
      Timeseries[RepI,1:n_t,'S_t','True'] = S_t
      Timeseries[RepI,1:n_t,'R_t','True'] = R_t
      Timeseries[RepI,1:n_t,'F_t','True'] = F_t
      Timeseries[RepI,1:n_t,'D_t','True'] = S_t/sum_S0
    }
    # Save estimates
    Timeseries[RepI,1:n_t,'N_t',ModelSet[ModelI]] = Report$N_t_hat
    Timeseries[RepI,1:n_t,'S_t',ModelSet[ModelI]] = Report$S_t_hat
    Timeseries[RepI,1:n_t,'R_t',ModelSet[ModelI]] = Report$R_t_hat
    Timeseries[RepI,1:n_t,'F_t',ModelSet[ModelI]] = Report$F_t
    Timeseries[RepI,1:n_t,'D_t',ModelSet[ModelI]] = Report$S_t/Report$sum_S0
    
    #### Random field values
    # true values
    # Save true values
    if(ModelI==1){
      Omega[RepI,,'True'] = Omega_s[1:n_s] 
      Epsilon[RepI,1:n_t,,'True'] = t(Epsilon_s[1:n_s,])
    }
    # Save estimates
    Omega[RepI,,ModelSet[ModelI]] = Report$Omega[1:n_s]
    Epsilon[RepI,1:n_t,,ModelSet[ModelI]] = t(Report$Epsilon[1:n_s,])
    
    #### Parameters and derived values
    Param[RepI,"R2_omega",ModelSet[ModelI]] = 1 - sum( (Omega[RepI,,'True']-Omega[RepI,,ModelSet[ModelI]])^2 )/sum( Omega[RepI,,'True']^2 )
    Param[RepI,"R2_epsilon",ModelSet[ModelI]] = 1 - base::sum( (Epsilon[RepI,1:n_t,,'True']-Epsilon[RepI,1:n_t,,ModelSet[ModelI]])^2 )/base::sum( Epsilon[RepI,1:n_t,,'True']^2 )
    Param[RepI,"Cor_omega",ModelSet[ModelI]] = cor( Omega[RepI,,'True'], Omega[RepI,,ModelSet[ModelI]] ) 
    Param[RepI,"Cor_epsilon",ModelSet[ModelI]] = cor( as.vector(Epsilon[RepI,1:n_t,,'True']), as.vector(Epsilon[RepI,1:n_t,,ModelSet[ModelI]]) )
  }} # End ModelI, RepI loops
  
  RE = (Timeseries[,,,ModelSet] - outer(Timeseries[,,,'True'],rep(1,3))) / outer(Timeseries[,,,'True'],rep(1,3))
  Bias = apply( RE[,,,,drop=FALSE], MARGIN=3:4, FUN=median, na.rm=TRUE)
  MARE = apply( abs(RE[,,,,drop=FALSE]), MARGIN=3:4, FUN=median, na.rm=TRUE)
  RMSE = sqrt(apply( RE[,,,,drop=FALSE]^2, MARGIN=3:4, FUN=median, na.rm=TRUE))

  par(mfrow=c(dim(RE)[3],dim(RE)[4]), mar=c(3,3,2,0), mgp=c(2,0.5,0))
  for(StatI in 1:dim(RE)[3]){
  for(ModelI in 1:dim(RE)[4]){
    matplot( t(RE[,,StatI,ModelI]), type="l", col="black", lty="solid", ylim=range(RE[,,StatI,],na.rm=TRUE))
  }}  

  # Correlations between Omega and Epsilon                                    # col=rgb(1,0,0,0.2), 
  png(file=paste(Dir_save,"Random_field_summaries.png",sep=""), width=5, height=5, res=200, units="in")
    par(mfrow=c(2,2), mar=c(0,2,0,0), mgp=c(2,0.5,0), tck=-0.02, oma=c(3,3,2,0.5), xaxs="i", yaxs="i")
    hist(Param[,"R2_omega","Spatial"], breaks=seq(0,1,length=25), xlim=c(0,1), prob=TRUE, main="", xlab="", ylab="", xaxt="n", col="darkgrey" ); box()
    mtext( side=3, text="Variance explained")
    mtext( side=2, text=expression(Omega), line=2)
    hist(Param[,"Cor_omega","Spatial"], breaks=seq(0,1,length=25), xlim=c(0,1), prob=TRUE, main="", xlab="", ylab="", xaxt="n", col="darkgrey" ); box()
    mtext( side=3, text="Correlation")
    hist(Param[,"R2_epsilon","Spatial"], breaks=seq(0,1,length=25), xlim=c(0,1), prob=TRUE, main="", xlab="", ylab="", xaxt="n", col="darkgrey" ); box()
    mtext( side=2, text=expression(Epsilon), line=2)
    axis(1)
    hist(Param[,"Cor_epsilon","Spatial"], breaks=seq(0,1,length=25), xlim=c(0,1), prob=TRUE, main="", xlab="", ylab="", xaxt="n", col="darkgrey" ); box()
    axis(1)
  dev.off()
}
       