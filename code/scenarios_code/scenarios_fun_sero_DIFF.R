

vars.plus.Vx.eff <-  c("S",
                       "I_detect_new",
                       "I",
                       "Idetectcum",
                       "Itot",
                       "Itotcum",
                       "H_new",
                       "Htot",
                       "Htotcum",
                       "Q",
                       "Qcum",
                       "V",
                       "Vcum",
                       "D_new",
                       "D",
                       "R_out",
                       #"Vx",
                       #"Vx_new",
                       "R_vx",
                       "R0_t",
                       "Vx.given",
                       "Vx.given.tot",
                       "Vx.given.65",
                       "Vx.given.65.tot",
                       "sero_vx"
)

######################################################################################################
## Logic check to keep the number vaccinated less than the coverage ratio and or population of LAC
######################################################################################################

nvx_logic_check <- function(nvx_t,nvx_y,vx.delay,vx.coverage,delay.start.65.vx){
  
  ## NOTE: nvx_t.feas times based on input vx rate from fn_t_readin
  
  nvx_t <- nvx_t + vx.delay
  
  ## Here inputing logic check to keep the number vaccinated less than the coverage ratio and or population of LAC
  
  # Expand the nvx(t,y) function over each time step
  interpolate.out <- approx(x=nvx_t, y=nvx_y, method="constant", n = max(nvx_t)-min(nvx_t)+1)
  nvx_t.feas <- interpolate.out$x
  nvx_y.feas <- interpolate.out$y
  nvx_y.feas_cum <- cumsum(nvx_y.feas)
  
  ##################################################
  ## Implementing logic check: nvx_y.feas    
  for (time in 1:length(nvx_y.feas)){
    if (nvx_y.feas_cum[time] > vx.coverage*1e7-100000 ) {
      nvx_y.feas[time] <- 0
    }
  }
  
  ##################################################
  ## Implementing logic check: nvx_65_y.feas    
  #  start.vx <- round(as.numeric(as.Date(start.date.vx) - as.Date("2020-03-01")) + start_time) + vx.delay - 1  # This relates start time to the one in vx.scenarios.mat, not the one in fn_t_readin
  
  # This start date for vaccines overall is tied to the vax rate and date scenarios read-in via fn_t_readin
  if (sum(nvx_y.feas)==0){
    nvx_65_y.feas <- nvx_y.feas
  }
  
  if (sum(nvx_y.feas)!=0){
    
    start.vx <- min(which(nvx_y.feas!=0))
    
    delay.start.65.vx <- as.numeric(delay.start.65.vx)
    
    nvx_65_y.feas <- nvx_y.feas
    
    for (time in (1:(start.vx + delay.start.65.vx))){
      nvx_65_y.feas[time] <- 0
    }
    
    for (time in (1 + start.vx + delay.start.65.vx):length(nvx_65_y.feas)){
      
      # This is the offset factor: the total # vax given before 65+ start
      n.vx.when.65.start <- nvx_y.feas_cum[start.vx+delay.start.65.vx]
      
      # This next line says: When the # of vaccines given to 65+ (which is the total # of vaccines given offset by the # given in Phase 1a before 65+)
      #                       > the vax coverage rate for 65+
      #                       then stop vaccinating 65+      
      if ((nvx_y.feas_cum[time] - n.vx.when.65.start) > vx.coverage*1.2e6 ) {   # 1.2e6: Population of 65+ in LA
        nvx_65_y.feas[time] <- 0
      }
    }
  }
  
  # Matrix::nnzero(as.matrix(nvx_65_y.feas))
  
  return(cbind(nvx_t.feas,nvx_y.feas,nvx_65_y.feas))
}

#nvx_ty_out <- nvx_logic_check(nvx_t=nvx_t,nvx_y=nvx_y,vx.delay=vx.delay,vx.coverage=vx.coverage) 

# max_nvx <- max(cumsum(nvx_y.feas))


########################################################################################
## Get Alpha, Kappa, Delta with populations protected: for vaccination scenarios
########################################################################################

# weighted.avg.scenarios.0 <- weighted.avg.scenario.65(X.mat=X.mat, freq.PREV.q=freq.PREV.q, freq.LAC.obs.age=freq.LAC.obs.age, logit.SEIR.est=logit.SEIR.est, psi.mat=psi.mat, factor.to.remove=factor.to.remove,
#                                                      n.days.65.vx=1, vx.coverage.65=.1)

weighted.avg.scenario.65 <- function(X.mat, freq.PREV.q, freq.LAC.obs.age, logit.SEIR.est, psi.mat, factor.to.remove,
                                     n.days.65.vx, vx.coverage.65){
  if (n.days.65.vx==0){
    n.days.65.vx=1
  }
  vx.coverage.65 <- as.numeric(vx.coverage.65)
  n.days.65.vx <- as.numeric(n.days.65.vx)
  percent.to.remove <- round(seq(0,vx.coverage.65, length = n.days.65.vx),3)
  n.protect.scenario <- length(percent.to.remove)
  factor.to.remove <- 4  # Age.65.
  
  weighted.avg.scenarios.overall <- 
    weighted.avg.protect.SCENARIOS(X.mat=X.mat, freq.PREV.q = freq.PREV.q, 
                                   freq.LAC.obs.age=freq.LAC.obs.age, logit.SEIR.est=logit.SEIR.est,
                                   psi.mat=psi.mat, percent.to.remove=percent.to.remove,
                                   factor.to.remove=factor.to.remove)
  return(weighted.avg.scenarios.overall)
}

########################################################################################
## Get Alpha, Kappa, Delta with populations protected
########################################################################################

weighted.avg.protect.SCENARIOS <- function(X.mat, freq.PREV.q, freq.LAC.obs.age, logit.SEIR.est, psi.mat, percent.to.remove, factor.to.remove){
  
  n.dates <- ncol(freq.LAC.obs.age)
  
  #percent.to.remove <- c(1,0.5) # Remove 100%
  #factor.to.remove <- 4  # Age.65.
  weighted.avg <- vector("list", n.dates)
  weighted.avg.scenarios <- vector("list",length(percent.to.remove))
  weighted.avg.scenarios.overall <- NULL
  
  for (remove in 1:length(percent.to.remove)){
    for (t in 1:n.dates){
      percent.remove <- percent.to.remove[remove]
      time = t
      name.date <- colnames(freq.LAC.obs.age)[t]
      
      freq.I <- freq.ILL(X.mat, freq.PREV.q, freq.LAC.obs.age, time = time)
      
      Pr.H.I.q <- get.Pr.q(model=1, time=time, logit.SEIR.est, X.mat, freq.q.IN = freq.I[[2]],  psi.mat)
      freq.H <- get.population.freq(prob.q.vec.IN = Pr.H.I.q, freq.q.IN = freq.I[[1]], X.mat)
      
      Pr.Q.H.q <- get.Pr.q(model=2, time=time, logit.SEIR.est, X.mat, freq.q.IN = freq.H[[2]],  psi.mat)
      freq.Q <- get.population.freq(prob.q.vec.IN = Pr.Q.H.q, freq.q.IN = freq.H[[1]], X.mat)
      
      Pr.D.Q.q <- get.Pr.q(model=3, time=time, logit.SEIR.est, X.mat, freq.q.IN = freq.Q[[2]],  psi.mat)
      freq.D <- get.population.freq(prob.q.vec.IN = Pr.D.Q.q, freq.q.IN = freq.Q[[1]], X.mat)
      
      ##
      freq.I.filter <- freq.I[[1]]*(1-percent.remove*X.mat[,factor.to.remove])
      weighted.avg.H.filter <- t(Pr.H.I.q) %*% freq.I.filter
      
      freq.H.filter <- (get.population.freq(Pr.H.I.q, freq.I.filter, X.mat))[[1]]
      weighted.avg.Q.filter <- t(Pr.Q.H.q) %*% freq.H.filter
      
      freq.Q.filter <- (get.population.freq(Pr.Q.H.q, freq.H.filter, X.mat))[[1]]
      weighted.avg.D.filter <- t(Pr.D.Q.q) %*% freq.Q.filter
      
      weighted.avg.t <- cbind( weighted.avg.H.filter, weighted.avg.Q.filter, weighted.avg.D.filter )
      colnames(weighted.avg.t) <- c("Alpha.","Kappa.","Delta.")
      colnames(weighted.avg.t) <- paste0( colnames(weighted.avg.t), name.date)
      weighted.avg[[t]] <- weighted.avg.t
      ##
    }
    weighted.avg.scenarios[[remove]] <- do.call(cbind, weighted.avg)
  }
  weighted.avg.scenarios.overall <- do.call(rbind,weighted.avg.scenarios)
  #weighted.avg.scenarios.overall <- as.data.frame(weighted.avg.scenarios)
  rownames(weighted.avg.scenarios.overall) <- paste0("Protect.",percent.to.remove*100)
  
  return(weighted.avg.scenarios.overall)
  
}


########################################################################################
## FUNCTION TO GET CI FOR EACH SCENARIO (embedded in correlated.param.SCENARIOS)
########################################################################################



get.CI.SCENARIOS <- function(TEST.out) {
  
  library(data.table)
  init.date.data="2020-03-01"
  
  ### Add CFR and IFR to the list (EXTRA STEP NOW THAT THIS IS BEING USED ALSO FOR summary_table)
  traj <- dplyr::mutate(TEST.out, Itot=I+A, CFRobs=(D/Idetectcum), CFRactual=(D/(Itotcum)) )
  traj <-  dplyr::select(traj,c(1:9,CFRobs,CFRactual,vars.plus.Vx.eff))
  ###
  
  ## TO SAVE MEMORY
  rm(TEST.out)
  
  print("Starting CI calc")
  
  ### MELTING AND APPLYING SUMMARY STAT FUNCTIONS
  df.traj <- reshape2::melt(traj, measure.vars = c(10:ncol(traj)), variable.name = "state.name")
  df.traj_dt <- as.data.table(df.traj)
  
  traj.CI <- df.traj_dt[, list(
    N=.N,
    mean = mean(value),
    median = quantile(value, c(.5),na.rm=TRUE),
    low_95 = quantile(value, c(.025),na.rm=TRUE),
    up_95 = quantile(value, c(.975),na.rm=TRUE),
    up_50 = quantile(value,.75,na.rm=TRUE),
    low_50 = quantile(value,.25,na.rm=TRUE)),
    #by = c("date","scenario.id",  "protect.id", "NPI.id")]#,"vx.id", "vx.coverage.id", "vx.start.date.id", "protect.id","NPI.id")]
    by = c("date", "state.name","scenario.id","vx.id", "vx.start.date.id", "delay.start.65.vx.id", "sero.id", "R0.id")]
  traj.CI <- as.data.frame(traj.CI)
  
  ## TO ALIGN DATES: MODEL
  init.date = init.date.data #"2020-01-03"
  init.date <- as.Date(init.date)
  traj.CI[["date"]] <- traj.CI[["date"]] + init.date
  
  return(traj.CI)
  
}


########################################################################################
## SPECIFYING EPIDEMIC MODEL TO BE SIMULATED LOOPING OVER SCENARIOS
## NOTE: NPI scenarios are called in from data "fn_readin.csv"
## NOTE: Protect scenarios are called in from weighted.avg.protect.mat
########################################################################################

correlated.param.SCENARIOS.diff <- function(fn_t_readin_input=fn_t_readin_input, ABC.out.mat, iter, time.steps, vx.scenarios.mat, sero.scenario.in, mu.scenario.in,
                                       X.mat=X.mat, freq.PREV.q=freq.PREV.q, freq.LAC.obs.age=freq.LAC.obs.age, logit.SEIR.est=logit.SEIR.est, psi.mat=psi.mat,
                                       vx.delay=vx.delay){
  ## Number of scenarios
  
  n.scenarios <- nrow(vx.scenarios.mat) * length(mu.scenario.in)
  
  ##########################################################################################
  ## Read in csv files with Beta_t, Alpha_t, Kappa_t, Delta_t
  
  start_time = 45
  
  fn_t_readin_path <- path(data.dir, fn_t_readin_input)
  fn_t_readin = as.data.frame(read.csv(fn_t_readin_path, sep=",",stringsAsFactors = FALSE))
  
  ## Get r_t
  r_t_dates <- as.Date(fn_t_readin$r_t)
  r_t_dates <- na.omit(r_t_dates)
  r_t_dates[1] <- r_t_dates[1]-start_time
  r_t <- round(as.numeric(r_t_dates - as.Date("2020-03-01")) + start_time)
  
  ## Get Beta_t
  Beta_t_dates <- as.Date(fn_t_readin$Beta_t)
  Beta_t_dates[1] <- Beta_t_dates[1]-start_time
  Beta_t <- round(as.numeric(Beta_t_dates - as.Date("2020-03-01")) + start_time)
  
  ## Get Alpha_t
  alpha_t_readin_path <- path(data.dir, "alpha_t_readin.csv")
  alpha_t_readin = as.data.frame(read.csv(alpha_t_readin_path, sep=",",stringsAsFactors = FALSE))
  
  Alpha_t_dates <- as.Date(alpha_t_readin$Alpha_t)
  Alpha_t_dates[1] <- Alpha_t_dates[1]-start_time
  Alpha_t_dates <- Alpha_t_dates[1:length(Alpha_t_dates)-1]
  Alpha_t <- round(as.numeric(Alpha_t_dates - as.Date("2020-03-01")) + start_time)
  
  ## Get Alpha Kappa Delta until last estimated date (i.e. before forward scenario starts)
  # This "last date" can be modified in the alpha_t_readin file
  # The first value in the weighted.avg.scenarios[1,""] will always be the remove.0 scenario so is not scenario-specific
  weighted.avg.scenarios.0 <- weighted.avg.scenario.65(X.mat=X.mat, freq.PREV.q=freq.PREV.q, freq.LAC.obs.age=freq.LAC.obs.age, logit.SEIR.est=logit.SEIR.est, psi.mat=psi.mat, factor.to.remove=factor.to.remove,
                                                       n.days.65.vx=1, vx.coverage.65=.1)
  # Alpha1 <- weighted.avg.scenarios.0[1,"Alpha.t1"]
  # Alpha2 <- weighted.avg.scenarios.0[1,"Alpha.t2"]
  # Kappa1 <- weighted.avg.scenarios.0[1,"Kappa.t1"]
  # Kappa2 <- weighted.avg.scenarios.0[1,"Kappa.t2"]
  # Delta1 <- weighted.avg.scenarios.0[1,"Delta.t1"]
  # Delta2 <- weighted.avg.scenarios.0[1,"Delta.t2"]
  
  Alpha1 <- 0.14#.155
  Kappa1 <- 0.6 #.65
  Delta1 <- .56
  
  Alpha2 <- 0.05 #.06
  Kappa2 <- .55 #.25
  Delta2 <- .52 #.75
  
  ## Alpha_t
  
  Alpha_y_chr <- alpha_t_readin$Alpha_y
  assign("Alpha1",Alpha1)
  assign("Alpha2", Alpha2)
  
  Alpha_y <- as.vector(length(Alpha_t))
  for (z in 1:length(Alpha_t)){
    Alpha_y[z] = get(Alpha_y_chr[z])
  }
  
  ## Kappa_t
  
  Kappa_y_chr <- alpha_t_readin$Kappa_y
  assign("Kappa1",Kappa1)
  assign("Kappa2", Kappa2)
  
  Kappa_y <- as.vector(length(Alpha_t))
  for (z in 1:length(Alpha_t)){
    Kappa_y[z] = get(Kappa_y_chr[z])
  }
  
  ## Delta_t
  
  Delta_y_chr <- alpha_t_readin$Delta_y
  assign("Delta1",Delta1)
  assign("Delta2", Delta2)
  
  Delta_y <- as.vector(length(Alpha_t))
  for (z in 1:length(Alpha_t)){
    Delta_y[z] = get(Delta_y_chr[z])
  }
  
  ## vx(t)
  nvx_t_dates <- as.Date(fn_t_readin$nvx_t)
  nvx_t_dates <- na.omit(nvx_t_dates)
  nvx_t_dates[1] <- nvx_t_dates[1]-start_time-vx.delay
  nvx_t <- round(as.numeric(nvx_t_dates - as.Date("2020-03-01")) + start_time)
  
  start.date.vx <- nvx_t_dates[3]
  
  #nvx_t <- nvx_t + vx.delay
  
  # print(nvx_t_dates)
  # print(nvx_t)
  
  ##########################################################################################
  ## INITIALIZE SCENARIOS.out
  
  SCENARIOS.out <- vector("list", n.scenarios)
  SCENARIOS.out.diff <- vector("list", n.scenarios)
  
  scenario.idx <- 1
  
  for (sero.idx in 1:length(sero.scenario.in)){
    
    frac.sero.vax = sero.scenario.in[sero.idx]
    
    print(paste0("frac.sero.vax", frac.sero.vax))
    
    for (NPI.idx in 1:length(mu.scenario.in)){
      
      ########################################################################
      ## NPI SCENARIOS Beta_y
      
      NPI.scenario = NPI.idx
      mu.scenario <- mu.scenario.in[NPI.idx]
      #NPI.name <- NPI.scenarios.mat[NPI.idx, "NPI.name"]
      #Beta_y_vals <- NPI.scenarios.mat[NPI.idx, "Beta_y_vals"]
      
      for (protect.idx in 1:nrow(vx.scenarios.mat)){
        
        print(paste0("Starting scenario ", scenario.idx))
        
        ########################################################################
        ## VACCINATION OVERALL SCENARIOS nvx_y
        
        protect.scenario = protect.idx
        name.vx <- vx.scenarios.mat[protect.idx, "name.vx"]
        #start.date.vx <- as.Date(vx.scenarios.mat[protect.idx, "start.date.vx"])
        vx.coverage <- as.numeric(vx.scenarios.mat[protect.idx, "vx.coverage"])
        delay.start.65.vx <- as.numeric(vx.scenarios.mat[protect.idx, "delay.start.65.vx"])
        
        nvx_y_vals <- as.character(vx.scenarios.mat[protect.idx,"vx_y.vals"])
        nvx_y <- as.numeric(na.omit(fn_t_readin[,nvx_y_vals]))
        
        ## Put in here to get dates vaccination becomes effective and the date range during which vaccination can logically be happening (based on population size)
        nvx_ty_out <- nvx_logic_check(nvx_t=nvx_t,nvx_y=nvx_y,vx.delay=vx.delay,vx.coverage=vx.coverage, delay.start.65.vx = delay.start.65.vx) 
        
        nvx_t.feas <- nvx_ty_out[,1]
        nvx_y.feas <- nvx_ty_out[,2]
        
        ########################################################################
        ## VACCINATE 65+ (LINKED TO VACCINATION OVERALL) SCENARIOS
        ## Specified in vx.scenarios.mat
        
        # Coverage for 65+ is same as for overall population (can modify this)
        vx.coverage.65 <- vx.coverage
        
        # Get dates during which intervention is run for 65+ and accounting for vx.delay
        nvx_65_y.feas <- nvx_ty_out[,3]
        nvx_65_y.feas.cum <- cumsum(nvx_65_y.feas)
        
        # Get number of days 65+ are vaccinated under this scenario timeline
        n.days.65.vx <- Matrix::nnzero(as.matrix(nvx_65_y.feas))
        
        ## Adjust the number of days it takes to decrease alpha, kappa, delta by removing 65+
        ## By making it more efficient (faster) if only sero- are vaccinated
        no.65.total <- 1.2e6
        no.recovered.total <- 5e6 #4e6
        frac.inf.65 <- .11
        
        # First find the percent of the overall recovered population that is 65+, based on their % in the infected population
        percent.immune.65 <- frac.inf.65 * no.recovered.total / no.65.total
        
        # Second find the factor reduction in the number of days it takes to decrease alpha kappa delta to their min values
        days.65.vx.redux.by.sero <- 1 - (1-frac.sero.vax)*percent.immune.65
        
        # Third find the reduction in the number of days it takes to decrease alpha kappa delta to their min values
        n.days.65.vx.sero <- round(n.days.65.vx * days.65.vx.redux.by.sero)
        
        #if (n.days.65.vx!=0) {n.days.65.vx.sero=7}
        
        print(paste0("n.days.65.vx = ", n.days.65.vx))
        print(paste0("n.days.65.vx = ", n.days.65.vx.sero))
        
        ## Get the % 65+ removed by day according to n.days.65.vx
        weighted.avg.scenarios <- weighted.avg.scenario.65(X.mat=X.mat, freq.PREV.q=freq.PREV.q, freq.LAC.obs.age=freq.LAC.obs.age, logit.SEIR.est=logit.SEIR.est, psi.mat=psi.mat, factor.to.remove=factor.to.remove,
                                                           n.days.65.vx=n.days.65.vx.sero, vx.coverage.65=vx.coverage.65)
        
        weighted.avg.scenarios.MIN <- weighted.avg.scenarios
        for (prob.idx in 1:length(weighted.avg.scenarios[,1])){
          weighted.avg.scenarios.MIN[prob.idx,4] <- min(weighted.avg.scenarios[prob.idx,4], Alpha2)
          weighted.avg.scenarios.MIN[prob.idx,5] <- min(weighted.avg.scenarios[prob.idx,5], Kappa2)
          weighted.avg.scenarios.MIN[prob.idx,6] <- min(weighted.avg.scenarios[prob.idx,6], Delta2)
        }
        weighted.avg.scenarios <- weighted.avg.scenarios.MIN
        
        ########################################################################
        ## Alpha Kappa Delta reduction given the 65+ population vaccinated and removed
        
        # This start date is computed in nvx_logic_check, and is based on fn_t_readin + vax efficacy delay + start 65+ vax delay
        start.time.effective.65 <- min(which(nvx_65_y.feas!=0))
        
        if (n.days.65.vx==0){
          Alpha_t_add <- start.time.effective.65
        }
        if (n.days.65.vx!=0){
          Alpha_t_add <- c(start.time.effective.65 : (start.time.effective.65+n.days.65.vx.sero-1))
        }
        #Alpha_t_add <- c(start.time.effective.65 : (start.time.effective.65+n.days.65.vx-1))
        Alpha_t.scenario <- c(Alpha_t, Alpha_t_add, 900)
        
        # Get values to add to Alpha_y
        Alpha_y.add <- weighted.avg.scenarios[,"Alpha.t2"]
        Alpha_y.scenario <- c(Alpha_y, Alpha_y.add, last(Alpha_y.add))
        
        Kappa_y.add <- weighted.avg.scenarios[,"Kappa.t2"]
        Kappa_y.scenario <- c(Kappa_y, Kappa_y.add, last(Kappa_y.add))
        
        Delta_y.add <- weighted.avg.scenarios[,"Delta.t2"]
        Delta_y.scenario <- c(Delta_y, Delta_y.add, last(Delta_y.add))
        
        ########################################################################
        ## Assign parameter combinations from the fitted model to iterate over, adding to the overall scenarios
        
        TEST.out <- vector("list", nrow(ABC.out.mat))
        
        for (idx in 1:nrow(ABC.out.mat)){
          
          R0 <- ABC.out.mat[idx,1]
          r1 <- ABC.out.mat[idx,2]
          start_time <- round(ABC.out.mat[idx,3])
          R0_redux1 <- ABC.out.mat[idx,4]
          # Delta1 <- ABC.out.mat[idx,5]
          # Alpha1 <- ABC.out.mat[idx,6]
          # Kappa1 <- ABC.out.mat[idx,7]
          p_V <- ABC.out.mat[idx,8]
          R0_redux2 <- ABC.out.mat[idx,9]
          # Delta2 <- ABC.out.mat[idx,10]
          # Alpha2 <- ABC.out.mat[idx,11]
          # Kappa2 <- ABC.out.mat[idx,12]
          r2 <- ABC.out.mat[idx,13]
          
          ##############
          ## r_y
          r_y_chr <- fn_t_readin$r_y
          assign("r1",r1)
          assign("r2", r1)
          
          r_y <- as.vector(length(r_t))
          for (z in 1:length(r_t)){
            r_y[z] = get(r_y_chr[z])
          }
          
          ##############
          ## Beta_y
          
          #Beta_y_vals <- as.character(NPI.scenarios.mat[NPI.scenario,"Beta_y_vals"])
          
          mu_y_chr <- fn_t_readin[,"Beta_y"]
          assign("mu.0",1)
          assign("mu.1", R0_redux1)
          assign("mu.2", R0_redux2)
          assign("mu.3",R0_redux3)
          assign("mu.4",mu.scenario)
          
          mu_y <- as.vector(length(Beta_t))
          for (z in 1:length(Beta_t)){
            mu_y[z] = get(mu_y_chr[z])
          }
          R0_y <- R0*mu_y
          
          R0_scenario <- round(last(R0_y),2)
          
          ## Get Beta_y as a function of R0, R0_redux, r, and Alpha
          
          Br.function <- function(R0.in, r.in, Alpha.in){
            d_IH <- 10   #days between illness onset and hospitalization
            d_IR <- 7    #days between illness onset and recovery (hospitalization not required)
            Br <- R0.in * ( 1 / ( (r.in/ ((Alpha.in/d_IH) + ((1-Alpha.in)/d_IR)))  + (1-r.in)*d_IR ))
            return(Br)
          }
          
          Beta_y <- as.vector(length(Beta_t))
          Beta_y<- c(
            Br.function(R0.in<-R0_y[1], r.in<-r1, Alpha.in<-Alpha1) ,
            Br.function(R0.in<-R0_y[2], r.in<-r1, Alpha.in<-Alpha1),
            Br.function(R0.in<-R0_y[3], r.in<-r1, Alpha.in<-Alpha1),
            Br.function(R0.in<-R0_y[4], r.in<-r1, Alpha.in<-Alpha1),
            Br.function(R0.in<-R0_y[5], r.in<-r1, Alpha.in<-Alpha1),
            Br.function(R0.in<-R0_y[6], r.in<-r2, Alpha.in<-Alpha2),
            Br.function(R0.in<-R0_y[7], r.in<-r2, Alpha.in<-Alpha2),
            Br.function(R0.in<-R0_y[8], r.in<-r2, Alpha.in<-Alpha2),
            Br.function(R0.in<-R0_y[9], r.in<-r2, Alpha.in<-Alpha2),
            Br.function(R0.in<-R0_y[10], r.in<-r2, Alpha.in<-Alpha2),
            Br.function(R0.in<-R0_y[11], r.in<-r2, Alpha.in<-Alpha2),
            Br.function(R0.in<-R0_y[12], r.in<-r2, Alpha.in<-Alpha2),
            Br.function(R0.in<-R0_y[13], r.in<-r2, Alpha.in<-Alpha2),
            Br.function(R0.in<-R0_y[14], r.in<-r2, Alpha.in<-Alpha2),
            Br.function(R0.in<-R0_y[15], r.in<-r2, Alpha.in<-Alpha2),
            Br.function(R0.in<-R0_y[16], r.in<-r2, Alpha.in<-Alpha2),
            Br.function(R0.in<-R0_y[17], r.in<-r2, Alpha.in<-Alpha2),
            Br.function(R0.in<-R0_y[18], r.in<-r2, Alpha.in<-Alpha2),
            Br.function(R0.in<-R0_y[19], r.in<-r2, Alpha.in<-Alpha2),
            Br.function(R0.in<-R0_y[20], r.in<-r2, Alpha.in<-Alpha2),
            Br.function(R0.in<-R0_y[21], r.in<-r2, Alpha.in<-Alpha2),
            Br.function(R0.in<-R0_y[22], r.in<-r2, Alpha.in<-Alpha2),
            Br.function(R0.in<-R0_y[23], r.in<-r2, Alpha.in<-Alpha2)
          )
          
          ##########################################################################################
          ## RUN THE MODEL
          
          ## COMPILE 
          x <- Vx_generator(Alpha_t=Alpha_t.scenario, Alpha_y=Alpha_y.scenario, Kappa_t=Alpha_t.scenario, Kappa_y=Kappa_y.scenario, Delta_t=Alpha_t.scenario, Delta_y=Delta_y.scenario, Beta_t=Beta_t, Beta_y=Beta_y, r_t=r_t, r_y=r_y, 
                            nvx_t=nvx_t.feas, nvx_y=nvx_y.feas, S_ini=1e7, E_ini=10, p_QV=p_V, x=frac.sero.vax)
          
          ## Prepare to add scenario inputs to output DF
          interpolate.out <- approx(x=Beta_t, y=R0_y, method="constant", n = max(Beta_t)-min(Beta_t)+1)
          R0_add <- as.data.frame(interpolate.out)
          colnames(R0_add) <- c("step","R0_t")
          
          nvx_t.given <- nvx_t.feas - vx.delay
          vx.given.df <- as.data.frame(cbind(nvx_t.given, nvx_y.feas, cumsum(nvx_y.feas), nvx_65_y.feas, nvx_65_y.feas.cum))
          colnames(vx.given.df) <- c("step", "Vx.given", "Vx.given.tot", "Vx.given.65", "Vx.given.65.tot")
          
          # Simulate a number of time steps accounting for model start and vaccine delay
          time.steps.simulate <- time.steps + 45
          
          ## SIMULATE
          TEST<-as.data.frame(plyr::rdply(iter, x$run(0:time.steps.simulate),.id="iter"))
          
          ## ADD R0 PLOT
          TEST <- merge(TEST,R0_add,by="step")
          TEST <- merge(TEST,vx.given.df,by="step")
          
          ## BIND INCLUDING OFFSETING OBSERVED DATA BY START DATE
          #TEST.out[[idx]] <- cbind(data.frame(par.id = idx, date = -start_time+TEST$step), TEST)
          
          # print(paste0("R(t) = ", mu.scenario*3.5))
          # print(paste0("R(t) = ", round(mu.scenario*3.5,1)))
          
          R0_save <- round(mu.scenario*3.5,1)
          
          ## BIND INCLUDING OFFSETING OBSERVED DATA BY START DATE
          TEST.out[[idx]] <- cbind(data.frame(scenario.id = paste0("target.", round(1-frac.sero.vax,1), "_" , name.vx ) , vx.id = name.vx, vx.start.date.id=start.date.vx, delay.start.65.vx.id = delay.start.65.vx, sero.id= round(1-frac.sero.vax,1), R0.id = R0_save, par.id = idx, date = -start_time+TEST$step), TEST)
          #TEST.out[[idx]] <- cbind(data.frame(scenario.id = paste0(scenario.idx, "_", name.vx, "_Rt.", R0_save ) , vx.id = name.vx, vx.start.date.id=start.date.vx, vx.coverage.id = vx.coverage, delay.start.65.vx.id = delay.start.65.vx, R0.id = R0_save, par.id = idx, date = -start_time+TEST$step), TEST)
          
        }  # end over idx
        
        ## Add to a dataframe over all idx
        TEST.out <- do.call(rbind, TEST.out)
        
        if (name.vx=="0"){
          TEST.out.0 <- TEST.out
          #colnames(TEST.out.0) <- paste0(colnames(TEST.out.0),".0")
          #colnames(TEST.out.0[,c(11:32)]) <- paste0(colnames(TEST.out.0[,c(11:32)]),".0")
        }
        
        n.cols.TEST.out <- ncol(TEST.out)
        TEST.out.diff <- abs ( TEST.out.0[,c(11:n.cols.TEST.out)] - TEST.out[,c(11:n.cols.TEST.out)] )
        TEST.out.diff.joined <- cbind(TEST.out[,c(1:10)], TEST.out.diff)
        
        ## Get CI for scenario
        SCENARIO.CI <- get.CI.SCENARIOS(TEST.out = TEST.out)
        SCENARIO.CI.diff <- get.CI.SCENARIOS(TEST.out = TEST.out.diff.joined)
        
        ## Put TEST.out dataframe into SCENARIOS.out
        # SCENARIOS.out[[scenario.idx]] <- TEST.out#SCENARIO.CI
        SCENARIOS.out[[scenario.idx]] <- SCENARIO.CI
        SCENARIOS.out.diff[[scenario.idx]] <- SCENARIO.CI.diff
        
        ##
        scenario.idx <- scenario.idx + 1 
        
        rm(TEST.out)
        rm(TEST.out.diff)
        
      } # end over y NPI.scenario
    } # end over x protect.scenario protect.idx
  } # end over sero.scenario.in
  
  # Output
  all.scenarios <- do.call(rbind,SCENARIOS.out)
  all.scenarios.diff <- do.call(rbind,SCENARIOS.out.diff)
  
  output <- vector("list", 2)
  output[[1]] <- all.scenarios
  output[[2]] <- all.scenarios.diff
  
  return(output)
}

# scenarios.out.test <- correlated.param.SCENARIOS(ABC.out.mat=ABC.out.mat[c(1:20),],iter=iter,time.steps=time.steps,protect.scenarios.mat=protect.scenarios.mat, NPI.scenarios.mat=NPI.scenarios.mat,
#                                        X.mat=X.mat, freq.PREV.q=freq.PREV.q, freq.LAC.obs.age=freq.LAC.obs.age, logit.SEIR.est=logit.SEIR.est, psi.mat=psi.mat,
#                                        vx.delay=vx.delay)


########################################################################################
## FUNCTION TO GET CI FOR EACH SCENARIO (embedded in correlated.param.SCENARIOS)
########################################################################################

plot.SCENARIOS <- function(traj.CI, data.in, startDatePlot, endDatePlot, vars.to.plot, filter.scenarios) {
  
  ## Filter only to variables of interest
  traj.CI <- traj.CI %>%  dplyr::filter(state.name %in% vars.to.plot) 
  
  if (!is.null(filter.scenarios)){
    traj.CI <- traj.CI %>% dplyr::filter(scenario.id %in% levels(traj.CI$scenario.id)[filter.scenarios])
  }
  
  ## Select only more recent dates
  init.date <- as.Date("2020-03-01")
  startDatePlot <- as.Date(startDatePlot) #init.date #- date.offset.4plot -1 #15
  endDatePlot <- as.Date(endDatePlot) #startDatePlot + time.steps.4plot #- 40  # the constant 40 because the traj are not aligned to start date
  traj.CI <- traj.CI %>% dplyr::filter(date >= startDatePlot) %>% dplyr::filter(date < endDatePlot)
  
  if(!is.null(data.in)){
    ## ALIGN DATES: DATA
    no_obs <- nrow(data.in)
    step <- 0:(no_obs-1)
    date <- init.date + step
    data.date <- cbind(date,data.in)
    rownames(data.date) <- step
    
    ## Select only more recent dates
    data.date <- data.date %>% dplyr::filter(date > startDatePlot) %>% dplyr::filter(date < endDatePlot)
    data <- reshape2::melt(data.date, measure.vars = c(2:ncol(data.date)), variable.name = "state.name")
  }
  
  ## PLOTTING
  traj.CI.line <- reshape2::melt(traj.CI[c("date","scenario.id",  "vx.id", "sero.id", "R0.id", "state.name", "median")], id.vars = c("date", "scenario.id",  "vx.id","sero.id", "R0.id","state.name"))
  traj.CI.area <- reshape2::melt(traj.CI[c("date","scenario.id",  "vx.id", "sero.id", "R0.id", "state.name", "low_95", "low_50", "up_50", "up_95")], id.vars = c("date", "scenario.id",  "vx.id", "sero.id", "R0.id","state.name"))
  traj.CI.area$type <- sapply(traj.CI.area$variable, function(x) {str_split(x, "_")[[1]][1]})
  traj.CI.area$CI <- sapply(traj.CI.area$variable, function(x) {str_split(x, "_")[[1]][2]})
  traj.CI.area$variable <- NULL
  traj.CI.area <- reshape2::dcast(traj.CI.area, "date+scenario.id+state.name+CI~type")
  
  p <- ggplot(transform(traj.CI.area, state.name = factor(state.name, levels=vars.to.plot)))
  
  #####################
  ### colors and names
  #####################
  
  longnames <- c("Suscept.",
                 "New Obs. Infected",
                 "Current Obs. Infected",
                 "Cum. Obs. Infected",
                 "Current Tot. Infected",
                 "Cum. Tot. Infected",
                 "New in Hospital",
                 "in Hospital",
                 "Cum. in Hospital",
                 "Current in ICU",
                 "Cum. in ICU",
                 "Current Ventilation",
                 "Cum. Ventilation",
                 "New Deaths",
                 "Cum. Deaths",
                 "Recovered",
                 #"Total Vax",
                 #"New Vax",
                 "Tot.Vx Recov.",
                 "CFR",
                 "IFR",
                 "R(t)",
                 "New Vx",
                 "Tot. Vax",
                 "New Vx 65+",
                 "Tot Vx 65+",
                 "Vx.Frac.Recov")
  
  names(longnames) <- c("S",
                        "I_detect_new",
                        "I",
                        "Idetectcum",
                        "Itot",
                        "Itotcum",
                        "H_new",
                        "Htot",
                        "Htotcum",
                        "Q",
                        "Qcum",
                        "V",
                        "Vcum",
                        "D_new",
                        "D",
                        "R_out",
                        #"Vx",
                        #"Vx_new",
                        "R_vx",
                        "CFRobs",
                        "CFRactual",
                        "R0_t",
                        "Vx.given",
                        "Vx.given.tot",
                        "Vx.given.65",
                        "Vx.given.65.tot",
                        "sero_vx"
  )
  
  ## Colors
  
  cols.list <- c(
    "salmon",
    "sandybrown",
    "navajowhite3",
    "olivedrab4",
    "olivedrab2",
    "mediumseagreen",
    "mediumaquamarine",
    "mediumturquoise",
    "cyan2",
    "lightskyblue",
    "steelblue2",
    "mediumpurple",
    "mediumorchid",
    "plum1",
    "violetred1",
    "deeppink3",
    "deeppink4",
    #"grey50",
    #"grey50",
    "grey50",
    "grey50",
    "grey50",
    "grey50",
    "grey50",
    "grey50",
    "grey50",
    "grey50"
  )
  
  names(cols.list) <- names(longnames)
  color.this.var <- as.character(cols.list[vars.to.plot])
  
  ## CAPACITY DATA FRAME
  capacity.vals <- as.data.frame(matrix(NA, nrow=length(levels(traj.CI$state.name)), ncol=2))
  capacity.vals[,1] <- levels(traj.CI$state.name)
  rownames(capacity.vals) <- levels(traj.CI$state.name)
  capacity.vals["Htot",2] <- 4000 
  capacity.vals["Q",2] <- 2245
  capacity.vals["V",2] <-1000
  colnames(capacity.vals) <- c("state.name","capacity")
  
  ## PLOT OPTIONS
  #  p <- p + facet_grid(state.name ~ scenario.id, labeller=labeller(state.name=longnames, scenario.id=longnames.scenarios), scales='free')
  p <- p + facet_grid(state.name ~ scenario.id, labeller=labeller(state.name=longnames), scales='free')
  p <- p + geom_ribbon(data = traj.CI.area, aes_string(x = "date", ymin = "low", ymax = "up", alpha = "CI", fill = "state.name"),show.legend = c(fill=FALSE))
  p <- p + geom_line(data = traj.CI.line, aes_string(x = "date", y = "value", linetype = "variable", colour = "state.name"), size = 1, show.legend = c(colour=FALSE))
  
  p <- p + scale_alpha_manual("Percentile", values = c("95" = 0.20, "50" = 0.50), labels = c("95" = "95th", "50" = "50th"))
  p <- p + scale_linetype("Stats")
  p <- p + guides(linetype = guide_legend(order = 1))
  
  
  ## ADD CAPACITY
  capacity.vals <- capacity.vals %>% filter(state.name %in% vars.to.plot)
  p <- p + geom_hline(data= capacity.vals, aes(yintercept=capacity),linetype = "dashed")
  
  ## ADD DATA
  if(!is.null(data.in)){
    p <- p + geom_point(data = data, aes_string(x = "date", y = "value"), size = .5, colour = "black")
  }
  
  p <- p + theme_bw() + theme(legend.position = "top", legend.box = "horizontal")
  p <- p + scale_x_date(limits = as.Date(c(startDatePlot,endDatePlot)), date_breaks = "1 month" , date_labels = "%d-%b-%y")
  p <- p + theme(axis.text.x = element_text(angle = 90),
                 strip.text.x = element_text(size = 8, face = "bold"))
  p <- p + ylab("Numbers in Compartments") + xlab(NULL)
  p <- p + scale_y_continuous(labels = scales::comma)
  p <- p + theme(strip.background = element_rect(colour="black", fill="white", 
                                                 size=1, linetype="solid"))
  p
}



