


########################################################################################
## SPECIFYING EPIDEMIC MODEL TO BE SIMULATED LOOPING OVER SCENARIOS
## NOTE: NPI scenarios are called in from data "fn_readin.csv"
## NOTE: Protect scenarios are called in from weighted.avg.protect.mat
########################################################################################

correlated.param.SCENARIOS.ratio <- function(fn_t_readin_input=fn_t_readin_input, ABC.out.mat, iter, time.steps, vx.scenarios.mat, sero.scenario.in, mu.scenario.in,
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
  SCENARIOS.out.ratio <- vector("list", n.scenarios)
  
  scenario.idx <- 1
  ratio.idx <- 1
  
  for (sero.idx in 1:length(sero.scenario.in)){
    
    frac.sero.vax = sero.scenario.in[sero.idx]
    
    print(paste0("frac.sero.vax=", frac.sero.vax))
    
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
          
          # print(paste0("R(t) = ", mu.scenario*3.5))
          # print(paste0("R(t) = ", round(mu.scenario*3.5,1)))
          
          R0_save <- round(mu.scenario*3.5,1)
          
          ## BIND INCLUDING OFFSETING OBSERVED DATA BY START DATE
          TEST.out[[idx]] <- cbind(data.frame(scenario.id = paste0("target.", round(1-frac.sero.vax,1), "_" , name.vx ) , vx.id = name.vx, vx.start.date.id=start.date.vx, delay.start.65.vx.id = delay.start.65.vx, sero.id= round(1-frac.sero.vax,1), R0.id = R0_save, par.id = idx, date = -start_time+TEST$step), TEST)
          #TEST.out[[idx]] <- cbind(data.frame(scenario.id = paste0(scenario.idx, "_", name.vx, "_Rt.", R0_save ) , vx.id = name.vx, vx.start.date.id=start.date.vx, vx.coverage.id = vx.coverage, delay.start.65.vx.id = delay.start.65.vx, R0.id = R0_save, par.id = idx, date = -start_time+TEST$step), TEST)
        }  # end over idx
        
        ## Add to a dataframe over all idx
        TEST.out <- do.call(rbind, TEST.out)
        
        ###################################
        ## Get CI for TEST.out (raw values)
        CI.TEST.out <- get.CI.SCENARIOS(TEST.out = TEST.out)
        
        ###################################
        ## Get CI for difference from 0
        scenario.id = paste0("target.", round(1-frac.sero.vax,1), "_" , name.vx )
        
        print(paste0("scenario.id = ", scenario.id))
        
        if (name.vx=="0"){
          TEST.out.0 <- TEST.out  # Get TEST.out.0
        }
        
        n.cols.TEST.out <- ncol(TEST.out)
        diff0 <- abs ( TEST.out.0[,c(11:n.cols.TEST.out)] - TEST.out[,c(11:n.cols.TEST.out)] )  # Get difference TEST.out.0 - TEST.out
        TEST.out.diff0 <- cbind(TEST.out[,c(1:10)], diff0)
        CI.TEST.out.diff0 <- get.CI.SCENARIOS(TEST.out = TEST.out.diff0)  # Get CI for TEST.out.diff0

        ###################################
        ## Get CI for ratio TEST.out.diff.0 / BASELINE.diff.0
        if (scenario.id == "target.0_10k"){
          TEST.out.BASELINE <- TEST.out   # Get TEST.out.Baseline
          BASELINE.diff0 <- abs ( TEST.out.0[,c(11:n.cols.TEST.out)] - TEST.out.BASELINE[,c(11:n.cols.TEST.out)] )  # Get difference TEST.out.0 - TEST.out.BASELINE
          TEST.out.BASELINE.diff0 <- cbind(TEST.out.BASELINE[,c(1:10)], BASELINE.diff0)
        }
        
        if (name.vx!="0"){
          #if (scenario.id!="target.0_10k"){
            ratio.0.BASELINE <- TEST.out.diff0[,c(11:n.cols.TEST.out)] / TEST.out.BASELINE.diff0[,c(11:n.cols.TEST.out)]  # Get ratio TEST.out.diff.0 - TEST.out.BASELINE.diff0
            TEST.out.ratio.0.BASELINE <- cbind(TEST.out[,c(1:10)], ratio.0.BASELINE)
            CI.ratio.0.BASELINE <- get.CI.SCENARIOS(TEST.out = TEST.out.ratio.0.BASELINE)  # Get CI for TEST.out.ratio.0.BASELINE
            SCENARIOS.out.ratio[[ratio.idx]] <- CI.ratio.0.BASELINE
            ratio.idx = ratio.idx+1
          #}
        }
        
        print(paste0("ratio.idx=", ratio.idx))
        
        ## Put TEST.out dataframe into SCENARIOS.out
        SCENARIOS.out[[scenario.idx]] <- CI.TEST.out
        SCENARIOS.out.diff[[scenario.idx]] <- CI.TEST.out.diff0
        
        ##
        scenario.idx <- scenario.idx + 1 
        
        rm(TEST.out)
        rm(TEST.out.diff0)
        rm(TEST.out.ratio.0.BASELINE)
        
      } # end over protect.idx
    } # end over x protect.scenario NPI.idx
  } # end over sero.idx
  
  # Output
  all.scenarios <- do.call(rbind,SCENARIOS.out)
  all.scenarios.diff <- do.call(rbind,SCENARIOS.out.diff)
  all.scenarios.ratio <- do.call(rbind,SCENARIOS.out.ratio)
  
  output <- vector("list", 2)
  output[[1]] <- all.scenarios
  output[[2]] <- all.scenarios.diff
  output[[3]] <- all.scenarios.ratio
  
  return(output)
}
