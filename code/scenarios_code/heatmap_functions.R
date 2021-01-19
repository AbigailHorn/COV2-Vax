
# scenario.id
# vx.id
# vx.coverage.id
# vx.start.date.id
# delay.start.65.vx.id
# NPI.id


######################################################
######################################################
## R(t) vs. coverage %
## Requires filter.65.delay
# data <- traj.CI
# var.to.plot <- "D"
# y.var <- "R0.id"
# x.var <- "vx.coverage.id"
# plot.title <- "100% Vaccine Uptake"
# filter.coverage <- 0
# filter.65.delay <- 10
# # Come back to modify labels for the x.axis based on what I want to visualize (if it is always customized, make it customized)
# 
# out <- plot_heatmap(data=data,x.var=x.var,y.var=y.var, var.to.plot = var.to.plot, filter.coverage = filter.coverage, filter.65.delay = filter.65.delay, plot.title = plot.title)
# out


######################################################
## R(t) vs. delay start 65+
## Requires filter.coverage

# data <- traj.CI
# var.to.plot <- "D"
# y.var <- "R0.id"
# x.var <- "delay.start.65.vx.id"
# plot.title <- "100% Vaccine Uptake"
# filter.coverage <- 0.8
# filter.65.delay <- 0
# # Come back to modify labels for the x.axis based on what I want to visualize (if it is always customized, make it customized)
# 
# out <- plot_heatmap(data=data,x.var=x.var,y.var=y.var, var.to.plot = var.to.plot, filter.coverage = filter.coverage, filter.65.delay = filter.65.delay, plot.title = plot.title)
# out


######################################################
## R(t) vs. vx.id
## DOES NOT require filter.coverage

#testing <- traj.CI %>% filter(state.name=="Htot") %>% filter(date==max(date))

# data <- traj.CI
# var.to.plot <- "Htot"
# y.var <- "R0.id"
# x.var <- "vx.id"
# plot.title <- "100% Vaccine Uptake"
# filter.coverage <- 0
# filter.65.delay <- 0
# #scale.x.labels <- c(vx.FAST="Fast", vx.EXPECTED="Expected", vx.SLOW="Slow" )
# scale.x.labels <- apply(expand.grid(c("Fast ", "Expected ","Slow "), c("80%","60%")) , 1, paste, collapse="")
# names(scale.x.labels) <- unique(data[,x.var])
# #scale.x.lables <- as.character(scale.x.labels)
# 
# out <- plot_heatmap(data=data,x.var=x.var,y.var=y.var, var.to.plot = var.to.plot, filter.coverage = filter.coverage, filter.65.delay = filter.65.delay, plot.title = plot.title, scale.x.labels=scale.x.labels)
# out



#heatmap.coverage.1 + heatmap.coverage.025 + plot_layout(guides = "collect")

plot_heatmap <- function(data, x.var, y.var, var.to.plot, filter.coverage, filter.65.delay, plot.title, scale.x.labels=NULL, scale.y.labels=NULL){
  
  if (!is.null(scale.x.labels)){
    scale.x.labels.in <- scale.x.labels
    #c(vx.FAST="Fast", vx.EXPECTED="Expected", vx.SLOW="Slow" )
  }
  if (!is.null(scale.y.labels)){
    scale.y.labels.in <- scale.y.labels
  }
  
  longnames <- c("Peak Obs. Infected",
                 "Cum. Obs. Infected",
                 "Peak Tot. Infected",
                 "Cum. Tot. Infected",
                 "New in Hospital",
                 "Peak in Hospital",
                 "Cum. in Hospital",
                 "New Deaths",
                 "Cum. Deaths",
                 "CFR",
                 "IFR")
  
  names(longnames) <- c("I",
                        "Idetectcum",
                        "Itot",
                        "Itotcum",
                        "H_new",
                        "Htot",
                        "Htotcum",
                        "D_new",
                        "D",
                        "CFRobs",
                        "CFRactual"
  )
  
  name.var.to.plot <- longnames[names(longnames)==var.to.plot]
  
  data <- data %>% filter(state.name==var.to.plot) #%>% filter(date==max(date)) %>% filter(vx.coverage.id==1)
  
  if (var.to.plot %in% c("D", "Idetectcum","Itotcum","Htotcum","CFRactual")){
    data <- data %>% filter(date==max(date)) }
  if (var.to.plot == "Htot" || state.name=="Itot" || state.name=="I"){
    data <- data %>% 
      group_by(scenario.id) %>%
      filter(median == max(median)) %>% as.data.frame()
  }
  
  max.median <- max(data$median)
  min.median <- min(data$median)
  
  if (filter.coverage!=0){
  data <- data %>% filter(vx.coverage.id==filter.coverage)}
  
  if (filter.65.delay!=0){
    data <- data %>% filter(delay.start.65.vx.id==filter.65.delay)}
  
  #data$NPI.id <- factor(data$NPI.id, c("NPI.Good","NPI.Obs","NPI.Bad"))
  
  #  try: label=coverage.id in aes()
  
  # ggplot(data, aes(x = factor(vx.id), y = factor(R0.id), fill = median))+
  #   geom_tile() + 
  #   geom_text(aes(label=round( median,-2 )  ), color="white" ) 
  
  p <-   
    ggplot(data, aes(x = factor(data[,x.var]), y = factor(data[,y.var]), fill = median))+
    geom_tile() + 
    geom_text(aes(label=round( median,-2 )  ), color="white" ) +
    scale_fill_gradient(name=paste0(name.var.to.plot, "\n (Median)"), low = "blue", high = "red", limits = range(min.median,max.median)) + theme_minimal() +
    # scale_y_discrete(name = "R(t) Scenario", labels = c(NPI.Bad = "3", 
    #                                                     NPI.Obs = "2", NPI.Good = "1.3")) +
    # scale_x_discrete(name = "Vaccination rate", labels=c(vx.FAST="Fast", vx.EXPECTED="Expected", vx.SLOW="Slow" ))+
    theme(plot.title = element_text(size=20, hjust = 0.5), axis.text.y = element_text(size=15), axis.text.x = element_text(size=15)) +
    labs(title = plot.title)
  if (y.var == "R0.id"){
    p <- p + scale_y_discrete(name = "R(t) Scenario") #, labels = c(NPI.Bad = "3", NPI.Obs = "2", NPI.Good = "1.3"))
  }
  if (y.var == "sero.id"){
    p <- p + scale_y_discrete(name = "% Recovered (Sero+) Vaccinated", labels = percent(scale.y.labels.in))
  }
  if (x.var == "vx.id"){
    p <- p +  scale_x_discrete(name = "Vaccination rate")
  }
  if (x.var == "vx.coverage.id"){
    p <- p + scale_x_discrete( name = "Vax % Uptake", labels= percent(data$vx.coverage.id))
  }
  if (x.var == "delay.start.65.vx.id"){
    p <- p + scale_x_discrete( name = "Phase II Start Date")
  }
  
  p <- p + geom_vline(xintercept = 1:length(unique(data[,x.var]))+0.5, color = "white", size = 1)
  p <- p + geom_hline(yintercept = 1:length(unique(data[,y.var]))+0.5, color = "white", size = 1)
  
}
  
 



########################################################################################################################
########################################################################################################################

# 
# vars.to.plot <- "Htot"
# y.max.in <- 14000
# startDatePlot <- "2020-11-01"
# chart.title <- "Hospitalizations: Scenarios with R(t) from 1 to 2.5"
# plot.scenarios.single.state <- plot.together.noribbon(traj.CI=data.plot.grid , data.in=NULL, endDatePlot=endDatePlot, vars.to.plot=vars.to.plot, y.lab.in=vars.to.plot, y.max.in=y.max.in, chart.title=chart.title, plot.capacity=NULL, plot.annotations=NULL, startDatePlot=startDatePlot)
# plot.scenarios.single.state

plot.together.noribbon <- function(traj.CI=traj.CI, data.in=data.in, endDatePlot=endDatePlot, vars.to.plot, y.lab.in, y.max.in, chart.title, plot.capacity, plot.annotations,startDatePlot) {
  
  ###########
  ### traj.CI
  ###########
  
  ## Filter only to variable of interest
  traj.CI <- traj.CI %>%  dplyr::filter(state.name %in% vars.to.plot)
  
  ## Select only more recent dates
  init.date <- as.Date("2020-03-01")
  startDatePlot <- as.Date(startDatePlot) #init.date #
  endDatePlot <- as.Date(endDatePlot) #startDatePlot + time.steps.4plot #- 40  # the constant 40 because the traj are not aligned to start date
  traj.CI <- traj.CI %>% dplyr::filter(date >= startDatePlot) %>% dplyr::filter(date < endDatePlot)
  
  ## Add title
  traj.CI$title <- chart.title
  
  ###########
  ### data.in
  ###########
  
  ## Data in -- plot only for the selected variable
  if(!is.null(data.in)){
    
    if(any(vars.to.plot %in% colnames(data.in))) {  # FIX LATER -- REMOVE TMP
      
      ## Filter only to variable of interest
      vars.to.extract <- vars.to.plot[vars.to.plot %in% colnames(data.in) ]
      
      data.in<- data.in %>% dplyr::select(vars.to.extract)
      
      ## ALIGN DATES: DATA
      no_obs <- nrow(data.in)
      step <- 0:(no_obs-1)
      date <- init.date + step
      data.date <- cbind(date,data.in)
      rownames(data.date) <- date
      #data.date$date <- NULL
      
      ## Select only more recent dates
      data.date <- data.date %>% dplyr::filter(date > startDatePlot)
      #      data.processed <- reshape2::melt(data.date, measure.vars = c(2:ncol(data.date)), variable.name = "state.name")
      data.processed <- reshape2::melt(data.date, measure.vars = c(2:ncol(data.date)), variable.name = "scenario.id")
      
    }
    
    else {data.processed = NULL}
  }
  
  #####################
  ### colors and names
  #####################
  
  longnames <- c("Susceptible",
                 "New Obs. Infected",
                 "Current Obs. Infected",
                 "Cum. Obs. Infected",
                 "Current Tot. Infected",
                 "Cum. Tot. Infected",
                 "New in Hospital",
                 "Current in Hospital",
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
                 "Vax Recovered",
                 "CFR",
                 "IFR",
                 "R(t) Scenario",
                 "New Vax",
                 "Total Vax",
                 "New Vax 65+",
                 "Total Vax 65+")
  
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
                        "Vx.given.65.tot"
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
    "grey50"
  )
  
  
  names(cols.list) <- names(longnames)
  color.this.var <- as.character(cols.list[vars.to.plot])
  
  ##################
  ### CREATE PLOT
  ##################
  
  p <- ggplot(data = traj.CI,
              aes(x = date,
                  y = median, ymin = low_95, ymax = up_95,
                  color = scenario.id,
                  fill = scenario.id,
                  group = scenario.id))
  
  # p <- p +  geom_ribbon(data = traj.CI,
  #                       aes(x = date,
  #                           y = median, ymin = low_50, ymax = up_50,
  #                           color = scenario.id,
  #                           fill = scenario.id,
  #                           group = scenario.id),alpha = .5, inherit.aes=TRUE, color=FALSE)  
  # p <- ggplot(data = traj.CI,
  #             aes(x = date,
  #                 y = median, ymin = low_95, ymax = up_95,
  #                 color = state.name,
  #                 fill = state.name,
  #                 group = state.name))
  # 
  # p <- p +  geom_ribbon(data = traj.CI,
  #                       aes(x = date,
  #                           y = median, ymin = low_50, ymax = up_50,
  #                           color = state.name,
  #                           fill = state.name,
  #                           group = state.name),alpha = .5, inherit.aes=TRUE, color=FALSE)
  
  #p <- p +  scale_fill_manual(values = c(color.this.var),labels = longnames) + scale_color_manual(values = c(color.this.var), labels = longnames)
  p <- p + geom_line() #+ geom_ribbon(alpha = 0.2, color = FALSE)
  
  if(!is.null(data.in)){
    # p <- p + geom_point(data = data.processed,
    #                     aes(x = date, y = value,
    #                         color = state.name),
    #                     alpha = 0.7,
    #                     inherit.aes = FALSE)
    p <- p + geom_point(data = data.processed,
                        aes(x = date, y = value,
                            color = "black"),
                        alpha = 0.7,
                        inherit.aes = FALSE)
  }
  
  ##################
  ## ADD CAPACITY
  if (!is.null(plot.capacity)){
    ##################
    ### CREATE CAPACITY DATA FRAME
    capacity.vals <- as.data.frame(matrix(NA, nrow=length(levels(traj.CI$state.name)), ncol=2))
    capacity.vals[,1] <- levels(traj.CI$state.name)
    rownames(capacity.vals) <- levels(traj.CI$state.name)
    capacity.vals["Htot",2] <- 4000
    capacity.vals["Q",2] <- 2245
    capacity.vals["V",2] <-1000
    colnames(capacity.vals) <- c("state.name","capacity")
    ##################
    ### ADD CAPACITY LINES
    capacity.vals <- capacity.vals %>% filter(state.name %in% vars.to.plot)
    p <- p + geom_hline(data = capacity.vals, aes(yintercept=capacity),linetype = "dashed", colour="azure4")
  }
  
  #################
  ## ADD DATE ANNOTATIONS
  if (!is.null(plot.annotations)){
    ######### Create data frame with annotations
    traj.CI.date <- as.data.frame(matrix(NA, nrow=7, ncol=3))
    colnames(traj.CI.date) <- c("date","date.label","y.place")
    traj.CI.date$date <- c(as.Date("2020-03-19"),as.Date("2020-05-08"),as.Date("2020-06-12"),as.Date("2020-07-01"),as.Date("2020-08-18"),as.Date("2020-10-31"),as.Date("2020-11-26"))
    traj.CI.date$date.label <- c("Stage I", "Stage II", "Stage III", "Modifications", "School Year", "Halloween", "Thanksgiving")
    traj.CI.date$y.place <- c(1:7)
    ######### Add data frame with annotations
    p <- p + geom_vline(data=traj.CI.date, aes(xintercept=as.Date(date)), linetype="dashed",colour="azure4", size=.35) +
      # annotate("text", label = traj.CI.date$date.label, x = traj.CI.date$date, y = (y.max.in/2)+(y.max.in/20)*traj.CI.date$y.place, size = 3.5, colour = "black")
      annotate("text", label = traj.CI.date$date.label, x = traj.CI.date$date, y = (y.max.in)-(y.max.in/25)*traj.CI.date$y.place, size = 3.5, colour = "black")
  }
  
  ##################
  ## FINISH PLOT
  p <- p + theme_bw() + theme(legend.title = element_blank())
  p <- p + scale_x_date(limits = as.Date(c(startDatePlot,endDatePlot)), date_breaks = "1 month" , date_labels = "%d-%b-%y")
  p <- p + scale_y_continuous(limits = c(0,y.max.in), breaks = seq(from = 0, to = y.max.in, by = y.max.in/10))
  p <- p + theme(axis.text.x = element_text(angle = 90),
                 strip.text.x = element_text(size = 12, face = "bold"))
  #  p <- p + ylab(paste0("Number  ", as.character(longnames[var.to.plot]))) + xlab(NULL)
  #p <- p + ylab("Probability") + xlab(NULL)
  p <- p + ylab(y.lab.in) + xlab(NULL)
  #p <- p + labs(title = title.input)
  #p<-p+theme(plot.title = element_text(size = 12, hjust = 0.5, face="bold"))
  
  p <- p + scale_color_viridis(discrete=TRUE)#scale_color_brewer()
  
  p <- p + facet_grid(. ~ title)
  
  
  p
  
  
  
}


