
scenario.id
vx.id
vx.coverage.id
vx.start.date.id
delay.start.65.vx.id
NPI.id

# data %>% 
#   group_by(NPI.id, vx.id) %>%
#   filter(median == max(median)) 
# 
# data %>% 
#   group_by(scenario.id) %>%
#   filter(median == max(median)) 



data <- traj.CI %>% filter(state.name=="D") %>% filter(date==max(date)) # %>% filter(vx.coverage.id==1)

data$NPI.id <- factor(data$NPI.id, c("NPI.Good","NPI.Obs","NPI.Bad"))

#plot_heatmap = function(data, num, param){
#  p <- 

min.median <- min(traj.CI$median)
max.median <- max(traj.CI$median)

#coverage1 <-   
  ggplot(data, aes(x = factor(delay.start.65.vx.id), y = factor(NPI.id), fill = (median)))+
  geom_tile() + 
  geom_text(aes(label=round( median,-2 )  ), color="white" ) +
  scale_fill_gradient(name="Deaths (Median)", low = "blue", high = "red", limits = range(min.median,max.median)) + theme_minimal() +
  scale_y_discrete(name = "R(t) Scenario", labels = c(NPI.Bad = "3", 
                              NPI.Obs = "2", NPI.Good = "1.3")) +
  scale_x_discrete(name = "Vaccination rate", labels=c(vx.FAST.coverage1="Fast", vx.EXPECTED.coverage1="Expected", vx.SLOW.coverage1="Slow" ))+
  theme(plot.title = element_text(size=20, hjust = 0.5), axis.text.y = element_text(size=15), axis.text.x = element_text(size=15)) +
  labs(title = "100% Vaccine Uptake")
#plot.title = element_text(size = (15))
  
  #coverage1 + coverage.25 + plot_layout(guides = "collect")
#scale_fill_continuous(breaks = c(2000, 4000), labels = c("2k", "4k"))    

#plot_heatmap = function(data, num, param){
#  p <- 



######################################################
######################################################

#before function:
data <- traj.CI
#data <- data %>% filter(vx.coverage.id==1)

#inputs:
var.to.plot <- "Htot"
y.var <- "NPI.id"
x.var <- "vx.id"
plot.title <- "100% Vaccine Uptake"
filter.coverage <- 0.25

heatmap.coverage.025 <- plot_heatmap(data=data,x.var=x.var,y.var=y.var, var.to.plot = var.to.plot, filter.coverage = filter.coverage, plot.title = plot.title)

heatmap.coverage.1 + heatmap.coverage.025 + plot_layout(guides = "collect")

plot_heatmap <- function(data, x.var, y.var, var.to.plot, filter.coverage, plot.title){
  
  longnames <- c("Current Obs. Infected",
                 "Cum. Obs. Infected",
                 "Current Tot. Infected",
                 "Cum. Tot. Infected",
                 "New in Hospital",
                 "Current in Hospital",
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
  
  if (state.name %in% c("D", "Idetectcum","Itotcum","Htotcum","CFRactual")){
    data <- data %>% filter(date==max(date))
  }
  if (state.name == "Htot" || state.name=="Itot" || state.name=="I"){
    data <- data %>% 
      group_by(scenario.id) %>%
      filter(median == max(median)) %>% as.data.frame()
  }
  
  max.median <- max(data$median)
  min.median <- min(data$median)
  
  data <- data %>% filter(vx.coverage.id==filter.coverage)
  
  data$NPI.id <- factor(data$NPI.id, c("NPI.Good","NPI.Obs","NPI.Bad"))
  
  #  try: label=coverage.id in aes()
  
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
  if (y.var == "NPI.id"){
    p <- p + scale_y_discrete(name = "R(t) Scenario", labels = c(NPI.Bad = "3", NPI.Obs = "2", NPI.Good = "1.3"))
  }
  if (x.var == "vx.id"){
    p <- p +  scale_x_discrete(name = "Vaccination rate", labels=c(vx.FAST="Fast", vx.EXPECTED="Expected", vx.SLOW="Slow" ))
  }
  if (x.var == "vx.coverage.id"){
    p <- p + scale_x_discrete( name = "Vax % Uptake", labels= percent(vx.coverage.id))
  }
  if (x.var == "delay.start.65.vx.id"){
    p <- p + scale_x_discrete( name = "Phase II Start Date")
  }
  
  p <- p + geom_vline(xintercept = 1:length(unique(data[,x.var]))+0.5, color = "white", size = 1)
  p <- p + geom_hline(yintercept = 1:length(unique(data[,y.var]))+0.5, color = "white", size = 1)
  
}
  
  

  #   #scale_fill_hue(h=c(200,400))+theme_minimal()
  #   
  # 
  # 
  # 
  #   
  #   scale_fill_manual(breaks = seq(0, 5, by = 1),  
  #                     values = c("#E6E6E6", col_youngadults, col_adults, col_elderly, col_kids, col_all)) +
  #   theme(legend.position = "none",
  #         axis.title.x = element_blank(),
  #         axis.text = element_text(size = 10)) +
  #   scale_x_discrete(expand = c(0,0), breaks = seq(0, 50, by = 10)) + 
  #   coord_fixed(50*4/(6*num))
  # 
  # if (param == "R0"){
  #   p <- p + scale_y_discrete(expand = c(0,0), breaks = seq(1, num, by = 2),
  #                             labels = seq(1.3, 2.6, by = 0.2))
  # } else if (param == "ve"){
  #   p <- p + scale_y_discrete(expand = c(0,0), breaks = seq(1, num, by = 2),
  #                             labels = c(" 30", " 50", " 70", " 90"))
  # } else if (param == "country") {
  #   p <- p + scale_y_discrete(expand = c(0,0), breaks = seq(1, num, by = 1),
  #                             labels = c("BEL", "USA", "IND", "ESP", "ZWE", "BRA", "CHN", "ZAF", "POL"))
  # } else if (param == "rollout") {
  #   p <- p + scale_y_discrete(expand = c(0,0), breaks = seq(1, num, by = 1),
  #                             labels = seq(0.25, 2, by = 0.25))
  # } else if (param == "ve_I"){
  #   p <- p + scale_y_discrete(expand = c(0,0), breaks = seq(1, num, by = 2),
  #                             labels = c(" 0"," 20", " 40", " 60", " 80", "100"))
  # }
  # 
  # 
  # for (i in seq(1.5, num + 0.5, by = 1)){
  #   p <- p + geom_hline(yintercept = i, col = "white")
  # }
  # 
  # for (i in seq(10.5, 40.5, by = 10)){
  #   p <- p + geom_vline(xintercept = i, col = "white")
  # }
  # 
  # p


plot.together.ribbon <- function(traj.CI=traj.CI, data.in=data.in, endDatePlot=endDatePlot, vars.to.plot, y.lab.in, y.max.in, chart.title, plot.capacity, plot.annotations) {
  
  ###########
  ### traj.CI
  ###########
  
  ## Filter only to variable of interest
  traj.CI <- traj.CI %>%  dplyr::filter(state.name %in% vars.to.plot)
  
  ## Select only more recent dates
  init.date <- as.Date("2020-03-01")
  startDatePlot <- init.date #
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
  
  p <- p +  geom_ribbon(data = traj.CI,
                        aes(x = date,
                            y = median, ymin = low_50, ymax = up_50,
                            color = scenario.id,
                            fill = scenario.id,
                            group = scenario.id),alpha = .5, inherit.aes=TRUE, color=FALSE)  
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
  p <- p + geom_line() + geom_ribbon(alpha = 0.2, color = FALSE)
  
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
  p <- p + facet_grid(. ~ title)
  
  
  p
  
  
  
}



########################################################################################################################
########################################################################################################################

traj.CI.save <- traj.CI
vars.to.plot <- "Htot"
y.max.in <- 25000
startDatePlot <- "2020-11-01"
plot.scenarios.single.state <- plot.together.noribbon(traj.CI=traj.CI, data.in=NULL, endDatePlot=endDatePlot, vars.to.plot=vars.to.plot, y.lab.in=vars.to.plot, y.max.in=y.max.in, chart.title="Scenario output", plot.capacity=NULL, plot.annotations=NULL, startDatePlot=startDatePlot)
plot.scenarios.single.state

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


