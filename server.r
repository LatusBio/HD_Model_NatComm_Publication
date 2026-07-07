options(gganimate.renderer = gifski_renderer())
server <- function(input, output, session) {
  
  safe_min_age <- function(df, fallback) {
    if (nrow(df) == 0) return(fallback)
    return(min(df$Age))
  }
  
  # Instantiate progress bar
  updateProgressBar(session, id = "pbRun", value = 0, title = "Awaiting input...", status = "default")
  
  
  ########################################
  # Somatic Repeat Expansion Model Parameters
  ########################################
  
  ##### HD degeneration variables ##### 
  # The below length change mutation variable encodes the ratio of expansion to contraction events (recommend using probability ratio of 1 contraction event per 6 expansion events)
  length_change_mutation <- reactive({
    c(rep(-1, input$n_contractions), rep(1, input$n_expansions))
  })
  rate1 <- .035 # 36to55_CAGs Handsaker et al (supp Note 4) 3.51% +/- 0.83%.
  rate2 <- .576 # >55 CAGs Handsaker et al (supp Note 4) 57.6% +/- 12.76%.
  
  # Two phase linear model (Handsaker et al)
  T1 <- 33.5
  T2 <- 55
  r1 <- .035
  r2 <- .576
  
  # Create the function
  twoPhaseLinearModel <- function(x){ r1*max(x-T1,0) + r2*max(x-T2,0)}
  
  # Create plot renderer tracker
  renderState <- reactiveValues(animationDone = FALSE, linePlotDone = FALSE, healthyDone = FALSE)
  renderState$animationDone <- FALSE
  renderState$plotDone <- FALSE
  renderState$healthyDone <- FALSE
  
  stored_plots <- reactiveValues(linePlot = NULL, healthyPlot = NULL)
  
  # Build descriptive filename prefix from current input parameters
  get_file_prefix <- reactive({
    paste0("HD_sim",
           "_CAG", input$germlineHDRepeatLength,
           "_intAge", input$age_at_intervention,
           "_life", input$maxLifetime,
           "_MSH3kd", input$therapeuticModifier * 100, "pct",
           "_transd", input$transductionRate * 100, "pct",
           "_nMSN", input$nMSNs_CaudatePutamen,
           "_thresh", input$HealthyCellCAG_Threshold,
           "_contr", input$n_contractions,
           "_exp", input$n_expansions)
  })
  
  ########################################
  # Define object classes
  ########################################
  
  
  # Create object class representing an individual cell
  Cell <- R6Class("Cell",
                  public = list(
                    repeat_length = NULL,
                    therapeutic_status = 0,
                    year_counter = 0,
                    
                    initialize = function(repeat_length = 0) {
                      self$repeat_length <- repeat_length
                      self$therapeutic_status <- 0
                      self$year_counter <- 0
                    },
                    
                    get_repeat_length = function() {
                      return(self$repeat_length)
                    },
                    
                    therapeutic_transduction_event = function() {
                      transductionOutcome <- rbinom(1,1, input$transductionRate)
                      if (transductionOutcome == 1) {
                        self$therapeutic_status <- "transduced"
                      } else {
                        self$therapeutic_status <- self$therapeutic_status
                      }
                    },
                    
                    get_therapeutic_transduction_status = function() {
                      return(self$therapeutic_status)
                    },
                    
                    mutate_repeat = function() {
                      # Compute probability of mutation event
                      AnnualMutationRate <- twoPhaseLinearModel(self$repeat_length)
                      
                      # Function to mutate repeat
                      mutate_length <- function(times) {
                        self$repeat_length <- self$repeat_length + sum(sample(length_change_mutation() , times, replace = TRUE))
                      }
                      
                      if(AnnualMutationRate <= 1){
                        # single mutation event
                        if (rbinom(1, 1, AnnualMutationRate) == 1) {
                          mutate_length(1)
                        }

                      } else {
                        # Multiple mutation events
                        cycles <- min(as.integer(AnnualMutationRate), 100) # Cap cycles at 100
                        remainder <- AnnualMutationRate %% 1
                        
                        # Process full cycles
                        mutate_length(cycles)
                        
                        # Process remainder
                        if (rbinom(1, 1, remainder) == 1) {
                          mutate_length(1)
                        }
                      }
                    },
                    
                    mutate_repeat_with_therapy = function() {
                      # Compute probability of mutation event
                      # Here we use therapeuticModifier value to reduce the annual mutation rate  
                      AnnualMutationRate <- twoPhaseLinearModel(self$repeat_length) * (1 - input$therapeuticModifier)
                      
                      # Function to mutate repeat
                      mutate_length <- function(times) {
                        self$repeat_length <- self$repeat_length + sum(sample(length_change_mutation() , times, replace = TRUE))
                      }
                      
                      if(AnnualMutationRate <= 1){
                        # single mutation event
                        if (rbinom(1, 1, AnnualMutationRate) == 1) {
                          mutate_length(1)
                        }
                        
                      } else {
                        # multiple mutation events
                        cycles <- min(as.integer(AnnualMutationRate), 100) # Cap cycles at 100
                        remainder <- AnnualMutationRate %% 1 
                        
                        # Process full cycles
                        mutate_length(cycles)
                        
                        # Process remainder
                        if (rbinom(1, 1, remainder) == 1) {
                          mutate_length(1)
                        }
                      }
                    },
                    
                    set_repeat_length = function(new_length) {
                      self$repeat_length <- new_length
                    },
                    
                    describe = function() {
                      cat("This cell has a repeat length of", self$repeat_length, "\n")
                    }
                    
                  )
  )
  
  # Create object class representing a population of cells
  PopulationManager <- R6Class("PopulationManager",
                               public = list(
                                 cells = NULL,
                                 
                                 initialize = function(cell_collection) {
                                   self$cells <- cell_collection
                                 },
                                 
                                 mutate_all_cells = function() {
                                   for (cell in self$cells) {
                                     if (cell$get_repeat_length() <= 300) {
                                       cell$mutate_repeat()
                                     }
                                   }
                                 },
                                 
                                 mutate_transduced_cells_with_therapy = function() {
                                   for (cell in self$cells) {
                                     if (cell$get_repeat_length() <= 300) {
                                       if (cell$get_therapeutic_transduction_status() == "transduced") {
                                         cell$mutate_repeat_with_therapy()
                                       } else {
                                         cell$mutate_repeat()
                                       }
                                     }
                                   }
                                 },
                                 
                                 administer_gene_therapy = function() {
                                   for (cell in self$cells) {
                                     cell$therapeutic_transduction_event()
                                   }
                                 },
                                 
                                 run_lifecycle = function(years) {
                                   data <- data.frame()
                                   for (year in 0:years) {
                                     repeat_lengths <- sapply(self$cells, function(cell) cell$get_repeat_length())
                                     transduction_status <- 0
                                     data <- rbind(data, data.frame(Year = year, RepeatLength = repeat_lengths, CellID = c(1:length(repeat_lengths)), Transduction = transduction_status))
                                     self$mutate_all_cells()
                                   }
                                   return(data)
                                 },
                                 
                                 run_lifecycle_with_therapy = function(years, untreated) {
                                   if (untreated != TRUE) {
                                     data <- data.frame()
                                     # run repeat expansion before gene therapy intervention
                                     for (year in 0:(input$age_at_intervention-1)) {
                                       repeat_lengths <- sapply(self$cells, function(cell) cell$get_repeat_length())
                                       transduction_status <- sapply(self$cells, function(cell) cell$get_therapeutic_transduction_status())
                                       data <- rbind(data, data.frame(Year = year, RepeatLength = repeat_lengths, CellID = c(1:length(repeat_lengths)), Transduction = transduction_status))
                                       self$mutate_all_cells()
                                     }
                                     
                                     # Administer gene therapy
                                     self$administer_gene_therapy()
                                     
                                     # run repeat expansion after gene therapy intervention
                                     for (year in input$age_at_intervention:years) {
                                       repeat_lengths <- sapply(self$cells, function(cell) cell$get_repeat_length())
                                       transduction_status <- sapply(self$cells, function(cell) cell$get_therapeutic_transduction_status())
                                       data <- rbind(data, data.frame(Year = year, RepeatLength = repeat_lengths, CellID = c(1:length(repeat_lengths)), Transduction = transduction_status))
                                       self$mutate_transduced_cells_with_therapy()
                                     }
                                     return(data)
                                   } else {
                                     data <- data.frame()
                                     for (year in 0:years) {
                                       repeat_lengths <- sapply(self$cells, function(cell) cell$get_repeat_length())
                                       transduction_status <- sapply(self$cells, function(cell) cell$get_therapeutic_transduction_status())
                                       data <- rbind(data, data.frame(Year = year, RepeatLength = repeat_lengths, CellID = c(1:length(repeat_lengths)), Transduction = transduction_status))
                                       self$mutate_transduced_cells_with_therapy()
                                     }
                                     return(data)
                                   }
                                 }
                               )
  )
  
  
  observe({
    if (renderState$animationDone && renderState$linePlotDone && renderState$healthyDone) {
      updateProgressBar(session, id = "pbRun", value = 100, title = "Complete!", status = "success")
      
      # Re-enable Run button
      session$sendCustomMessage(type = "enableButton", message = list(id = "run"))
      
      # Reset flags for next run
      renderState$animationDone <- FALSE
      renderState$linePlotDone <- FALSE
      renderState$healthyDone <- FALSE
    }
  })
  
  
  ########################################
  # Run Simulation
  ########################################

  # Simulation must happen within this reactive block
  simulation_data <- eventReactive(input$run, {
    validate(
      need(input$maxLifetime > input$age_at_intervention,
           "Error: Max Lifetime must be greater than Age at Intervention."),
    )
    
    # Disable run simulation button after it is pressed
    session$sendCustomMessage(type = "disableButton", message = list(id = "run"))
    
    # Reset progress bar to 0
    updateProgressBar(session, id = "pbRun", value = 0, status = "info", title = "Running...")
    
    # Update progress bar
    updateProgressBar(session, id = "pbRun", value = 10, status = "info", title = "Starting simulation...")
    Sys.sleep(0.3)
    
    ######## CAP SCORE #############
    # Create a CAP score function
    # CAP = Age x (CAG-L)+ / K
    # (CAG-L)+ is a function that equals x when x >= 0 and 0 otherwise.
    
    CAG_L_fun <- function(x) {
      if (x >= 0) {
        return(x)
      } else {
        return(0)
      }
    }
    
    CAP_Score_Fun <- function(CurrentAge){
      L <- 30 # an estimate of the lower limit of CAG expansion values for which toxicity occurs
      K <- 6.49
      CAP_score <- CurrentAge * CAG_L_fun(input$germlineHDRepeatLength - L) / K
      return(list(CAP_score, K))
    }
    
    CAP_Score_Lookup_Table <- sapply(c(1:100), function(x){
      CAP_Score_Fun(x)[[1]]
    })
    
    Predicted_AgeOfOnset <- length(CAP_Score_Lookup_Table[CAP_Score_Lookup_Table <= 100])
    
    # # Compute stage data frames
    # stage1_df <- barplot_data_pct_dead_cells_merged[grep("1", barplot_data_pct_dead_cells_merged$HD_ISS_Stage), ]
    # stage2_df <- barplot_data_pct_dead_cells_merged[grep("2", barplot_data_pct_dead_cells_merged$HD_ISS_Stage), ]
    # stage3_df <- barplot_data_pct_dead_cells_merged[grep("3", barplot_data_pct_dead_cells_merged$HD_ISS_Stage), ]
    
    ### Instantiate Cells ###
    # Update progress bar
    updateProgressBar(session, id = "pbRun", value = 20, title = "Generating cell populations...")
    Sys.sleep(0.3)
    
    # Create cell collection
    cell_collection <- list()
    cell_collection_untreated <- list()

    #Instantiate cell objects with germline repeat length
    for (MSN in 1:input$nMSNs_CaudatePutamen) {
      MSN <- Cell$new(input$germlineHDRepeatLength)
      cell_collection[[length(cell_collection) + 1]] <- MSN
    }

    # Instantiate cell objects with germline repeat length
    for (MSN in 1:input$nMSNs_CaudatePutamen) {
      MSN <- Cell$new(input$germlineHDRepeatLength)
      cell_collection_untreated[[length(cell_collection_untreated) + 1]] <- MSN
    }
    

    ### Create and Manipulate Populations of Cells ###
    # Create PopulationManager
    population <- PopulationManager$new(cell_collection)
    population_untreated <- PopulationManager$new(cell_collection_untreated)
    
    # Update progress bar
    updateProgressBar(session, id = "pbRun", value = 30, title = "Running cell mutation simulation...")
    Sys.sleep(0.3)
    
    # Run simulation
    #if (input$therapy_on == TRUE) {
      # Get data for each year in lifetime (therapy on)
      plot_data <- population$run_lifecycle_with_therapy(input$maxLifetime, untreated = FALSE)
      plot_data_untreated <- population_untreated$run_lifecycle_with_therapy(input$maxLifetime, untreated = TRUE)
    #} else {
    #  plot_data <- population$run_lifecycle(input$maxLifetime)
    #  plot_data_untreated <- data.frame()
    #}
    
    ### Calculate Effect Size ###
    # Get pct of cells <= 150 
    get_healthy_cells_fun <- function(z) {
      levels <- levels(as.factor(z$Year))
      barplot_data_pct_healthy_cells <- sapply(levels, function(x) {
        y <- as.numeric(x)
        year_subset <- z[z$Year == y,]
        pcthealthy <- (sum(year_subset$RepeatLength <= input$HealthyCellCAG_Threshold, na.rm = TRUE) / length(year_subset$RepeatLength)) * 100
        return(c(y, pcthealthy))
      } )
      barplot_data_pct_healthy_cells <- barplot_data_pct_healthy_cells
      barplot_data_pct_healthy_cells <- t(as.data.frame(barplot_data_pct_healthy_cells))
      barplot_data_pct_healthy_cells <- as.data.frame(barplot_data_pct_healthy_cells)
      colnames(barplot_data_pct_healthy_cells) <- c("Age","pct_healthy_Cells")
      return(barplot_data_pct_healthy_cells)
    }
    
    barplot_data_pct_healthy_cells <- get_healthy_cells_fun(plot_data)
    barplot_data_pct_healthy_cells_untreated <- get_healthy_cells_fun(plot_data_untreated)
    barplot_data_pct_healthy_cells_merged <- merge(barplot_data_pct_healthy_cells,barplot_data_pct_healthy_cells_untreated, by = "Age")
    colnames(barplot_data_pct_healthy_cells_merged) <- c("Age","HealthyTreated","HealthyUntreated")
    
    # Calculate intersect of CAP-100 with untreated cells
    PctHealthCells_at_Predicted_AgeOfOnset <- barplot_data_pct_healthy_cells_untreated[Predicted_AgeOfOnset == barplot_data_pct_healthy_cells_untreated$Age,][1,2]
    
    # # Get Therapeutic benefit (years of delayed onset)
    shiftedAge <- min(barplot_data_pct_healthy_cells[PctHealthCells_at_Predicted_AgeOfOnset >= barplot_data_pct_healthy_cells$pct_healthy_Cells,][,1])
    therapeuticBenefitYears <- shiftedAge - Predicted_AgeOfOnset
    
    # Update progress bar
    updateProgressBar(session, id = "pbRun", value = 40, title = "Running cell mutation simulation...")
    Sys.sleep(0.3)
    
    # Calculate therapeutic benefit for treatment across all ages
    therapeuticBenefitYears_fun <- function(ages){
      shiftedAge <- min(barplot_data_pct_healthy_cells[PctHealthCells_at_Predicted_AgeOfOnset >= barplot_data_pct_healthy_cells$pct_healthy_Cells,][,1])
      return(shiftedAge)
    }
    
    # Calculate dead cell percentages for treated cells
    barplot_data_pct_dead_cells <- get_dead_cells_fun(plot_data)
    # Calculate dead cell percentages for untreated cells
    barplot_data_pct_dead_cells_untreated <- get_dead_cells_fun(plot_data_untreated)
    # Merge the two dead cell data frames by Age
    barplot_data_pct_dead_cells_merged <- merge(barplot_data_pct_dead_cells,
                                                barplot_data_pct_dead_cells_untreated,
                                                by = "Age")
    colnames(barplot_data_pct_dead_cells_merged) <- c("Age", "DeadTreated", "DeadUntreated")
    
    # Define HD-ISS Stage function
    Run_HDISS_score_fun_production <- function(pct_cells_alive, age, cag){
      age <- as.numeric(age)
      if (age %in% colnames(TotalFunctionalCapacity_Score)) {
        #print(paste0("age is ", age, " class is ", class(age)))
        
        # Extract functional capacity and ensure it's a single value
        FunctionalCapacityVal <- as.numeric(TotalFunctionalCapacity_Score[grep(paste0('\\b', cag, '\\b'), TotalFunctionalCapacity_Score$CAG), grep(paste0('\\b', age, '\\b'), colnames(TotalFunctionalCapacity_Score))])
        if (length(FunctionalCapacityVal) != 1) {
          stop("Functional Capacity should be a single value.")
          #FunctionalCapacityVal <- 1
        } else if (is.na(FunctionalCapacityVal)) {
           FunctionalCapacityVal <- 1
        } else if (FunctionalCapacityVal == 0) {
          FunctionalCapacityVal <- 1
        } else {
           FunctionalCapacityVal <- as.numeric(FunctionalCapacityVal)
        }
        #print(paste0("functional capacity is ", FunctionalCapacityVal, " class is ", class(FunctionalCapacityVal)))
        
        # Extract motor score and ensure it's a single value
        MotorScoreVal <- as.numeric(Motor_Score[grep(paste0('\\b', cag, '\\b'), Motor_Score$CAG), grep(paste0('\\b', age, '\\b'), colnames(Motor_Score))])
        if (length(as.numeric(MotorScoreVal)) != 1) {
          stop("Motor Score should be a single value.")
          #MotorScoreVal <- 0
        } else if (is.na(MotorScoreVal)) {
          MotorScoreVal <- 80
        } else {
          MotorScoreVal <- as.numeric(MotorScoreVal)
        }
        #print(paste0("Motor score is ", MotorScoreVal, " class                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         s is ", class(MotorScoreVal)))
        
        # Call the function with validated inputs
        HD_ISS_Stage <- HD_ISS_Stage_Calc_FUN(
          CurrentAge = as.numeric(age),
          cag = as.numeric(input$germlineHDRepeatLength), 
          brainVolume = as.numeric(pct_cells_alive), 
          motorScore = as.numeric(MotorScoreVal), 
          functionalCapacity = as.numeric(FunctionalCapacityVal)
        )
        
        # Check to confirm the function is working
        print(paste0("age is ", age, " ISS stage is ", HD_ISS_Stage, " Motor Score is ", MotorScoreVal, " functionalCapacity is ", FunctionalCapacityVal))
        
        # output only stage
        return(HD_ISS_Stage)
      } else {
        return(NA)
      }
    }
    
    # Update progress bar
    updateProgressBar(session, id = "pbRun", value = 60, title = "Running cell mutation simulation...")
    Sys.sleep(0.3)
    
    # Calculate HD-ISS stage across all ages
    HDISS_Stages <- apply(barplot_data_pct_dead_cells_merged, 1, function(x){
      #x <- as.character(x) 
      #print(x)
      #return(x)
      age <- x[1]
      #print(paste0("age is ", age))
      pctAliveCells <- (100 - x[3]) / 100
      #print(paste0("pctAliveCells is ", pctAliveCells))
      HDISS_Stage <- Run_HDISS_score_fun_production(pctAliveCells, age, input$germlineHDRepeatLength)
      return(HDISS_Stage)
    })
    barplot_data_pct_dead_cells_merged["HD_ISS_Stage"] <- HDISS_Stages
    
    stage1_df <- barplot_data_pct_dead_cells_merged[grep("1", barplot_data_pct_dead_cells_merged$HD_ISS_Stage),]
    stage2_df <- barplot_data_pct_dead_cells_merged[grep("2", barplot_data_pct_dead_cells_merged$HD_ISS_Stage),]
    stage3_df <- barplot_data_pct_dead_cells_merged[grep("3", barplot_data_pct_dead_cells_merged$HD_ISS_Stage),]
    
    s1_start <- safe_min_age(stage1_df, input$maxLifetime)
    s2_start <- safe_min_age(stage2_df, input$maxLifetime)
    s3_start <- safe_min_age(stage3_df, input$maxLifetime)
    
    background_data <- data.frame(
      xmin = as.numeric(c(-Inf, s1_start, s2_start, s3_start)),
      xmax = as.numeric(c(s1_start, s2_start, s3_start, input$maxLifetime)),
      ymin = rep(-Inf, 4),
      ymax = rep(Inf, 4),
      fill = factor(c("white", "lightgrey", "grey", "darkgrey"), 
                    levels = c("white", "lightgrey", "grey", "darkgrey"))
    )
    

    # background_data <- data.frame(
    #   xmin = as.numeric(c(-Inf, min(stage1_df$Age), min(stage2_df$Age), min(stage3_df$Age))),
    #   xmax = as.numeric(c(min(stage1_df$Age), min(stage2_df$Age), min(stage3_df$Age), input$maxLifetime)),
    #   ymin = rep(-Inf, 4),
    #   ymax = rep(Inf, 4),
    #   fill = factor(c("white", "lightgrey", "grey", "darkgrey"), levels = c("white", "lightgrey", "grey", "darkgrey"))
    # )
    
    # Run Function on healthy cells
    barplot_data_pct_healthy_cells <- get_healthy_cells_fun(plot_data)
    barplot_data_pct_healthy_cells_untreated <- get_healthy_cells_fun(plot_data_untreated)
    barplot_data_pct_healthy_cells_merged <- merge(barplot_data_pct_healthy_cells,barplot_data_pct_healthy_cells_untreated, by = "Age")
    colnames(barplot_data_pct_healthy_cells_merged) <- c("Age","HealthyTreated","HealthyUntreated")
    
    linePlot_data <- merge(barplot_data_pct_healthy_cells_merged, barplot_data_pct_dead_cells_merged, by = "Age")
    
    # Combine stages into a data frame with new colors
    convertToTreated_fun <- function(x){
      selected_line <- linePlot_data[linePlot_data$Age == x,]
      pctHealthUntreatedCells <- selected_line[3]
      selected_line_2 <- linePlot_data[linePlot_data$HealthyTreated >= pctHealthUntreatedCells[1,1],]
      age <- max(selected_line_2$Age) + 1
      return(age)
    }
    
    background_data_treated <- data.frame(
      xmin = as.numeric(c(-Inf, convertToTreated_fun(s1_start), convertToTreated_fun(s2_start), convertToTreated_fun(s3_start))),
      xmax = as.numeric(c(convertToTreated_fun(s1_start), convertToTreated_fun(s2_start), convertToTreated_fun(s3_start), input$maxLifetime)),
      ymin = rep(-Inf, 4),
      ymax = rep(Inf, 4),
      fill = factor(c("white", "lightgrey", "grey", "darkgrey"), levels = c("white", "lightgrey", "grey", "darkgrey"))
    )
    
    # Update progress bar
    updateProgressBar(session, id = "pbRun", value = 70, title = "Rendering output plots...")
    Sys.sleep(0.3)
    
    #plot_data_list <- list(plot_data,plot_data_untreated,background_data,background_data_treated)
    # Return all needed objects as a list
    plot_data_list <- list(
      plot_data = plot_data,
      plot_data_untreated = plot_data_untreated,
      background_data = background_data,
      background_data_treated = background_data_treated,
      stage1_df = stage1_df,
      stage2_df = stage2_df,
      stage3_df = stage3_df,
      barplot_data_pct_dead_cells_merged = barplot_data_pct_dead_cells_merged,
      barplot_data_pct_healthy_cells_merged = barplot_data_pct_healthy_cells_merged,
      therapeuticBenefitYears = therapeuticBenefitYears,
      linePlot_data = linePlot_data)
  
    return(plot_data_list)
  })
  

  #######################################
  # Create and Render Output Plots
  #######################################
  
  
  # Create a CAP score function
  # CAP = Age x (CAG-L)+ / K
  # (CAG-L)+ is a function that equals x when x >= 0 and 0 otherwise.
  
  CAG_L_fun <- function(x) {
    if (x >= 0) {
      return(x)
    } else {
      return(0)
    }
  }
  
  CAP_Score_Fun <- function(CurrentAge, repLen){
    L <- 30 # an estimate of the lower limit of CAG expansion values for which toxicity occurs
    K <- 6.49
    CAP_score <- CurrentAge * CAG_L_fun(repLen - L) / K
    return(list(CAP_score, K))
  }
  
  GetPredAgeOfOnset <- function(repLen){
    CAP_Score_Lookup_Table <- sapply(c(1:100), function(x){
      CAP_Score_Fun(x, repLen)[[1]]
    })
    
    Predicted_AgeOfOnset <- length(CAP_Score_Lookup_Table[CAP_Score_Lookup_Table <= 100])
    return(Predicted_AgeOfOnset)
  }
  

  ### Create animated plot of cell loss ###
  output$animationPlot <- renderImage({
    plot_data <- simulation_data()[[1]]
    plot_data_untreated <- simulation_data()[[2]]
    #barplot_data_pct_dead_cells_merged <- simulation_data()[[8]]
    #background_data <- simulation_data()[[3]]
    #background_data_treated <- simulation_data()[[4]]
    
    # Isolate inputs to prevent gif from reloading with any change
    repeat_len <- isolate(input$germlineHDRepeatLength)
    highlight <- isolate(input$highlightTransducedCells)
    max_life <- isolate(input$maxLifetime)
    
    if (highlight != TRUE) {
      # Create the plot
      p <- ggplot(plot_data, aes(x = RepeatLength, y = CellID)) +
        geom_point(aes(color = RepeatLength), alpha = 0.7) +
        geom_vline(xintercept = 150, color = "red", linetype = "dashed") +
        scale_color_gradient2(low = "forestgreen", mid = "red", high = "red", midpoint = 200) +
        labs(title = "Distribution of Repeat Lengths", 
             subtitle = "Year: {frame_time}. CAP Score: {round(CAP_Score_Fun(frame_time, repeat_len)[[1]], 2)}",          x = "Repeat Length", 
          y = "") +
        theme_minimal() +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        coord_cartesian(xlim = c(0, 285)) +
        transition_time(Year) +
        ease_aes('linear')
    } else {
      p <- ggplot(plot_data, aes(x = RepeatLength, y = CellID)) +
        geom_point(aes(color = RepeatLength), alpha = 0.7) +
        geom_point(data = subset(plot_data, Transduction == "transduced"),
                   aes(x = RepeatLength, y = CellID),
                   color = "violet", shape = 21, size = 3, fill = NA) +
        geom_vline(xintercept = 150, color = "red", linetype = "dashed") +
        scale_color_gradient2(low = "forestgreen", mid = "red", high = "red", midpoint = 200) +
        labs(title = "Distribution of Repeat Lengths", 
             subtitle = "Year: {frame_time}. CAP Score: {round(CAP_Score_Fun(frame_time, repeat_len)[[1]], 2)}",  
             x = "Repeat Length", 
             y = "") +
        theme_minimal() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        transition_time(Year) +
        ease_aes('linear') +
        coord_cartesian(xlim = c(0, 285))
    }
    
    
    # Save the cell animation as a GIF
    anim_file <- tempfile(fileext = '.gif')
    anim_save(anim_file, animate(p, nframes = max_life, fps = 4))
    
    isolate({ renderState$animationDone <- TRUE })
    
    # Return a list containing the filename
    list(src = anim_file,
         contentType = 'image/gif',
         width = 800,
         height = 600,
         alt = "Animated plot")
  }, deleteFile = TRUE)


  ### Collect information about dead cells (>300 CAGs) ###
  # Define function
  get_dead_cells_fun <- function(z){
    levels <- levels(as.factor(z$Year))
    barplot_data_pct_dead_cells <- sapply(levels, function(x) {
      y <- as.numeric(x)
      year_subset <- z[z$Year == y,]
      pctDead <- (sum(year_subset$RepeatLength >= 300, na.rm = TRUE) / length(year_subset$RepeatLength)) * 100
      return(c(y, pctDead))
    } )
    barplot_data_pct_dead_cells <- barplot_data_pct_dead_cells
    barplot_data_pct_dead_cells <- t(as.data.frame(barplot_data_pct_dead_cells))
    barplot_data_pct_dead_cells <- as.data.frame(barplot_data_pct_dead_cells)
    colnames(barplot_data_pct_dead_cells) <- c("Age","pct_Dead_Cells")
    return(barplot_data_pct_dead_cells)
  }

  # Render plot
  output$deadCellsPlot <- renderPlot({

    # Make plot_data dataframes generated in simulation dataframes availible ###
    plot_data <- simulation_data()[[1]]
    plot_data_untreated <- simulation_data()[[2]]

    # Run Function on dead cells
    barplot_data_pct_dead_cells <- get_dead_cells_fun(plot_data)
    #barplot_data_pct_dead_cells_untreated <- get_dead_cells_fun(plot_data_untreated)
    #barplot_data_pct_dead_cells_merged <- merge(barplot_data_pct_dead_cells,barplot_data_pct_dead_cells_untreated, by = "Age")
    #colnames(barplot_data_pct_dead_cells_merged) <- c("Age","DeadTreated","DeadUntreated")

 
    dead_plot <- ggplot(barplot_data_pct_dead_cells, aes(x = Age, y = pct_Dead_Cells)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_cowplot() +
      ylim(0, 100) +
      labs(title = "Percentage of Dead Cells", x = "Age", y = "% Dead Cells")
    
    return(dead_plot)
  })
  

  ### Collect information about healthy cells (<150 CAGs) ###
  # Define function
  get_healthy_cells_fun <- function(z) {
    levels <- levels(as.factor(z$Year))
    barplot_data_pct_healthy_cells <- sapply(levels, function(x) {
      y <- as.numeric(x)
      year_subset <- z[z$Year == y,]
      pcthealthy <- (sum(year_subset$RepeatLength <= input$HealthyCellCAG_Threshold, na.rm = TRUE) / length(year_subset$RepeatLength)) * 100
      return(c(y, pcthealthy))
    } )
    barplot_data_pct_healthy_cells <- barplot_data_pct_healthy_cells
    barplot_data_pct_healthy_cells <- t(as.data.frame(barplot_data_pct_healthy_cells))
    barplot_data_pct_healthy_cells <- as.data.frame(barplot_data_pct_healthy_cells)
    colnames(barplot_data_pct_healthy_cells) <- c("Age","pct_healthy_Cells")
    return(barplot_data_pct_healthy_cells)
  }


  ### Create histogram of healthy cells ###
  output$healthyCellsPlot <- renderPlot({
    # Make plot_data dataframes generated in simulation dataframes availible ###
    plot_data <- simulation_data()[[1]]
    plot_data_untreated <- simulation_data()[[2]]
    #barplot_data_pct_dead_cells_merged <- simulation_data()[[8]]

    barplot_data_pct_healthy_cells <- get_healthy_cells_fun(plot_data)
    barplot_data_pct_healthy_cells_untreated <- get_healthy_cells_fun(plot_data_untreated)
    barplot_data_pct_healthy_cells_merged <- merge(barplot_data_pct_healthy_cells,barplot_data_pct_healthy_cells_untreated, by = "Age")
    colnames(barplot_data_pct_healthy_cells_merged) <- c("Age","HealthyTreated","HealthyUntreated")
    
    healthy_plot <- ggplot(barplot_data_pct_healthy_cells, aes(x = Age, y = pct_healthy_Cells)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_cowplot() +
      ylim(0, 100) +
      labs(title = paste0("Percentage of Cells Under ", input$HealthyCellCAG_Threshold, " Repeats"), 
           x = "Age", y = paste0("% Cells Under ", input$HealthyCellCAG_Threshold, " Repeats"))
    

    
    isolate({ renderState$healthyDone <- TRUE })
    
    # Update progress bar
    updateProgressBar(session, id = "pbRun", value = 80, title = "Rendering output plots...")
    
    # Store the plot before returning it
    stored_plots$healthyPlot <- healthy_plot
    
    return(healthy_plot)
  })
  
  ### Create a line plot of from both treated and untreated datasets ###
  output$linePlot <- renderPlot({
    repeat_len <- isolate(input$germlineHDRepeatLength)
    intervention_age <- isolate(input$age_at_intervention)
    max_life <- isolate(input$maxLifetime)
    
    plot_data <- simulation_data()[[1]]
    plot_data_untreated <- simulation_data()[[2]]
    background_data <- simulation_data()[[3]]
    background_data_treated <- simulation_data()[[4]]
    stage1_df <- simulation_data()[[5]]
    stage2_df <- simulation_data()[[6]]
    stage3_df <- simulation_data()[[7]]
    #barplot_data_pct_dead_cells_merged <- simulation_data()[[8]]
    #barplot_data_pct_healthy_cells_merged <- simulation_data()[[9]]
    therapeuticBenefitYears <- simulation_data()[[10]]

    # # Run Function on dead cells
    barplot_data_pct_dead_cells <- get_dead_cells_fun(plot_data)
    barplot_data_pct_dead_cells_untreated <- get_dead_cells_fun(plot_data_untreated)
    barplot_data_pct_dead_cells_merged <- merge(barplot_data_pct_dead_cells,barplot_data_pct_dead_cells_untreated, by = "Age")
    colnames(barplot_data_pct_dead_cells_merged) <- c("Age","DeadTreated","DeadUntreated")
     
    # Run Function on healthy cells
    barplot_data_pct_healthy_cells <- get_healthy_cells_fun(plot_data)
    barplot_data_pct_healthy_cells_untreated <- get_healthy_cells_fun(plot_data_untreated)
    barplot_data_pct_healthy_cells_merged <- merge(barplot_data_pct_healthy_cells,barplot_data_pct_healthy_cells_untreated, by = "Age")
    colnames(barplot_data_pct_healthy_cells_merged) <- c("Age","HealthyTreated","HealthyUntreated")

    # Create a line plot with two lines 
    linePlot_data <- merge(barplot_data_pct_healthy_cells_merged, barplot_data_pct_dead_cells_merged, by = "Age")
    
    # Calculate predicted age of onset
    Predicted_AgeOfOnset <- GetPredAgeOfOnset(repeat_len)

    # Update progress bar
    updateProgressBar(session, id = "pbRun", value = 90, title = "Rendering output plots...")
    
    # Combine stages into a data frame with new colors
    convertToTreated_fun <- function(x){
      if (!is.finite(x) || x >= max_life) return(max_life)
      selected_line <- linePlot_data[linePlot_data$Age == x,]
      if (nrow(selected_line) == 0) return(max_life)
      pctHealthUntreatedCells <- selected_line[3]
      selected_line_2 <- linePlot_data[linePlot_data$HealthyTreated >= pctHealthUntreatedCells[1,1],]
      if (nrow(selected_line_2) == 0) return(max_life)
      age <- max(selected_line_2$Age) + 1
      return(age)
    }
    
    s1_start <- safe_min_age(stage1_df, max_life)
    s2_start <- safe_min_age(stage2_df, max_life)
    s3_start <- safe_min_age(stage3_df, max_life)
    
    background_data_treated <- data.frame(
      xmin = as.numeric(c(-Inf, convertToTreated_fun(s1_start), convertToTreated_fun(s2_start), convertToTreated_fun(s3_start))),
      xmax = as.numeric(c(convertToTreated_fun(s1_start), convertToTreated_fun(s2_start), convertToTreated_fun(s3_start), max_life)),
      ymin = rep(-Inf, 4),
      ymax = rep(Inf, 4),
      fill = factor(c("white", "lightgrey", "grey", "darkgrey"), levels = c("white", "lightgrey", "grey", "darkgrey"))
    )
    
    line_plot <- ggplot(linePlot_data, aes(x = Age)) +
      geom_rect(data = background_data_treated, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), alpha = 0.2, inherit.aes = FALSE) +
      geom_line(aes(y = HealthyTreated, color = "healthy_treated"), linewidth = 1.2) +
      geom_line(aes(y = HealthyUntreated, color = "healthy_untreated"), linewidth = 1.2) +
      geom_vline(xintercept = intervention_age, color = "purple", linetype = "dashed") +
      annotate("text", x = intervention_age, y = 50, label = "Gene Therapy Administered", 
               vjust = -1, hjust = 1, color = "purple", angle = 90) +
      geom_vline(xintercept = convertToTreated_fun(Predicted_AgeOfOnset), color = "purple", linetype = "dashed") +
      annotate("text", x = convertToTreated_fun(Predicted_AgeOfOnset), y = 50, label = "Motor Symptom Onset (CAP=100)", 
               vjust = -1, hjust = 0.3, color = "purple", angle = 90) +
      geom_vline(xintercept = convertToTreated_fun(Predicted_AgeOfOnset-14), color = "purple", linetype = "dashed") +
      annotate("text", x = convertToTreated_fun(Predicted_AgeOfOnset-14), y = 50, label = "Typical Volumetric Change Detection", 
               vjust = -1, hjust = 1, color = "purple", angle = 90) +
      annotate("text", x=background_data_treated$xmin[2], y=-10, label=paste0(round(background_data_treated$xmin[2],1)), 
               vjust=1.5, hjust=0.5, color="red", size=3) +
      annotate("text", x=background_data_treated$xmin[3], y=-10, label=paste0(round(background_data_treated$xmin[3],1)), 
               vjust=1.5, hjust=0.5, color="red", size=3) +
      scale_color_manual(
        values = c("healthy_treated" = "purple", 
                   "healthy_untreated" = "grey50"),
        labels = c("healthy_treated" = paste0("Healthy Cells (<", input$HealthyCellCAG_Threshold, ") treated"),
                   "healthy_untreated" = paste0("Healthy Cells (<", input$HealthyCellCAG_Threshold, ") untreated"))
      ) +
      scale_fill_manual(name="HD-ISS Stages", values=c("white"="white", "lightgrey"="lightpink", 
                                                       "grey"="violet", "darkgrey"="purple"),
                        labels=c("Stage 0", "Stage 1", "Stage 2", "Stage 3")) +
      labs(title="Percentage of Cells Over Time",
           subtitle=paste0("Therapeutic Benefit (Years Delayed Onset) = ", therapeuticBenefitYears),
           x="Age",
           y="Percentage of cells",
           color="Cell Type") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))

    isolate({ renderState$linePlotDone <- TRUE })
    
    # Store the line_plot before returning it
    stored_plots$linePlot <- line_plot
    
    return(line_plot)
  })
  
  # ============================================================================
  # Download Handlers
  # ============================================================================
  
  # --- Animation GIF ---
  output$download_animation_gif <- downloadHandler(
    filename = function() { paste0(get_file_prefix(), "_animation.gif") },
    content = function(file) {
      plot_data <- simulation_data()$plot_data
      repeat_len <- input$germlineHDRepeatLength
      highlight <- input$highlightTransducedCells
      max_life <- input$maxLifetime
      
      if (highlight != TRUE) {
        p <- ggplot(plot_data, aes(x = RepeatLength, y = CellID)) +
          geom_point(aes(color = RepeatLength), alpha = 0.7) +
          geom_vline(xintercept = input$HealthyCellCAG_Threshold, color = "red", linetype = "dashed") +
          scale_color_gradient2(low = "forestgreen", mid = "red", high = "red", midpoint = 200) +
          labs(title = "Distribution of Repeat Lengths",
               subtitle = "Year: {frame_time}",
               x = "Repeat Length", y = "") +
          theme_minimal() +
          theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
          coord_cartesian(xlim = c(0, 285)) +
          transition_time(Year) + ease_aes('linear')
      } else {
        p <- ggplot(plot_data, aes(x = RepeatLength, y = CellID)) +
          geom_point(aes(color = RepeatLength), alpha = 0.7) +
          geom_point(data = subset(plot_data, Transduction == "transduced"),
                     aes(x = RepeatLength, y = CellID),
                     color = "violet", shape = 21, size = 3, fill = NA) +
          geom_vline(xintercept = input$HealthyCellCAG_Threshold, color = "red", linetype = "dashed") +
          scale_color_gradient2(low = "forestgreen", mid = "red", high = "red", midpoint = 200) +
          labs(title = "Distribution of Repeat Lengths",
               subtitle = "Year: {frame_time}",
               x = "Repeat Length", y = "") +
          theme_minimal() +
          theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
          coord_cartesian(xlim = c(0, 285)) +
          transition_time(Year) + ease_aes('linear')
      }
      
      anim_save(file, animate(p, nframes = max_life, fps = 4))
    }
  )
  
  # --- Animation underlying data (CSV) ---
  output$download_animation_csv <- downloadHandler(
    filename = function() { paste0(get_file_prefix(), "_cell_data.csv") },
    content = function(file) {
      write.csv(simulation_data()$plot_data, file, row.names = FALSE)
    }
  )
  
  # --- Line Plot PNG ---
  output$download_linePlot_png <- downloadHandler(
    filename = function() { paste0(get_file_prefix(), "_linePlot.png") },
    content = function(file) {
      ggsave(file, stored_plots$linePlot, width = 10, height = 7, dpi = 300)
    }
  )
  
  # --- Line Plot data (CSV) ---
  output$download_linePlot_csv <- downloadHandler(
    filename = function() { paste0(get_file_prefix(), "_linePlot_data.csv") },
    content = function(file) {
      write.csv(simulation_data()$linePlot_data, file, row.names = FALSE)
    }
  )
  
  # --- Healthy Cells Plot PNG ---
  output$download_healthyPlot_png <- downloadHandler(
    filename = function() { paste0(get_file_prefix(), "_healthyCells.png") },
    content = function(file) {
      ggsave(file, stored_plots$healthyPlot, width = 10, height = 7, dpi = 300)
    }
  )
  
  # --- Healthy Cells data (CSV) ---
  output$download_healthyPlot_csv <- downloadHandler(
    filename = function() { paste0(get_file_prefix(), "_healthyCells_data.csv") },
    content = function(file) {
      write.csv(simulation_data()$barplot_data_pct_healthy_cells_merged, file, row.names = FALSE)
    }
  )
  
  # --- Run Configuration file ---
  output$download_config <- downloadHandler(
    filename = function() { paste0(get_file_prefix(), "_config.csv") },
    content = function(file) {
      config <- data.frame(
        Parameter = c(
          "Germline CAG Repeat Length",
          "Age at Intervention",
          "Max Lifetime",
          "MSH3 Knockdown Percentage",
          "Transduction Rate",
          "Number of MSNs Simulated",
          "Healthy Cell CAG Threshold",
          "Contraction Events per Cycle",
          "Expansion Events per Cycle",
          "Expansion:Contraction Ratio",
          "Phase 1 Threshold (T1)",
          "Phase 1 Rate (r1)",
          "Phase 2 Threshold (T2)",
          "Phase 2 Rate (r2)",
          "Cell Death Threshold (CAG)",
          "Timestamp"
        ),
        Value = c(
          input$germlineHDRepeatLength,
          input$age_at_intervention,
          input$maxLifetime,
          input$therapeuticModifier,
          input$transductionRate,
          input$nMSNs_CaudatePutamen,
          input$HealthyCellCAG_Threshold,
          input$n_contractions,
          input$n_expansions,
          paste0(input$n_expansions, ":", input$n_contractions),
          33.5,
          0.035,
          55,
          0.576,
          300,
          format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        )
      )
      write.csv(config, file, row.names = FALSE)
    }
  )
  # #### Line plot showing shift in volumetric change ####
  # output$line_plot_VolChangeShift <- renderPlot({
  #   plot_data <- simulation_data()[[1]]
  #   plot_data_untreated <- simulation_data()[[2]]
  #   background_data <- simulation_data()[[3]]
  #   background_data_treated <- simulation_data()[[4]]
  #   
  #   # Run Function on dead cells
  #   barplot_data_pct_dead_cells <- get_dead_cells_fun(plot_data)
  #   barplot_data_pct_dead_cells_untreated <- get_dead_cells_fun(plot_data_untreated)
  #   barplot_data_pct_dead_cells_merged <- merge(barplot_data_pct_dead_cells,barplot_data_pct_dead_cells_untreated, by = "Age")
  #   colnames(barplot_data_pct_dead_cells_merged) <- c("Age","DeadTreated","DeadUntreated")
  #   
  #   # Run Function on healthy cells
  #   barplot_data_pct_healthy_cells <- get_healthy_cells_fun(plot_data)
  #   barplot_data_pct_healthy_cells_untreated <- get_healthy_cells_fun(plot_data_untreated)
  #   barplot_data_pct_healthy_cells_merged <- merge(barplot_data_pct_healthy_cells,barplot_data_pct_healthy_cells_untreated, by = "Age")
  #   colnames(barplot_data_pct_healthy_cells_merged) <- c("Age","HealthyTreated","HealthyUntreated")
  #   
  #   # Create a line plot with two lines 
  #   linePlot_data <- merge(barplot_data_pct_healthy_cells_merged, barplot_data_pct_dead_cells_merged, by = "Age")
  #   
  #   # Calculate predicted age of onset
  #   Predicted_AgeOfOnset <- GetPredAgeOfOnset(input$germlineHDRepeatLength)
  #   
  #   line_plot_VolChangeShift <- ggplot(linePlot_data, aes(x = Age)) +
  #     #geom_rect(data = background_data_treated, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), alpha = 0.2, inherit.aes = FALSE) +
  #     geom_line(aes(y = HealthyTreated, color = "Healthy Cells (<150) treated")) +
  #     geom_line(aes(y = HealthyUntreated, color = "Healthy Cells (<150) untreated")) +
  #     #geom_line(aes(y = DeadTreated, color = "Dead Cells (>300) treated")) +
  #     #geom_line(aes(y = DeadUntreated, color = "Dead Cells (>300) untreated")) +
  #     geom_vline(xintercept = age_at_intervention, color = "purple", linetype = "dashed") +
  #     annotate("text", x = age_at_intervention, y = 50, label = "Gene Therapy Administered", 
  #              vjust = -1, hjust = 1, color = "purple", angle = 90) +
  #     #geom_vline(xintercept = convertToTreated_fun(Predicted_AgeOfOnset), color = "purple", linetype = "dashed") +
  #     #annotate("text", x = convertToTreated_fun(Predicted_AgeOfOnset), y = 50, label = "Motor Symptom Onset (CAP=100)", 
  #     #         vjust = -1, hjust = -0.1, color = "purple", angle = 90) +
  #     geom_vline(xintercept = Predicted_AgeOfOnset - 14, color = "grey", linetype = "dashed") +
  #     annotate("text", x = Predicted_AgeOfOnset - 14, y = 50, label = "Typical Volumetric Change Detection Untreated", 
  #              vjust = -1, hjust = 1, color = "grey", angle = 90) +
  #     geom_vline(xintercept = convertToTreated_fun(Predicted_AgeOfOnset-14), color = "purple", linetype = "dashed") +
  #     annotate("text", x = convertToTreated_fun(Predicted_AgeOfOnset-14), y = 50, label = "Typical Volumetric Change Detection Treated", 
  #              vjust = -1, hjust = 1, color = "purple", angle = 90) +
  #     scale_color_manual(values = c("Healthy Cells (<150) treated" = "purple", 
  #                                   "Healthy Cells (<150) untreated" = "grey",
  #                                   "Dead Cells (>300) treated" = "violet", 
  #                                   "Dead Cells (>300) untreated" = "lightgrey")) +
  #     #scale_fill_manual(name="HD-ISS Stages", values=c("white"="white", "lightgrey"="lightpink", 
  #     #                                                 "grey"="violet", "darkgrey"="purple"),
  #     #                  labels=c("Stage 0", 
  #     #                           "Stage 1",
  #     #                           "Stage 2",
  #     #                           "Stage 3")) +
  #     labs(title="Percentage of Cells Over Time",
  #          subtitle=paste0("Delay in onset of volumetric changes = ", as.numeric(convertToTreated_fun(Predicted_AgeOfOnset-14) - (Predicted_AgeOfOnset - 14)), " years"),
  #          x="Age",
  #          y="Percentage of cells",
  #          color="Cell Type") +
  #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #           panel.background = element_blank(), axis.line = element_line(colour = "black"))
  #   
  #   return(line_plot_VolChangeShift)
  # })
}