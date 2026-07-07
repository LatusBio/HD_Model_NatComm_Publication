required_packages <- c("shiny", "ggplot2", "cowplot", "viridis",
                       "R6", "reshape2", "gifski", "av",
                       "RCurl", "httr", "remotes", "shinyWidgets","gganimate")

# Install any missing packages
installed <- rownames(installed.packages())
to_install <- setdiff(required_packages, installed)
if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cran.rstudio.com/")
}

#library(remotes)
#remotes::install_github("thomasp85/gganimate")

# Now load libraries
library(shiny)
library(ggplot2)
library(cowplot)
library(viridis)
library(R6)
library(reshape2)
library(readr)
library(RCurl)
library(httr)
library(gganimate)
library(shinyWidgets)

# Disable SSL verification (for servers with bad CA config)
set_config(config(ssl_verifypeer = 0L))

ui <- fluidPage(
    tags$head(
    HTML("
      <!-- Google tag (gtag.js) -->
      <script async src='https://www.googletagmanager.com/gtag/js?id=G-L0DN986BNX'></script>
      <script>
        window.dataLayer = window.dataLayer || [];
        function gtag(){dataLayer.push(arguments);}
        gtag('js', new Date());

        gtag('config', 'G-L0DN986BNX');
      </script>
    "),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('disableButton', function(message) {
        document.getElementById(message.id).disabled = true;
      });
      Shiny.addCustomMessageHandler('enableButton', function(message) {
        document.getElementById(message.id).disabled = false;
      });
  "))
  ),
  titlePanel("Latus Bio - Huntingon's Disease Therapeutic Model"),
  
  # Description box
  fluidRow(
    column(12,
           div(style = "background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 8px; padding: 20px; margin-bottom: 20px;",
               p(style = "font-size: 25px; line-height: 1.6; color: #333:",
                 "Welcome:"),
               p(style = "font-size: 15px; line-height: 1.6; color: #333;",
                 "This application simulates somatic CAG repeat expansion in human medium spiny neurons (MSNs) ",
                 "of the caudate and putamen using a two-phase linear model of CAG repeat expansion reported by ",
                 tags$a(href = "https://doi.org/10.1016/j.cell.2024.11.038", "Handsaker et al. (2025)", target = "_blank"), ". ",
                 "Here we build on this foundational work by enabling users to simulate the impact of a therapeutic intervention targeting somatic expansion.",
                 "Users can use this web-app to model therapeutic impact on repeat expansion dynamics, ",
                 "predict the health of striatal medium spiny neurons over time, and predict the delay to onset of HD milestones like the onset of motor symptoms.",
                 "A complete description of this model can be found in our publication ", tags$a(href = "https://www.biorxiv.org/content/10.64898/2026.01.06.697909v1", "Simpson and Ranum et al. (2026)", target = "_blank"), ". "
               ),
               p(style = "font-size: 25px; line-height: 1.6; color: #333:",
                 "Instructions:"),
               p(style = "font-size: 14px; line-height: 1.6; color: #555;",
                 "Adjust the parameters in the toolbar below to set the germline CAG repeat length, age at therapeutic intervention, ",
                 "percentage of cells targeted, level of somatic instability lowering, and more, then click ", tags$b("Run Simulation"), " to generate results. ",
                 "Outputs include an animated visualization of repeat expansion, treated vs. untreated cell health trajectories, ",
                 "and HD-ISS stage progression. All plots and underlying data are downloadable."
               ),
               p(style = "font-size: 10px; line-height: 1.0; color: #333:",
                 "This web app was developed and is maintained by Latus Bio. Address correspondance to paul.ranum@latusbio.com. Code is available on", tags$a(href = "https://github.com/LatusBio/HD_Model_NatComm_Publication/", "Github", target = "_blank"),". "
               ),
           )
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      
      sliderInput("germlineHDRepeatLength", "Germline CAG Repeat Length:", value = 40, min = 36, max = 56, step = 1),
      helpText("The inherited CAG repeat length on the mutant HTT allele. Typical HD range is 36–55."),
      
      br(),
      
      sliderInput("age_at_intervention", "Age at Intervention:", value = 30, min = 0, max = 100),
      helpText("The age (in years) at which the gene therapy is administered."),
      
      br(),
      
      sliderInput("maxLifetime", "Max Lifetime:", value = 85, min = 50, max = 85),
      helpText("The maximum age simulated. The model tracks cell expansion dynamics from birth to this age."),
      
      br(),
      
      sliderInput("therapeuticModifier", "Percentage of MSH3 Knockdown:", value = 0.5, min = 0.1, max = 1, step = 0.1),
      helpText("The fraction by which MSH3-driven expansion rate is reduced in transduced cells. A value of 0.5 means a 50% reduction in expansion rate."),
      
      br(),
      
      sliderInput("transductionRate", "Percentage of Cells Targeted:", value = 0.5, min = 0.1, max = 1, step = 0.1),
      helpText("The probability that each neuron is successfully targeted by the therapeutic intervention."),
      
      br(),
      
      sliderInput("nMSNs_CaudatePutamen", "Number of MSNs to simulate:", value = 1000, min = 100, max = 3000, step = 100),
      helpText("The number of neurons to be simulated, increasing this value will increase wait time."),
      
      br(),
      
      sliderInput("HealthyCellCAG_Threshold", "Maximum CAG threshold for 'healthy' cells", value = 150, min = 50, max = 250, step = 5),
      helpText("Cells with repeat lengths below this threshold are counted as 'healthy'. Based on Handsaker et al., transcriptional dysregulation begins around 150 CAG."),
      
      br(),
      
      sliderInput("n_contractions", "Contraction events per cycle (-1 CAG):", 
                  value = 1, min = 0, max = 10, step = 1),
      sliderInput("n_expansions", "Expansion events per cycle (+1 CAG):", 
                  value = 6, min = 1, max = 10, step = 1),
      helpText("Controls the expansion-to-contraction bias. Default 6:1 means each mutation event has a 6/7 chance of adding +1 CAG and a 1/7 chance of -1 CAG."),
      
      br(),
      
      #checkboxInput("therapy_on", "Therapy On", value = TRUE),
      checkboxInput("highlightTransducedCells", "Highlight Transduced Cells", value = TRUE),
      helpText("When checked, transduced cells are outlined in violet in the animation."),
      
      
      actionButton("run", "Run Simulation"),
      br(),
      br(),
      div(style = "margin-bottom: -10px; margin-top: 5px;",
          progressBar(id = "pbRun", value = 0, display_pct = TRUE, status = "info")
      )
    ),
    
    mainPanel(
      #plotOutput("animationPlot"),
      imageOutput("animationPlot", height = "600px"),
      fluidRow(
        column(3, downloadButton("download_animation_gif", "Download GIF")),
        column(3, downloadButton("download_animation_csv", "Download Data (CSV)"))
      ),
      br(),
      br(),
      br(),
      br(),
      br(),
      br(),
      plotOutput("linePlot"),
      fluidRow(
        column(3, downloadButton("download_linePlot_png", "Download PNG")),
        column(3, downloadButton("download_linePlot_csv", "Download Data (CSV)"))
      ),
      br(),
      #plotOutput("line_plot_VolChangeShift"),
      br(),
      plotOutput("healthyCellsPlot"),
      fluidRow(
        column(3, downloadButton("download_healthyPlot_png", "Download PNG")),
        column(3, downloadButton("download_healthyPlot_csv", "Download Data (CSV)"))
      ),
      br(),
      hr(),
      downloadButton("download_config", "Download Run Configuration")
    )
  )
)
