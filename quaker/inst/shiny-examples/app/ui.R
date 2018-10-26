library(leaflet)
library(shiny)

# Choices for drop-downs
time_period <- c(
  "In The last Hour" = "PAST_HOUR",
  "In The last Day" = "PAST_DAY",
  "In The last Week" = "PAST_WEEK",
  "In The last Month" = "PAST_MONTH"
)

magnitude <- c(
  "Greater than 0" = "all",
  "Greater than 1" = "1.0",
  "Greater than 2.5" = "2.5",
  "Greater than 4.5" = "4.5"
)

navbarPage("Quaker - Seimic Activity Monitor", id="nav",

  tabPanel("Interactive map",
    div(class="outer",
      tags$head(
        # Include  CSS
        includeCSS("styles.css")
      ),

      leafletOutput("map", width = "100%", height = "600px"),

      absolutePanel(id = "controls", class = "panel panel-default", fixed = TRUE,
                    draggable = TRUE, top = 60, left = "auto", right = 20, bottom = "auto",
                    width = 330, height = "auto",

                    h2("Filters"),

                    selectInput("timePeriodId", "Time Period", time_period, selected = "PAST_WEEK"),
                    selectInput("magnitude", "Magnitude", magnitude, selected = "2.5")
                    )
      )
    ,tags$div(id="cite",
               'Data sourced from ', tags$em('United States Geological Survey')
      )
  )
  ,

  tabPanel("Raw Data explorer",
    DT::dataTableOutput("rawData")
  ),
  conditionalPanel("false", icon("crosshair"))
)
