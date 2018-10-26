library(leaflet)
library(shiny)
library(dplyr)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

  geoJson <- reactive({
    x <- get_seismic_data(timeFrame = input$timePeriodId, minMagnitude = input$magnitude)
    x
  })

  output$map <- renderLeaflet({
      x <- geoJson()
      plot(x)
  })

  output$rawData <- DT::renderDataTable({

    x <- flatten_to_table(geoJson())

    #get dataframe in correct format for rendering
    df = x %>% mutate(
      time = as.POSIXct(time/1000, origin="1970-01-01"),
      place  = paste('<a href="',url,'">',place,'</a>', sep=""),
      `tsunami risk` = ifelse(tsunami == 0,'No','Yes')
    ) %>% select(place, time, longtitude, latitude, depth, magnitude, `tsunami risk`)

    DT::datatable(df, options = list(pageLength=20), escape = FALSE)
  })

})
