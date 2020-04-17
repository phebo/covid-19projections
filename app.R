# Copyright (C) 2020, Phebo Wibbens, Wesley Koo, and Anita McGahan

#   This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# For this app to work, you first need to create projection estimates using the script run-model.R

library(shiny)
library(tidyverse)
library(RColorBrewer)

df <- read_csv("output/model-out.csv") %>%
  mutate(Var = factor(Var, levels  = c("Death", "Reported case", "Infection")))
vGeo2 <- read_csv("input/geo-sel.csv") %>% pull(Geo)

Date1 <- as.Date("2020-01-15")
Date2 <- max(df$Date)

ui <- fluidPage(
  titlePanel("Covid projections and reported events"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("dates", "Dates", Date1, Date2, value = c(Date1+17, Date2 - 30), timeFormat = "%d %b"),
      checkboxInput("cum", "Cumulative"),
      checkboxInput("log", "Logarithmic scale", value = T),
      selectInput("geo", "Select geography", choices = c("All", vGeo2))
    ),
    mainPanel(plotOutput("plot", height = 800))
  )
)

server <- function(input, output) {
  output$plot <- renderPlot({
    if(input$geo == "All") {
      df2 <- df
      dotsize <- 0.25; linesize <- 0.25
    } else {
      df2 <- df %>% filter(Geo == input$geo) 
      dotsize <- 1.5; linesize <- 0.5
    }
    if(input$cum == T) {
      ggOut <- df2 %>%
        filter(between(Date, input$dates[1], input$dates[2]), Geo %in% vGeo2) %>% mutate(Geo = factor(Geo, levels = vGeo2)) %>%
        ggplot(aes(x=Date, y=CumEst, color=Var)) +
        geom_line(aes(y=CumLow), linetype = 2, size = linesize) + geom_line(aes(y=CumHigh), linetype = 2, size = linesize) +
        geom_point(aes(y=Cum), size=dotsize) + geom_line(size = linesize) 
    } else {
      ggOut <- df2 %>% 
        filter(between(Date, input$dates[1], input$dates[2]), Geo %in% vGeo2) %>% mutate(Geo = factor(Geo, levels = vGeo2)) %>%
        ggplot(aes(x=Date, y=NewEst, color=Var)) +
        geom_line(aes(y=NewLow), linetype = 2, size = linesize) + geom_line(aes(y=NewHigh), linetype = 2, size = linesize) +
        geom_point(aes(y=New), size=dotsize) + geom_line(size = linesize) 
    }
    ggOut <- ggOut + facet_wrap(~Geo, ncol = 5) +
      scale_x_date(labels = function(x) substr(format(x, format="%d %b"),2,10), date_breaks = "1 month") +
      scale_color_manual(values = brewer.pal(3,"Set1"), guide = guide_legend(reverse = TRUE)) +
      ggtitle("Dot = Reported data, Line = Model estimate, Dash = 95% interval") + 
      xlab(element_blank()) + ylab(element_blank()) +
      theme(legend.title=element_blank(), legend.position="top", axis.text.x = element_text(angle = 45, hjust = 1))
      
    if(input$log) {
      ggOut <- ggOut + scale_y_continuous(labels = scales::comma, trans="log10", limits = if(input$cum) c(1,1e7) else c(1,3e5))
    } else {
      ggOut <- ggOut + ylim(if(input$cum) c(0,3e6) else c(0,2e5))
    }
    ggOut
  })}

shinyApp(ui = ui, server = server)
