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

load("fit-model-data.RData")

Date1 <- as.Date("2020-01-15")
Date2 <- max(dfOut2$Date)
Date3 <- max(dfOut2 %>% filter(!is.na(New)) %>% pull(Date))
dfOut3 <- dfOut2 %>% mutate(New = ifelse(Date <= Date3 & is.na(New) & Var != "Infection", 0, New))
  
ui <- fluidPage(
  titlePanel("Covid projections and reported events"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("dates", "Dates", Date1, Date2, value = c(Date3 - 30, Date3), timeFormat = "%d %b"),
      selectInput("geo", "Select geography", choices = vGeo2)
    ),
    mainPanel(tableOutput("table"))
  )
)

input <- list(geo = "Austria", dates = c(Date3-30,Date3))

server <- function(input, output) {
  output$table <- renderTable({
    dfOutSum <- dfOut3 %>% filter(Geo == input$geo, Date >= input$dates[1], Date <= input$dates[2]) %>%
      group_by(Var) %>% summarize_at(c("New", "NewEst", "NewLow", "NewHigh"), sum) %>%
      rename(Reported = New, Estimate = NewEst, `Lower bound` = NewLow, `Upper bound` = NewHigh)
    dfOutSum2 <- dfOut3 %>% filter(Geo == input$geo, Date == input$dates[2]) %>%
      select(Var, Reported = New, Estimate = NewEst, `Lower bound` = NewLow, `Upper bound` = NewHigh)
    bind_rows(
      bind_cols(
        tibble(Statistic = paste0("Total infections (", format(input$dates[1], format="%b %d"),"-", format(input$dates[2], format="%b %d"),")" )),
        dfOutSum %>% filter(Var == "Reported case") %>% select(Reported),
        dfOutSum %>% filter(Var == "Infection") %>% select(Estimate:`Upper bound`)),
      dfOutSum %>% filter(Var == "Death") %>% select(-Var) %>%
        mutate(Statistic = paste0("Total deaths (", format(input$dates[1], format="%b %d"),"-", format(input$dates[2], format="%b %d"),")" )),
      bind_cols(
        tibble(Statistic = paste0("New infections per day (", format(input$dates[2], format="%b %d"),")" )),
        dfOutSum2 %>% filter(Var == "Reported case") %>% select(Reported),
        dfOutSum2 %>% filter(Var == "Infection") %>% select(Estimate:`Upper bound`)),
      dfOutSum2 %>% filter(Var == "Death") %>% select(-Var) %>%
        mutate(Statistic = paste0("New deaths per day (", format(input$dates[2], format="%b %d"),")" )),
      dfGeo %>% filter(Geo == input$geo, Var == "g") %>% select(Estimate, `Lower bound` = Low, `Upper bound` = High) %>%
        mutate_all(function(x) 100 * x) %>%
        mutate(Statistic = paste0("Current percentage growth rate of new infections (", format(Date3, format="%b %d"),")")),
      dfGeo2 %>% filter(Geo == input$geo) %>% mutate(Statistic = ifelse(Var == "DaysTo100", "Number of days to 100 new infections / day", "Number of days to 1 new infection / day")) %>%
        select(Statistic, Estimate, `Lower bound` = Low, `Upper bound` = High)
    ) %>%
      mutate_at(vars(Reported), function(x) format(x, big.mark = ',', digits = 0, nsmall = 0, scientific = F)) %>%
      mutate_at(vars(Estimate:`Upper bound`), function(x) format(x, big.mark = ',', digits = 0, nsmall = 1, scientific = F))
  }, 
  align = "r"
  )}

shinyApp(ui = ui, server = server)
