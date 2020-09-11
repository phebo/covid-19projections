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

library(shiny)
library(tidyverse)

#setwd("pol-app")
dfDg <- read_csv("dg-data.csv")
dfGbase <- read_csv("gbase-data.csv")
dfDesc <- read_csv("policy-descriptions.csv")

dfDg <- dfDg %>% mutate(polCode = substr(pol, 1, 2)) %>% arrange(iter, pol, level) %>%
  group_by(iter, pol) %>% mutate(cum = cumsum(value)) %>% ungroup()

dfPol <- dfDg %>% select(polCode, pol, level) %>% distinct() %>% left_join(dfDesc) %>%
  mutate(polDesc = paste(level, polDesc, sep = ". "))

dfPol2 <- dfPol %>% select(pol, polCode, polDesc) %>%
  bind_rows(tibble(pol = unique(dfPol$pol), polCode = unique(dfPol$polCode), polDesc = "0. No/limited policies")) %>%
  arrange(pol, polDesc) %>% group_by(pol, polCode) %>% nest() %>%
  mutate(buttons = pmap(list(pol, polCode, data), function(pol, polCode, data)
    selectInput(polCode, pol, choices = data$polDesc, selected = "0. No/limited policies") ))

ui <- fluidPage(
  titlePanel("Covid spread"),
  sidebarLayout(
    sidebarPanel(dfPol2$buttons),
    mainPanel(
      plotOutput("plot", height = 800))
  )
)

server <- function(input, output) {
  output$plot <- renderPlot({
    dfPol3 <- dfPol2 %>%
      mutate(
        input = map_chr(polCode, ~ input[[.]]),
        level = as.numeric(substr(input,1,1))) %>%
      select(polCode, level)
    dfG <- inner_join(dfDg, dfPol3) %>% bind_rows(tibble(iter = unique(dfDg$iter), cum = 0)) %>%
      group_by(iter) %>% summarize(dg = sum(cum)) %>% ungroup() %>%
      full_join(dfGbase %>% rename(gbase = value, base = name)) %>%
      mutate(value = gbase - dg) %>%
      group_by(base) %>% summarize(estimate = median(value), low = quantile(value, probs=0.025),
                                   high = quantile(value, probs=0.975), prob = sum(value < 0) / n()) %>%
      mutate(base = factor(base, levels = c("low", "medium", "high"), labels = c("Best 10%", "Median", "Worst 10%")))
    dfG %>% ggplot(aes(x = base, y = estimate, ymin = low, ymax = high,
                       fill = ifelse(estimate > 0, "red", "green"))) +
      geom_col() + geom_errorbar(size=1) + ylim(c(-1.2,1.8)) + scale_fill_manual(values = c(red = "red3", green = "green4")) +
      xlab("Geography-specific policy effectiveness") + ylab("Weekly growth rate of infections (estimate and 95% interval)") +
      theme(text = element_text(size=20), legend.position = "none")
  })}

shinyApp(ui = ui, server = server)
