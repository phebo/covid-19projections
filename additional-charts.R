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

library(tidyverse)
library(RColorBrewer)

df <- read_csv("model-out.csv") %>%
  mutate(Var = factor(Var, levels  = c("Death", "Reported case", "Infection")))
vGeo2 <- read_csv("input/geo-sel.csv") %>% pull(Geo)

dfLastDate <- df %>% arrange(Geo, Var, Date) %>%
  filter(!is.na(New), !is.na(NewEst), Var == "Reported case", substr(Geo,1,5) != "China",
         substr(Geo, 1, 2) != "US" | Geo == "US - 12 states total") %>%
  group_by(Geo) %>% summarize(Date = last(Date))
dfCo <- bind_rows(
  dfLastDate %>% inner_join(df %>% filter(Var == "Reported case") %>% select(Geo, Date, Var, Cum)),
  dfLastDate %>% inner_join(df %>% filter(Var == "Infection") %>% select(Geo, Date, Var, Cum = CumEst))
  ) %>%
  mutate(Var = fct_recode(Var, `Reported cases` = "Reported case", `Estimated infections` = "Infection"))
dfCo <- dfCo %>% filter(Geo %in% (dfCo %>% filter(Var == "Estimated infections") %>% top_n(8, Cum) %>% pull(Geo)))

dfCo %>% ggplot(aes(x = Var, y = Cum, fill = reorder(Geo, Cum))) + geom_col() +
  scale_y_continuous(labels = scales::comma) + 
  scale_fill_manual(values = brewer.pal(8,"Set1")) +
  ylab(element_blank()) + xlab(element_blank()) +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/fig-est-rep.png", width=4, height=5)
dfCo %>% group_by(Var) %>% summarize(sum(Cum))
