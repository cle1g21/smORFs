install.packages("tidyverse")
library(tidyverse)
install.packages("ggplot2")
library(ggplot2)

data("midwest", package = "ggplot2") 
ggplot(midwest, aes(x=area, y=poptotal))

