

library(lavaan)
library(devtools)
library(tabuSearch)
library(GA)


setwd("/Users/RJacobucci/Documents/Github/autoSEM")
setwd("C:/Users/jacobucc/Documents/Github/autoSEM")


#setwd("./regsem")
document()

setwd("..")
install("autoSEM")
#document()
#document("C:/Users/jacobucc/Documents/GitHub/regsem")
detach("package:autoSEM", unload=TRUE)
library(autoSEM)
