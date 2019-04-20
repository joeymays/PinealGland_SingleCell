setwd("D:/Data Analysis/PinealGland_SingleCell")
library(knitr)
source("./Scripts/MetadataInit.R")
source("./Scripts/SetupDependencies.R")

purl(input = "./Day1.Rmd", output = "./notebook_extracts/Day1.R")
purl(input = "./Day2.Rmd", output = "./notebook_extracts/Day2.R")
purl(input = "./Day1n2_NoDoublets.Rmd", output = "./notebook_extracts/Day1n2_nd.R")

source("./notebook_extracts/Day1.R")
source("./notebook_extracts/Day2.R")
source("./notebook_extracts/Day1n2_nd.R")
