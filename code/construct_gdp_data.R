if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table)

# GDP data from worldbank https://databank.worldbank.org/source/world-development-indicators#
raw_data <- fread("../data/real_gdp.csv")
raw_data[, `Series Name`:=NULL]
raw_data[, `Series Code`:=NULL]
raw_data <- raw_data[`Country Code` != ""]
setnames(raw_data, old = c("Country Name", "Country Code"), 
         new = c("country", "code"))

# transform to long table
data_long <- melt(raw_data, id.vars = c("country", "code"), 
                  variable.name = "year", value.name = "gdp")

# clean the data
## year
data_long[, year:=as.integer(substr(year, 1, 4))]
## GDP as billion dollors 
data_long[, gdp:=as.numeric(gdp)/10^9]
## sort by country code and year
setorderv(data_long, cols = c("code", "year"))

# write data
fwrite(data_long, "../data/real_gdp_long.csv")