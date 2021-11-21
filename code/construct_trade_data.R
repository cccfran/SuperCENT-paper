if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table)

# Trade data from https://correlatesofwar.org/data-sets/bilateral-trade
cow_data <- fread("../data/COW_Trade_4.0/Dyadic_COW_4.0.csv")
cow_data <- cow_data[year >= 1984]

# transform to long table: importer1 imports from importer2
trade_1 <- cow_data[,.(year, to = importer1, from = importer2, import = smoothflow1)]
trade_2 <- cow_data[,.(year, from = importer1, to = importer2, import = smoothflow2)]
trade <- rbind(trade_1, trade_2)

# merge with gdp data to get iso3 code
gdp_data <- fread("../data/real_gdp_long.csv")
country_list <- unique(gdp_data[,.(country, code)])
country_list[grep("United States", country), country := "United States"]
country_list[grep("KOR", code), country := "South Korea"]

trade <- merge(trade, country_list, by.x = "to", by.y = "country")
setnames(trade, old = "code", new = "code1")
trade <- merge(trade, country_list, by.x = "from", by.y = "country")
setnames(trade, old = "code", new = "code2")
trade <- trade[,.(year, from = code2, to = code1, value = import)]

# HKG from the International Monetary Fund (IMF) Direction of Trade Statistics
# https://data.imf.org/?sk=9D6028D4-F14A-464C-A2F2-59B2CD424B85
hk_trade <- fread("../data/hkg_trade.csv")
trade <- rbindlist(list(trade, hk_trade))

# merge with FX_sub to get a subset
fx <- fread("../data/FX.csv")
fx_country <- unique(fx[,.(code = iso3)])
trade <- merge(trade, fx_country, by.x = "to", by.y = "code")
trade <- merge(trade, fx_country, by.x = "from", by.y = "code")

# write trade data
fwrite(trade, "../data/trade_data_sub.csv")


