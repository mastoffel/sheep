library(RJDBC)
library(dplyr)
library(dbplyr)

dbname <- "data/db/StKilda_Data.accdb"
driver <- "net.ucanaccess.jdbc.UcanloadDriver"
driverpath <- "data/db/UCanAccess/loader/ucanload.jar"
options <- paste0("jdbc:ucanaccess://", dbname, ";memory=false")

con <- DBI::dbConnect(JDBC(driver, driverpath), options)
# src <- src_dbi(con)

tbls <- dbGetTables(con)
flds <- dbGetFields(con, "Sheep")
Sheep <- dbGetQuery(con, "Select * from Sheep")
Census <- dbGetQuery(con, "Select * from CensusData")
Capture <- dbGetQuery(con, "Select * from CaptureData")

dbDisconnect(con)

library(ggplot2)
Capture %>% 
    sample_frac(1) %>% 
    ggplot(aes(Measurer, Foreleg)) +
    geom_boxplot()


