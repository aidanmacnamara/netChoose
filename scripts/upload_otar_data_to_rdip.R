
# example of writing data to rdip, using custom function

# Load libraries
require(tidyverse)
require(futile.logger)
require(glue)
require(sparklyr)
require(implyr)

# load function from genedomainr repo
library(devtools)
source("R/upload_tmptable.R")

# define table name and target location
file_name = "v2g/otar_v2g.txt"
table_name = "otar_v2g"
db = "am673712"

# initialize connections
# for uploading large tables, is it better to use spark, which is based on hive and more slow but more resiliant to errors
sc = spark_connect(master="yarn-client", spark_home="/opt/cloudera/parcels/SPARK2/lib/spark2", version="2.3.1")

# you will also need to initialize an impala connection - the upload_table function
# will need to run some manual commands (INVALIDATE METADATA) in impala 
imp = src_impala(odbc::odbc(), dsn="impaladsn")

upload_table(file_name, table_name, db, sc=sc, imp=imp)


# MERGE TABLES ------------------------------------------------------------

# connect to the tables in rdip

sc = spark_connect(master="yarn-client", spark_home="/opt/cloudera/parcels/SPARK2/lib/spark2", version="2.3.1")

sc %>% tbl_change_db("am673712")
src_tbls(sc)
v2g <- tbl(sc, "otar_v2g")

system.time({
  v2g_head = v2g %>% head() %>% collect()
})


