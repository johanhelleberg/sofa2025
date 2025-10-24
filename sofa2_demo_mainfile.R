#2025-10-23
#SOFA/SOFA2 demo file
library(data.table)
library(Rcpp)
library(pryr)
library(glue)
rebuild_from_source = TRUE
if(rebuild_from_source){
  source("./algorithm/R/sofa_intervals.R")
  source("./algorithm/R/options_construction.R")
  source("./algorithm/R/sofa_utility_functions.R")
  source("./algorithm/R/sofa_weights_karda.R")
  sourceCpp("./algorithm/src/tp_basic_outliers.cpp")
  sourceCpp("./algorithm/src/tp_basic_interval_match.cpp")
}


#load the master file
#spell_master_key_path = "P:\\GLUC-ICU\\users\\johanh\\sofa2025\\datasets\\spell_master_key_2501023.csv"
spell_master_key_path = "P:\\GLUC-ICU\\users\\johanh\\sofa2025\\datasets\\spell_master_dataset_from_jm_stat_file_with_spell_id_251023.csv"
spell_master_key = fread(spell_master_key_path, 
                         select = c("StudyID", 
                                    "spell_id", 
                                    "PatientID", 
                                    "VtfHuvudId",
                                    "spell_sex",
                                    "icuadm_time", 
                                    "icudis_time", 
                                    "spell_start", 
                                    "spell_stop"))[!is.na(StudyID)]

#convert time stamps to target timezone for the run
shift_dt_times(spell_master_key, target_tz = options$time$ptz)

intdt = unique(copy(spell_master_key[, .(spell_id, spell_start, spell_stop)]))

intdt = intdt[spell_id %in% sample(spell_id, 5)]
setkey(intdt, spell_id, spell_start, spell_stop)
setkey(intdt, spell_id)


intervals = make_intervals(intdt, 
                           options,
                           "sliding", 
                           interval_length = 24*3600, 
                           slide_duration = 3600, 
                           max_aggregate_before = 24*3600)

library(RSQLite)
db = dbConnect(SQLite(), options$db$db_path)


