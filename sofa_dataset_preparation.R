#2025-10-06
#Dataset preparation and imputation
library(data.table)
score_file_path = "P:\\GLUC-ICU\\data\\processed\\2025-09\\sofa_scores\\sofa_highres_15_min_2025_09_08_08_55_41.csv"
#sofadt = fread(score_file_path, nrows = 1000)
sofadt = fread(score_file_path, select = c(1,2,6,64,66,67,69,70,72,73,74))

#use the Ellis' inversion when applicable
sofadt[, `:=` (SOFA_RESP = fcoalesce(SOFA_RESP, SOFA_RESP_SpO2_imputed),
               SOFA_RESP_SpO2_imputed = NULL)]

#locf first, then assume normal
sofacols = c("SOFA_RESP","SOFA_GCS","SOFA_CIRC","SOFA_RENAL","SOFA_LIVER","SOFA_COAG")
for (cn in sofacols){
  sofadt[, (cn) := data.table::nafill(.SD[[cn]], "locf"), PatientID]
  sofadt[, (cn) := data.table::nafill(.SD[[cn]], "const", 0), PatientID]
}

#compute the sofa score
sofadt[, SOFA_SCORE := SOFA_RESP+SOFA_GCS+SOFA_CIRC+SOFA_RENAL+SOFA_LIVER+SOFA_COAG]

#clean up a bit
{
  rm(score_file_path)
  rm(sofacols)
  rm(cn)
}

#shift the colnames to indicate non-raw variables from the dataset, i.e. post-imputations
{
  colnames(sofadt)[colnames(sofadt) == "AdmissionTime"] = "adm_time"
  colnames(sofadt)[colnames(sofadt) == "SOFA_end"] = "sofa_time"
  colnames(sofadt)[colnames(sofadt) == "SOFA_RESP"] = "resp"
  colnames(sofadt)[colnames(sofadt) == "SOFA_GCS"] = "cns"
  colnames(sofadt)[colnames(sofadt) == "SOFA_CIRC"] = "cv"
  colnames(sofadt)[colnames(sofadt) == "SOFA_RENAL"] = "renal"
  colnames(sofadt)[colnames(sofadt) == "SOFA_LIVER"] = "liver"
  colnames(sofadt)[colnames(sofadt) == "SOFA_COAG"] = "coag"
  colnames(sofadt)[colnames(sofadt) == "SOFA_SCORE"] = "sofa"
}

#find the hourly scores
sofadt[, dt := as.numeric(sofa_time - adm_time, units = "hours")]
sofadt[, icu_hour := floor(dt)]
sofadt[, dtdiff := dt - icu_hour]
sofadt = sofadt[, .SD[which.min(dtdiff)], .(PatientID, icu_hour)]
sofadt[, `:=`(dt = NULL,
              dtdiff = NULL)]


#get the outcomes too
master_datapath = "P:\\GLUC-ICU\\users\\navid\\HMT data\\data\\csv\\Demographics_January2025.csv"
master = fread(master_datapath)

sofadt[master[,.(PatientID = patientid, deathdate)], deathdate := deathdate, on = "PatientID"]

#clean up a bit
{
  rm(master_datapath)
  rm(master)
}

sofadt[, time_to_death := as.numeric(deathdate - adm_time, units = "days")]
sofadt[, death_at_30 := 0L]
sofadt[time_to_death <=30, death_at_30 := 1L]
sofadt[, deathdate := NULL]
