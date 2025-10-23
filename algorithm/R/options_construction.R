library(jsonlite)

#set up the list of options so we don't clutter things up too much
{
  options = list()
  options$time = list()
  options$file = list()
  options$db = list()
  options$intervals = list()
  options$sofa = list()
  
}

#the time options
{
  options$time$ptz = "CET"
  options$time$po = "1970-01-01"
  options$time$format = "%Y-%m-%d %H:%M:%OS"
  options$time$non_existent = "roll-forward"
  options$time$ambiguous = "earliest" 
}

#the sofa options
{
  options$sofa$global = list()
  {
    options$sofa$global$locf_time = 24*3600
    options$sofa$global$nocb_time = 0
    options$sofa$global$use_saps_labs = TRUE
    options$sofa$global$infinitesimal = 1e-6
    
  }
  options$sofa$limits = list()
  {
    options$sofa$limits$values = "./algorithm/limits/values.csv"
    options$sofa$limits$codebook = "./algorithm/limits/limits_241014.csv"
  }
  
  options$sofa$antro = list()
  
  {
    options$sofa$antro$weight_switch_cutoff = 120
    options$sofa$antro$heigth_switch_cutoff = 150
    
    options$sofa$antro$weight_decimal_shift_lower_bound = 20
    options$sofa$antro$weight_decimal_shift_upper_bound = 400
    
    options$sofa$antro$heigth_decimal_shift_lower_bound = 100
    options$sofa$antro$heigth_decimal_shift_upper_bound = 270
    
    options$sofa$antro$weight_outlier_function = "hampel"
    options$sofa$antro$weight_outlier_before = 3600*24*7
    options$sofa$antro$weight_outlier_after = 3600*24*7
    options$sofa$antro$weight_outlier_k = 5
    
    
    options$sofa$antro$sofa_weight_aggregation_function = "first"
    #241015 : probably better to have a very old own weight than to impute from population...
    options$sofa$antro$sofa_weight_locf = 24*3600*365*10
    options$sofa$antro$sofa_weight_nocb = 24*3600*365*10
    
    options$sofa$antro$sofa_weight_round = 1L
    
    options$sofa$antro$popdata = data.table(SexMale = c(NA, 0L, 1L), mean_height = c(172.8236, 163.9311, 177.5733), mean_weight = c(79.66054, 71.18388, 84.16502), mean_bmi = c(26.67086, 26.48859, 26.69170))
  }
  
  
  
  options$sofa$renal = list()
  options$sofa$renal$assume_NA_is_zero = TRUE
  options$sofa$renal$min_obs_time_for_NA_zero = 12*3600
  options$sofa$renal$new_catheter_remove_outliers = FALSE
  options$sofa$renal$max_break_until_new_catheter = 7*24*3600
  options$sofa$renal$new_catheter_outlier_timerange = 6*3600
  options$sofa$renal$new_catheter_outlier_method = "tukey"
  options$sofa$renal$new_catheter_outlier_k = 1.5
  options$sofa$renal$urine_per_kg_per_hr_1st_percentile = 0.04988096
  options$sofa$renal$urine_per_kg_per_hr_2_5th_percentile = 0.16306
  
  options$sofa$gcs = list()
  options$sofa$gcs$maximum_component_locf_time = 24*3600
  options$sofa$gcs$maximum_imv_break = 6*3600
  options$sofa$gcs$mminimum_imv_duration = 3600
  options$sofa$gcs$sed_outliers_timerange = 9900
  options$sofa$gcs$sed_outliers_method = "tukey"
  options$sofa$gcs$sed_outliers_k = 1.5
  options$sofa$gcs$sedative_gcs_duration = 24*3600
  options$sofa$gcs$minimum_sedation_duration = 3600
  options$sofa$gcs$maximum_preintubation_gcs_locf_time = -6*3600
  options$sofa$gcs$maximum_post_intubation_nocb_time = 24*3600
  options$sofa$gcs$maximum_imvsed_locf_time = 24*3600
  options$sofa$gcs$carry_forward_sedated_gcs_values = FALSE
  
  options$sofa$resp = list()
  options$sofa$resp$max_gap_before_new_interval = 3900
  options$sofa$resp$added_time_from_last_mode_update = 0
  options$sofa$resp$pfi_rollback_time = 300
  options$sofa$resp$impute_pao2_from_sao2 = TRUE
  options$sofa$resp$impute_normal_abg_values = TRUE
  options$sofa$resp$pfi_outliers_spo2_window = 1800
  options$sofa$resp$pfi_outliers_fio2_after_time = 600
  options$sofa$resp$pfi_model_path = "./algorithm/xgboost/bg_xgb_9feats.xgboost"
  options$sofa$resp$pfi_outlier_cutoff_prob = 0.7279894 #found by maximizing the f0.5 of the pr curve of the model
  options$sofa$resp$use_spo2_ratio = TRUE
  options$sofa$resp$spo2_outliers_timerange = 1800
  options$sofa$resp$spo2_outliers_method = "tukey"
  options$sofa$resp$spo2_outliers_k = 3
  options$sofa$resp$maximum_rs_break = 900
  options$sofa$resp$minimum_rs_duration = 900
  options$sofa$resp$spo2_matched_time_window = 600
  options$sofa$resp$time_weight_duration = 24*3600
  options$sofa$resp$bgvars = data.table(VariableID = c(20000300, 20001200, 20000200, 20001300, 24000714, 20004200, 24000290, 20000800),
                                        Parameter = c("pH", "pCO2", "pO2", "BE", "SBE", "HCO3", "stHCO3", "SO2"),
                                        SampleType = "a",
                                        Min = c(6, 0, 0, -40, -40, 0, 0, 0),
                                        Max = c(9, 30, 101, 40, 40, 80, 80, 100))
  
  options$sofa$resp$labactions = rbind(
    data.table(Statement = c(210,292,722,875,1010), action = "Delete"),
    data.table(Statement = c(83:86, 580, 855), action = "Ok"),
    data.table(Statement = c(211, 452, 476, 581, 664), action = "null")
  )
  
  
  options$sofa$circ = list()
  options$sofa$circ$map_model_path = "./algorithm/xgboost/xgb_10_feats_220518.xgboost"
  options$sofa$circ$map_model_features_path = "./algorithm/xgboost/features_xgb_10_feats_220518.csv"
  options$sofa$circ$map_outlier_cutoff_prob = 0.26359290
  options$sofa$circ$nibp_outliers_timerange = 1800
  options$sofa$circ$nibp_outliers_method = "hampel"
  options$sofa$circ$nibp_outliers_k = 1.4826*2
  options$sofa$circ$vd_instance_mergetime = 300
  options$sofa$circ$vd_min_inf_time = 3600
  options$sofa$circ$vd_cof_coefficient = 2.25
}

#the file options
{
  options$file$chunk_size = 5e6
  options$file$labres_chunk_size = 1e6
  options$file$max_chunks_to_read = 200
  options$pharma$concentration_signif = 2
  
  options$file$sep = ";"
  options$file$dec = ","
  options$file$encoding = "Latin-1"
  
  options$file$paths = list()
  #options$file$paths$sir_master = file.choose()
  #options$file$paths$sir_mortality = file.choose()
  #options$file$paths$sir_diagnosis = file.choose()
  #options$file$paths$labres = file.choose()
  options$file$paths$sir_master = "P:\\GLUC-ICU\\data\\raw\\2022-11\\SIR_Masterfile_221108.xlsx"
  options$file$paths$sir_mortality = "P:\\GLUC-ICU\\data\\raw\\2022-11\\SIR_Mortalitet_221108.xlsx"
  options$file$paths$sir_diagnosis = "P:\\GLUC-ICU\\data\\raw\\2022-05\\SIR_Diagnoser_210607.xlsx"
  options$file$paths$sir_saps = "P:\\GLUC-ICU\\data\\raw\\2022-05\\SIR_Saps_210607.xlsx"
  options$file$paths$sir_atgarder = "P:\\GLUC-ICU\\data\\raw\\2022-05\\SIR_?tg?rder_210607.xlsx"
  options$file$paths$sir_komplikationer1 = "P:\\GLUC-ICU\\data\\raw\\2022-05\\SIR_Komplikationer_2010_2011.xlsx"
  options$file$paths$sir_komplikationer2 = "P:\\GLUC-ICU\\data\\raw\\2022-05\\SIR_Komplikationer_2012_2021.xlsx"
  
  #options$file$paths$labres = "P:\\GLUC-ICU\\data\\processed\\2022-10\\Gluc_Icu_CCC_LabRes_221122_semicolonfix.csv"
  options$file$paths$labres = "P:\\GLUC-ICU\\data\\processed\\2024-09\\Gluc_Icu_CCC_LabRes_240905_semicolonfix.csv"
  
  options$file$paths$karda_tpk = "P:\\GLUC-ICU\\data\\raw\\2024-10\\Gluc-Icu TC Tpk 241024.csv"
  options$file$paths$karda_bili = "P:\\GLUC-ICU\\data\\raw\\2024-10\\Gluc-Icu TC Bilirubin 241024.csv"
  options$file$paths$karda_krea = "P:\\GLUC-ICU\\data\\raw\\2024-10\\Gluc-Icu TC Kreatinin 241024.csv"
  
  
  options$file$paths$obsrec = "P:\\GLUC-ICU\\data\\raw\\2023-05\\Gluc_Icu_CCC_Observrec_SOFA_230522.csv"
  options$file$paths$karda_weight =  "P:\\GLUC-ICU\\data\\processed\\2024-06\\Gluc_Icu_TC_Langd_Vikt_240605.csv"
  
  options$file$paths$obsrec2 = "P:\\GLUC-ICU\\data\\processed\\2024-01\\Gluc_Icu_CCC_Observrec_Linn_H_1_240102.csv"
  options$file$paths$obsrec3 = "P:\\GLUC-ICU\\data\\raw\\2024-01\\Gluc_Icu_CCC_Observrec_Linn_H_2_240102.csv"
  options$file$paths$obsrec4 = "P:\\GLUC-ICU\\data\\raw\\2024-01\\Gluc_Icu_CCC_Observrec_Linn_H_3_240102.csv"
  
  options$file$paths$monvals2 = "P:\\GLUC-ICU\\data\\raw\\2024-01\\Gluc_Icu_CCC_MonVals_Linn_H_240102.csv"
  options$file$paths$monvals3 = "P:\\GLUC-ICU\\data\\processed\\2025-09\\Gluc-Icu - MonVals Hjartfrekvens 250915 - rowfix.csv"
  
  options$file$paths$pharma_with_pharmaid = "P:\\GLUC-ICU\\data\\raw\\2023-03\\Gluc_Icu_CCC_Pharma_221108.csv"
  options$file$paths$pharma_comp_key = "P:\\GLUC-ICU\\data\\raw\\2022-06\\Info om PharmaID och CopmID 20220607.xlsx"
  options$file$paths$pharma = "P:\\GLUC-ICU\\data\\raw\\2022-05\\Gluc_Icu_CCC_Pharma.csv"
  options$file$paths$daily_orders  = "P:\\GLUC-ICU\\data\\raw\\2023-05\\GluC_Icu_CCC_DailyOrderSchedule_230508.csv"
  
  options$file$redefinitions_directory = "./algorithm/variable_redefinitions/"
}

#interval options
{
  options$intervals$non_overlap_duration = 1e-5
  options$intervals$index_column_name = "SOFA_N"
  options$intervals$interval_start_name = "SOFA_start"
  options$intervals$interval_end_name = "SOFA_end"
}


#db options
{
  options$db$remove_duplicates = TRUE
  options$db$vacuum_db = TRUE
  
  options$db$sql_directory = "./algorithm/sql/"
  
  options$db$db_path = "P:\\GLUC-ICU\\users\\johanh\\sql\\db_v9_updated_pharma_with_pharmaID.sqlite"
  
}
{
  options_json = jsonlite::toJSON(options, pretty = TRUE, auto_unbox = TRUE)
  con = file ("options.json")
  writeLines(options_json, con)
  close(con)
  rm(con)
  rm(options_json)
}
