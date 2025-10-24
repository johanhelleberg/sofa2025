#the function that actually calculates SOFA scores.
calculate_sofa = function(interval_data, db, master, options, verbose = TRUE, key_column = NULL){
  if(TRUE){
    loglist = list()
    
    #dummy data
    if (TRUE){
      interval_data = copy(intervals)
      master = copy(spell_master_key)
      key_column = NULL
      #interval_data = intervallist
    }
    
    multiple_intervals = NULL
    
    
    
    #check what type of data was passed on - list of data.table or list of data.table
    if (any(grepl("data.table", attr(interval_data, "class")))){
      multiple_intervals = FALSE
    } else {
      if (all(unlist(lapply(intervallist, function(x) any(grepl("data.table", attr(x, "class"))))))){
        multiple_intervals = TRUE
      } else {
        return("interval_data must either be a single data.table or a list of data.tables")
      }
    }
    
    #prepare the intervals and determine which patientIDs we need to query
    pids = integer(0)
    
    if(!multiple_intervals){
      if(is.null(key_column)) key_column = colnames(interval_data)[1]
      interval_data = copy(interval_data)
      setkeyv(interval_data, c(key_column, 
                               options$interval$interval_start_name, 
                               options$interval$interval_end_name))
      
      
      pids = c(pids, master[spell_id %in% interval_data$spell_id, unique(PatientID)])
      
    } else {
      for (n in names(inteval_data)){
        if(is.null(key_column)) key_column = colnames(inteval_data[[i]])[1]
        interval_data[[n]] = copy(inteval_data[[n]])
        setkeyv(interval_data[[n]], c(key_column, 
                                      options$interval$interval_start_name, 
                                      options$interval$interval_end_name))
        
        pids = c(pids, master[spell_id %in% interval_data[[n]]$spell_id, unique(PatientID)])
      }
    }
    
    pids = sort(unique(pids))
    intkey = master[PatientID %in% pids, .(StudyID, spell_id, PatientID)]
    
    loglist[[length(loglist)+1]] = data.table(event = "prepared intervals and pid selection",
                                              time = Sys.time(),
                                              memory_usage = pryr::mem_used())
    
    
    
    

    
    
    
    if(verbose){
      endtime = Sys.time()
      print(sofa_echo_time("Ticking away, the moments that make up a dull day", starttime, prevtime, endtime))
      prevtime = endtime
    }
    
    
    #Initialize intervals
    if(TRUE){
      intervals = copy(interval_data)
      
      setkey(intervals, key(interval_data))
      urval = unique(intervals$PatientID)#[1:100]
      #for pasting into various queries
      urval_sql = collapse_to_sql(urval)
      #store the names for future use
      vartypes = dtReadTable(db, "vartypes")
      limits = fread(options$sofa$limits$codebook, sep = options$file$sep, dec = options$file$dec)[, .(VariableID, Min, Max, plausible_min_sign, plausible_min_numeric = as.numeric(plausible_min_numeric), plausible_max_sign, plausible_max_numeric= as.numeric(plausible_max_numeric))]
      
      #>= is straight-forward
      limits[!is.na(plausible_min_numeric) & plausible_min_sign == ">=", Min := plausible_min_numeric]
      limits[!is.na(plausible_min_numeric) & plausible_min_sign == ">", Min := plausible_min_numeric + options$sofa$global$infinitesimal]
      
      
      limits[!is.na(plausible_max_numeric) & plausible_max_sign == "<=", Max := plausible_max_numeric]
      limits[!is.na(plausible_max_numeric) & plausible_max_sign == "<", Max := plausible_max_numeric - options$sofa$global$infinitesimal]
      
      limits = limits[, .(VariableID, Min, Max)]
      
    }
    
    logdt = rbind(logdt, 
                  data.table(event = "initialization", time = Sys.time(), memory_usage = mem_used()))
    
    #master
    #sql = glue_sql("SELECT * FROM hv_master WHERE PatientID IN ({urval*})", .con = db)
    
    
    #get the weights
    antrodata = get_weights(db, options, urval, master, vartypes, verbose)
    logdt = rbind(logdt, 
                  data.table(event = "weights", time = Sys.time(), memory_usage = mem_used()))
    
    #get the resp settings
    ventilator_periods = calculate_ventilator_periods(intervals, db, options, limits, urval_sql, verbose)
    logdt = rbind(logdt, 
                  data.table(event = "ventilator_periods", time = Sys.time(), memory_usage = mem_used()))
    
    
    #the MAPs will be used in both circulation and respiration
    
    sql = glue_sql("SELECT PatientID, DateTime, NumericValue FROM monvals WHERE VariableID = 110 AND PatientID IN({urval*})", .con = db)
    map_values = dtGetQuery(db, sql)
    
    # map_values = dtGetQuery(db, paste0("SELECT PatientID, DateTime, NumericValue FROM monvals WHERE VariableID IN ",
    #                                    collapse_to_sql(110),
    #                                    " AND PatientID IN ",
    #                                    urval_sql,
    #                                    ""))
    #map_values[, DateTime := as.POSIXct(DateTime, tz = options$time$ptz, origin = options$time$po)]
    setkey(map_values, PatientID, DateTime)
    
    logdt = rbind(logdt, 
                  data.table(event = "map_values", time = Sys.time(), memory_usage = mem_used()))
    
    #Extract the SAPS-labs from the SIR_SAPS-table
    if (TRUE){
      prevtime = Sys.time()
      #get the saps for these patientID
      
      sql = glue_sql("SELECT PatientID, StudyID, FirstAdmissionDateTime as AdmissionDateTime FROM 'patients'  WHERE PatientID IN ({urval*})", .con = db)
      patients = dtGetQuery(db, sql)
      setkey(patients, StudyID, PatientID, AdmissionDateTime)
      
      sql = glue_sql("SELECT PatientID, StudyID, VtfHuvudID, AdmissionDateTime FROM 'sir_master'  WHERE PatientID IN ({urval*})", .con = db)
      sir_master = dtGetQuery(db, sql)
      setkey(sir_master, StudyID, PatientID, AdmissionDateTime)
      #merge
      sir_key = sir_master[patients]
      
      #always fetch GCS, no questions asked
      saps_cols = c("VtfHuvudId", "GCS", "GCSEye", "GCSVerbal")
      #if saps-labs-requested, also load those
      if (options$sofa$global$use_saps_labs){
        saps_cols = c(saps_cols, 
                      "TemperatureMax",
                      "BilirubinMax ", 
                      "CreatinineMax", 
                      "ThrombocytesMin", 
                      "SystolicBloodPressureMin", 
                      "FiO2",
                      "PaO2")
      }
      
      saps_cols_sql = paste(saps_cols, collapse = ",")
      
      sql = paste0("SELECT ", saps_cols_sql, " FROM 'sir_saps'")
      rm(saps_cols)
      
      saps_vars = dtGetQuery(db, sql)
      
      saps_vars = saps_vars[sir_key, on = "VtfHuvudId"]
      
      saps_vars[, `:=` (DateTime = AdmissionDateTime,
                        StudyID = NULL,
                        AdmissionDateTime = NULL,
                        VtfHuvudId = NULL,
                        GCSMotor = GCS - GCSEye - GCSVerbal)]
      
      rm(sir_master)
      rm(sir_key)
      rm(patients)
      
      #convert all to double before the melt
      for (cn in colnames(saps_vars)[!colnames(saps_vars) %in% c("PatientID", "DateTime")]){
        saps_vars[, (cn) := as.numeric(.SD[[cn]])]
      }
      #saps_vars2 = copy(saps_vars)
      
      #prepare saps to be joined with the rest
      saps_vars[, DateTime := as.POSIXct(DateTime, tz = options$time$ptz, origin = options$time$po)]
      saps_vars = melt(saps_vars, id.vars = c("PatientID", "DateTime"), variable.name = "VarName", value.name = "Value")
      saps_vars[, EnterTime := DateTime]
      saps_vars[VarName == "GCS", VariableID := 30000300]
      saps_vars[VarName == "GCSVerbal", VariableID := 10000100]
      saps_vars[VarName == "GCSMotor", VariableID := 10000200]
      saps_vars[VarName == "GCSEye", VariableID := 10000300]
      saps_vars[VarName =="TemperatureMax",VariableID := 0]
      saps_vars[VarName =="BilirubinMax", VariableID := 20004300]
      saps_vars[VarName =="CreatinineMax", VariableID := 20000600]
      saps_vars[VarName =="ThrombocytesMin", VariableID := 20000110]
      saps_vars[VarName =="SystolicBloodPressureMin", VariableID := 100]
      saps_vars[VarName =="FiO2",VariableID := 2000]
      saps_vars[VarName =="PaO2",VariableID := 20000200]
      
      saps_vars[,VarName := NULL]
      
      setcolorder(saps_vars, c("PatientID", "VariableID", "DateTime", "EnterTime", "Value"))
      
      if(verbose){
        endtime = Sys.time()
        print(sofa_echo_time("SAPS-vars fetched from db", starttime, prevtime, endtime))
        prevtime = Sys.time()
      }
    }
    
    logdt = rbind(logdt, 
                  data.table(event = "saps", time = Sys.time(), memory_usage = mem_used()))
    
    
    #coagulation
    sofa_coagulation(intervals, db, options, limits, saps_vars, urval_sql, verbose)
    logdt = rbind(logdt, 
                  data.table(event = "coagulation", time = Sys.time(), memory_usage = mem_used()))
    
    #liver
    sofa_liver(intervals, db, options, limits, saps_vars, urval_sql, verbose)
    logdt = rbind(logdt, 
                  data.table(event = "liver", time = Sys.time(), memory_usage = mem_used()))
    
    #renal
    sofa_renal(intervals, db, options, limits, saps_vars, urval_sql, verbose)
    logdt = rbind(logdt, 
                  data.table(event = "renal", time = Sys.time(), memory_usage = mem_used()))
    
    #resp
    sofa_resp(intervals, db, options, limits, saps_vars, ventilator_periods, map_values, urval_sql, verbose)
    logdt = rbind(logdt, 
                  data.table(event = "respiration", time = Sys.time(), memory_usage = mem_used()))
    
    #gcs
    sofa_cns(intervals, db, options, limits, saps_vars, ventilator_periods, urval_sql, verbose)
    logdt = rbind(logdt, 
                  data.table(event = "cns", time = Sys.time(), memory_usage = mem_used()))
    
    #remove auxillary tables
    rm(saps_vars)
    rm(ventilator_periods)
    
    #circulation
    sofa_circ(intervals, db, options, limits, antrodata, map_values, urval_sql, verbose)
    logdt = rbind(logdt, 
                  data.table(event = "circulation", time = Sys.time(), memory_usage = mem_used()))
    
    if(multiple_intervals){
      for (interval in intervallist){
        # #Liver
        # interval[Bilirubin<20, SOFA_LIVER := 0]
        # interval[Bilirubin %between% c(20,32), SOFA_LIVER := 1]
        # interval[Bilirubin %between% c(33,101), SOFA_LIVER := 2]
        # interval[Bilirubin %between% c(102,204), SOFA_LIVER := 3]
        # interval[Bilirubin>204, SOFA_LIVER := 4]
        # 
        # #Coagulation
        # interval[TPK>=150, SOFA_COAG := 0]
        # interval[TPK<150, SOFA_COAG := 1]
        # interval[TPK<100, SOFA_COAG := 2]
        # interval[TPK<50, SOFA_COAG := 3]
        # interval[TPK<20, SOFA_COAG := 4]
        # 
        
        # if(sofa_global_icuregswe){
        #   intervals = rbind(intervals, lastinterval)
        #   intervals = intervals[order(PatientID, IntervalN)]
        #   intervals = intervals[,IntervalN := fifelse(IntervalN==Inf, shift(IntervalN)+1, IntervalN)]
        #   if(verbose){
        #     endtime = Sys.time()
        #     print(sofa_echo_time("Intervals merged", starttime, prevtime, endtime))
        #     prevtime = endtime
        #   }
        # }
        
        
        #interval[,UrinePer24 := (Urine * 24*3600)/ as.numeric(IntervalEnd - IntervalStart, units = "secs")]
        
        
        # }
        if(TRUE){
          #RESP
          interval[PFI > (400/7.5), SOFA_RESP := 0]
          interval[PFI %between% c(300/7.5,400/7.5), SOFA_RESP := 1]
          interval[PFI %between% c(0*200/7.5,300/7.5), SOFA_RESP := 2]
          # if (sofa_global_icuregswe){
          #   intervals[PFI %between% c(100/7.5,200/7.5), SOFA_RESP := 3]
          #   intervals[PFI < 100/7.5, SOFA_RESP := 4]
          #   
          # } else {
          interval[PFI_rs %between% c(100/7.5, 200/7.5), SOFA_RESP := 3]
          interval[PFI_rs < 100/7.5, SOFA_RESP := 4]
          # }
          #Cutoffs for SpO2/FiO2:
          #1 : < 512
          #2 : < 357
          #3 : < 214
          #4 : < 89
          #From Pandharipande et al. 2009
          if(options$sofa$resp$use_spo2_ratio){
            interval[SpO2FiO2ratio >= 512, SOFA_RESP_SpO2 := 0]
            interval[SpO2FiO2ratio %between% c(357,512), SOFA_RESP_SpO2 := 1]
            interval[SpO2FiO2ratio %between% c(214*0,357), SOFA_RESP_SpO2 := 2]
            # if (sofa_global_icuregswe){
            #   intervals[SpO2FiO2ratio %between% c(89,214), SOFA_RESP_SpO2 := 3]
            #   intervals[SpO2FiO2ratio < 89 , SOFA_RESP_SpO2 := 4]
            # } else {
            interval[SpO2FiO2ratio_rs %between% c(89,214), SOFA_RESP_SpO2 := 3]
            interval[SpO2FiO2ratio_rs < 89 , SOFA_RESP_SpO2 := 4]
            
            # }
            
            interval[PFI_spo2_imputed  > (400/7.5), SOFA_RESP_SpO2_imputed := 0]
            interval[PFI_spo2_imputed  %between% c(300/7.5,400/7.5), SOFA_RESP_SpO2_imputed := 1]
            interval[PFI_spo2_imputed  %between% c(0*200/7.5,300/7.5), SOFA_RESP_SpO2_imputed := 2]
            # if (sofa_global_icuregswe){
            #   intervals[PFI_spo2_imputed  %between% c(100/7.5,200/7.5), SOFA_RESP_SpO2_imputed := 3]
            #   intervals[PFI_spo2_imputed  < 100/7.5, SOFA_RESP_SpO2_imputed := 4]
            #
            # } else {
            interval[PFI_spo2_imputed_rs  %between% c(100/7.5, 200/7.5), SOFA_RESP_SpO2_imputed := 3]
            interval[PFI_spo2_imputed_rs  < 100/7.5, SOFA_RESP_SpO2_imputed := 4]
          }
          
          interval[anyGCS == 15, SOFA_GCS := 0]
          interval[anyGCS %between% c(13,14), SOFA_GCS := 1]
          interval[anyGCS %between% c(10,12), SOFA_GCS := 2]
          interval[anyGCS %between% c(6,9), SOFA_GCS := 3]
          interval[anyGCS < 6, SOFA_GCS := 4]
          
          interval[clinicianGCS == 15, SOFA_GCS_CLIN := 0]
          interval[clinicianGCS %between% c(13,14), SOFA_GCS_CLIN := 1]
          interval[clinicianGCS %between% c(10,12), SOFA_GCS_CLIN := 2]
          interval[clinicianGCS %between% c(6,9), SOFA_GCS_CLIN := 3]
          interval[clinicianGCS < 6, SOFA_GCS_CLIN := 4]
          
          #Circulation
          interval[MAP>70, SOFA_CIRC := 0]
          interval[MAP<=70, SOFA_CIRC := 1]
          interval[DobutaminmcgPerKgPerMin>0, SOFA_CIRC := 2]
          interval[DopaminmcgPerKgPerMin>0 & DopaminmcgPerKgPerMin<=5, SOFA_CIRC := 2]
          
          interval[LambdenNorEpiEquiv>0 & 
                     LambdenNorEpiEquiv<=0.1, SOFA_CIRC := 3]
          interval[DopaminmcgPerKgPerMin>5 & DopaminmcgPerKgPerMin<=15, SOFA_CIRC := 3]
          interval[LambdenNorEpiEquiv>0.1, SOFA_CIRC := 4]
          interval[DopaminmcgPerKgPerMin>15, SOFA_CIRC := 4]
          
          #Renal
          interval[,UrinePer24 := (Urine * 24*3600)/ as.numeric(SOFA_end - SOFA_start, units = "secs")]
          interval[Kreatinin <110, SOFA_RENAL := 0]
          interval[Kreatinin %between% c(110,170), SOFA_RENAL := 1]
          interval[Kreatinin %between% c(171,299), SOFA_RENAL := 2]
          interval[Kreatinin %between% c(300,440), SOFA_RENAL := 3]
          interval[UrinePer24<500, SOFA_RENAL := 3]
          interval[Kreatinin > 440, SOFA_RENAL := 4]
          interval[UrinePer24<200, SOFA_RENAL := 4]
          
          #KDIGO
          interval[krea_baseline_ratio < 1.5, KDIGO := 0]
          interval[, krea_increase := Kreatinin - (Kreatinin / krea_baseline_ratio)]
          interval[krea_increase >= 26.5, KDIGO := 1]
          interval[, krea_increase := NULL]
          
          interval[krea_baseline_ratio >= 1.5, KDIGO := 1]
          interval[perkgperhr6 <= 0.5, KDIGO := 1]
          interval[krea_baseline_ratio >= 2.0, KDIGO := 2]
          interval[perkgperhr12 <= 0.5, KDIGO := 2]
          interval[krea_baseline_ratio >= 3.0, KDIGO := 3]
          interval[perkgperhr24 <= 0.3, KDIGO := 3]
          interval[anuria_duration >= 12*3600, KDIGO := 3]
          interval[Kreatinin >= 353.6, KDIGO := 3]
          interval[crrt_duration > 0, KDIGO := 3]
          
          
          #Liver
          interval[Bilirubin<20, SOFA_LIVER := 0]
          interval[Bilirubin %between% c(20,32), SOFA_LIVER := 1]
          interval[Bilirubin %between% c(33,101), SOFA_LIVER := 2]
          interval[Bilirubin %between% c(102,204), SOFA_LIVER := 3]
          interval[Bilirubin>204, SOFA_LIVER := 4]
          
          #Coagulation
          interval[TPK>=150, SOFA_COAG := 0]
          interval[TPK<150, SOFA_COAG := 1]
          interval[TPK<100, SOFA_COAG := 2]
          interval[TPK<50, SOFA_COAG := 3]
          interval[TPK<20, SOFA_COAG := 4]
          
          
          
          
        }
        
      }
    } else {
      # if(sofa_global_icuregswe){
      #   intervals = rbind(intervals, lastinterval)
      #   intervals = intervals[order(PatientID, IntervalN)]
      #   intervals = intervals[,IntervalN := fifelse(IntervalN==Inf, shift(IntervalN)+1, IntervalN)]
      #   if(verbose){
      #     endtime = Sys.time()
      #     print(sofa_echo_time("Intervals merged", starttime, prevtime, endtime))
      #     prevtime = endtime
      #   }
      # }
      
      
      #RESP
      intervals[PFI > (400/7.5), SOFA_RESP := 0]
      intervals[PFI %between% c(300/7.5,400/7.5), SOFA_RESP := 1]
      intervals[PFI %between% c(0*200/7.5,300/7.5), SOFA_RESP := 2]
      # if (sofa_global_icuregswe){
      #   intervals[PFI %between% c(100/7.5,200/7.5), SOFA_RESP := 3]
      #   intervals[PFI < 100/7.5, SOFA_RESP := 4]
      #   
      # } else {
      intervals[PFI_rs %between% c(100/7.5, 200/7.5), SOFA_RESP := 3]
      intervals[PFI_rs < 100/7.5, SOFA_RESP := 4]
      # }
      #Cutoffs for SpO2/FiO2:
      #1 : < 512
      #2 : < 357
      #3 : < 214
      #4 : < 89
      #From Pandharipande et al. 2009
      if(options$sofa$resp$use_spo2_ratio){
        intervals[SpO2FiO2ratio >= 512, SOFA_RESP_SpO2 := 0]
        intervals[SpO2FiO2ratio %between% c(357,512), SOFA_RESP_SpO2 := 1]
        intervals[SpO2FiO2ratio %between% c(214*0,357), SOFA_RESP_SpO2 := 2]
        # if (sofa_global_icuregswe){
        #   intervals[SpO2FiO2ratio %between% c(89,214), SOFA_RESP_SpO2 := 3]
        #   intervals[SpO2FiO2ratio < 89 , SOFA_RESP_SpO2 := 4]
        # } else {
        intervals[SpO2FiO2ratio_rs %between% c(89,214), SOFA_RESP_SpO2 := 3]
        intervals[SpO2FiO2ratio_rs < 89 , SOFA_RESP_SpO2 := 4]
        
        # }
        
        intervals[PFI_spo2_imputed  > (400/7.5), SOFA_RESP_SpO2_imputed := 0]
        intervals[PFI_spo2_imputed  %between% c(300/7.5,400/7.5), SOFA_RESP_SpO2_imputed := 1]
        intervals[PFI_spo2_imputed  %between% c(0*200/7.5,300/7.5), SOFA_RESP_SpO2_imputed := 2]
        # if (sofa_global_icuregswe){
        #   intervals[PFI_spo2_imputed  %between% c(100/7.5,200/7.5), SOFA_RESP_SpO2_imputed := 3]
        #   intervals[PFI_spo2_imputed  < 100/7.5, SOFA_RESP_SpO2_imputed := 4]
        #
        # } else {
        intervals[PFI_spo2_imputed_rs  %between% c(100/7.5, 200/7.5), SOFA_RESP_SpO2_imputed := 3]
        intervals[PFI_spo2_imputed_rs  < 100/7.5, SOFA_RESP_SpO2_imputed := 4]
      }
      
      intervals[anyGCS == 15, SOFA_GCS := 0]
      intervals[anyGCS %between% c(13,14), SOFA_GCS := 1]
      intervals[anyGCS %between% c(10,12), SOFA_GCS := 2]
      intervals[anyGCS %between% c(6,9), SOFA_GCS := 3]
      intervals[anyGCS < 6, SOFA_GCS := 4]
      
      intervals[clinicianGCS == 15, SOFA_GCS_CLIN := 0]
      intervals[clinicianGCS %between% c(13,14), SOFA_GCS_CLIN := 1]
      intervals[clinicianGCS %between% c(10,12), SOFA_GCS_CLIN := 2]
      intervals[clinicianGCS %between% c(6,9), SOFA_GCS_CLIN := 3]
      intervals[clinicianGCS < 6, SOFA_GCS_CLIN := 4]
      
      #Circulation
      intervals[MAP>70, SOFA_CIRC := 0]
      intervals[MAP<=70, SOFA_CIRC := 1]
      intervals[DobutaminmcgPerKgPerMin>0, SOFA_CIRC := 2]
      intervals[DopaminmcgPerKgPerMin>0 & DopaminmcgPerKgPerMin<=5, SOFA_CIRC := 2]
      
      intervals[LambdenNorEpiEquiv>0 & 
                  LambdenNorEpiEquiv<=0.1, SOFA_CIRC := 3]
      intervals[DopaminmcgPerKgPerMin>5 & DopaminmcgPerKgPerMin<=15, SOFA_CIRC := 3]
      intervals[LambdenNorEpiEquiv>0.1, SOFA_CIRC := 4]
      intervals[DopaminmcgPerKgPerMin>15, SOFA_CIRC := 4]
      
      #Renal
      intervals[,UrinePer24 := (Urine * 24*3600)/ as.numeric(SOFA_end - SOFA_start, units = "secs")]
      intervals[Kreatinin <110, SOFA_RENAL := 0]
      intervals[Kreatinin %between% c(110,170), SOFA_RENAL := 1]
      intervals[Kreatinin %between% c(171,299), SOFA_RENAL := 2]
      intervals[Kreatinin %between% c(300,440), SOFA_RENAL := 3]
      intervals[UrinePer24<500, SOFA_RENAL := 3]
      intervals[Kreatinin > 440, SOFA_RENAL := 4]
      intervals[UrinePer24<200, SOFA_RENAL := 4]
      
      #KDIGO
      intervals[krea_baseline_ratio < 1.5, KDIGO := 0]
      intervals[, krea_increase := Kreatinin - (Kreatinin / krea_baseline_ratio)]
      intervals[krea_increase >= 26.5, KDIGO := 1]
      intervals[, krea_increase := NULL]
      
      intervals[krea_baseline_ratio >= 1.5, KDIGO := 1]
      intervals[perkgperhr6 <= 0.5, KDIGO := 1]
      intervals[krea_baseline_ratio >= 2.0, KDIGO := 2]
      intervals[perkgperhr12 <= 0.5, KDIGO := 2]
      intervals[krea_baseline_ratio >= 3.0, KDIGO := 3]
      intervals[perkgperhr24 <= 0.3, KDIGO := 3]
      intervals[anuria_duration >= 12*3600, KDIGO := 3]
      intervals[Kreatinin >= 353.6, KDIGO := 3]
      intervals[crrt_duration > 0, KDIGO := 3]
      
      
      #Liver
      intervals[Bilirubin<20, SOFA_LIVER := 0]
      intervals[Bilirubin %between% c(20,32), SOFA_LIVER := 1]
      intervals[Bilirubin %between% c(33,101), SOFA_LIVER := 2]
      intervals[Bilirubin %between% c(102,204), SOFA_LIVER := 3]
      intervals[Bilirubin>204, SOFA_LIVER := 4]
      
      #Coagulation
      intervals[TPK>=150, SOFA_COAG := 0]
      intervals[TPK<150, SOFA_COAG := 1]
      intervals[TPK<100, SOFA_COAG := 2]
      intervals[TPK<50, SOFA_COAG := 3]
      intervals[TPK<20, SOFA_COAG := 4]
      
    }
    
    logdt = rbind(logdt, 
                  data.table(event = "sofa_calculation", time = Sys.time(), memory_usage = mem_used()))
  }
}


for (interval in intervallist){
  interval[, SOFA_SCORE := fcoalesce(SOFA_RESP, SOFA_RESP_SpO2_imputed) + SOFA_GCS + SOFA_CIRC+SOFA_RENAL + SOFA_LIVER + SOFA_COAG]
}

for (interval in intervallist){
  setcolorder(interval, rightcolorder)
}


t1 = Sys.time()
identical_rows_highres = detect_identical_rows(intervallist[[3]]$PatientID, intervallist[[3]][, 7:63])
print(Sys.time() - t1)


last_row = intervallist[[3]][, SOFA_N == max(SOFA_N), PatientID][[2]]

#intervallist[[3]][(!identical_rows_highres) | last_row]


tmpdir = "P:\\GLUC-ICU\\data\\processed\\2025-09\\sofa_scores\\"
for (intervalname in names(intervallist)){
  print(intervalname)
  print(Sys.time())
  fn = paste0(tmpdir, "sofa_", intervalname, "_",gsub("[- :]", "_", as.character(round(Sys.time()))) ,".csv")
  fwrite(intervallist[[intervalname]], file = fn, encoding = "UTF-8", bom = TRUE)
  #fwrite(intervallist[[3]][(!identical_rows_highres) | last_row], file = fn, encoding = "UTF-8", bom = TRUE)
  #fwrite(intervallist[[3]][, .SD[(SOFA_N - 1) %% 3 == 0 | SOFA_N == max(SOFA_N)], PatientID], file = fn, encoding = "UTF-8", bom = TRUE)
}


library(pROC)
fd = copy(intervallist[["highres"]][, .(PatientID, AdmissionTime, DischargeTime, SOFA_end, SOFA_N, SOFA_GCS, SOFA_RESP = fcoalesce(SOFA_RESP, SOFA_RESP_SpO2_imputed), SOFA_RENAL, KDIGO, SOFA_COAG, SOFA_LIVER)])
fd[, dt := as.numeric(SOFA_end - AdmissionTime)]

add_cols(fd, 
         unique(master[!is.na(DeathDateTime), .(PatientID, DeathDateTime = as.POSIXct(DeathDateTime, tz = "CET"))]),
         by = "PatientID")
#fd[unique(master[!is.na(DeathDateTime), .(PatientID, DeathDateTime = as.POSIXct(DeathDateTime, tz = "CET"))]), DeathDateTime := DeathDateTime, on = "PatientID"]
fd[, time_to_death := as.numeric(DeathDateTime - DischargeTime)]
fd[time_to_death < 31*24*3600, deathat30 := 1L]
fd[is.na(deathat30), deathat30 := 0L]


# #kdigoset
# kdset = copy(intervallist$icudays)
# add_cols(kdset, 
#          unique(master[!is.na(DeathDateTime), .(PatientID, DeathDateTime = as.POSIXct(DeathDateTime, tz = "CET"))]),
#          by = "PatientID")
# setcolorder(kdset, c("PatientID", "AdmissionTime", "DischargeTime", "DeathDateTime", "SOFA_N", "SOFA_start", "SOFA_end", "SOFA_RENAL", "KDIGO","SOFA_RESP", "SOFA_RESP_SpO2_imputed", "SOFA_GCS", "SOFA_LIVER" , "SOFA_COAG"))
# fwrite(kdset, "sofa_except_circulation_and_kdigo_241107.csv")


# fd[, `:=` (AdmissionTime = NULL,
#            SOFA_end = NULL)]
# fd
setkey(fd, SOFA_N, deathat30, PatientID)


restablelist = list()
nmax = 2*1344L
for (param in c("SOFA_RESP", "SOFA_COAG", "SOFA_LIVER", "SOFA_GCS", "SOFA_RENAL", "KDIGO")){
  print(param)
  print(Sys.time())
  restable = data.table(N = integer(nmax), auc = numeric(nmax), cilow = numeric(nmax), cihigh = numeric(nmax))
  for (i in 1L:nmax){
    
    thisdata = na.omit(fd[SOFA_N == i, c(param, "deathat30"), with = FALSE])
    thisroc = roc(response = thisdata$deathat30, predictor = thisdata[[param]], levels = c(0,1), direction ="<")
    thisci = ci.auc(thisroc)
    set(restable, i, 1L, i)
    set(restable, i, 2L, thisci[[2]])
    set(restable, i, 3L, thisci[[1]])
    set(restable, i, 4L, thisci[[3]])
  }
  restable[, parameter := param]
  restablelist[[param]] = restable
}
restablefull = rbindlist(restablelist)


restablefull[, ICU_day := (4+N*0.25)/24.0]
splineformula = as.formula('y ~ s(x, bs = "tp", m = 3)')
ggplot(data = restablefull[parameter %in% c("KDIGO", "SOFA_RENAL")], aes(x = ICU_day, y = auc, color = parameter))+
  geom_smooth(method = "gam", formula = splineformula)+
  geom_smooth(aes(y = cilow), method = "gam", formula = splineformula, linetype = "dashed")+
  geom_smooth(aes(y = cihigh), method = "gam", formula = splineformula, linetype = "dashed")+
  coord_cartesian(xlim = c(0,14))+
  theme_light()

rt2 = melt(restable, id.vars = c("N", "ICU_day", "parameter"))
ggplot(data = rt2[ICU_day < 24], aes(x = ICU_day, y = value, color = variable))+geom_smooth()+labs(x = "ICU day", y = "AUROC")



