#2023-04-18
#shifting the weight script into its own function to reduce the amound of lines

#2024-10-14
#including the weights from TC

get_weights = function(db, options, pids, intkey, vartypes, verbose){
  if (FALSE){
    pids = copy(pids)
    intkey = copy(intkey)
  }
  
  if (verbose){
    print("Weight data not supplied to function, generating complete set of IdealWeights... query sent to db")
  }
  
  #load the karda-weights
  kwh = fread(options$file$paths$karda_weight, encoding = options$file$encoding, colClasses=c("character"))
  kwh[, DateTime := as.POSIXct(handelsetidpunkt, tz = options$time$ptz, origin = options$time$po, format = options$time$format)]
  kwh = kwh[, .(StudyID = as.integer(StudyID), DateTime, termnamn, Value = as.numeric(gsub(",", ".", varde)))]
  
  
  
  #only subset the relevant patients
  kwh = kwh[StudyID %in% intkey$StudyID]
  
  #fix patientID == 34641 whichi has two studyIDs
  
  #kwh[StudyID %in% master2[PatientID == 34641, StudyID]]
  
  #Harmonize column names with CCC names
  kwh[, EnterTime := DateTime]
  setcolorder(kwh, c(1,2,5,3,4))
  colnames(kwh)[4] = "Abbreviation"
  kwh[Abbreviation == "Längd", Abbreviation := "Längd_TC"]
  kwh[Abbreviation == "Vikt", Abbreviation := "Vikt_TC"]
  setkey(kwh, StudyID, DateTime, EnterTime, Abbreviation)
  
  if(verbose){
    endtime = Sys.time()
    print(sofa_echo_time("Karda weights read", starttime, prevtime, endtime))
    prevtime = Sys.time()
  }
  
  #Fetch the CCC weights
  prevtime = Sys.time()
  sql = glue_sql("SELECT PatientID, VariableID, DateTime, EnterTime, NumericValue FROM obsrec WHERE VariableID IN (10000400, 10000450, 15002666,15003298) AND PatientID IN ({pids*})", .con = db)

  wh = dtGetQuery(db, sql)[, .(PatientID,
                               VariableID,
                               DateTime = as.POSIXct(DateTime, tz = options$time$ptz, origin = options$time$po),
                               EnterTime = as.POSIXct(EnterTime, tz = options$time$ptz, origin = options$time$po),
                               Value = NumericValue)]
  
  endtime = Sys.time()
  
  setkey(wh, PatientID, VariableID, DateTime)
  
  #add the units
  ref_left_join(wh, vartypes[, .(VariableID, Abbreviation, Unit)], by =  "VariableID")
  
  #merge with the key to StudyIDs
  wh2 = merge(wh, 
        unique(intkey[PatientID %in% pids, .(PatientID, StudyID)]),
        "PatientID",
        all.x = TRUE,
        allow.cartesian = TRUE)
  
  #ditch columns that aren't in the karda data
  wh2[, `:=` (PatientID = NULL,
              Unit = NULL,
              VariableID = NULL)]
  setcolorder(wh2, c("StudyID", "DateTime", "EnterTime", "Abbreviation", "Value"))
  
  setkey(wh2, StudyID, DateTime, EnterTime, Abbreviation)
  
  
  if(verbose){
    endtime = Sys.time()
    print(sofa_echo_time("Weights fetched from db", starttime, prevtime, endtime))
    prevtime = Sys.time()
  }
  
  
  
  

  #Generate the 'raw' dataset for others to enjoy
  if(FALSE){
    
    #full weights and heihgs, i.e. KARDA + CCC
    fwh = unique(rbind(kwh, wh2))
    #sort
    setkey(fwh, StudyID, DateTime, EnterTime, Abbreviation)
    #if there are multiple records per time point, track them
    fwh[, RowNumber := 1:.N, .(StudyID, DateTime, EnterTime, Abbreviation)]
    setkey(fwh, StudyID, DateTime, EnterTime, Abbreviation, RowNumber)
    
    fwh = dcast(fwh, StudyID+DateTime+EnterTime+RowNumber~Abbreviation, value.var = "Value", fun.aggregate = first)
    
    setcolorder(fwh, c("StudyID", "DateTime", "EnterTime", "RowNumber", "L?ngd", "L?ngd_TC", "Ankomstvikt", "Daglig vikt", "Vikt", "Vikt_TC"))
    
    #the 'raw' merged dataset
    fwh = fwh[, lapply(.SD, as.character)]
    
    #change NA to "" for easy viewing in excel
    fwh[is.na(Ankomstvikt), Ankomstvikt := ""]
    fwh[is.na(Vikt), Vikt := ""]
    fwh[is.na(Vikt_TC), Vikt_TC := ""]
    fwh[is.na(`Daglig vikt`), `Daglig vikt` := ""]
    fwh[is.na(`L?ngd`), `L?ngd` := ""]
    fwh[is.na(`L?ngd_TC`), `L?ngd_TC` := ""]
    
    #a cleaned raw dataset for whoever wants it
    fwrite(fwh, "gluc_icu_weights_heights_tc_ccc_raw_241104.csv", sep = ";", encoding = "UTF-8", bom = TRUE)
    
  }

  #now redo the whole thing but clean the data
  #only keep the first registration of 'ankomstvikt'
  wh3 = rbind(wh[VariableID != 15002666],
             wh[VariableID == 15002666, lapply(.SD, first), by = PatientID]
  )

  #merge with the key to StudyIDs
  wh3 = merge(wh3, 
              unique(intkey[PatientID %in% pids, .(PatientID, StudyID)]),
              "PatientID",
              all.x = TRUE,
              allow.cartesian = TRUE)
  
  #ditch columns that aren't in the karda data
  wh3[, `:=` (PatientID = NULL,
              Unit = NULL,
              VariableID = NULL)]
  setcolorder(wh3, c("StudyID", "DateTime", "EnterTime", "Abbreviation", "Value"))
  
  setkey(wh3, StudyID, DateTime, EnterTime, Abbreviation)
  
  #distinct values per time point
  wh3 = unique(wh3[, .(StudyID, DateTime, EnterTime, Abbreviation, Value)])
  fwh = unique(rbind(kwh, wh3))
  
  setkey(fwh, StudyID, DateTime, EnterTime, Abbreviation)
  #add the row number to track
  fwh[, RowNumber := 1:.N, .(StudyID, DateTime, EnterTime, Abbreviation)]
  
  #Spread over two columns for future comparisons
  fwh2 = dcast(fwh, StudyID+DateTime+EnterTime+RowNumber~Abbreviation, value.var = "Value", fun.aggregate = dplyr::first)
  if(verbose){
    endtime = Sys.time()
    print(sofa_echo_time("Weights harmonized between sources", starttime, prevtime, endtime))
    prevtime = Sys.time()
  }
  
  #ditch measurements from children
  sql = glue_sql("SELECT PatientID, StudyID, AdmissionYear, Age FROM sir_master WHERE PatientID IN ({pids*})", .con = db)
  myob = dtGetQuery(db, sql)
  #whenever a query has studyid and patientid, force this conversion!
  #myob[PatientID == 34641, StudyID := 23092]
  
  myob = unique(myob[, .(StudyID, YoB = AdmissionYear - Age)])
  myob = myob[, .(YoB = floor(mean(YoB))), StudyID]
  fwh2[myob, YoB := YoB, on = "StudyID"]
  fwh2[, age_at_measurement := as.integer(year(DateTime)) - YoB]
  fwh2 = fwh2[age_at_measurement >=18]
  fwh2[, `:=` (YoB = NULL,
               age_at_measurement = NULL)]
  rm(myob)
  
  #
  fwh2[, Weight := fcoalesce(Vikt, Vikt_TC, `Daglig vikt`, `Ankomstvikt`)]
  fwh2[, Height := fcoalesce(`Längd`, `Längd_TC`)]
  fwh2 = unique(fwh2[, .(StudyID, DateTime, EnterTime, Height, Weight)])
  
  if(verbose){
    endtime = Sys.time()
    print(sofa_echo_time("Weights from non-adult period removed", starttime, prevtime, endtime))
    prevtime = Sys.time()
  }
  
  #This is a very simple outlier detection algorithm but it seems to work for weights
  fwh2[Weight==0, Weight := NA]
  fwh2[Height==0, Height := NA]
  
  #this could be written as !is.na(Height) | !is.na(Weight) too ;)
  fwh2 = fwh2[!(is.na(Height) & is.na(Weight))]

  #try to decimal shift the inputs
  #reasonable range of human height : 50 - 270 cm
  #Sorry, Robert Waldow
  #fwh2c = copy(fwh2)
  #For height, use population distribution
  fwh2[, Height := decimal_shift_input(Height, 
                                       options$sofa$antro$heigth_decimal_shift_lower_bound, 
                                       options$sofa$antro$heigth_decimal_shift_upper_bound, 
                                       SDs = 2)]
  #For weights, use patient distribution
  fwh2[, Weight := decimal_shift_input(Weight, 
                                       options$sofa$antro$weight_decimal_shift_lower_bound, 
                                       options$sofa$antro$weight_decimal_shift_upper_bound, 
                                       SDs = 5), 
       by = StudyID]
  
  fwh2 = fwh2[!(is.na(Height) & is.na(Weight))]
  
  
  
  
  #First, find errors where the weight is switched with the height
  fwh2[Weight > Height & 
         Weight > options$sofa$antro$weight_switch_cutoff & 
         Height < options$sofa$antro$heigth_switch_cutoff, 
       `:=` (Weight = Height, Height = Weight)]
  #Second, try to find errors where the weight is likely to be temperature
  #Skip this, will be caught by Hampel Filter later
  #fwh2[Weight %between% c(34,42) & StudyID %in% fwh2[Weight>60]$StudyID, Weight := NA]
  
  #make sure the table is sorted and without NAs before the c++ function call
  fwh2 = fwh2[!(is.na(Height) & is.na(Weight))]
  setkey(fwh2, StudyID, DateTime)
  #use the Hampel filter to filter the weights
  fwh2[!is.na(Weight), Weight2 := tp_outliers(StudyID,
                                              DateTime, 
                                              Weight,
                                              type = options$sofa$antro$weight_outlier_function,
                                              before = options$sofa$antro$weight_outlier_before,
                                              after = options$sofa$antro$weight_outlier_after,
                                              k = options$sofa$antro$weight_outlier_k)]
  
  if(verbose){
    endtime = Sys.time()
    print(sofa_echo_time("Weights outliers removed", starttime, prevtime, endtime))
    prevtime = Sys.time()
  }
  
  
  if(FALSE){
    #cleaned dataset
    fwh2clean = copy(fwh2)
    fwh2clean[Weight != Weight2, Weight := NA]
    fwh2clean = fwh2clean[!is.na(Height) | !is.na(Weight), 1:5, with = FALSE]
    fwh2cc = copy(fwh2clean)
    setkey(fwh2cc, StudyID, DateTime)
    fwh2cc = fwh2cc[, lapply(.SD, as.character)]
    fwh2cc[is.na(Height), Height := ""]
    fwh2cc[is.na(Weight), Weight := ""]
    fwh2cc[, Height := gsub("\\.", ",", Height)]
    fwh2cc[, Weight := gsub("\\.", ",", Weight)]
    fwrite(fwh2cc, "gluc_icu_weight_height_ccc_tc_cleaned_241104.csv", sep = ";", dec = ",", encoding = "UTF-8", bom = TRUE)
    
    #generate the population and gender means
    if(FALSE){
      fwh2clean[unique(master[, .(StudyID, SexMale)]), SexMale := SexMale, on = "StudyID"]
      
      popdata = rbind(
        fwh2clean[!is.na(Weight), .(StudyID, DateTime, Measurement = "Weight", Value = Weight, SexMale)],
        fwh2clean[!is.na(Height), .(StudyID, DateTime, Measurement = "Height", Value = Height, SexMale)]
      )
      popdata = rbind(
        popdata[, .(meanval = mean(Value)), .(StudyID, SexMale, Measurement)][, .(mean = mean(meanval)), .(Measurement, SexMale)],
        popdata[, .(meanval = mean(Value)), .(StudyID, Measurement)][, .(SexMale = NA_integer_, mean = mean(meanval)), .(Measurement)]
      )
      
      fwh2clean[,SexMale := NULL]
      setkey(popdata, Measurement, SexMale)
      popdata = dcast(popdata, SexMale~Measurement, value.var = "mean")
      popdata[, bmi := Weight / (Height*Height*1e-04)]
    }
    
  }
  
  #Ditch the outliers
  fwh2[Weight != Weight2, Weight := NA]
  fwh2[, Weight2 := NULL]
  fwh2 = fwh2[!is.na(Height) | !is.na(Weight)]
  setkey(fwh2, StudyID, DateTime, EnterTime)
  
  
  
  
  #Genarate the 'SOFA weights'
  #sql = paste0("SELECT PatientID, StudyID, SexMale, FirstAdmissionDateTime as AdmissionTime, LastDischargeDateTime as DischargeTime FROM patients WHERE PatientID IN ", urval_sql)
  antrodata = copy(intkey)[PatientID %in% pids]
  
  ref_left_join(antrodata, unique(spell_master_key[, .(spell_id, spell_sex)]), by = "spell_id")
  
  antrodata = antrodata[, .(StudyID, spell_id, SexMale = fifelse(spell_sex == "Male",1,0),spell_start, spell_stop)]
  antrodata = unique(antrodata)
  
  # antrodata = copy(master)[PatientID %in% urval, .(StudyID, 
  #                                                  PatientID, 
  #                                                  SexMale, 
  #                                                  AdmissionTime = as.POSIXct(AdmissionDateTime, tz = options$time$ptz, origin = options$time$po), 
  #                                                  DischargeTime = as.POSIXct(DischargeDateTime, tz = options$time$ptz, origin = options$time$po))]
  # 
  # #antrodata = dtGetQuery(db, sql)
  # # antrodata[, `:=` (AdmissionTime = as.POSIXct(AdmissionDateTime, tz = options$time$ptz, origin = options$time$po),
  # #                   DischargeTime = as.POSIXct(DischargeDateTime, tz = options$time$ptz, origin = options$time$po))]
  # antrodata = antrodata[, .(SexMale = first(SexMale),
  #                           AdmissionTime = min(AdmissionTime),
  #                           DischargeTime = max(DischargeTime)),
  #                       .(StudyID, PatientID)]
  
  
  
  
  #we're going to join on spell_id eventually...
  spellkey = unique(antrodata[, .(StudyID, spell_id)])
  
  #from now on it's spell_ids all the way
  antrodata[, StudyID := NULL]
  setkey(antrodata, spell_id, spell_start, spell_stop)
  
  #first weights, merge to allow spell_ids
  wts = copy(fwh2)[!is.na(Weight), .(StudyID, DateTime, Weight)]
  wts = merge(spellkey, wts, "StudyID", allow.cartesian = TRUE)
  wts = unique(wts[, .(spell_id, DateTime, Weight)])
  setkey(wts, spell_id, DateTime)
  
  #make sure they're keyed first before this c++ call...
  interval_aggregate(antrodata, 
                     wts, 
                     "Weight",
                     options, 
                     fun.aggregate = options$sofa$antro$sofa_weight_aggregation_function, 
                     max_locf_time = options$sofa$antro$sofa_weight_locf, 
                     max_nocb_time = options$sofa$antro$sofa_weight_nocb)
  
  
  hts = copy(fwh2)[!is.na(Height), .(StudyID, DateTime, Height)]
  hts = merge(spellkey, hts, "StudyID", allow.cartesian = TRUE)
  hts = unique(hts[, .(spell_id, DateTime, Height)])
  setkey(hts, spell_id, DateTime)
  
  #make sure they're keyed first before this c++ call...
  interval_aggregate(antrodata, 
                     hts, 
                     "Height",
                     options, 
                     fun.aggregate = options$sofa$antro$sofa_weight_aggregation_function, 
                     max_locf_time = options$sofa$antro$sofa_weight_locf, 
                     max_nocb_time = options$sofa$antro$sofa_weight_nocb)
  
  if(verbose){
    endtime = Sys.time()
    print(sofa_echo_time("Weights aggregated per patient", starttime, prevtime, endtime))
    prevtime = Sys.time()
  }
  
  
  
  #fetch the population averages
  popdata = data.table(copy(options$sofa$antro$popdata))
  ref_left_join(antrodata, popdata, by = "SexMale")
  rm(popdata)
  
  #estimate weight and height backwards from average BMI if one is present
  antrodata[, `:=` (weight_calc = mean_bmi * Height_first * Height_first * 1e-4,
                    height_calc = sqrt(1e4 * Weight_first / mean_bmi))]
  
  #Priority: first any measured weight, then estimated from height, third population averages
  antrodata[, `:=` (Weight = round(fcoalesce(Weight_first, weight_calc, mean_weight), options$sofa$antro$sofa_weight_round),
                    Height = round(fcoalesce(Height_first, height_calc, mean_height), options$sofa$antro$sofa_weight_round),
                    Weight_first = NULL,
                    Height_first = NULL,
                    mean_weight = NULL,
                    mean_height = NULL,
                    weight_calc = NULL,
                    height_calc = NULL,
                    mean_bmi = NULL,
                    SexMale = NULL)]
  
  setkey(antrodata, spell_id, spell_start)
  #generate the 'sofa weights'-table
  if(FALSE){
    anc = copy(antrodata)
    setcolorder(anc, c(1,2,3,4,5,7,6,8))
    colnames(anc)[5] = "SOFA_weight_DateTime"
    colnames(anc)[6] = "SOFA_weight"
    colnames(anc)[7] = "SOFA_height_DateTime"
    colnames(anc)[8] = "SOFA_height"
    setkey(anc, StudyID, PatientID, AdmissionTime)
    anc = anc[, lapply(.SD, as.character)]
    anc[is.na(SOFA_weight_DateTime), SOFA_weight_DateTime := ""]
    anc[is.na(SOFA_height_DateTime), SOFA_height_DateTime := ""]
    anc[, SOFA_weight := gsub("\\.", ",", SOFA_weight)]
    anc[, SOFA_height := gsub("\\.", ",", SOFA_height)]
    
    fwrite(anc, "sofa_weight_height_241104.csv", sep = ";", encoding = "UTF-8", bom = TRUE)
    
    ad2 = copy(antrodata)
    ad2[Weight_first_DateTime >= AdmissionTime & Weight_first_DateTime <= DischargeTime, WeightDataSource := "Measured in ICU"]
    ad2[Weight_first_DateTime < AdmissionTime, WeightDataSource := "Last observation carried forward"]
    ad2[Weight_first_DateTime > DischargeTime, WeightDataSource := "Next observation carried backward"]
    ad2[is.na(Weight_first_DateTime) & !is.na(Height_first_DateTime), WeightDataSource := "Imputed from height and population BMI"]
    ad2[is.na(Weight_first_DateTime) & is.na(Height_first_DateTime), WeightDataSource := "Imputed from population"]
    ad2[, .N, WeightDataSource][order(N, decreasing = TRUE)]
  }
  
  #for the rest of the code, only the weight and datetime is needed
  antrodata = antrodata[, .(spell_id,
                            WeightDateTime = Weight_first_DateTime, 
                            Weight,
                            HeightDateTime = Height_first_DateTime,
                            Height)]
  antrodata = unique(antrodata)
  
  setkey(antrodata, spell_id)
  
  return(antrodata)
}
