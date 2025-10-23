#antibiotic intensity estimator script
library(Rcpp)

#helper function for db reads
dtGetQuery = function(db, sql, silent = TRUE, print_time = FALSE){
  if (!silent){
    print(sql)
  }
  
  QStime = Sys.time()
  res = dbGetQuery(db, sql)
  QEtime = Sys.time()
  if (!print_time){
    print(QEtime - QStime)
  }
  setDT(res)
  return(res)
}

abx_path = "P:\\GLUC-ICU\\users\\johanh\\gluc_icu_general\\src\\ab_intensity.cpp"
sourceCpp(abx_path, rebuild = TRUE)



library(data.table)
antibiotic_intensity_key_path = "P:\\GLUC-ICU\\users\\anna\\sepsis\\antibiotika_key\\se_nl_linked.csv"
ab_intkey = fread(antibiotic_intensity_key_path)
antibiotic_intensity_key_path2 = "P:\\GLUC-ICU\\users\\anna\\sepsis\\antibiotika_key\\unmatched_ranking_compounds_251015AS.csv"
ab_intkey2 = fread(antibiotic_intensity_key_path2)

library(readxl)
library(glue)
pharma_compid_key = "P:\\GLUC-ICU\\data\\raw\\2022-06\\Info om PharmaID och CopmID 20220607.xlsx"
pharmacomp = readxl::read_excel(pharma_compid_key)
setDT(pharmacomp)


library(RSQLite)
db_path = "P:\\GLUC-ICU\\users\\johanh\\sql\\db_v9_updated_pharma_with_pharmaID.sqlite"
db = dbConnect(SQLite(), db_path)

colnames(ab_intkey)[1] = "PharmaID"
colnames(ab_intkey)[2] = "PharmaName"
colnames(ab_intkey)[3] = "CompName"
ab_intkey[pharmacomp, CompID := CompID, on = "CompName"]

setcolorder(ab_intkey, c(8,3,1,2,4,6,5,7))
colnames(ab_intkey)

colnames(ab_intkey)[6] = "ab_intensity"
ab_intkey[, CompID := as.integer(CompID)]

ab_intkey_comps = unique(ab_intkey[!is.na(ab_intensity), .(CompID, ab_intensity = as.integer(ab_intensity))])
ab_intkey_comps2 = unique(ab_intkey2[, .(CompID, ab_intensity = as.integer(rank))][!is.na(ab_intensity)])

ab_intkey = funion(ab_intkey_comps, ab_intkey_comps2)

ab_intkey[unique(pharmacomp[, .(CompID = as.integer(CompID), CompName)]), CompName := CompName, on = "CompID"]

setcolorder(ab_intkey, c(1,3,2))

setkey(ab_intkey, CompID)

#fwrite(ab_intkey, file.choose())

#in the new database, the sql does the heavy lifting - i.e. joining the concentrations and so on
navid_cohort = unique(demo$PatientID)

#let's use the original master file (here translated to english but anyway)
master = dbReadTable(db, "hv_master")
setDT(master)
#no children
master = master[Age >= 18]

#get all antibtioics we know of for these patients
ids = unique(master$PatientID)
ab_comps = unique(ab_intkey$CompID)
#query the db
sql = glue_sql("SELECT * FROM hv_pharmacomp_fast WHERE PatientID IN ({ids*}) AND CompID IN ({ab_comps*});", .con = db)
{
  t1 = Sys.time()
  ab_pharma = dtGetQuery(db, sql)
  t2 = Sys.time()
  print(t2 - t1)
}
setkey(ab_pharma, PatientID, OrderNumber, PharmaID, CompID, DateTime, RowNumber)

#ab_pharma2 = copy(ab_pharma)

#ab_pharma = copy(ab_pharma2)


ab_pharma[ab_intkey, ab_intensity := ab_intensity, on = "CompID"]

ab_pharma = ab_pharma[GivenDose > 0 & !is.na(ab_intensity)]
ab_pharma[,DateTime := as.POSIXct(DateTime, tz = "CET")]

#fetch the most common dosage intervals for the antibiotics in question

setkey(ab_pharma, PatientID, CompID, OrderNumber, DateTime)
abx_doseinterval = ab_pharma[, .(median_interval = as.numeric(median(DateTime - shift(DateTime, 1L), na.rm = TRUE), units = "secs")), .(PatientID, CompID, OrderNumber)]

#determine the common intervals for these antibiotics
abx_doseinterval = abx_doseinterval[!is.na(median_interval), .N, .(CompID, median_interval)][order(CompID, N, decreasing = TRUE)]
abx_doseinterval = abx_doseinterval[order(CompID, N, decreasing = TRUE)]
setkey(abx_doseinterval, CompID)
abx_doseinterval[ab_intkey, CompName := CompName, on = "CompID"]
abx_doseinterval = abx_doseinterval[, .(median_interval, N = N, percent = N/sum(N)), .(CompID, CompName)]
abx_doseinterval = abx_doseinterval[, .SD[1], CompID]


#join with the intensity
abx_doseinterval[ab_intkey, ab_intensity := ab_intensity, on = "CompID"]

#sort
setkey(abx_doseinterval, CompID)

library(LaF)

#now, add the antibiotics from TC to make sure we REALLY capture all escalations
abx_tc_path = "P:\\GLUC-ICU\\data\\processed\\2025-04\\Gluc_Icu_TC_Lakemedel_ATC_J0_Ordination_Adm_250327.csv"

abx_tc = fread(abx_tc_path, encoding = "Latin-1")
colnames(abx_tc) = gsub("~~", "", colnames(abx_tc))

abx_codes = abx_tc[, .(lakemedel_preparat, ATC_kod)]
abx_codes = abx_codes[, .(n = .N), c("lakemedel_preparat", "ATC_kod")]
abx_codes = abx_codes[, lapply(.SD, function(x) gsub("~~", "", x)), n]
setcolorder(abx_codes, c("ATC_kod", "lakemedel_preparat"))
setkey(abx_codes, ATC_kod, n)

abx_tc[,ATC_kod := gsub("~~", "", ATC_kod)]


#use this filter from John to make sure it's cleaned
abx_tc = abx_tc[grepl("^J0", ATC_kod) & #
                  !(ATC_kod %in% c("J01CA08", "J01EA01", "J01XE01", "J01XX05") | #hard coded to exclude UTI abx
                      grepl("^J05|^J06|^J07", ATC_kod) #hard code removal of antiviral, immunoglobulins, and vaccine
                  )
]



#match the final sepsis antimicrobials with the ATC codes
final_abx_path = "P:\\GLUC-ICU\\users\\anna\\sepsis\\documents eSepsisscript\\Sepsisantimicrobials FINAL 250206.xlsx"
as_abx_final = readxl::read_excel(final_abx_path)
setDT(as_abx_final)
#melt(as_abx_final, id.vars = c("PharmaID", "CompID"), measure.vars = c("CompName", "trade_name"))


#set up this loop to match the trade and compound names with ATC codes
#first join to get both compname and compid for these...
sepsis_abx = merge(as_abx_final, pharmacomp[, .(PharmaID, CompID, CompName)], by = "PharmaID", all.x = TRUE)
sepsis_abx[, CompID := as.integer(CompID)]
sepsis_abx[, trade_name := str_split_fixed(PharmaName, " |\\(", 2)[,1]]

#get the long list of string values to match with the drugs...
sepsis_abx_long = melt(sepsis_abx, id.vars = c("PharmaID", "CompID"), measure.vars = c("CompName", "trade_name"), variable.name = "string_type", value.name = "string_value")
#remove energy, glucose, na/cl
sepsis_abx_long = sepsis_abx_long[!string_value %in% c("Na+", "Cl-", "NULL", "Energi tot", "Glukos", "Energi(kolhydr)")]
#remove mianserin
sepsis_abx_long = sepsis_abx_long[!grepl("ianserin", string_value)]

#we're going to match on the string_value later on...
name_match_dt = data.table(string_value = sepsis_abx_long[, unique(string_value)])
name_match_dt[, lowercase := tolower(string_value)]


#string matching to make an educated guess on ATC code for the compnames
library(stringdist)
matchlist = list()
lowernames = tolower(abx_codes$lakemedel_preparat)

#loop through strings to match
for (i in 1:nrow(name_match_dt)){
  this_lower_case = name_match_dt$lowercase[i]
  this_string_value = name_match_dt$string_value[i]
  #exact matches
  exact_matches = unique(abx_codes[grep(this_lower_case, lowernames), .(ATC_kod, lakemedel_preparat, similarity = 1.0)])
  #partial matches
  partial_best_match = abx_codes[, .(ATC_kod, lakemedel_preparat, similarity = stringsim(this_lower_case, lowernames))][order(similarity, decreasing = TRUE)]
  #keep best matches and top 3 best partial matches
  thisres = rbind(exact_matches, partial_best_match[1:3])
  #add the long name to allow the join back later on
  thisres[, string_value := this_string_value]
  matchlist[[length(matchlist)+1]] = thisres
  
  
}
matchdt = rbindlist(matchlist)

#for manual checking
#fwrite(matchdt[, .SD[1], .(string_value, ATC_kod, similarity)], file.choose())

#after manual check, there will be a xlsx-file that links these
atc_key_path = "P:\\GLUC-ICU\\users\\anna\\sepsis\\antibiotika_key\\sepsis_atc_codes_checked_251022.xlsx"
atc_key = read_excel(atc_key_path)
setDT(atc_key)

setkey(sepsis_abx_long, string_value)
setkey(atc_key, string_value)

comp_atc_pharma_key = merge(sepsis_abx_long, atc_key, allow.cartesian = TRUE, all.x = TRUE)

hardcoded_atcs = data.table(PharmaID = c(1001219, #pentacarinat
                                         1002290, #noxafil
                                         1002386), #noxafil
                            ATC_kod = c("P01CX01",
                                        "J02AC04",
                                        "J02AC04"))
hardcoded_atcs = merge(hardcoded_atcs, pharmacomp[, .(PharmaID, PharmaName)])

matched_atcs = unique(comp_atc_pharma_key[, .(PharmaID, ATC_kod)])
matched_atcs = merge(matched_atcs, pharmacomp[, .(PharmaID, PharmaName)])
pharma_atcs = rbind(hardcoded_atcs, matched_atcs)


#ab_pharma[!PharmaID %in% pharma_atcs$PharmaID, unique(PharmaID)]
unmatched = copy(pharmacomp)[PharmaID %in% ab_pharma[!PharmaID %in% pharma_atcs$PharmaID, unique(PharmaID)]][, .(PharmaID, PharmaName, CompID, CompName)]
unmatched[, CompID := as.numeric(CompID)]

unmatched = merge(unmatched, unique(comp_atc_pharma_key[, .(CompID, ATC_kod, string_value)])[!is.na(CompID)], "CompID", allow.cartesian = TRUE)

#match them on name again
unmatched[, pharma_part1 := str_split_fixed(PharmaName, " |\\(", 2)[,1]]


unmatched = unique(unmatched[pharma_part1 == string_value, .(PharmaID, ATC_kod)])
unmatched = merge(unmatched, unique(pharmacomp[, .(PharmaID, PharmaName)]))
unmatched = unmatched[!PharmaName %in% c("Fucidin kräm 2 %",
                            "Trimetoprim tabl 160mg",
                            "Trimetoprim mixt 10mg/ml")]

pharma_atc_key = unique(rbind(pharma_atcs, unmatched))




ab_pharma[!PharmaID %in% pharma_atc_key$PharmaID, unique(PharmaID)]








pharma_atc_key = unique(pharma)


comp_atc_key = unique(comp_atc_pharma_key[, .(CompID, ATC_kod)])[!is.na(CompID)]
setkey(comp_atc_key)

pcomp2 = pharmacomp[, .(CompID = as.numeric(CompID),PharmaID, CompName)][!is.na(CompID)]
setkey(pcomp2, CompID)
comp_pharma_atc_key = merge(comp_atc_key, pcomp2, "CompID")
comp_pharma_atc_key = comp_pharma_atc_key[!CompName %in% c("Na+", 
                                                           "Energi tot",
                                                           "Glukos",
                                                           "Energi(kolhydr)",
                                                           "Cl-")]

comp_atc_key = unique(comp_pharma_atc_key[, .(CompID, ATC_kod)])

abp2 = copy(ab_pharma)
abp2 = merge(abp2, comp_atc_key, "CompID", all.x = TRUE, allow.cartesian = TRUE)

unique(pharmacomp[, .(PharmaID, PharmaName)])
abp2 = merge(abp2, unique(pharmacomp[, .(PharmaID, PharmaName)]), "PharmaID")

abp2 = abp2[!is.na(ATC_kod)]






comp_atc_pharma_key = merge(comp_atc_pharma_key, 
                            pharmacomp[,.(PharmaID, CompID = as.numeric(CompID), CompName)], 
                            c("PharmaID", "CompID"),
                            all.x = TRUE)

comp_atc_pharma_key = merge(comp_atc_pharma_key, 
                            unique(pharmacomp[, .(PharmaID, PharmaName)]),
                            "PharmaID")
comp_atc_pharma_key[, PharmaID := as.integer(PharmaID)]

setkey(comp_atc_pharma_key, PharmaID, CompID, ATC_kod)
#fwrite(comp_atc_pharma_key, file.choose())

pkey = unique(comp_atc_pharma_key[, .(PharmaID, ATC_kod, PharmaName)])
setkey(pkey, PharmaID)

ab_pharma_atc = merge(ab_pharma, pkey, "PharmaID", all.x = TRUE, allow.cartesian = TRUE)






#load cultures
bl_path = "P:/GLUC-ICU/users/Logan/RawData/Gluc-Icu TC Blododling 241022.csv"
#to match how Logan read this file
bl = read.csv2(bl_path, encoding = "latin1")
setDT(bl)

#bl[, .N, undersokning][order(N, decreasing = TRUE)]
#bl[, .N, .(provmaterial)][order(N, decreasing = TRUE)]

#fix this error
#bl[StudyID %in% c(13921, 25581)]
#bl[, .(nchar = nchar(provtagning_datum))][, .(n = .N), nchar]
#bl[nchar(provtagning_datum) != 10]

bl[nchar(provtagning_tid) == 10, `:=`(provtagning_datum = provtagning_tid,
                                      provtagning_tid = analysnamn,
                                      analysnamn = ersattningsvarde)]

#only include where a date is known
bl = bl[nchar(provtagning_datum) == 10]
#force sample time to 12:00 if not known
bl[nchar(provtagning_tid) == 0, provtagning_tid := "12:00:00"]

#bl[, .(chars = nchar(provtagning_tid))][, .(n = .N), chars]
bl[, DateTimeString := paste0(provtagning_datum, " ", provtagning_tid)]
bl[, DateTime := as.POSIXct(DateTimeString, tz = "CET")]

#unique time points for blood sampling
blood_sample_times = unique(bl[, .(StudyID = as.numeric(StudyID), DateTime)])
setkey(blood_sample_times, StudyID, DateTime)

#unique time points for initiating antibiotics
ab_starttimes = ab_pharma[, .(PatientID, OrderNumber, DateTime)]

ab_starttimes[unique(master[, .(PatientID, StudyID)]), StudyID := StudyID, on = "PatientID"]
setcolorder(ab_starttimes, c("StudyID", "PatientID"))
setkey(ab_starttimes, StudyID, DateTime)

#this will need the ATC code later on but let's go with order number for now
ab_starttimes[, n := .N, .(PatientID, OrderNumber)]

#remove single administrations
ab_starttimes = ab_starttimes[n>1]
setkey(ab_starttimes, StudyID, OrderNumber, DateTime)

#time for first administration of antibiotics per order number
ab_starttimes = ab_starttimes[, .SD[1], .(StudyID, PatientID, OrderNumber)]

#full merge, i.e. many-to-many matching
sidata = merge(blood_sample_times[, .(StudyID, culture_time = DateTime)],
               unique(ab_starttimes[, .(StudyID, OrderNumber, antibiotic_initiation_time = DateTime)]),
               "StudyID", 
               allow.cartesian = TRUE)

sidata[, culture_to_antibiotic := as.numeric(antibiotic_initiation_time - culture_time, units = "hours")]

sidata[culture_to_antibiotic %between% list(0,72), option_1 := 1L]
sidata[culture_to_antibiotic %between% list(-24,0), option_2 := 1L]
sidata[is.na(option_1), option_1 := 0L]
sidata[is.na(option_2), option_2 := 0L]
sidata[, suspected_infection := fifelse(option_1 == 1 | option_2 == 1, 1L, 0L)]
sidata[suspected_infection == 1, suspected_infection_time := pmin(culture_time, antibiotic_initiation_time)]

sidata = unique(sidata[suspected_infection == 1, .(StudyID, suspected_infection_time)])

#sort and count time between suspected infection times
setkey(sidata, StudyID, suspected_infection_time)
sidata[, dt := as.numeric(suspected_infection_time - shift(suspected_infection_time, 1L), units = "hours"), StudyID]






#lockout period : no new suspected infection within 72 hours
sidata = sidata[is.na(dt) | !is.na(dt) & dt > 72]

fwrite(sidata[, .(StudyID, suspected_infection_time)], file.choose())



fwrite(sepsis_ab_key)
sepsis_ab_key = merge(ab_intkey, atc_key, by = "CompName", all.x = TRUE)
sepsis_ab_key[is.na(ATC_kod)]

setkey(ab_pharma, PatientID, DateTime)
#the c++ call to determine the intensity at any given moment
ab_intensity = level_at_time_rcpp(ab_pharma$PatientID,
                                  ab_pharma$DateTime,
                                  ab_pharma$CompID,
                                  abx_doseinterval$CompID,
                                  abx_doseinterval$median_interval*3,
                                  abx_doseinterval$ab_intensity)
setDT(ab_intensity)
ab_pharma = cbind(ab_pharma, ab_intensity)

#setting up the wide dataset for the compounds for later on...
ab_pharma[ab_intkey, CompName := CompName, on = "CompID"]
unique(ab_pharma$PharmaUnit)

#shift all doses to mg when applicable
ab_pharma[PharmaUnit == "g", `:=` (GivenDose = GivenDose * 1000,
                                   PharmaUnit = "mg")]
ab_pharma[, CompNameUnit := paste0(CompName, " (", PharmaUnit, ")")]
#wide version
ab_wide = dcast(ab_pharma, PatientID + DateTime ~ CompNameUnit, value.var = "GivenDose", fun.aggregate = sum)
#add the antibiotic level and count
ab_wide[unique(ab_pharma[, .(PatientID, DateTime, RowNumber, level, count)]), `:=` (level = level, count = count), on = c("PatientID", "DateTime")]

#save for the future
#fwrite(ab_wide, file.choose())
rm(ab_wide)

ab_pharma[, dt := as.numeric(DateTime - shift(DateTime, 1L), units = "secs"), PatientID]


#compute escalations
ab_pharma[, `:=` (delta_level = level - shift(level, 1L),
                  delta_count = count - shift(count, 1L)), PatientID]

#N/A means no previous, so in that case deltas are just the current maximum
ab_pharma[is.na(delta_level), delta_level := level]
ab_pharma[is.na(delta_count), delta_count := count]

#if more than 3 days gone since last antibiotic - assume this is new antibiotics regardless
ab_pharma[dt > 24*3*3600, `:=` (delta_level = level,
                                delta_count = count)]

#now, compute escalations
ab_pharma[, escalation_level := 0L]
ab_pharma[delta_level > 0, escalation_level := 1L]

ab_pharma[, escalation_count := 0L]
ab_pharma[delta_count > 0, escalation_count := 1L]

ab_pharma[, escalation := fifelse(escalation_count == 1 | escalation_level == 1, 1,0)]

#setup the escalation dataset
ab_esc = ab_pharma[escalation == 1, .(PatientID, DateTime, level, count, escalation_level, escalation_count)]

#now, lockout the escalations by 72 hours
setkey(ab_esc, PatientID, DateTime)

#here's a set of escalations!
ab_esc[, dt := as.numeric(DateTime - shift(DateTime, 1L), units = "secs"), PatientID]


