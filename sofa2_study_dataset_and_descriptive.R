#2025-10-21
#code for a table one for the trial, and for general dataset preparation
if (FALSE){
  
  
  #2025-10-23
  #Use Navid's clean table 1-dataset as master as it's the same spells that JM generated
  navid_t1_path = "P:\\GLUC-ICU\\users\\navid\\TillJohanH\\SIR_spells_221108_20251021_ALLPATIENTID.rds"
  navidt1 = readRDS(navid_t1_path)
  setDT(navidt1)
  navidt1[, RowNo := 1:.N]
  #navidt1 = navidt1[1:5, 1:16, with = FALSE]
  
  
}

#master datapath
master_datapath_old = "P:\\GLUC-ICU\\data\\raw\\2022-11\\SIR_Masterfile_221108.xlsx"
master_old = read_excel(master_datapath_old)
setDT(master_old)

#2025-10-23
#Use the spells Johan generated to be consistent with previous
master_datapath = "P:\\GLUC-ICU\\users\\johanh\\sofa2025\\datasets\\spell_master_dataset_from_jm_stat_file_with_spell_id_251023.csv"
master = fread(master_datapath)

if(FALSE){
  #convert to CET
  timecols = master[, names(which(sapply(.SD, is_posixct)))]
  print(master[, ..timecols])
  is_posixct = function(x) {
    inherits(x, "POSIXct")
  }
  
  master[, (timecols) := lapply(.SD, function(x) {
    as.POSIXct(format(x, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
               tz = "CET")
  }), .SDcols = timecols]
  
  print(master[, ..timecols])
}


#helper function to generate a row for a consort diagram
make_consort_row = function(data, name, patname = "StudyID"){
  data.table(name = name,
             nadm = nrow(data),
             npat = length(unique(data[[patname]])))
}

if (FALSE){
  #a smaller master set for now, we don't need 300 columns
  master = copy(navidt1)[, .(StudyID = studyid,
                             PatientID = patientid,
                             VtfHuvudId = vtfhuvudid,
                             RowNo,
                             AdmissionTime = start,
                             DischargeTime = stop,
                             age
  )]
  
  master[, rowno := 1:.N, .(StudyID, AdmissionTime)]
  master[rowno == 1, spellno := 1:.N, StudyID]
  master[, spellno := nafill(spellno, "locf"), StudyID]
  master[, rowno := NULL]
  
  unique(master[, .(StudyID, AdmissionTime, DischargeTime)])
}


reload_data_files = TRUE
#first, compute the 'spells' - this will be useful later on

m2 = copy(master)
setkey(m2, spell_id, icuadm_time)
#first row per spell
m2 = m2[, .SD[1], spell_id]

if (FALSE){
  m2 = copy(master)[, .(StudyID,
                        PatientID,
                        VtfHuvudId,
                        Avd,
                        Alder,
                        Kon,
                        Akutinlaggning,
                        Opererad,
                        Moderklinik,
                        InskrTidpunkt,
                        UtskrTidpunkt,
                        AnkomstVag,
                        UtskrivenTill,
                        VardResultat,
                        AvlidenTid)]
  m2[, InskrTidpunkt := as.POSIXct(InskrTidpunkt, tz = "CET")]
  m2[, UtskrTidpunkt := as.POSIXct(UtskrTidpunkt, tz = "CET")]
  
  #sort on time
  setkey(m2, StudyID, InskrTidpunkt, UtskrTidpunkt)
  m2[, time_since_discharge := as.numeric(InskrTidpunkt - shift(UtskrTidpunkt, 1L), units = "hours"), StudyID]
  
  #new spell = more than 24 hours since last discharge or no last discharge recorded
  m2[ time_since_discharge > 24 | is.na(time_since_discharge), SpellNumber := 1:.N, StudyID]
  m2[, SpellNumber := nafill(SpellNumber, "locf"), StudyID]
  
  setcolorder(m2, c("StudyID", "SpellNumber"))
  setkey(m2, StudyID, SpellNumber, InskrTidpunkt, UtskrTidpunkt)
  
  #now, m2 is the master dataset that contains the spells and patientIDs, this will be needed later on
  #in the SOFA compute to link PatientID data to spells
  #generate a unique ID for the spell, can be useful later on
  spellmaster_key = copy(m2)[, .(StudyID, SpellNumber, PatientID, VtfHuvudId, SpellID = as.integer(StudyID * 1000 + SpellNumber))][, lapply(.SD, as.integer)]
  #add the times per part of the spell
  spellmaster_key = cbind(spellmaster_key, m2[, .(AdmissionTime = InskrTidpunkt, DischargeTime = UtskrTidpunkt)])
  fwrite(spellmaster_key[, lapply(.SD, as.character)], file.choose())
  
  #it's a rare opportunity when 'spell_master' is a suitable variable name
  #it simply can't go to waste
  spell_master = m2[, .(AdmissionTime = min(InskrTidpunkt),
                        DischargeTime = max(UtskrTidpunkt),
                        AdmissionWard = first(Avd),
                        DischargeWard = last(Avd),
                        VtfHuvudId = first(VtfHuvudId),
                        Akutinlaggning = first(Akutinlaggning),
                        Opererad = first(Opererad),
                        Moderklinik = first(Moderklinik),
                        AdmissionAge = first(Alder),
                        Sex = first(Kon),
                        AnkomstVag = first(AnkomstVag),
                        UtskrivenTill = last(UtskrivenTill),
                        VardResultat = last(VardResultat),
                        AvlidenTid = last(AvlidenTid)), .(StudyID, SpellNumber)]
}

#it's a rare opportunity when 'spell_master' is a suitable variable name
#it simply can't go to waste
spell_master = copy(m2)[!is.na(StudyID)]

#fetch more ages that might be known
spell_master = merge(spell_master, 
                     master_old[StudyID %in% spell_master[is.na(Alder), StudyID] & !is.na(Alder), .(StudyID, Age2 = Alder)],
                     by = "StudyID", all.x = TRUE)

spell_master[, `:=` (Alder = fcoalesce(Alder, as.integer(Age2)),
                     Age2 = NULL)]




#generate data to make a consort diagram later on
consortlist = list()
consortlist[[length(consortlist)+1]] = make_consort_row(spell_master, "total")

#remove the kids
spell_master = spell_master[Alder >= 18]
consortlist[[length(consortlist)+1]] = make_consort_row(spell_master, "children_removed")

#load the mortality datafile
if(reload_data_files){
  mortality_path = "P:\\GLUC-ICU\\data\\raw\\2022-11\\SIR_Mortalitet_221108.xlsx"
  mortality = read_excel(mortality_path, col_types = c("numeric",
                                                       "numeric",
                                                       "text",
                                                       "numeric",
                                                       "date",
                                                       "numeric",
                                                       "numeric"))
  setDT(mortality)
}


#only keeping one death entry per patient
mortality = unique(mortality[, .(StudyID, SIR_AvlidenTid, Avliden, AvlidenDatum)])
setkey(mortality, StudyID)

#helper function to find the first non-na-entry in a vector
first_non_na = function(x){
  if(all(is.na(x))) return(x[1]) else return(x[!is.na(x)][1])
}

#for multiple rows, fetch the one with the true time for death
mortality[SIR_AvlidenTid == "NULL", SIR_AvlidenTid := NA_character_]

#go with the first stated death date if there are multiple
#and also remove the duplicate rows of SIR_Avliden_Tid when relevant
setkey(mortality, StudyID, AvlidenDatum)
mortality = mortality[,lapply(.SD, first_non_na), StudyID]
setkey(mortality, StudyID)
mortality[, has_known_outcome := 1L]

#add the outcomes
spell_master = merge(spell_master, mortality, all.x = TRUE)
spell_master = spell_master[has_known_outcome == 1]

#check valid personal numbers
#valid personal ids
if(reload_data_files){
  valid_id_path =  "P:\\GLUC-ICU\\data\\raw\\2023-09\\GLUC-ICU StudyID som har valida pnr 20230822.xlsx"
  valid_id_dt = read_excel(valid_id_path)
  setDT(valid_id_dt)
}


#keep only those where outcome is known
spell_master = spell_master[has_known_outcome == 1 & StudyID %in% valid_id_dt$StudyID]
#consortlist[[length(consortlist)+1]] = make_consort_row(spell_master, "outcome_known")
#consorttable = rbindlist(consortlist)

consortlist[[length(consortlist)+1]] = make_consort_row(spell_master, "has_known_outcome")



#flag admissions from other ICUs
#don't remove completely as it may be useful for sensitivity analyses later on
spell_master[adm_source == "Other ICU", from_other_icu := 1L]
spell_master[is.na(from_other_icu), from_other_icu := 0L]
consortlist[[length(consortlist)+1]] = make_consort_row(spell_master[from_other_icu == 0], "transfers_from_other_icu_removed")

#rank admissions per patient
setkey(spell_master, StudyID, spell_start)
spell_master[, adult_spell_number := 1:.N, StudyID]

spell_master[from_other_icu == 0, adult_not_from_icu_spell_number := 1:.N, StudyID]

#spell_master = spell_master[, .SD[1], StudyID]
consortlist[[length(consortlist)+1]] = make_consort_row(spell_master[from_other_icu == 0, .SD[which.min(adult_spell_number)], StudyID], "only_first_admissions")
consorttable = rbindlist(consortlist)



#determine the true death date from all known facts

#spell_master[spell_dis_dest == "Avliden", AvlidenTid := spell_]
#spell_master[AvlidenTid == "NULL", AvlidenTid := NA_character_]

#first, flag the ICU-mortality
#spell_master[, .N, VardResultat]
spell_master[, icu_mortality := fifelse(spell_dis_dest == "Avliden",1L,0L)]
spell_master[icu_mortality == 1, AvlidenTid := spell_end]

#next determine death date with maximum precision
spell_master[, SIR_deathDateTime := as.POSIXct(SIR_AvlidenTid, tz = "UTC")]
#spell_master[is.na(SIR_deathDateTime) & !is.na(AvlidenDatum)]


#if we only know the date, we can only be certain that they were dead at the end of that day
spell_master[!is.na(AvlidenDatum), KnownDeathTime := as.POSIXct(paste0(as.character(AvlidenDatum), " 23:59:59"), tz = "CET")]
#perfer the SIR timestamps when they are known
spell_master[, DeathDateTime := fcoalesce(SIR_deathDateTime, KnownDeathTime)]

#spell_master[adult_not_from_icu_spell_number == 1]


#clean up the dataset a bit
# spell_master = spell_master[, .(StudyID,
#                                 SpellNumber,
#                                 VtfHuvudId,
#                                 AdmissionTime,
#                                 DischargeTime,
#                                 Akutinlaggning,
#                                 AdmissionWard,
#                                 DischargeWard,
#                                 AdmissionAge,
#                                 Opererad,
#                                 Moderklinik,
#                                 Sex,
#                                 AnkomstVag,
#                                 UtskrivenTill,
#                                 has_known_outcome,
#                                 from_other_icu,
#                                 adult_spell_number,
#                                 adult_not_from_icu_spell_number,
#                                 icu_mortality,
#                                 DeathDateTime)]

spell_master[is.na(has_known_outcome), has_known_outcome := 0L]
#these known dead/alive flag
spell_master[has_known_outcome == 1, dead := fifelse(!is.na(DeathDateTime),1L,0L)]
spell_master[, days_to_death := as.numeric(DeathDateTime - spell_start, units = "days")]

#compute the outcomes
spell_master[has_known_outcome == 1, `:=` (dead_at_28 = 0L,
                                           dead_at_30 = 0L,
                                           dead_at_90 = 0L,
                                           dead_at_360 = 0L)]
spell_master[has_known_outcome == 1 & !is.na(days_to_death), `:=` (dead_at_28 = fifelse(days_to_death <= 28,1L,0L),
                                                                   dead_at_30 = fifelse(days_to_death <= 30,1L,0L),
                                                                   dead_at_90 = fifelse(days_to_death <= 90,1L,0L),
                                                                   dead_at_360 = fifelse(days_to_death <= 360,1L,0L))]

#compute length of stay too
#spell_master[, length_of_stay := as.numeric(DischargeTime - AdmissionTime, units = "days")]



if(FALSE){
  #add the trauma subgrouping and SAPS points
  if(reload_data_files){
    sir_saps_path = "P:\\GLUC-ICU\\data\\raw\\2022-05\\SIR_Saps_210607.xlsx"
    sirsaps = read_excel(sir_saps_path)
    setDT(sirsaps)
  }
  
  
  spell_master = merge(spell_master, sirsaps[, .(VtfHuvudId,
                                                 PoangTot,
                                                 Trauma,
                                                 Operationstyp)], by = "VtfHuvudId", all.x = TRUE)
  
  #fetch the COVID-19 patients
  if(reload_data_files){
    siri_path = "P:\\GLUC-ICU\\data\\raw\\2022-05\\SIR_Siri_210607.xlsx"
    sirsiri = read_excel(siri_path, col_types = "text")
    setDT(sirsiri)
    #explicit conversion to posixCT with CET time zone!
    sirsiri[, AdmissionTime := datetimeutils::convert_date(as.numeric(VtfStart), type = "Excel", fraction = TRUE, tz = "CET")]
    sirsiri[, VtfHuvudId := as.numeric(VtfHuvudId)]
  }
  
  
  
  #determine the covid patients
  covid = sirsiri[Coronavirus == "COVID-19", .(VtfHuvudId, AdmissionTime, Coronavirus, AdmissionAge = as.numeric(Alder), Sex = Kon)]
  
  spell_master[covid, Coronavirus1 := Coronavirus, on = "VtfHuvudId"]
  spell_master[covid, Coronavirus2 := Coronavirus, on = c("AdmissionTime", "AdmissionAge", "Sex")]
  spell_master[, covid19 := fcoalesce(Coronavirus1, Coronavirus2)]
  spell_master[, `:=` (Coronavirus1 = NULL,
                       Coronavirus2 = NULL)]
}


#Finally, we also have to setup some sepsis3 definitions
suspected_infection_path = "P:\\GLUC-ICU\\users\\johanh\\sofa2025\\datasets\\unique_suspected_infections_251024.csv"

sidata = fread(suspected_infection_path)
#sidata[, suspected_infection_time2 := as.POSIXct(suspected_infection_time, tz = "UTC")]
#sidata[, suspected_infection_time := NULL]
sidata[, dt := NULL]
colnames(sidata)[2] = "suspected_infection_time"
sidata[, suspected_infection_time2 := suspected_infection_time]
setkey(sidata, StudyID, suspected_infection_time, suspected_infection_time2)

#setup the foverlaps
spell_master[, `:=` (minTime = spell_start - 24*3600,
                     maxTime = spell_start + 24*3600)]
setkey(spell_master, StudyID, minTime, maxTime)
spell_master = foverlaps(spell_master, sidata)
spell_master[!is.na(suspected_infection_time), suspected_infection_datetime := suspected_infection_time]
spell_master[, `:=` (suspected_infection_time = NULL,
                     suspected_infection_time2 = NULL,
                     minTime = NULL,
                     maxTime = NULL)]





#SOFA at admission will be needed to be certain of sepsis
sofa_1_path = "P:\\GLUC-ICU\\data\\processed\\2025-10\\sofa_hourly_with_imputations_and_studyID_251022.csv"
sofa1 = fread(sofa_1_path, colClasses = c("character"))
sofa1 = sofa1[, .SD[1], PatientID]
sofa1[, `:=` (adm_time = as.POSIXct(gsub("T", " ", gsub("Z", "", adm_time)), tz = "UTC"),
              sofa_time = as.POSIXct(gsub("T", " ", gsub("Z", "", sofa_time)), tz = "UTC"))]

cns = colnames(sofa1)[!grepl("time", colnames(sofa1))]
for (cn in cns){
  sofa1[, (cn) := as.numeric(.SD[[cn]])]
}
rm(cns)

setkey(sofa1, StudyID, adm_time)
spell_master = 
  merge(spell_master, 
        sofa1[, .(StudyID, spell_start = adm_time, sofa)],
        c("StudyID", "spell_start"), all.x = TRUE)

#true sepsis!
spell_master[!is.na(suspected_infection_datetime) & sofa >= 2, sepsis_at_admission := 1L]
spell_master[is.na(sepsis_at_admission), sepsis_at_admission := 0L]
if(FALSE){
  #get ccc from navid's file
  navid_t1_path = "P:\\GLUC-ICU\\users\\navid\\Study4\\data\\Raw\\SIR_spells_221108_20251021.rds"
  #navid_t1_path = "P:\\GLUC-ICU\\users\\navid\\Study4\\data\\Raw\\SIR_spells_221108_20250130.rds"
  navidt1 = readRDS(navid_t1_path)
  setDT(navidt1)
  navidt1[, .N, admsource]
  cciw = navidt1[, .(StudyID = studyid, AdmissionTime = as.POSIXct(as.character(start), tz = "CET"), PreICU_CCIunw)]
  
  colnames(navidt1)[grepl("c", colnames(navidt1))]
}


spell_comorb_path = "P:\\GLUC-ICU\\users\\johanh\\sofa2025\\datasets\\spells_with_comorbidities_and_charlston_age_251023.csv"
comorbs = fread(spell_comorb_path)
comorbs[, `:=` (StudyID = NULL,
                Alder = NULL,
                spell = NULL,
                spell_start = NULL,
                spell_stop = NULL,
                spell_hospital_discharge = NULL, 
                spell_admsource = NULL)]


spell_master = merge(spell_master,
                     comorbs,
                     "spell_id",
                     all.x = TRUE)




#fwrite(spell_master, file.choose())

#now, setup the table 1
library(table1)
tdt = copy(spell_master)#[adult_not_from_icu_spell_number == 1]


#pool the admission ways
{
  # tdt[AnkomstVag %in% c("Operation (inkluderar interventionell radiologi, skopier mm)",
  #                       "Postoperativ vård (på annan uppvakningsenhet)",
  #                       "Konvertering från vårdtyp 'Postop' på samma IVA"),
  #     `Admitted from` := "Operation, Intervention, or Postoperative Care"]
  # tdt[AnkomstVag %in% c("Vårdavdelning", "Intermediärvård"),
  #     `Admitted from` := "Ward"]
  # tdt[AnkomstVag %in% c("Annat sjukhus (ej IVA)"),
  #     `Admitted from` := "Other Hospital"]
  # tdt[AnkomstVag %in% c("Akutmottagning"),
  #     `Admitted from` := "Emergency Department"]
  
  tdt[spell_admsource %in% c("Operating theatre",
                         "Recovery"#,
                         #"Konvertering från vårdtyp 'Postop' på samma IVA"
                         ),
      `Admitted from` := "Operation, Intervention, or Postoperative Care"]
  tdt[spell_admsource %in% c("Ward"),
      `Admitted from` := "Ward"]
  tdt[spell_admsource %in% c("Other hospital"),
      `Admitted from` := "Other Hospital"]
  tdt[spell_admsource %in% c("Emergency department"),
      `Admitted from` := "Emergency Department"]
  
}
#setup types of surgery
{
  #tdt[, .N, s3_surgerytype]
  
  #tdt[, .N, Opererad]
  #tdt[, .N, Moderklinik]
}

#the subgroups
{
  tdt[, .N, spell_adm_trauma]
  tdt[, Trauma := as.factor(fifelse(is.na(spell_adm_trauma), "No", fifelse(spell_adm_trauma == 1, "Yes", "No")))]
  
  
  tdt[, `COVID-19` := as.factor(fifelse(is.na(spell_siri_covid19), "No", fifelse(spell_siri_covid19 == 1, "Yes", "No")))]
  
  tdt[, .N, spell_surgery]
  
  tdt[, `Elective Admission` := as.factor(fifelse((spell_surgery == "Ja-elektivt") & 
                                                    (spell_emergencyadm == 0),
                                                  "Elective Surgical Admission",
                                                  "Other admission"))]
  tdt[, `Sepsis at Admission` := factor(sepsis_at_admission, labels = c("No", "Yes"))]
}



#pretty names for outcomes and other basic variables
{
  tdt[, Sex := as.factor(fifelse(sex == "Male", "Male", "Female"))]
  tdt[, `Death in ICU` := factor(icu_mortality, labels = c("No", "Yes"))]
  tdt[, `Death at day 30` := factor(dead_at_30, labels = c("No", "Yes"))]
  tdt[, `Death at day 90` := factor(dead_at_90, labels = c("No", "Yes"))]
  tdt[, `Length of stay` := spell_los]
  tdt[, Age := spell_age]
  tdt[, `CCI (weighted with age)` := CCIage]
  
}





t1 = table1(~Age +
              Sex + 
              #Myocardial_infarction+
              #Cerebrovascular_disease+
              #Renal_disease+
              #Severe_liver_disease+
              #Malignancy+
              `CCI (weighted with age)`+
              `Admitted from` +
              `Trauma` +
              `COVID-19` +
              `Sepsis at Admission` +
              `Length of stay`+
              `Death in ICU` + 
              `Death at day 30`+
              `Death at day 90` 
            | `Elective Admission`
            , data = tdt[adult_not_from_icu_spell_number == 1]
)

t1

t1_sepsis = table1(~Age +
                     Sex + 
                     `Admitted from` +
                     `Length of stay`+
                     `Trauma` +
                     `COVID-19` +
                     `Elective Admission` +
                     `Death in ICU` + 
                     `Death at day 30`+
                     `Death at day 90` 
                   | `Sepsis at Admission`
                   , data = tdt#[adult_not_from_icu_spell_number == 1]
)

t1_sepsis
fwrite(spell_master, file.choose())
