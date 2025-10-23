#2025-10-21
#code for a table one for the trial, and for general dataset preparation

#master datapath
master_datapath = "P:\\GLUC-ICU\\data\\raw\\2022-11\\SIR_Masterfile_221108.xlsx"
master = read_excel(master_datapath)
setDT(master)

#helper function to generate a row for a consort diagram
make_consort_row = function(data, name, patname = "StudyID"){
  data.table(name = name,
             nadm = nrow(data),
             npat = length(unique(data[[patname]])))
}


reload_data_files = FALSE
#first, compute the 'spells' - this will be useful later on
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

#generate data to make a consort diagram later on
consortlist = list()
consortlist[[length(consortlist)+1]] = make_consort_row(spell_master, "total")

#remove the kids
spell_master = spell_master[AdmissionAge >= 18]
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
spell_master[AnkomstVag == "Annan IVA", from_other_icu := 1L]
spell_master[is.na(from_other_icu), from_other_icu := 0L]
consortlist[[length(consortlist)+1]] = make_consort_row(spell_master[from_other_icu == 0], "transfers_from_other_icu_removed")

#rank admissions per patient
setkey(spell_master, StudyID, AdmissionTime)
spell_master[, adult_spell_number := 1:.N, StudyID]

spell_master[from_other_icu == 0, adult_not_from_icu_spell_number := 1:.N, StudyID]

#spell_master = spell_master[, .SD[1], StudyID]
consortlist[[length(consortlist)+1]] = make_consort_row(spell_master[from_other_icu == 0, .SD[which.min(adult_spell_number)], StudyID], "only_first_admissions")
consorttable = rbindlist(consortlist)



#determine the true death date from all known facts
spell_master[AvlidenTid == "NULL", AvlidenTid := NA_character_]

#first, flag the ICU-mortality
spell_master[, .N, VardResultat]
spell_master[, icu_mortality := fifelse(VardResultat == "Avliden",1L,0L)]

#next determine death date with maximum precision
spell_master[, SIR_deathDateTime := as.POSIXct(SIR_AvlidenTid, tz = "CET")]
spell_master[is.na(SIR_deathDateTime) & !is.na(AvlidenDatum)]


#if we only know the date, we can only be certain that they were dead at the end of that day
spell_master[!is.na(AvlidenDatum), KnownDeathTime := as.POSIXct(paste0(as.character(AvlidenDatum), " 23:59:59"), tz = "CET")]
#perfer the SIR timestamps when they are known
spell_master[, DeathDateTime := fcoalesce(SIR_deathDateTime, KnownDeathTime)]

spell_master[adult_not_from_icu_spell_number == 1]


#clean up the dataset a bit
spell_master = spell_master[, .(StudyID,
                                SpellNumber,
                                VtfHuvudId,
                                AdmissionTime,
                                DischargeTime,
                                Akutinlaggning,
                                AdmissionWard,
                                DischargeWard,
                                AdmissionAge,
                                Opererad,
                                Moderklinik,
                                Sex,
                                AnkomstVag,
                                UtskrivenTill,
                                has_known_outcome,
                                from_other_icu,
                                adult_spell_number,
                                adult_not_from_icu_spell_number,
                                icu_mortality,
                                DeathDateTime)]

spell_master[is.na(has_known_outcome), has_known_outcome := 0L]
#these known dead/alive flag
spell_master[has_known_outcome == 1, dead := fifelse(!is.na(DeathDateTime),1L,0L)]
spell_master[, days_to_death := as.numeric(DeathDateTime - AdmissionTime, units = "days")]

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
spell_master[, length_of_stay := as.numeric(DischargeTime - AdmissionTime, units = "days")]



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

#Finally, we also have to setup some sepsis3 definitions
suspected_infection_path = "P:\\GLUC-ICU\\users\\anna\\sepsis\\suspected_infection\\suspected_infection_times_from_cultures_and_pharma_abs_251022.csv"
suspected_infection_path = "C:\\Users\\johhell\\OneDrive - KI.SE\\Dokument\\tte_prototyping\\sidata_unique_times_251022.csv"

sidata = fread(suspected_infection_path)
sidata[, suspected_infection_time2 := as.POSIXct(suspected_infection_time, tz = "CET")]
sidata[, suspected_infection_time := NULL]
sidata[, dt := NULL]
colnames(sidata)[2] = "suspected_infection_time"
sidata[, suspected_infection_time2 := suspected_infection_time]
setkey(sidata, StudyID, suspected_infection_time, suspected_infection_time2)

#setup the foverlaps
spell_master[, `:=` (minTime = AdmissionTime - 24*3600,
                     maxTime = AdmissionTime + 24*3600)]
setkey(spell_master, StudyID, minTime, maxTime)
spell_master = foverlaps(spell_master, sidata)
spell_master[!is.na(suspected_infection_time), suspected_infection_datetime := suspected_infection_time]
spell_master[, `:=` (suspected_infection_time = NULL,
                     suspected_infection_time2 = NULL,
                     minTime = NULL,
                     maxTime = NULL)]





#SOFA at admission will be needed to be certain of sepsis
sofa_1_path = "P:\\GLUC-ICU\\data\\processed\\2025-09\\sofa_scores\\sofa_hourly_imputed.csv"
sofa_1_path = "P:\\GLUC-ICU\\data\\processed\\2025-09\\sofa_scores\\sofa_highres_15_min_2025_09_08_08_55_41.csv"
sofa_1_path = "C:\\Users\\johhell\\OneDrive - KI.SE\\Dokument\\sofa2025\\sofa_hourly_with_imputations_and_studyID_251022.csv"
sofa1 = fread(sofa_1_path, colClasses = c("character"))
sofa1 = sofa1[, .SD[1], PatientID]
sofa1[, `:=` (adm_time = as.POSIXct(gsub("T", " ", gsub("Z", "", adm_time)), tz = "CET"),
              sofa_time = as.POSIXct(gsub("T", " ", gsub("Z", "", sofa_time)), tz = "CET"))]

cns = colnames(sofa1)[!grepl("time", colnames(sofa1))]
for (cn in cns){
  sofa1[, (cn) := as.numeric(.SD[[cn]])]
}
rm(cns)

setkey(sofa1, StudyID, adm_time)
spell_master = 
  merge(spell_master, 
        sofa1[, .(StudyID, AdmissionTime = adm_time, sofa)],
        c("StudyID", "AdmissionTime"), all.x = TRUE)

#true sepsis!
spell_master[!is.na(suspected_infection_datetime) & sofa >= 2, sepsis_at_admission := 1L]
spell_master[is.na(sepsis_at_admission), sepsis_at_admission := 0L]






#now, setup the table 1
library(table1)
tdt = copy(spell_master)#[adult_not_from_icu_spell_number == 1]


#pool the admission ways
{
  tdt[AnkomstVag %in% c("Operation (inkluderar interventionell radiologi, skopier mm)",
                        "Postoperativ vård (på annan uppvakningsenhet)",
                        "Konvertering från vårdtyp 'Postop' på samma IVA"),
      `Admitted from` := "Operation, Intervention, or Postoperative Care"]
  tdt[AnkomstVag %in% c("Vårdavdelning", "Intermediärvård"),
      `Admitted from` := "Ward"]
  tdt[AnkomstVag %in% c("Annat sjukhus (ej IVA)"),
      `Admitted from` := "Other Hospital"]
  tdt[AnkomstVag %in% c("Akutmottagning"),
      `Admitted from` := "Emergency Department"]
  
}
#setup types of surgery
{
  tdt[, .N, Opererad]
  tdt[, .N, Moderklinik]
}

#the subgroups
{
  tdt[, Trauma := as.factor(fifelse(is.na(Trauma), "No", fifelse(Trauma == "Trauma", "Yes", "No")))]
  tdt[, `COVID-19` := as.factor(fifelse(is.na(covid19), "No", fifelse(covid19 == "COVID-19", "Yes", "No")))]
  tdt[, `Elective Admission` := as.factor(fifelse((Opererad == "Ja-elektivt") & 
                                                    (Akutinlaggning == "Nej"),
                                                  "Elective Surgical Admission",
                                                  "Other admission"))]
  tdt[, `Sepsis at Admission` := factor(sepsis_at_admission, labels = c("No", "Yes"))]
}



#pretty names for outcomes and other basic variables
{
  tdt[, Sex := as.factor(fifelse(Sex == "M", "Male", "Female"))]
  tdt[, `Death in ICU` := factor(icu_mortality, labels = c("No", "Yes"))]
  tdt[, `Death at day 30` := factor(dead_at_30, labels = c("No", "Yes"))]
  tdt[, `Death at day 90` := factor(dead_at_90, labels = c("No", "Yes"))]
  tdt[, `Length of stay` := length_of_stay]
  tdt[, Age := AdmissionAge]
}





t1 = table1(~Age +
              Sex + 
              `Admitted from` +
              `Length of stay`+
              `Trauma` +
              `COVID-19` +
              `Sepsis at Admission` + 
              `Death in ICU` + 
              `Death at day 30`+
              `Death at day 90` 
            | `Elective Admission`
            , data = tdt[adult_not_from_icu_spell_number == 1]
            )

t1

