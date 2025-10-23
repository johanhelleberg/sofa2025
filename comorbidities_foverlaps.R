library(readstata13)
library(data.table)

starttime = Sys.time()
spell_path = "P:\\GLUC-ICU\\users\\navid\\Study2\\data\\Raw\\SIR_spells_221108.dta"
df = readstata13::read.dta13(spell_path, 
                             convert.dates = T, 
                             convert.factors = T, 
                             nonint.factors = T)

setDT(df)
dt = copy(df)

dt = dt[, .(StudyID, Alder, spell, spell_admsource, spell_start, spell_stop)]
dt = dt[!is.na(StudyID)]
dt = unique(dt)

#helper to find datetime columns
is_posixct = function(x) {
  inherits(x, "POSIXct")
}
timecols = dt[, names(which(sapply(.SD, is_posixct)))]

#force attribute to UTC for foverlaps later on
dt[, (timecols) := lapply(.SD, function(x) {
  as.POSIXct(format(x, "%Y-%m-%d %H:%M:%S", tz = "GMT"),
             tz = "UTC")
}), .SDcols = timecols]


#ICU discharge date + 30 days becomes the hospital discharge date after consensus with JM 30/5-24
dt[, spell_hospital_discharge := spell_stop + 30*24*3600]

#generate a spell_id for the joins later on
dt[, spell_id := as.integer(StudyID * 1000 + spell)]

#for this purpose, it's basically how far to look back
dt[is.na(Alder), Alder := 100]
dt = dt[, .SD[which.max(Alder)], spell_id]

comorb_path = "P:\\GLUC-ICU\\users\\navid\\Study2\\data\\Raw\\Gluc_Icu_TC_Diagnoser_240520.csv"
dt_comorb = fread(comorb_path, encoding = "Latin-1")


#starttime = Sys.time()
#comorbidity DT
dtc = copy(dt_comorb)[, .(StudyID,
                          DiagnosKod,
                          Inskrivningstidpunkt,
                          Utskrivningstidpunkt,
                          Dodsorsak)]
#like Navid did
dtc = dtc[Dodsorsak == FALSE]


#store
dt_backup = copy(dt)
dtc_backup = copy(dtc)


#prototyping, small subset
#subsetselect = sample(unique(df[spell_count>1]$StudyID), 1000, replace = FALSE)

starttime = Sys.time()
dt = unique(copy(dt_backup))#[StudyID %in% subsetselect]
dtc = unique(copy(dtc_backup))#[StudyID %in% subsetselect]


#let's go with a massive foverlap here...

dt = dt[, .(StudyID, spell_id, spell_hospital_discharge, Alder)]

#time for overlap = lifetime of patient + one year
dt[, `:=` (mintime = spell_hospital_discharge - 3600*24*365*(fifelse(is.na(Alder), 100, Alder + 1)),
           maxtime = spell_hospital_discharge)]
setkey(dt, StudyID, mintime, maxtime)
dtc = dtc[, .(StudyID, Utskrivningstidpunkt, DiagnosKod)]
dtc[, `:=` (mintime = Utskrivningstidpunkt,
            maxtime = Utskrivningstidpunkt)]
dtc = na.omit(dtc)
setkey(dtc, StudyID, mintime, maxtime)

#dtf will contain comorbidities per spell
dtf = foverlaps(dtc, dt, mult = "all")[!is.na(spell_hospital_discharge), .(StudyID, spell_id, spell_hospital_discharge, Utskrivningstidpunkt, DiagnosKod)]

#some more logical names
colnames(dtf)[4] = "datum"
colnames(dtf)[5] = "code"

#helper function to determine which icd era to use based on dates
icd_era = function(x, 
                   cutoffs = c("icd8" = as.POSIXct("1969-01-01 00:00:00", tz = "UTC"),
                               "icd9" = as.POSIXct("1987-01-01 00:00:00", tz = "UTC"),
                               "icd10" = as.POSIXct("1997-01-01 00:00:00", tz = "UTC"))
                   , default_era = "icd7")
{
  out = rep(default_era, length(x))
  for (cutoff in names(cutoffs)){
    out[x>cutoffs[[cutoff]]] = cutoff
  }
  
  return(out)
}

#add the correct era to these comorbidities
dtf[, era := icd_era(datum)]

#sort codes by era and code
setkey(dtf, era, code)
#setup list of regex matching strings per era
{
  dlist = list()
  dlist[["Myocardial_infarction"]] = c(icd7  = "\\<4201",
                                       icd8  = "\\<410|\\<411|\\<41201|\\<41291",
                                       icd9  = "\\<410|\\<412",
                                       icd10 = "\\<I21|\\<I22|\\<I252")
  
  dlist[["Congestive_heart_failure"]] = c(
    icd7  = "\\<42221|\\<42222|\\<4341|\\<4342",
    icd8  = "\\<42508|\\<42509|\\<4270|\\<4271|\\<428",
    icd9  = paste(c("\\<402A", "402B", "402X", "404A","404B","404X","425E","425F","425H","425W","425X","428"),collapse="|\\<"),
    icd10 = "\\<I110|\\<I130|\\<I132|\\<I255|\\<I420|\\<I426|\\<I427|\\<I428|\\<I429|\\<I43|\\<I50"
  )
  
  dlist[["Peripheral_vascular_disease"]] = c(
    icd7  = "\\<4501|\\<451|\\<453",
    icd8  = "\\<440|\\<441|\\<4431|\\<4439",
    icd9  = "\\<440|\\<441|\\<443B|\\<443X|\\<447B|\\<557",
    icd10 = "\\<I70|\\<I71|\\<I731|\\<I738|\\<I739|\\<I771|\\<I790|\\<I792|\\<K55"
  )
  
  dlist[["Cerebrovascular_disease"]] = c(
    icd7  = paste(c("\\<330",331:334),collapse="|\\<"),
    icd8  = "\\<430|\\<431|\\<432|\\<433|\\<434|\\<435|\\<436|\\<437|\\<438",
    icd9  = "\\<430|\\<431|\\<432|\\<433|\\<434|\\<435|\\<436|\\<437|\\<438",
    icd10 = "\\<G45|\\<I60|\\<I61|\\<I62|\\<I63|\\<I64|\\<I67|\\<I69"
  )
  
  dlist[["Chronic_obstructive_pulmonary_disease"]] = c(
    icd7  = "\\<502|\\<5271",
    icd8  = "\\<491|\\<492",
    icd9  = "\\<491|\\<492|\\<496",
    icd10 = "\\<J43|\\<J44"
  )
  
  dlist[["Chronic_other_pulmonary_disease"]] = c(
    icd7  = paste(c("\\<241",501,523:526),collapse="|\\<"),
    icd8  = paste(c("\\<490",493,515:518),collapse="|\\<"),
    icd9  = paste(c("\\<490",493:495,500:508,516,517),collapse="|\\<"),
    icd10 = paste(c("\\<J45",41,42,46,47,60:70),collapse="|\\<J")
  )
  
  dlist[["Rheumatic_disease"]] = c(
    icd7  = paste(c("\\<72200","72201","72210","72220","72223","4560","4561","4562","4563"),collapse="|\\<"),
    icd8  = paste(c("\\<446",696,"7120","7121","7122","7123","7125", 716, "7340", "7341", "734,9"),collapse="|\\<"),
    icd9  = paste(c("\\<446","696A","710A","710B","710C","710D","710E",714,"719D",720,725),collapse="|\\<"),
    icd10 = paste(c("\\<M05","06",123,"070","071","072","073","08",13,30,313:316,32:34,350:351,353,45:46),collapse="|\\<M")
  )
  
  dlist[["Dementia"]] = c(
    icd7 = "\\<304|\\<305",
    icd8  = "\\<290",
    icd9  = "\\<290|\\<294B|\\<331A|\\<331B|\\<331C|\\<331X",
    icd10 = "\\<F00|\\<F01|\\<F02|\\<F03|\\<F051|\\<G30|\\<G311|\\<G319"
  )
  
  dlist[["Hemiplegia"]] = c(
    icd7 = "\\<351|\\<352|\\<35700",
    icd8 = "\\<343|\\<344",
    icd9 = "\\<342|\\<343|\\<344A|\\<344B|\\<344C|\\<344D|\\<344E|\\<344F",
    icd10 = "\\<G114|\\<G80|\\<G81|\\<G82|\\<G830|\\<G831|\\<G832|\\<G833|\\<G838"
  )
  
  dlist[["Diabetes_without_chronic_complication"]] = c(
    icd7 = "\\<26009",
    icd8 =  "\\<25000|\\<25007|\\<25008",
    icd9 = "\\<250A|\\<250B|\\<250C",
    icd10 = "\\<E100|\\<E101|\\<E110|\\<E111|\\<E119|\\<E120|\\<E121|\\<E130|\\<E131|\\<E140|\\<E141"
  )
  
  dlist[["Diabetes_with_chronic_complication"]] = c(
    icd7 = "\\<2602|\\<26021|\\<26029|\\<2603|\\<2604|\\<26049|\\<26099",
    icd8 = "\\<25001|\\<25002|\\<25003|\\<25004|\\<25005",
    icd9 = "\\<250D|\\<250E|\\<250F|\\<250G",
    icd10 = "\\<E102|\\<E103|\\<E104|\\<E105|\\<E107|\\<E112|\\<E113|\\<E114|\\<E115|\\<E116|\\<E117|\\<E122|\\<E123|\\<E124|\\<E125|\\<E126|\\<E127|\\<E132|\\<E133|\\<E134|\\<E135|\\<E136|\\<E137|\\<E142|\\<E143|\\<E144|\\<E145|\\<E146|\\<E147"
  )
  
  dlist[["Renal_disease"]] = c(
    icd7 = "\\<592|\\<593|\\<792",
    icd8 = "\\<582|\\<583|\\<584|\\<792|\\<593|\\<40399|\\<40499|\\<79299|\\<Y2901",
    icd9 = "\\<403A|\\<403B|\\<403X|\\<582|\\<583|\\<585|\\<586|\\<588A|\\<V42A|\\<V45B|\\<V56",
    icd10 = "\\<I120|\\<I131|\\<N032|\\<N033|\\<N034|\\<N035|\\<N036|\\<N037|\\<N052|\\<N053|\\<N054|\\<N055|\\<N056|\\<N057|\\<N11|\\<N18|\\<N19|\\<N250|\\<Q611|\\<Q612|\\<Q613|\\<Q614|\\<Z49|\\<Z940|\\<Z992"
  )
  
  dlist[["Mild_liver_disease"]] = c(
    icd7  = "\\<581",
    icd8  = "\\<070|\\<571|\\<573",
    icd9 =  "\\<070|\\<571C|\\<571E|\\<571F|\\<573",
    icd10 = "\\<B15|\\<B16|\\<B17|\\<B18|\\<B19|\\<K703|\\<K709|\\<K73|\\<K746|\\<K754"
  )
  
  dlist[["Liver_special"]] = c(
    icd8  = "\\<7853",
    icd9 = "\\<789F",
    icd10 = "\\<R18"
  )
  
  #special logic will be needed but let's just do the regex matching first
  #JH
  dlist[["Severe_liver_disease"]] = c(
    icd7 = "\\<4621",
    icd8 = "\\<4560|\\<5719|\\<57302",
    icd9 = "\\<456A|\\<456B|\\<456C|\\<572C|\\<572D|\\<572E",
    icd10 =  "\\<I850|\\<I859|\\<I982|\\<I983"
  )
  
  dlist[["Peptic_ulcer_disease"]] = c(
    icd7  = "\\<540|\\<541|\\<542",
    icd8  = "\\<531|\\<532|\\<533|\\<534",
    icd9 = "\\<531|\\<532|\\<533|\\<534",
    icd10 ="\\<K25|\\<K26|\\<K27|\\<K28"
  )
  
  dlist[["Malignancy"]] = c(
    icd7  = paste(paste("\\<",paste(140:190,collapse = "|\\<"),sep=""), paste("|\\<",paste(192:197,collapse = "|\\<"),sep=""), paste("|\\<",paste(200:204,collapse = "|\\<"),sep=""),sep=""),
    icd8  = paste(paste("\\<",paste(c(140:172,174),collapse = "|\\<"),sep=""), paste("|\\<",paste(c(180:207,209),collapse = "|\\<"),sep=""),sep=""),
    icd9  = paste(paste("\\<",paste(140:172,collapse = "|\\<"),sep=""), paste("|\\<",paste(174:208,collapse = "|\\<"),sep=""),sep=""),
    icd10 = paste("\\<C00|\\<C0",paste(1:9,collapse = "|\\<C0",sep=""),paste("|\\<C",paste(c(10:41,43,45:58,60:76,81:86,88:97),collapse = "|\\<C"),sep=""),sep="")
  )
  
  dlist[["Metastatic_solid_tumor"]] = c(
    icd7 = "\\<15691|\\<198|\\<199",
    icd8 = "\\<196|\\<197|\\<198|\\<199",
    icd9 = "\\<196|\\<197|\\<198|\\<199A|\\<199B",
    icd10 = "\\<C77|\\<C78|\\<C79|\\<C80"
  )
  
  dlist[["Aids"]] = c(
    icd9  = "\\<079J|\\<279K",
    icd10 = "\\<B20|\\<B21|\\<B22|\\<B23|\\<B24|\\<F024|\\<O987|\\<R75|\\<Z219|\\<Z717"
  )
}
setkey(dt, spell_id)
comorbmatrix = copy(dt[, 1:3])

#prepare these for fast joins...

#first up, a unique table of icd eras and codes
dtfu = unique(dtf[, .(era, code)])
setkey(dtfu, era, code)
#second, the list of all codes per spell, but now keyed on era and code
dtf2 = copy(dtf)
setkey(dtf2, era, code, spell_id)

for (diagnosis in names(dlist)){
  print(diagnosis)
  print(Sys.time())
  #fetch the regex vector for this diagnosis
  dvec = dlist[[diagnosis]]
  #this will be a list of tables of matching regexes per era
  eralist = list()
  for (thisera in names(dvec)){
    #regex in the unique set of codes per era - much faster than the full list!
    eralist[[length(eralist)+1]] = dtfu[era == thisera & grepl(dvec[[thisera]], code)]
  }
  eradt = rbindlist(eralist)
  setkey(eradt, era, code)
  #select the matching rows from the list of codes per spell
  thisdata = unique(dtf2[eradt, .(spell_id, this_diag = 1L)])
  colnames(thisdata)[2] = diagnosis
  
  setkey(thisdata, spell_id)
  #join onto the comorbmatrix
  comorbmatrix= merge(comorbmatrix, thisdata, "spell_id", all.x = TRUE)
}

#fill with 0s when missing
comorbmatrix = cbind(
  comorbmatrix[, 1:3],
  comorbmatrix[, lapply(.SD, function(x) nafill(x, "const", fill = 0)), .SDcols = 4:ncol(comorbmatrix)]
)

endtime = Sys.time()
print(endtime - starttime)
comorbmatrix[, lapply(.SD, sum), .SDcols = 4:ncol(comorbmatrix)]


#version 1 - a bit slow
# for(diagnosis in names(dlist)){
#   print(diagnosis)
#   print(Sys.time())
#   dvec = dlist[[diagnosis]]
#   eralist = list()
#   for (icd in names(dvec)){
#     thisdata = dtf[era == icd][grepl(dvec[[icd]], code), .(spell_id
#                                                            #, era, code
#     )]
#     thisdata[, (diagnosis) := 1L]
#     eralist[[length(eralist)+1]] = thisdata
#     
#     
#   }
#   eradt = unique(rbindlist(eralist))
#   setkey(eradt, spell_id)
#   comorbmatrix = merge(comorbmatrix, eradt, "spell_id", all.x = TRUE)
#   
# }
# 
# #nafill to 0 if not found with grepl
# comorbmatrix = cbind(
#   comorbmatrix[, 1:3],
#   comorbmatrix[, lapply(.SD, function(x) nafill(x, "const", fill = 0)), .SDcols = 4:ncol(comorbmatrix)]
# )

#special liver logic
#if mild + special -> severe
comorbmatrix[Mild_liver_disease == 1 & Liver_special == 1, Severe_liver_disease := 1L]
#if severe - not mild
comorbmatrix[Severe_liver_disease == 1, Mild_liver_disease := 0L]


#Navid's weighting code
Matrix = copy(comorbmatrix)
# Calculate the unweighted comorbidity index
Matrix$CCIunw <- Matrix$Myocardial_infarction + Matrix$Congestive_heart_failure + Matrix$Peripheral_vascular_disease +
  Matrix$Cerebrovascular_disease + Matrix$Chronic_obstructive_pulmonary_disease + Matrix$Chronic_other_pulmonary_disease +
  Matrix$Rheumatic_disease + Matrix$Dementia + Matrix$Hemiplegia + Matrix$Diabetes_without_chronic_complication +
  Matrix$Diabetes_with_chronic_complication + Matrix$Renal_disease + Matrix$Mild_liver_disease + Matrix$Severe_liver_disease +
  Matrix$Peptic_ulcer_disease + Matrix$Malignancy + Matrix$Metastatic_solid_tumor + Matrix$Aids

# Calculate the weighted comorbidity index
Matrix$CCIw <- Matrix$Myocardial_infarction + Matrix$Congestive_heart_failure + Matrix$Peripheral_vascular_disease +
  Matrix$Cerebrovascular_disease + Matrix$Chronic_obstructive_pulmonary_disease + Matrix$Chronic_other_pulmonary_disease +
  Matrix$Rheumatic_disease + Matrix$Dementia + 2*Matrix$Hemiplegia + Matrix$Diabetes_without_chronic_complication +
  2*Matrix$Diabetes_with_chronic_complication + 2*Matrix$Renal_disease + Matrix$Mild_liver_disease + 3*Matrix$Severe_liver_disease +
  Matrix$Peptic_ulcer_disease + 2*Matrix$Malignancy + 6*Matrix$Metastatic_solid_tumor + 6*Matrix$Aids



#for join, not really needed
setkey(Matrix, spell_id)

spells_with_comorb = merge(dt_backup, Matrix[, -c(2:3), with = FALSE], by = "spell_id", all.x = TRUE)

{
  spells_with_comorb[, CCIage := CCIw]
  spells_with_comorb[Alder >= 50, CCIage := CCIage+1]
  spells_with_comorb[Alder >= 60, CCIage := CCIage+1]
  spells_with_comorb[Alder >= 70, CCIage := CCIage+1]
  spells_with_comorb[Alder >= 80, CCIage := CCIage+1]
}

fwrite(spells_with_comorb, file.choose())

#generate a master csv with spell_id too
master = df
setDT(master)
master[, spell_id := StudyID * 1000 + spell]
setcolorder(master, c("StudyID","spell_id"))

timecols = master[, names(which(sapply(.SD, is_posixct)))]

#force attribute to UTC for foverlaps later on
master[, (timecols) := lapply(.SD, function(x) {
  as.POSIXct(format(x, "%Y-%m-%d %H:%M:%S", tz = "GMT"),
             tz = "UTC")
}), .SDcols = timecols]

fwrite(master, file.choose())
