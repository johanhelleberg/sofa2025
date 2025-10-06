#2025-10-06
#Time to non-inferiority analysis for SOFA
library(data.table)
library(pROC)

file_path = "P:\\GLUC-ICU\\users\\johanh\\sofa2025\\datasets\\sofa_test_dataset_hourly_with_outcome_251006.csv"
sofadt = fread(file_path)
rm(file_path)

#compute the max
sofadt[, sofa_max := max(sofa), PatientID]

#subset only what we'll need for faster bootstraps
sdt = sofadt[, .(pid = PatientID, icu_hour, sofa)]

#setup the code for non-inferiority bootstrapping
delta = 0.02
nboots = 1000
pids = unique(sdt$pid)

#setup the meta table
pid_meta = unique(sofadt[, .(pid = PatientID, sofa_max, death_at_30)])

#helper function for auc computations
compute_aucdiff_delong = function(sdt_hour, pid_meta) {
  
  dt = copy(sdt_hour)
  #join with the data
  dt[pid_meta, `:=` (sofa_max = sofa_max,
                     death_at_30 = death_at_30), on = "pid"]
  
  #drop NAs
  dt = dt[!is.na(sofa) & !is.na(sofa_max) & !is.na(death_at_30)]
  if (nrow(dt) < 2L || length(unique(dt$death_at_30)) < 2L)
    return(list(delta=NA_real_, se=NA_real_, auc_t=NA_real_, auc_max=NA_real_, p_value = NA_real_))
  
  
  #setup the roc curves
  roc_t   = roc(response = dt$death_at_30, predictor = dt$sofa, quiet = TRUE)
  roc_max = roc(response = dt$death_at_30, predictor = dt$sofa_max, quiet = TRUE)
  
  #fetch the aucs
  auc_t   = as.numeric(auc(roc_t))
  auc_max = as.numeric(auc(roc_max))
  delta   = auc_t - auc_max
  
  #DeLong paired test, Z statistic = SE of the difference
  dltest = roc.test(roc_t, roc_max, paired = TRUE, method = "delong")
  z = unname(as.numeric(dltest$statistic)) 
  
  se = if (is.finite(z) && z != 0) abs(delta / z) else NA_real_
  
  list(delta = delta, se = se, auc_t = auc_t, auc_max = auc_max, p_value = dltest$p.value)
}

#helper for bootstrap resampling (stratified by survivors or not)
bootstrap_resample = function(sdt, pid_meta) {
  pos_ids = pid_meta[death_at_30 == 1, pid]
  neg_ids = pid_meta[death_at_30 == 0, pid]
  
  boot_ids = c(sample(pos_ids, length(pos_ids), replace = TRUE),
               sample(neg_ids, length(neg_ids), replace = TRUE))
  
  #duplicate the patients whole trajectory, note : quite slow
  dt_list = lapply(boot_ids, function(id) sdt[pid == id])
  sdt_b = rbindlist(dt_list, use.names = TRUE, fill = FALSE)
  
  # Build matching pid_meta for the resample (with duplicates)
  pm_list = lapply(boot_ids, function(id) pid_meta[pid == id])
  pid_meta_b = rbindlist(pm_list, use.names = TRUE, fill = FALSE)
  
  list(sdt = sdt_b, pid_meta = pid_meta_b)
}