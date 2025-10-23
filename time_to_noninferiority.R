#2025-10-06
#Time to non-inferiority analysis for SOFA
library(data.table)
library(pROC)

library(Rcpp)
library(readxl)

cpp_boot_path = "P:\\GLUC-ICU\\users\\johanh\\gluc_icu_general\\src\\fast_R_bootstrap.cpp"
sourceCpp(cpp_boot_path)

file_path = "P:\\GLUC-ICU\\users\\johanh\\sofa2025\\datasets\\sofa_test_dataset_hourly_with_outcome_251006.csv"
sofadt = fread(file_path)
rm(file_path)

#get the outcomes too
master_datapath = "P:\\GLUC-ICU\\users\\navid\\HMT data\\data\\csv\\Demographics_January2025.csv"
master = fread(master_datapath)
setkey(master, studyid, patientid, start)

#master datapath
master_datapath = "P:\\GLUC-ICU\\data\\raw\\2022-11\\SIR_Masterfile_221108.xlsx"
master = read_excel(master_datapath)
setDT(master)

master = master[Alder >= 18 & AnkomstVag != "Annan IVA", 
                .(StudyID, 
                  PatientID, 
                  Alder, 
                  InskrTidpunkt, 
                  UtskrTidpunkt, 
                  AnkomstVag, 
                  VardResultat, 
                  AvlidenTid)]

master[AvlidenTid != "NULL", `:=` (DeathDateTime = as.POSIXct(AvlidenTid, tz = "CET"))]
master[, InskrTidpunkt := as.POSIXct(InskrTidpunkt, tz = "CET")]
master[, UtskrTidpunkt := as.POSIXct(UtskrTidpunkt, tz = "CET")]


setkey(master, StudyID, InskrTidpunkt)
master = master[, .SD[1], StudyID]

mortality_path = "P:\\GLUC-ICU\\data\\raw\\2022-11\\SIR_Mortalitet_221108.xlsx"
mortality = read_excel(mortality_path, col_types = c("numeric",
                                                     "numeric",
                                                     "text",
                                                     "numeric",
                                                     "date",
                                                     "numeric",
                                                     "numeric"))
setDT(mortality)




library(RSQLite)
db_path = "P:\\GLUC-ICU\\users\\johanh\\sql\\db_v9_updated_pharma_with_pharmaID.sqlite"
db = dbConnect(SQLite(), db_path)



master2 = data.table(dbReadTable(db, "sir_master"))
setkey(master2, StudyID, PatientID, AdmissionDateTime)
pids = master2[AdmittedFrom != "Other ICU"][, .SD[1], StudyID]$PatientID



#setup the code for non-inferiority bootstrapping
delta = 0.02
nboots = 1000
#pids = unique(sofadt$PatientID)
#for faster prototyping!
sample_fraction = 1
pids = sample(pids, length(pids)*sample_fraction, replace = FALSE)

#subset only what we'll need for faster bootstraps
sdt = sofadt[PatientID %in% pids, .(pid = PatientID, icu_hour, sofa)]

#setup the meta table
pid_meta = unique(sofadt[PatientID %in% pids, .(pid = PatientID, sofa_max, death_at_30)])

#helper function for auc computations
compute_aucdiff_delong = function(sdt_hour, pid_meta) {
  
  dt = copy(sdt_hour)
  #join with the data
  dt[unique(pid_meta), `:=` (sofa_max = sofa_max,
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
  
  list(aucdelta = delta, se = se, auc_t = auc_t, auc_max = auc_max, p_value = dltest$p.value)
}

#helper for bootstrap resampling (stratified by survivors or not)
bootstrap_resample = function(sdt, pid_meta) {
  pos_ids = pid_meta[death_at_30 == 1, pid]
  neg_ids = pid_meta[death_at_30 == 0, pid]
  
  boot_ids = c(sample(pos_ids, length(pos_ids), replace = TRUE),
               sample(neg_ids, length(neg_ids), replace = TRUE))
  
  #duplicate the patients whole trajectory, note : quite slow
  #dt_list = lapply(boot_ids, function(id) sdt[pid == id])
  #sdt_b = rbindlist(dt_list, use.names = TRUE, fill = FALSE)
  sdt_b = cpp_integer_bootstrap_set(sdt, boot_ids, n_threads = 8)
  setDT(sdt_b)
  
  # Build matching pid_meta for the resample (with duplicates)
  #pm_list = lapply(boot_ids, function(id) pid_meta[pid == id])
  #pid_meta_b = rbindlist(pm_list, use.names = TRUE, fill = FALSE)
  pid_meta_b = cpp_integer_bootstrap_set(pid_meta, boot_ids)
  setDT(pid_meta_b)
  
  list(sdt = sdt_b, pid_meta = pid_meta_b)
}


#main function, compute a band for non-inferiority through bootstrapping
sofa_noninferiority_band = function(sdt,
                                    pid_meta,
                                    delta_cutoff = 0.02,
                                    nboots = 1000,
                                    alpha = 0.05,
                                    min_pos = 10L,
                                    min_neg = 10L,
                                    min_hour = NULL,
                                    max_hour = NULL,
                                    seed = 1L,
                                    verbose = TRUE) {
  #sanity checks
  stopifnot(is.data.table(sdt))
  stopifnot(is.data.table(pid_meta))
  req_cols = c("pid", "icu_hour", "sofa")
  if (!all(req_cols %in% names(sdt))) {
    stop("sdt must have columns: ", paste(req_cols, collapse = ", "))
  }
  
  req_cols = c("pid", "sofa_max", "death_at_30")
  if (!all(req_cols %in% names(pid_meta))) {
    stop("pid_meta must have columns: ", paste(req_cols, collapse = ", "))
  }
  
  #allow limits for hours
  if (is.null(min_hour)) min_hour = min(sdt$icu_hour)
  if (is.null(max_hour)) max_hour = max(sdt$icu_hour)
  
  #Ensure integer types, might be useful for later shift to C++
  dt = copy(sdt)[icu_hour %between% list(min_hour, max_hour), lapply(.SD, as.integer)]
  mdt = copy(pid_meta)[, lapply(.SD, as.integer)]
  setkey(dt, pid, icu_hour)
  setkey(mdt, pid)
  
  #Candidate hours with enough cases/controls
  by_hour_dt = merge(dt, mdt, all = TRUE)[, .(
    n = .N,
    pos = sum(death_at_30 == 1L, na.rm = TRUE),
    neg = sum(death_at_30 == 0L, na.rm = TRUE)), by = icu_hour][
      pos >= min_pos & neg >= min_neg][
        order(icu_hour)]
  #stop if no time points remain for analysis :)
  if (nrow(by_hour_dt) == 0) stop("No hours meet min_pos/min_neg criteria.")
  
  #Original-sample estimates (delta_hat, se_hat) for all hours
  hours = unique(by_hour_dt$icu_hour)
  
  delonglist = lapply(hours, function(h) {
    dt_h = dt[icu_hour == h, .(pid, sofa)]
    aucdiff = compute_aucdiff_delong(dt_h, mdt)
    c(hour = h, aucdiff)
  })
  base_dt = rbindlist(lapply(delonglist, as.data.table), fill = TRUE)
  
  #Failsafe to make sure the delta and se estimates are usable
  base_dt = base_dt[is.finite(aucdelta) & is.finite(se)]
  if (!nrow(base_dt)) stop("No usable hours after computing DeLong SEs.")
  hours = base_dt$hour
  
  #Bootstrap to get max-T critical value c_alpha
  if (verbose) message("Bootstrap (stratified by outcome), grab a coffee")
  set.seed(seed)
  M = numeric(nboots)
  M[] = NA_real_
  
  for (b in 1:nboots) {
    #resample (stratified by outcome), across all hours
    bs = bootstrap_resample(dt, mdt)
    dt_b = bs$sdt
    mdt_b = bs$pid_meta
    
    #For each hour, compute Z_b(t) = (delta_b(t) - delta_hat(t)) / SE_b(t)
    Zb = numeric(length(hours)); 
    Zb[] = NA_real_
    for (i in seq_along(hours)) {
      h = hours[i]
      dt_h_b = dt_b[icu_hour == h, .(pid, sofa)]
      bootdiff = compute_aucdiff_delong(dt_h_b, mdt_b)
      #fail-check to see that the equation will be numeric
      if (is.finite(bootdiff$aucdelta) && is.finite(bootdiff$se) && bootdiff$se > 0) {
        Zb[i] = (bootdiff$aucdelta - base_dt[hour == h]$aucdelta) / bootdiff$se
      }
    }
    
    #Max over available times for this bootstrap replicate
    if (any(is.finite(Zb))) {
      M[b] = max(Zb[is.finite(Zb)])
    }
    #prompt what bootstrap we're on
    if (verbose && b %% max(1L, floor(nboots/10)) == 0L) message("  Boot ", b, "/", nboots)
  }
  
  M =  M[is.finite(M)]
  if (!length(M)) stop("All bootstrap replicates were degenerate; consider reducing min_pos/min_neg or using larger data.")
  c_alpha = as.numeric(quantile(M, probs = 1 - alpha, type = 8, names = FALSE))
  
  
  #Uniform lower band and NI decision boundary
  base_dt[, lower_band := aucdelta - c_alpha * se]
  base_dt[, NI := lower_band >= (-delta_cutoff)]
  
  # earliest NI hour (if any)
  earliest_NI = base_dt[NI == TRUE, min(hour, na.rm = TRUE)]
  if (!is.finite(earliest_NI)) earliest_NI = NA_integer_
  
  out = list(
    grid        = base_dt[order(hour)],
    c_alpha     = c_alpha,
    alpha       = alpha,
    delta_margin= delta_cutoff,
    earliest_NI = earliest_NI,
    bootstrap_M = M
  )
  
  out
}


pdata = copy(out$grid)



library(scales)
library(patchwork)
library(ggplot2)

#95% Wald CI for ΔAUC
pdata[, `:=`(
  delta_ci_lo = aucdelta - qnorm(0.975) * se,
  delta_ci_hi = aucdelta + qnorm(0.975) * se
)]

p_top = ggplot(pdata, aes(x = hour)) +
  geom_ribbon(aes(x = hour, ymin = pmin(auc_t, auc_max), ymax = pmax(auc_t, auc_max)),
              fill = "grey95", inherit.aes = FALSE) +
  geom_line(aes(y = auc_t), linewidth = 1) +
  geom_point(aes(y = auc_t), size = 1.8) +
  geom_line(aes(y = auc_max), linetype = 2, linewidth = 1) +
  labs(y = "AUC", x = NULL,
       title = "SOFA AUC over time",
       subtitle = "Solid = AUC(t), Dashed = AUC(max over 0–7d)") +
  scale_x_continuous(breaks = pretty_breaks(8)) +
  coord_cartesian(ylim = c(0.68, 0.85)) +
  theme_minimal(base_size = 12) +
  theme(plot.title.position = "plot",
        panel.grid.minor = element_blank())

p_top

p_bot = ggplot(pdata, aes(x = hour)) +
  # CI for ΔAUC
  geom_ribbon(aes(x = hour, ymin = delta_ci_lo, ymax = delta_ci_hi),
              alpha = 0.2) +
  # ΔAUC curve
  geom_line(aes(y = aucdelta), linewidth = 1) +
  geom_point(aes(y = aucdelta), size = 1.6) +
  # Uniform lower band L(t)
  geom_line(aes(y = lower_band), linewidth = 1, linetype = 3) +
  # -δ threshold
  geom_hline(yintercept = -delta_cutoff, linetype = 2) +
  # Highlight hours meeting NI (L(t) >= -δ)
  geom_point(data = pdata[NI == TRUE],
             aes(y = aucdelta), size = 2.8, shape = 21, stroke = 0.8, fill = "white") +
  labs(x = "ICU hour", y = expression(Delta*" AUC (AUC(t) - AUC"[max]*")"),
       title = "ΔAUC with 95% CI and non-inferiority band",
       subtitle = expression(paste("Ribbon: 95% CI for ", Delta, "; dotted: -", delta, 
                                   "; dash-dot: uniform lower band L(t)"))) +
  scale_x_continuous(breaks = pretty_breaks(8)) +
  scale_y_continuous(labels = number_format(accuracy = 0.001)) +
  theme_minimal(base_size = 12) +
  theme(plot.title.position = "plot",
        panel.grid.minor = element_blank())

p_bot

p_top / p_bot + plot_layout(heights = c(2, 2.2))
