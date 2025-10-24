#2025-10-23
#helper function for POSIXct
is_posixct = function(x) {
  inherits(x, "POSIXct")
}

#helper function to shift time columns to another time zone
#while the apparent time is unchanged
shift_dt_times = function(dt, timecols = NULL, origin_tz ="UTC", target_tz = "CET"){
  #find timecols if not provided
  if(is.null(timecols)){
    timecols = dt[, names(which(sapply(.SD, is_posixct)))]
  }
  
  
  dt[, (timecols) := lapply(.SD, function(x) {
    as.POSIXct(format(x, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
               tz = options$time$ptz)
  }), .SDcols = timecols]
  
  dt
}

#Show the timestamps
sofa_echo_time = function(event_name, starttime, prevtime, endtime, decimals =3){
  paste0(event_name, ", time spent: ", 
         round(as.double(endtime-prevtime, units = "secs"),decimals), 
         " seconds. Total time elapsed: ", 
         round(as.double(endtime-starttime, units = "secs"),decimals), 
         " seconds. Memory used: ",
         as.character(round(as.numeric(pryr::mem_used() / 1000000000),2)),
         " GB.")
}


#get result of a query as data.table
dtGetQuery = function(db, sql, silent = TRUE, print_time = FALSE){
  if (!silent){
    print(sql)
  }
  
  QStime = Sys.time()
  res = dbGetQuery(db, sql)
  QEtime = Sys.time()
  if (print_time){
    print(QEtime - QStime)
  }
  setDT(res)
  return(res)
}

#read tables as data.tables
dtReadTable = function(db, ...){
  data.table(dbReadTable(db, ...))
}


#helper function for fast left joins by reference, common in this task
#2020-10-08
#this is a fast data.table way of adding columns, I think...
#2023-03-30
#now allows 'by' to be skipped if both x and y share a key
#2025-10-24
#rewrote to skip 'eval(parse)' and renamed from legacy 'add_cols()' to the more
#clear 'ref_left_join()'
ref_left_join = function(x, y, by = NULL, cols = NULL, newcolnames = NULL,
                         mult = "all") {
  stopifnot(is.data.table(x), is.data.table(y))
  
  # Infer 'by' only if both keyed and identical
  if (is.null(by)) {
    if (!is.null(key(x)) && !is.null(key(y)) && identical(key(x), key(y))) {
      by = key(x)
    } else {
      stop("'by' must be provided unless 'x' and 'y' are both keyed with the same key")
    }
  }
  
  #default - join all non-key columns from y
  if (is.null(cols)) cols = setdiff(names(y), if (is.character(by)) unname(by) else by)
  if (length(cols) == 0L) return(invisible(x))
  
  #new col names (allows the new column in x to have a different name than i y)
  if (is.null(newcolnames)) newcolnames = cols
  if (length(cols) != length(newcolnames)) {
    stop("'cols' and 'newcolnames' must have the same length")
  }
  
  #Only keep the columns needed for the join, i.e. don't join redundant columns from y
  z = unique(c(if (is.character(by)) unique(c(names(by), unname(by))) else by, cols))
  
  #Update-join: assign from i.y to x
  x[y[, ..z],
    on = by,
    (newcolnames) := mget(paste0("i.", cols)),
    #nomatch = nomatch,
    mult = mult]
  
  invisible(x)
}



#decimal shift method for correcting input errors
decimal_shift_input = function(x, min_range, max_range, SDs = 3, oom = 4){
  
  #split into 'within' limits and 'without' limits
  xin = x[x %between% c(min_range, max_range) & !is.na(x)]
  meanin = mean(xin)
  xout = x[!x %between% c(min_range, max_range) & !is.na(x)]
  if (length(xout) == 0){
    return(x)
  }
  
  #calculate the range of credible values, based on distribution of within-limits values
  limlow = max(min_range, mean(xin) - sd(xin) * SDs)
  limhigh = min(max_range, mean(xin) + sd(xin) * SDs)
  #calcluate multipliers for order of magnitude of the decimal shift
  mults = 10^c(-rev(1:oom), 0, 1:oom)
  
  #generaten a table of multiplied values
  shifts = as.data.table(lapply(mults, function(x) x * xout))
  #clear those that are outside the limits
  shifts = shifts[, lapply(.SD, function(x) fifelse(x %between% c(limlow, limhigh),abs(x - meanin),NA_real_))]
  shifts[, res := NA_real_]
  shifts[, no := 1:.N]
  allnas = shifts[, which(all(is.na(unlist(.SD)))), no][[1]]
  #shifts[!allnas, res := xout[!allnas] * mults[which.min(.SD)]]
  shifts[!allnas, bestmatch := mults[which.min(.SD)], no]
  
  #shifts[all(is.na(unlist(.SD[[.I]])))]
  
  #coalesce (this will mean 'use the lowest possible valid value')
  #find the closest absolut distance from group mean
  
  #shifts[,  fifelse(all(is.na(unlist(.SD))), NA_real_, which.min(unlist(.SD)[!is.na(unlist(.SD))])), .I]
  
  #shifts = xout * mults[shifts[,  fifelse(all(is.na(unlist(.SD))), NA_real_, which.min(unlist(.SD)[!is.na(unlist(.SD))])), .I][,2, with = FALSE][[1]]]  
  #coalesce (this will mean 'use the lowest possible valid value')
  #shifts = fcoalesce(shifts)
  #print(lenght(shifts))
  #print(length(xout))
  shifts = xout * shifts$bestmatch
  #alter the input vector
  res = x
  res[which(!res %between% c(min_range, max_range)& !is.na(res))] = shifts
  #return the vector
  return(res)
}


#a helper function to do the aggregation (really, value selection) by reference from the cpp function call
interval_aggregate = function(x, y, vcol, options, newcolnames = NULL, fun.aggregate = "min", max_locf_time = 0, max_nocb_time = 0){
  if (is.null(newcolnames)) newcolnames = c(paste0(vcol, "_", fun.aggregate), paste0(paste0(vcol, "_", fun.aggregate), "_",key(y)[2]))
  tempdt = cpp_interval_aggregate(groups = y[[key(y)[1]]],
                                  times = as.numeric(y[[key(y)[2]]], tz = options$time$ptz),
                                  values = y[[vcol]],
                                  int_groups = x[[key(x)[1]]],
                                  int_begin = as.numeric(x[[key(x)[2]]], tz = options$time$ptz),
                                  int_end = as.numeric(x[[key(x)[3]]], tz = options$time$ptz),
                                  aggregation_function = fun.aggregate,
                                  locf = max_locf_time,
                                  nocb = max_nocb_time)
  setDT(tempdt)
  
  #reverse the order as the function returns time first and value last, but the intervals want value first and time after
  newcolnames = rev(newcolnames)
  tempdt[, `:=` (time = as.POSIXct(as.integer(time), tz = options$time$ptz, origin = options$time$po))]
  for (cnum in 1:2){
    x[, (newcolnames[cnum]) := tempdt[[cnum]]]
    
  }
  x
}
