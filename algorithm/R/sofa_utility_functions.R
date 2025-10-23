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
