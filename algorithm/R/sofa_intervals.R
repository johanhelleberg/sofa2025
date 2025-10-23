

#generating intervals
library(data.table)

# options = list()
# options$time = list()
# options$time$ptz = "CET"
# options$time$po = "1970-01-01"
# options$time$format = "%Y-%m-%d %H:%M:%OS"
# 
# 
# data = data.table(subject = c(1,2,3,4,5),
#                   admissiontime = c("2014-03-05 16:14:03.4515",
#                                     "2018-08-24 02:30:46.1200",
#                                     "2015-12-13 13:00:21.8231",
#                                     "2016-07-25 22:06:05.0054",
#                                     "2021-01-10 08:33:51.6124"),
#                   dischargetime = c("2014-03-07 11:10:16.3256",
#                                     "2018-08-24 03:21:52.6528",
#                                     "2016-01-02 09:16:01.0249",
#                                     "2016-07-28 14:26:43.3568",
#                                     "2021-01-12 15:55:02.4563"))
# 
# data[, `:=` (admissiontime = as.POSIXct(admissiontime, tz = options$time$ptz, format = options$time$format),
#              dischargetime = as.POSIXct(dischargetime, tz = options$time$ptz, format = options$time$format))]
# 
# setkey(data, subject)

make_intervals = function(data,
                          options,
                          type,
                          interval_length = NULL,
                          interval_breakpoints = NULL,
                          time_columns = NULL,
                          slide_duration = NULL,
                          min_interval_duration = NULL,
                          start_interval_duration = 0){
  # if(FALSE){
  #   data = copy(demo)[PatientID == 18733, .(PatientID, AdmissionTime, DischargeTime)]
  #   type = "sliding"
  #   interval_length = 3600*24
  #   start_interval_duration = 0
  #   interval_breakpoints = NULL
  #   slide_duration = 3600*2
  #   min_interval_duration = 5*3600
  #   #interval_breakpoints = c("08:00:00", "12:00:00", "20:00:00")
  #   #interval_breakpoints = c("09:00:00")
  #   time_columns = NULL
  # }
  

  
  {
    non_overlap_duration = options$intervals$non_overlap_duration
    index_column_name = options$intervals$index_column_name
    interval_start_name = options$intervals$interval_start_name
    interval_end_name = options$intervals$interval_end_name
  }
  
  # rm(list = c("data",
  #             "type",
  #             "interval_length",
  #             "interval_breakpoints",
  #             "min_interval_duration",
  #             "time_columns",
  #             "slide_duration",
  #             "start_interval_duration",
  #             "non_overlap_duration",
  #             "index_column_name",
  #             "interval_start_name",
  #             "interval_end_name",
  #             "y"))
  
  #this will be our output table
  y = copy(data)
  
  
  #use first non-key columns if nothing else is specified
  if (is.null(time_columns)){
    time_columns = colnames(y)[!colnames(y) %in% key(y)][1:2]
    if (any(is.na(time_columns))){
      stop("No time columns specified and there are less than two non-key columns in the data frame")
    }
  }
  #sanity check that both columns are numeric
  if (!all(y[, lapply(.SD, typeof), .SDcols = time_columns] %in% c("numeric", "double", "integer"))){
    stop("Both time columns must be in some sort of numeric format")
  }
  
  #check that all time durations are >=0
  if (any(y[[time_columns[2]]] <= y[[time_columns[[1]]]])){
    stop(paste0("Some of '", time_columns[2], "' are less or equal than '", time_columns[2], "'"))
  }
  
  #only keep the relevant columns
  y = y[, c(key(data), time_columns), with = FALSE]
  
  #if a starting interval is requested, create those and also remove that time from the interval start column column
  startintervals = NULL
  if (start_interval_duration >0){
    startintervals = copy(y[, c(key(y), time_columns), with = FALSE])
    startintervals[, (index_column_name) := 0]
    startintervals[, (interval_start_name) := .SD[[time_columns[1]]] - start_interval_duration / 2]
    startintervals[, (interval_end_name) := pmin(.SD[[time_columns[1]]] + start_interval_duration / 2, .SD[[time_columns[2]]])]
    y[, (time_columns[1]) := startintervals[[interval_end_name]]]
  }
  
  #generate interval of fixed length
  if (type == "relative"){
    #calculate the duration
    y[, total_duration := as.numeric(.SD[[time_columns[2]]] - .SD[[time_columns[1]]], units = "secs")]
    y[, nintervals := ceiling(total_duration / interval_length)]
    #create index column
    y = y[, .(IntervalN  = 1:nintervals), c(key(y), time_columns)]
    #remove intervals of length 0
    y = y[y[[time_columns[1]]] < y[[time_columns[[2]]]]]
    y[, (interval_start_name) := .SD[[time_columns[1]]] + (IntervalN -1 ) * interval_length]
    y[, (interval_end_name) := pmin(.SD[[time_columns[1]]] + (IntervalN) * interval_length,.SD[[time_columns[2]]])]
    colnames(y)[colnames(y) == "IntervalN"] = index_column_name
  }
  
  if (type == "sliding"){
    #calculate the duration
    y[, total_duration := as.numeric(.SD[[time_columns[2]]] - .SD[[time_columns[1]]], units = "secs")]
    y[, nintervals := ceiling((total_duration)/slide_duration)]
    #create index column
    y = y[, .(IntervalN  = 1:nintervals), c(key(y), time_columns)]
    #remove intervals of length 0
    y = y[y[[time_columns[1]]] < y[[time_columns[[2]]]]]
    
    
    y[, (interval_start_name) := pmax(pmin(.SD[[time_columns[1]]] + (IntervalN) * slide_duration, .SD[[time_columns[2]]]) - interval_length,.SD[[time_columns[1]]])]
    y[, (interval_end_name) := pmin(.SD[[time_columns[1]]] + (IntervalN) * slide_duration,.SD[[time_columns[2]]])]
    if(!is.null(min_interval_duration)){
      y = y[which(y[[interval_end_name]] - y[[interval_start_name]] > as.difftime(min_interval_duration, units = "secs"))]
      y[, IntervalN := 1+IntervalN - min(IntervalN), c(key(data))]
    }
    
    #y[, (interval_start_name) := pmax(.SD[[time_columns[1]]] + (IntervalN) * slide_duration - interval_length],.SD[[time_columns[1]]])]   
    #y[, (interval_end_name) := pmin(.SD[[time_columns[1]]] + (IntervalN -1 ) * slide_duration + interval_length,.SD[[time_columns[2]]])]
    colnames(y)[colnames(y) == "IntervalN"] = index_column_name
  }
  
  if (type == "absolute"){
    #remove intervals of length 0
    y = y[y[[time_columns[1]]] < y[[time_columns[[2]]]]]
    #set up the data-table of potential intervals
    w = copy(y)
    #calculate number of days (+/- one day to make sure we overlap at all time points)
    w[, number_of_days := as.numeric(as.Date(.SD[[time_columns[2]]], origin = options$time$po))+2 - as.numeric(as.Date(.SD[[time_columns[1]]], origin = options$time$po))]
    #calculate the start of the intervals
    w = w[, .(day_number = (0:number_of_days)-1), c(key(w), time_columns)]
    w[, day_char := as.character(as.Date(trunc(.SD[[time_columns[1]]], "days"),origin = options$time$po) + day_number)]
    w[, day_number := NULL]
    w = w[, .(timepoint1 = paste0(day_char, " ", interval_breakpoints)), c(key(y), time_columns, "day_char")]
    
    
    #end of intervals = start of interval shifted by 1
    w[, (interval_start_name) := as.POSIXct(timepoint1, tz = options$time$ptz, origin = options$time$po)]
    w[, (interval_end_name) := shift(.SD[[interval_start_name]], -1), key(y)]
    
    #remove rows with NA
    w = na.omit(w[, c(key(y), interval_start_name, interval_end_name), with = FALSE])
    #set key and overlap join with y
    setkeyv(w, c(key(y), interval_start_name, interval_end_name))
    setkeyv(y, c(key(y), time_columns))
    y = na.omit(foverlaps(w,y))
    
    #make sure the intervals are capped at admissiontime/dischargetime
    y[, (interval_start_name) := pmax(.SD[[interval_start_name]], .SD[[time_columns[1]]])]
    y[, (interval_end_name) := pmin(.SD[[interval_end_name]], .SD[[time_columns[2]]])]
    
    #again, remove those that are of length 0
    y = y[y[[interval_start_name]] < y[[interval_end_name]]]
    
    #count the number of intervals
    y[, (index_column_name) := 1:.N, key(data)]
    #set up the correct column order
    setcolorder(y, c(key(data), time_columns, index_column_name, interval_start_name, interval_end_name))
  }
  
  #bind the start interval back if requested
  if (!is.null(startintervals)){
    y = rbind(startintervals, y)
    rm(startintervals)
  }
  
  #shorten the interval ends to make sure they don't overlap, if requested
  if(non_overlap_duration >0){
    y[, (interval_end_name) := fifelse(.SD[[index_column_name]] == max(.SD[[index_column_name]]), 
                                       .SD[[interval_end_name]], 
                                       .SD[[interval_end_name]]-non_overlap_duration),
      key(data)]
  }
  
  #sort y
  setkeyv(y, c(key(data), interval_start_name, interval_end_name))
  
  #restore the starting point if it had been modified with a start interval
  y[, (time_columns[1]) := min(.SD[[time_columns[1]]]), key(data)]
  
  return(y)
}

