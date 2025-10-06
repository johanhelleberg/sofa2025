#2025-10-06
#Time to non-inferiority analysis for SOFA
library(data.table)
file_path = "P:\\GLUC-ICU\\users\\johanh\\sofa2025\\datasets\\sofa_test_dataset_hourly_with_outcome_251006.csv"
sofadt = fread(file_path)
rm(file_path)
