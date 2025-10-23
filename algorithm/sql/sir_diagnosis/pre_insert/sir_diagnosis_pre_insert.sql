CREATE TABLE sir_diagnosis(			  
				  VtfHuvudId INTEGER,
				  PrimaryICD INTEGER,
				  ICD10 TEXT,
				  FOREIGN KEY (VtfHuvudId) REFERENCES sir_master(VtfHuvudId) ON UPDATE CASCADE ON DELETE SET NULL
                  );