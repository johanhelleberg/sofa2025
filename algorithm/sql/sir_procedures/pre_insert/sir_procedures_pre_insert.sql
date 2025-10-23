CREATE TABLE sir_procedures(			  
				  VtfHuvudId INTEGER,
				  KvaCode TEXT,
				  StartDateTime REAL,
				  EndDateTime REAL,
				  DurationMinutes REAL,
				  HasDuration INTEGER,
				  FOREIGN KEY (VtfHuvudId) REFERENCES sir_master(VtfHuvudId) ON UPDATE CASCADE ON DELETE SET NULL
                  );