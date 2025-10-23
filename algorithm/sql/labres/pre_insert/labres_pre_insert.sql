CREATE TABLE labres_temp(
                  PatientID INTEGER, 
				  VariableID INTEGER, 
				  EnterTime REAL, 
				  Finding TEXT,
				  ResultNo INT,				  
				  SampleTime REAL, 
				  LabStatement TEXT, 
				  StringValue TEXT,
				  VariableType CHAR,				  
				  NumericValue REAL,
				  ValueID CHAR
                  );
				  
CREATE TABLE vartypes(
						VariableID INTEGER PRIMARY KEY,
						Abbreviation TEXT,
						AbbreviationEng TEXT,
						Unit TEXT
						);
				  

CREATE TABLE labres(
                  PatientID INTEGER NOT NULL, 
				  VariableID INTEGER NOT NULL, 
				  EnterTime REAL NOT NULL, 
				  Finding TEXT,
				  ResultNo INT,
				  SampleTime REAL NOT NULL,
				  RowNumber INTEGER NOT NULL,
				  LabStatement TEXT, 
				  StringValue TEXT, 
				  VariableType CHAR,
				  NumericValue REAL,
				  ValueID CHAR,
				  PRIMARY KEY (VariableID, PatientID, SampleTime, EnterTime, RowNumber),
				  FOREIGN KEY (PatientID) REFERENCES patients(PatientID) ON UPDATE CASCADE ON DELETE SET NULL,
				  FOREIGN KEY (VariableID) REFERENCES vartypes(VariableID) ON UPDATE CASCADE ON DELETE SET NULL
                  ) WITHOUT ROWID;	