CREATE TABLE monvals(
	PatientID INTEGER NOT NULL, 
	VariableID INTEGER NOT NULL, 
	DateTime REAL NOT NULL,
	EnterTime REAL NOT NULL, 
	RowNumber INTEGER NOT NULL,
	Status REAL,
	NumericValue REAL,
	PRIMARY KEY (VariableID, PatientID, DateTime, EnterTime, RowNumber),
	FOREIGN KEY (PatientID) REFERENCES patients(PatientID) ON UPDATE CASCADE ON DELETE SET NULL,
	FOREIGN KEY (VariableID) REFERENCES vartypes(VariableID) ON UPDATE CASCADE ON DELETE SET NULL
	) WITHOUT ROWID;	
				  
CREATE TABLE monvals_temp(
	PatientID INTEGER, 
	VariableID INTEGER, 
	DateTime REAL,
	EnterTime REAL,
	Status REAL,
	NumericValue REAL,
	FOREIGN KEY (PatientID) REFERENCES patients(PatientID) ON UPDATE CASCADE ON DELETE SET NULL,
	FOREIGN KEY (VariableID) REFERENCES vartypes(VariableID) ON UPDATE CASCADE ON DELETE SET NULL
	);	
			  
