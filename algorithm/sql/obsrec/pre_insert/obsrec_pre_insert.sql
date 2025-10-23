CREATE TABLE obsrec(
	PatientID INTEGER NOT NULL, 
	VariableID INTEGER NOT NULL, 
	DateTime REAL NOT NULL,
	EnterTime REAL NOT NULL, 
	RowNumber INTEGER,
	OrderNumber INTEGER,
	Status INTEGER,
	NumericValue REAL,
	StringValue REAL,
	PRIMARY KEY (VariableID, PatientID, DateTime, EnterTime, RowNumber),
	FOREIGN KEY (PatientID) REFERENCES patients(PatientID) ON UPDATE CASCADE ON DELETE SET NULL,
	FOREIGN KEY (VariableID) REFERENCES vartypes(VariableID) ON UPDATE CASCADE ON DELETE SET NULL
	) WITHOUT ROWID;	
				  
CREATE TABLE factor_names(
	VariableID INTEGER,
	NumericValue REAL,
	StringValue TEXT,
	StringValueEng TEXT,
	PRIMARY KEY (VariableID, NumericValue),
	FOREIGN KEY (VariableID) REFERENCES vartypes(VariableID) ON UPDATE CASCADE ON DELETE SET NULL
	) WITHOUT ROWID;
				  
			  
