CREATE TABLE patients(
                  PatientID INTEGER,
                  StudyID INTEGER,
                  MinAge INTEGER,
                  SexMale INTEGER,
                  FirstAdmissionDateTime REAL,
                  LastDischargeDateTime REAL,
                  DeadInICU INTEGER,
				  PRIMARY KEY (StudyID, PatientID, FirstAdmissionDateTime)
                  FOREIGN KEY (StudyID) REFERENCES subjects(StudyID) ON UPDATE CASCADE ON DELETE SET NULL
                  ) WITHOUT ROWID;