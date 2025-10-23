PRAGMA table_info('labres_temp');				  
SELECT * FROM labres_temp LIMIT 10;
SELECT COUNT(1) FROM labres_temp;

PRAGMA table_info('vartypes');				  
SELECT * FROM vartypes LIMIT 10;
SELECT COUNT(1) FROM vartypes;

PRAGMA table_info('labres');				  
SELECT * FROM labres LIMIT 10;
select COUNT(1) FROM labres;
				 
				  
INSERT INTO labres SELECT 
    PatientID, 
	VariableID, 
	EnterTime, 
	Finding,
	ResultNo,
	SampleTime,
	ROW_NUMBER() 
		OVER(PARTITION BY 
				VariableID, 
				PatientID, 
				SampleTime, 
				EnterTime  
			ORDER BY 
				VariableID, 
				PatientID, 
				SampleTime, 
				EnterTime) 
		AS RowNumber,
	LabStatement, 
	StringValue, 
	VariableType,
	NumericValue,
	ValueID
	FROM (SELECT DISTINCT * FROM labres_temp);
				  
PRAGMA table_info('labres');				  
SELECT * FROM labres LIMIT 10;
SELECT COUNT(1) FROM labres;


--human view
--with timestamps as datetime and abbreviation joined in
CREATE VIEW 'hv_labres' AS
	SELECT 
		PatientID, 
		labres.VariableID, 
		Abbreviation,
		datetime(EnterTime, 'unixepoch', 'localtime') AS EnterTime,
		Finding, 
		ResultNo,
		datetime(SampleTime, 'unixepoch', 'localtime') AS SampleTime,
		LabStatement
		StringValue,
		VariableType,
		NumericValue,  
		ValueID,
		Unit
		 FROM labres 
			LEFT JOIN vartypes 
				ON labres.VariableID = vartypes.VariableID;
				
SELECT * FROM hv_labres LIMIT 10;

DROP TABLE 'labres_temp';
VACUUM;