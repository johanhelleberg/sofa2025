PRAGMA table_info('monvals_temp');				  
SELECT * FROM monvals_temp LIMIT 10;
--SELECT COUNT(1) FROM monvals_temp;

PRAGMA table_info('factor_names');				  
SELECT * FROM factor_names LIMIT 10;
SELECT COUNT(1) FROM factor_names;	

PRAGMA table_info('vartypes');				  
SELECT * FROM vartypes LIMIT 10;
SELECT COUNT(1) FROM vartypes;


INSERT INTO monvals SELECT
	PatientID,
	VariableID,
	DateTime,
	EnterTime,
	ROW_NUMBER() 
		OVER(PARTITION BY 
				VariableID, 
				PatientID, 
				DateTime,  
				EnterTime 
			ORDER BY 
				VariableID, 
				PatientID, 
				DateTime,
				EnterTime) 
		AS RowNumber,
	Status,
	NumericValue
	FROM (SELECT DISTINCT PatientID, VariableID, DateTime, EnterTime, Status, NumericValue FROM monvals_temp); 



			
-- human view
-- timestamps as datetime(), join with vartypes, and include the stringvalues from the factor_names
CREATE VIEW 'hv_monvals' AS
	SELECT 
		PatientID, 
		monvals.VariableID, 
		datetime(DateTime, 'unixepoch', 'localtime') AS DateTime,
		datetime(EnterTime, 'unixepoch', 'localtime') AS EnterTime,
		RowNumber,
		Abbreviation, 
		NumericValue, 
		Status
		FROM monvals 
			LEFT JOIN vartypes 
				ON monvals.VariableID = vartypes.VariableID;
				
SELECT * FROM hv_monvals LIMIT 10;

DROP TABLE 'monvals_temp';
VACUUM;	
PRAGMA optimize;
