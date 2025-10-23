

PRAGMA table_info('obsrec');				  
SELECT * FROM obsrec LIMIT 10;
SELECT COUNT(1) FROM obsrec;

PRAGMA table_info('factor_names');				  
SELECT * FROM factor_names LIMIT 10;
SELECT COUNT(1) FROM factor_names;		 
			
-- human view
-- timestamps as datetime(), join with vartypes, and include the stringvalues from the factor_names
CREATE VIEW 'hv_obsrec' AS
	SELECT 
		PatientID, 
		obsrec.VariableID, 
		datetime(DateTime, 'unixepoch', 'localtime') AS DateTime,
		datetime(EnterTime, 'unixepoch', 'localtime') AS EnterTime,
		RowNumber,
		OrderNumber,
		Status,
		Abbreviation, 
		obsrec.NumericValue, 
		COALESCE(obsrec.StringValue, factor_names.StringValue) as StringValue, 
		Unit
		FROM obsrec 
			LEFT JOIN vartypes 
				ON obsrec.VariableID = vartypes.VariableID
			LEFT JOIN factor_names
				ON obsrec.VariableID = factor_names.VariableID AND
				   obsrec.NumericValue = factor_names.NumericValue;
				
SELECT * FROM hv_obsrec LIMIT 10;
