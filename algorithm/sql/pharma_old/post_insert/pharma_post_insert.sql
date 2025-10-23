PRAGMA table_info('pharma_tmp');				  
SELECT * FROM pharma_tmp LIMIT 10;
SELECT COUNT(1) FROM pharma_tmp;

PRAGMA table_info('pharma_names');				  
SELECT * FROM pharma_names LIMIT 10;
SELECT COUNT(1) FROM pharma_names;

PRAGMA table_info('pharma_compounds');				  
SELECT * FROM pharma_compounds LIMIT 10;
SELECT COUNT(1) FROM pharma_compounds;

PRAGMA table_info('pharma');				  
SELECT * FROM pharma LIMIT 10;
select COUNT(1) FROM pharma;
				 
				  
INSERT INTO pharma SELECT
	PatientID, 
	CompID,
	PharmaID,
	OrderNumber,
	DateTime, 
	ROW_NUMBER() 
		OVER(PARTITION BY 
				CompID, 
				PatientID, 
				PharmaID, 
				OrderNumber, 
				DateTime 
			ORDER BY 
				CompID, 
				PatientID, 
				PharmaID, 
				OrderNumber, 
				DateTime) 
		AS RowNumber,
	GivenDose,
	Rate
	FROM (SELECT DISTINCT * FROM pharma_tmp);
				  
PRAGMA table_info('pharma');				  
SELECT * FROM pharma LIMIT 10;
SELECT COUNT(1) FROM pharma;


--human view
--with timestamps as datetime and abbreviation joined in
CREATE VIEW 'hv_pharma' AS
	SELECT 
		PatientID, 
		pharma.CompID,
		pharma_compounds.CompName,
		pharma.PharmaID,
		pharma_names.PharmaName,
		pharma_names.PharmaType,
		OrderNumber,
		datetime(DateTime, 'unixepoch', 'localtime') AS DateTime,
		RowNumber, 
		GivenDose, 
		Rate,
		pharma_names.Unit
		FROM pharma 
			LEFT JOIN pharma_compounds 
				ON pharma.CompID = pharma_compounds.CompID
			LEFT JOIN pharma_names
				ON pharma.PharmaID = pharma_names.PharmaID;
				
SELECT * FROM hv_pharma LIMIT 10;

--DROP TABLE 'pharma_tmp';
--VACUUM;