PRAGMA table_info('pharma_tmp');				  
SELECT * FROM pharma_tmp LIMIT 10;
--SELECT COUNT(1) FROM pharma_tmp;

--PRAGMA table_info('pharma_infusions_tmp');
--SELECT * FROM pharma_infusions_tmp LIMIT 10;
--SELECT COUNT(1) FROM pharma_infusions_tmp;
--compid index on the preparations
--CREATE INDEX pp_compid_index ON pharma_prep(CompID);
--CREATE INDEX po_patientid_index ON pharma_orders(PatientID);

PRAGMA table_info('pharma_names');				  
SELECT * FROM pharma_names LIMIT 10;
SELECT COUNT(1) FROM pharma_names;

PRAGMA table_info('pharma_compounds');				  
SELECT * FROM pharma_compounds LIMIT 10;
SELECT COUNT(1) FROM pharma_compounds;

PRAGMA table_info('pharma_prep');				  
SELECT * FROM pharma_prep LIMIT 10;
select COUNT(1) FROM pharma_prep;

PRAGMA table_info('daily_orders');				  
SELECT * FROM pharma_orders LIMIT 10;
select COUNT(1) FROM pharma_orders;

PRAGMA table_info('daily_orders_estimated');				  
SELECT * FROM pharma_orders LIMIT 10;
select COUNT(1) FROM pharma_orders;

CREATE TABLE pharma(
	OrderNumber INTEGER NOT NULL, 
	PharmaID INTEGER NOT NULL,
	CompID INTEGER NOT NULL, 	 
	DateTime REAL NOT NULL, 
	RowNumber INTEGER NOT NULL,
	Rate REAL,
	GivenDose REAL,
	PRIMARY KEY (OrderNumber, PharmaID, CompID, DateTime, RowNumber),
	FOREIGN KEY (PharmaID) REFERENCES pharma_names(PharmaID) ON UPDATE CASCADE ON DELETE SET NULL,
	FOREIGN KEY (CompID) REFERENCES pharma_compounds(CompID) ON UPDATE CASCADE ON DELETE SET NULL,
	FOREIGN KEY (OrderNumber) REFERENCES pharma_orders(OrderNumber) ON UPDATE CASCADE ON DELETE SET NULL
) WITHOUT ROWID;

INSERT INTO pharma SELECT
	OrderNumber,
	PharmaID,
	CompID,
	DateTime,
	ROW_NUMBER() 
		OVER(PARTITION BY 
				OrderNumber, 
				PharmaID, 
				CompID,  
				DateTime 
			ORDER BY 
				OrderNumber, 
				PharmaID, 
				CompID,
				DateTime) 
		AS RowNumber,
	Rate,
	GivenDose
	FROM (SELECT DISTINCT OrderNumber, PharmaID, CompID, DateTime, Rate, GivenDose FROM pharma_tmp);
	
-- CREATE TABLE pharma_infusions(
	-- OrderNumber INTEGER NOT NULL,
	-- DateTime REAL NOT NULL,
	-- RowNumber INTEGER NOT NULL,
	-- Rate REAL,
	-- PRIMARY KEY (OrderNumber, DateTime, RowNumber),
	-- FOREIGN KEY (OrderNumber) REFERENCES pharma_orders(OrderNumber) ON UPDATE CASCADE ON DELETE SET NULL
-- ) WITHOUT ROWID;

-- INSERT INTO pharma_infusions SELECT
	-- OrderNumber, 
	-- DateTime,
	-- ROW_NUMBER() 
		-- OVER(PARTITION BY 
				-- OrderNumber,  
				-- DateTime 
			-- ORDER BY 
				-- OrderNumber, 
				-- DateTime) 
		-- AS RowNumber,
	-- Rate
	-- FROM (SELECT DISTINCT OrderNumber, DateTime, Rate FROM pharma_infusions_tmp);


PRAGMA table_info('pharma');				  
SELECT * FROM pharma LIMIT 10;
--SELECT COUNT(1) FROM pharma_tmp;

-- PRAGMA table_info('pharma_infusions');
-- SELECT * FROM pharma_infusions LIMIT 10;

--human view
--with timestamps as datetime and abbreviation joined in
-- CREATE VIEW 'hv_pharma2' AS
	-- SELECT 
		-- PatientID, 
		-- pharma.CompID,
		-- CompName,
		-- pharma.PharmaID,
		-- PharmaName,
		-- PharmaType,
		-- pharma.OrderNumber,
		-- datetime(DateTime, 'unixepoch', 'localtime') AS DateTime,
		-- Rate,
		-- RowNumber, 
		-- GivenDose,
		-- Unit
		-- FROM pharma 
			-- LEFT JOIN pharma_orders
				-- ON pharma.OrderNumber = pharma_orders.OrderNumber
			-- LEFT JOIN pharma_compounds 
				-- ON pharma.CompID = pharma_compounds.CompID
			-- LEFT JOIN pharma_names
				-- ON pharma.PharmaID = pharma_names.PharmaID;
				
-- SELECT * FROM hv_pharma_noninfusion LIMIT 10;

-- CREATE VIEW 'hv_pharma_infusion' AS 
	-- SELECT
		-- PatientID,
		-- pharma_prep.CompID,
		-- CompName,
		-- pharma_prep.PharmaID,
		-- PharmaName,
		-- PharmaType,
		-- pharma_infusions.OrderNumber,
		-- datetime(DateTime, 'unixepoch', 'localtime') AS DateTime,
		-- RowNumber,
		-- Rate,
		-- Concentration,
		-- pharma_prep.Unit
		-- FROM (pharma_infusions
			-- JOIN pharma_prep ON	
				-- pharma_infusions.OrderNumber = pharma_prep.OrderNumber)
			-- LEFT JOIN pharma_orders
				-- ON pharma_infusions.OrderNumber = pharma_orders.OrderNumber		
			-- LEFT JOIN pharma_compounds ON
				-- pharma_prep.CompID = pharma_compounds.CompID
			-- LEFT JOIN pharma_names ON
				-- pharma_prep.PharmaID = pharma_names.PharmaID;
				
-- SELECT * FROM hv_pharma_infusion LIMIT 10;

-- CREATE VIEW 'hv_pharma_infusion_new' AS 
	-- SELECT
		-- po.PatientID,
		-- pi.OrderNumber,
		-- pi.RowNumber,
		-- pi.DateTime AS DateTimeNum,
		-- datetime(pi.DateTime, 'unixepoch', 'localtime') AS DateTime,
		--pi.DateTime,
		-- pi.Rate,
		-- pp.PharmaID,
		-- pn.PharmaName,
		-- pn.PharmaType,
		-- pp.CompID,
		-- pc.CompName,
		-- do.Concentration,
		-- do.ConcentrationUnit AS Unit
		-- FROM
			-- pharma_infusions AS pi
		-- JOIN
			-- pharma_orders AS po ON pi.OrderNumber = po.OrderNumber
		-- JOIN
			-- pharma_prep AS pp ON pi.OrderNumber = pp.OrderNumber
		-- JOIN 
			-- pharma_compounds AS pc ON pp.CompID = pc.CompID
		-- JOIN 
			-- pharma_names as pn ON pp.PharmaID = pn.PharmaID
		-- JOIN
			-- (
				-- SELECT
					-- do.OrderNumber,
					-- do.PharmaID,
					-- do.RowNumber,
					-- do.Concentration,
					-- do.ConcentrationUnit
				-- FROM
					-- daily_orders AS do
				-- JOIN
					-- (
						-- SELECT
							-- OrderNumber,
							-- PharmaID,
							-- MIN(RowNumber) AS MinRowNumber
						-- FROM
							-- daily_orders
						-- GROUP BY
							-- OrderNumber,
							-- PharmaID
					-- ) AS min_do ON do.OrderNumber = min_do.OrderNumber AND do.PharmaID = min_do.PharmaID AND do.RowNumber = min_do.MinRowNumber
			-- ) AS do ON pp.OrderNumber = do.OrderNumber AND pp.PharmaID = do.PharmaID ORDER BY
			  -- po.PatientID, pi.DateTime, pi.OrderNumber, pp.CompID, pp.PharmaID, pi.RowNumber;
			  

--a fast view of pharma_infusions, if queried on compid			  
CREATE VIEW 'hv_pharma2' AS 
	SELECT
		po.PatientID,
		pi.OrderNumber,
		pi.RowNumber,
		--pi.DateTime AS DateTimeNum,
		datetime(pi.DateTime, 'unixepoch', 'localtime') AS DateTime,
		--pi.DateTime,
		pi.Rate,
		pi.GivenDose,
		pi.PharmaID,
		pn.PharmaName,
		pn.PharmaType,
		pi.CompID,
		pc.CompName,
		do.Concentration,
		do.ConcentrationUnit AS Unit,
		pp.Unit AS PharmaUnit 
		FROM
			pharma2 AS pi
		JOIN
			pharma_orders AS po ON pi.OrderNumber = po.OrderNumber 
		JOIN
			pharma_prep AS pp ON pi.OrderNumber = pp.OrderNumber AND pi.CompID = pp.CompID AND pi.PharmaID = pp.PharmaID
		JOIN 
			pharma_compounds AS pc ON pp.CompID = pc.CompID
		JOIN 
			pharma_names as pn ON pp.PharmaID = pn.PharmaID			
		JOIN
			daily_orders_first AS do ON pp.OrderNumber = do.OrderNumber AND pp.PharmaID = do.PharmaID ORDER BY
			  po.PatientID, pi.DateTime, pi.OrderNumber, pp.CompID, pp.PharmaID, pi.RowNumber;	

CREATE VIEW 'hv_pharma3' AS 
	SELECT
		po.PatientID,
		pi.OrderNumber,
		pi.RowNumber,
		--pi.DateTime AS DateTimeNum,
		datetime(pi.DateTime, 'unixepoch', 'localtime') AS DateTime,
		--pi.DateTime,
		pi.Rate,
		pi.GivenDose,
		pi.PharmaID,
		--pn.PharmaName,
		--pn.PharmaType,
		pi.CompID,
		--pc.CompName,
		--do.Concentration,
		--do.ConcentrationUnit AS Unit,
		pp.Unit AS PharmaUnit 
		FROM
			pharma2 AS pi
		JOIN
			pharma_orders AS po ON pi.OrderNumber = po.OrderNumber 
		JOIN
			pharma_prep AS pp ON pi.OrderNumber = pp.OrderNumber AND pi.CompID = pp.CompID AND pi.PharmaID = pp.PharmaID
		-- JOIN 
			-- pharma_compounds AS pc ON pp.CompID = pc.CompID
		-- JOIN 
			-- pharma_names as pn ON pp.PharmaID = pn.PharmaID			
		-- JOIN
			-- daily_orders_first AS do ON pp.OrderNumber = do.OrderNumber AND pp.PharmaID = do.PharmaID ORDER BY
			  -- po.PatientID, pi.DateTime, pi.OrderNumber, pp.CompID, pp.PharmaID, pi.RowNumber;				  
			  

-- CREATE VIEW hv_pharma_infusions_full AS
    -- SELECT
        -- po.PatientID,
        -- pi.OrderNumber,
        -- pi.RowNumber,
        -- datetime(pi.DateTime, 'unixepoch', 'localtime') AS DateTime,
        -- pi.Rate,
        -- pp.PharmaID,
        -- pn.PharmaName,
        -- pp.CompID,
        -- pc.CompName,
        -- do.Concentration,
        -- do.ConcentrationUnit
    -- FROM
        -- pharma_infusions AS pi
    -- JOIN
        -- pharma_orders AS po ON pi.OrderNumber = po.OrderNumber
    -- JOIN
        -- pharma_prep AS pp ON pi.OrderNumber = pp.OrderNumber
		
		
    -- JOIN
        -- (
            -- SELECT
                -- do.OrderNumber,
                -- do.PharmaID,
                -- do.RowNumber,
                -- do.Concentration,
                -- do.ConcentrationUnit
            -- FROM
                -- daily_orders AS do
            -- JOIN
                -- (
                    -- SELECT
                        -- OrderNumber,
                        -- PharmaID,
                        -- MIN(RowNumber) AS MinRowNumber
                    -- FROM
                        -- daily_orders
                    -- GROUP BY
                        -- OrderNumber,
                        -- PharmaID
                -- ) AS min_do ON do.OrderNumber = min_do.OrderNumber AND do.PharmaID = min_do.PharmaID AND do.RowNumber = min_do.MinRowNumber
        -- ) AS do ON pp.OrderNumber = do.OrderNumber AND pp.PharmaID = do.PharmaID
    -- JOIN
        -- pharma_names AS pn ON pp.PharmaID = pn.PharmaID
    -- JOIN
        -- pharma_compounds AS pc ON pp.CompID = pc.CompID;
			  				
--DROP TABLE 'pharma2_tmp';
--DROP TABLE 'pharma_infusions_tmp';
--VACUUM;	
--PRAGMA optimize;

		 
/*				  
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
*/

--human view
--with timestamps as datetime and abbreviation joined in
/*
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
*/

--DROP TABLE 'pharma_tmp';
--VACUUM;