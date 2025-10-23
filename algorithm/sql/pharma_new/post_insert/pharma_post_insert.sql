--PRAGMA table_info('pharma_tmp');				  
--SELECT * FROM pharma_tmp LIMIT 10;
--SELECT COUNT(1) FROM pharma_tmp;

--PRAGMA table_info('pharma_infusions_tmp');
--SELECT * FROM pharma_infusions_tmp LIMIT 10;
--SELECT COUNT(1) FROM pharma_infusions_tmp;
--compid index on the preparations
--CREATE INDEX pp_compid_index ON pharma_prep(CompID);
--CREATE INDEX po_patientid_index ON pharma_orders(PatientID);

-- Index for fast queries on PharmaID
CREATE INDEX idx_comp_pharma_key_pharma ON comp_pharma_key(PharmaID, CompID);

-- Index for fast concentration joins on PharmaID
CREATE INDEX idx_doe_by_pharma   ON daily_orders_estimated(PharmaID, PatientID, OrderNumber);

--Index pharma on pharmaID
CREATE INDEX idx_pharma_by_pharma ON pharma(PharmaID);

--index for speeding up joins
CREATE INDEX idx_pharma_join ON pharma(PatientID, OrderNumber, PharmaID, DateTime);

CREATE INDEX IF NOT EXISTS idx_pharma_join
  ON pharma(PatientID, OrderNumber, PharmaID, DateTime);

PRAGMA table_info('pharma_names');				  
SELECT * FROM pharma_names LIMIT 10;
SELECT COUNT(1) FROM pharma_names;

PRAGMA table_info('pharma_compounds');				  
SELECT * FROM pharma_compounds LIMIT 10;
SELECT COUNT(1) FROM pharma_compounds;

--PRAGMA table_info('pharma_prep');				  
--SELECT * FROM pharma_prep LIMIT 10;
--select COUNT(1) FROM pharma_prep;

PRAGMA table_info('daily_orders');				  
SELECT * FROM pharma_orders LIMIT 10;
select COUNT(1) FROM pharma_orders;

PRAGMA table_info('daily_orders_estimated');				  
SELECT * FROM pharma_orders LIMIT 10;
select COUNT(1) FROM pharma_orders;

CREATE TABLE pharma(
	PatientID INTEGER NOT NULL,
	OrderNumber INTEGER NOT NULL, 
	PharmaID INTEGER NOT NULL,
	--CompID INTEGER NOT NULL, 	 
	DateTime REAL NOT NULL, 
	RowNumber INTEGER NOT NULL,
	Rate REAL,
	GivenDose REAL,
	PRIMARY KEY (PatientID, OrderNumber, PharmaID, DateTime, RowNumber),
	FOREIGN KEY (PharmaID) REFERENCES pharma_names(PharmaID) ON UPDATE CASCADE ON DELETE SET NULL,
	--FOREIGN KEY (CompID) REFERENCES pharma_compounds(CompID) ON UPDATE CASCADE ON DELETE SET NULL,
	FOREIGN KEY (OrderNumber) REFERENCES pharma_orders(OrderNumber) ON UPDATE CASCADE ON DELETE SET NULL
) WITHOUT ROWID;

INSERT INTO pharma SELECT
	OrderNumber,
	PharmaID,
	--CompID,
	DateTime2 as DateTime,
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
	FROM (SELECT DISTINCT OrderNumber, PharmaID,  DateTime, Rate, GivenDose FROM pharma_tmp);
	
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

--human legible view of pharma
CREATE VIEW hv_pharma AS
SELECT
  p.PatientID,
  p.OrderNumber,
  p.PharmaID,
  n.PharmaName,
  n.Unit            AS PharmaUnit,
  n.PharmaType,
  d.Concentration,
  d.ConcentrationUnit,
  datetime(p.DateTime, 'unixepoch', 'localtime') AS DateTime,
  p.RowNumber,
  p.GivenDose,
  p.Rate
FROM pharma AS p
JOIN pharma_names AS n
  ON n.PharmaID = p.PharmaID
LEFT JOIN daily_orders_estimated AS d
  ON d.PatientID = p.PatientID
 AND d.OrderNumber = p.OrderNumber
 AND d.PharmaID = p.PharmaID;
 --human legible view of pharma
CREATE VIEW hv_pharma2 AS
SELECT
  p.PatientID,
  p.OrderNumber,
  p.PharmaID,
  n.PharmaName,
  n.Unit AS PharmaUnit,
  n.PharmaType,
  d.Concentration,
  d.ConcentrationUnit,
  datetime(p.DateTime, 'unixepoch', 'localtime') AS DateTime,
  p.RowNumber,
  p.GivenDose,
  p.Rate
FROM daily_orders_estimated AS d
JOIN pharma AS p
  ON p.PatientID   = d.PatientID
 AND p.OrderNumber = d.OrderNumber
 AND p.PharmaID    = d.PharmaID
JOIN pharma_names AS n
  ON n.PharmaID    = d.PharmaID;
  
  
CREATE VIEW hv_pharmacomp AS
SELECT
  d.PatientID,
  d.OrderNumber,
  cpk.CompID,
  c.CompName,
  d.PharmaID,
  n.PharmaName,
  n.Unit                           AS PharmaUnit,
  n.PharmaType,
  cpk.ratio,
  d.Concentration,
  d.ConcentrationUnit,
  datetime(p.DateTime, 'unixepoch', 'localtime') AS DateTime,
  p.RowNumber,
  p.GivenDose,
  /* Doses scaled by compound ratio */
  CASE
    WHEN p.GivenDose IS NOT NULL AND cpk.ratio IS NOT NULL
      THEN p.GivenDose * cpk.ratio
    ELSE NULL
  END AS CompoundDose,
  p.Rate,
  /* DoseRate calculation, only makes sense if ConcentrationUnit = something/mL as rate is always mL/hour*/
  CASE
    WHEN p.Rate IS NOT NULL AND d.Concentration IS NOT NULL AND cpk.ratio IS NOT NULL
      THEN p.Rate * d.Concentration * cpk.ratio
    ELSE NULL
  END AS CompoundRatePerHour
FROM comp_pharma_key AS cpk
JOIN pharma_compounds AS c
  ON c.CompID   = cpk.CompID
JOIN pharma_names AS n
  ON n.PharmaID = cpk.PharmaID
JOIN daily_orders_estimated AS d
  ON d.PharmaID = cpk.PharmaID
JOIN pharma AS p
  ON p.PatientID   = d.PatientID
 AND p.OrderNumber = d.OrderNumber
 And p.PharmaID    = d.PharmaID;



CREATE VIEW hv_pharmacomp_fast AS
SELECT
  d.PatientID,
  d.OrderNumber,
  cpk.CompID,
  d.PharmaID,
  p.DateTime,
  p.RowNumber,
  p.GivenDose,
  p.Rate,
  n.Unit AS PharmaUnit,
  d.Concentration,
  d.ConcentrationUnit,
  cpk.ratio
FROM comp_pharma_key AS cpk
JOIN daily_orders_estimated AS d
  ON d.PharmaID = cpk.PharmaID
JOIN pharma AS p
  ON p.PatientID   = d.PatientID
 AND p.OrderNumber = d.OrderNumber
 AND p.PharmaID    = d.PharmaID
JOIN pharma_names AS n
  ON n.PharmaID    = d.PharmaID;  
--SELECT COUNT(1) FROM pharma_tmp;

 --CREATE VIEW 'hv_pharma2' AS
	 --SELECT 
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
/*	  
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
			  */

/*
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
			  
			  */

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