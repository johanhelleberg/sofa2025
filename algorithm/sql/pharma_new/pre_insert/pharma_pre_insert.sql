--compounds
CREATE TABLE pharma_compounds(
CompID INTEGER PRIMARY KEY NOT NULL,
CompName TEXT) WITHOUT ROWID;

--pharmaceuticals
CREATE TABLE pharma_names(
PharmaID INTEGER PRIMARY KEY NOT NULL, 
PharmaName TEXT,
Unit TEXT,
PharmaType INTEGER) WITHOUT ROWID;

-- linking table (many-to-many)
CREATE TABLE comp_pharma_key(
  CompID    INTEGER NOT NULL,
  PharmaID  INTEGER NOT NULL,
  ratio     REAL CHECK (ratio IS NULL OR ratio > 0),
  PRIMARY KEY (CompID, PharmaID),
  FOREIGN KEY (PharmaID) REFERENCES pharma_names(PharmaID)
    ON UPDATE CASCADE ON DELETE CASCADE,
  FOREIGN KEY (CompID)   REFERENCES pharma_compounds(CompID)
    ON UPDATE CASCADE ON DELETE CASCADE
) WITHOUT ROWID;


CREATE TABLE pharma_orders(
PatientID INTEGER NOT NULL,
OrderNumber INTEGER NOT NULL,
PRIMARY KEY (PatientID, OrderNumber),
FOREIGN KEY (PatientID) REFERENCES patients(PatientID) ON UPDATE CASCADE ON DELETE CASCADE
) WITHOUT ROWID;




--CREATE TABLE pharma_prep( 				   
--				  PharmaID INTEGER NOT NULL,
--				  CompID INTEGER NOT NULL,
--				  OrderNumber INTEGER NOT NULL,
--				  Abbreviation TEXT,
--				  EstimatedConcentration REAL,				  
--				  EstimatedConcentrationUnit TEXT,
--				  PRIMARY KEY (PharmaID, CompID, OrderNumber),					
--					FOREIGN KEY (PharmaID) REFERENCES pharma_names(PharmaID) ON UPDATE CASCADE ON DELETE SET NULL,
--					FOREIGN KEY (CompID) REFERENCES pharma_compounds(CompID) ON UPDATE CASCADE ON DELETE SET NULL,
--					FOREIGN KEY (OrderNumber) REFERENCES pharma_orders(OrderNumber) ON UPDATE CASCADE ON DELETE SET NULL
--                 ) WITHOUT ROWID;

CREATE TABLE daily_orders( 
				  PatientID INTEGER NOT NULL,
				  OrderNumber INTEGER NOT NULL, 
				  PharmaID INTEGER NOT NULL,
				  RowNumber INTEGER NOT NULL,
				  Concentration REAL,
				  ConcentrationUnit TEXT,
				  Dose REAL,
				  PlannedTime REAL,
				  PRIMARY KEY (PatientID, OrderNumber, PharmaID, RowNumber),
					FOREIGN KEY (PatientID) REFERENCES patients(PatientID) ON UPDATE CASCADE ON DELETE SET NULL,
					FOREIGN KEY (OrderNumber) REFERENCES pharma_orders(OrderNumber) ON UPDATE CASCADE ON DELETE SET NULL,
					FOREIGN KEY (PharmaID) REFERENCES pharma_names(PharmaID) ON UPDATE CASCADE ON DELETE SET NULL
                  ) WITHOUT ROWID;
				  
CREATE TABLE daily_orders_estimated( 
  PatientID          INTEGER NOT NULL,
  OrderNumber        INTEGER NOT NULL, 
  PharmaID           INTEGER NOT NULL,
  Concentration      REAL,              
  ConcentrationUnit  TEXT,              
  PRIMARY KEY (PatientID, OrderNumber, PharmaID),

  -- Composite FK: the order key is (PatientID, OrderNumber)
  FOREIGN KEY (PatientID, OrderNumber)
    REFERENCES pharma_orders(PatientID, OrderNumber)
    ON UPDATE CASCADE ON DELETE CASCADE,

  FOREIGN KEY (PharmaID)
    REFERENCES pharma_names(PharmaID)
    ON UPDATE CASCADE ON DELETE RESTRICT

) WITHOUT ROWID;
				  
				
CREATE TABLE pharma(
  PatientID   INTEGER NOT NULL,
  OrderNumber INTEGER NOT NULL,
  PharmaID    INTEGER NOT NULL,
  DateTime    REAL NOT NULL,     
  RowNumber   INTEGER NOT NULL,
  GivenDose   REAL,                   -- amount given
  Rate        REAL,                   -- infusion rate (mL/H)

  PRIMARY KEY (PatientID, OrderNumber, PharmaID, DateTime, RowNumber),

  FOREIGN KEY (PatientID)
    REFERENCES patients(PatientID)
    ON UPDATE CASCADE ON DELETE CASCADE,

  -- FK to pharma_names (usually RESTRICT to preserve history if a drug is removed from the catalog)
  FOREIGN KEY (PharmaID)
    REFERENCES pharma_names(PharmaID)
    ON UPDATE CASCADE ON DELETE RESTRICT,

  -- FK to orders must be composite
  FOREIGN KEY (PatientID, OrderNumber)
    REFERENCES pharma_orders(PatientID, OrderNumber)
    ON UPDATE CASCADE ON DELETE CASCADE

) WITHOUT ROWID;				

/*
CREATE TABLE pharma(
	PatientID INTEGER,
	OrderNumber INTEGER, 
	PharmaID INTEGER, 
	DateTime REAL, 
	RowNumber INTEGER,
	GivenDose REAL,
	Rate REAL,
	PRIMARY KEY (PatientID, OrderNumber, PharmaID, DateTime, RowNumber)--,
	--FOREIGN KEY (PatientID) REFERENCES patients(PatientID) ON UPDATE CASCADE ON DELETE SET NULL,
	--FOREIGN KEY (PharmaID) REFERENCES pharma_names(PharmaID) ON UPDATE CASCADE ON DELETE SET NULL,
	--FOREIGN KEY (OrderNumber) REFERENCES pharma_orders(OrderNumber) ON UPDATE CASCADE ON DELETE SET NULL
) WITHOUT ROWID;

*/
-- CREATE TABLE pharma_orders(
-- PatientID INTEGER NOT NULL,
-- OrderNumber INTEGER NOT NULL,
-- PRIMARY KEY (PatientID, OrderNumber),
-- FOREIGN KEY (PatientID) REFERENCES patients(PatientID) ON UPDATE CASCADE ON DELETE SET NULL
-- ) WITHOUT ROWID;

-- CREATE TABLE pharma_prep(
-- OrderNumber INTEGER NOT NULL, 
-- PharmaID INTEGER NOT NULL, 
-- CompID INTEGER NOT NULL, 
-- Concentration REAL NOT NULL, 
-- Unit TEXT NOT NULL, 
-- N_obs INTEGER NOT NULL,
-- PRIMARY KEY (OrderNumber, CompID, PharmaID),
-- FOREIGN KEY (OrderNumber) REFERENCES pharma_orders(OrderNumber) ON UPDATE CASCADE ON DELETE SET NULL,
-- FOREIGN KEY (PharmaID) REFERENCES pharma_names(PharmaID) ON UPDATE CASCADE ON DELETE SET NULL,
-- FOREIGN KEY (CompID) REFERENCES pharma_compounds(CompID) ON UPDATE CASCADE ON DELETE SET NULL
-- ) WITHOUT ROWID;

--CREATE TABLE pharma_tmp(
--	OrderNumber INTEGER, 
--	PharmaID INTEGER, 
--	CompID INTEGER, 
--	DateTime REAL,
--	Rate REAL,	
--	GivenDose REAL,
--	Abbreviation TEXT
--);

-- CREATE TABLE pharma2_infusions_tmp(
	-- OrderNumber INTEGER NOT NULL,
	-- CompID INTEGER,
	-- PharmaID INTEGER,
	-- DateTime REAL, 
	-- Rate REAL,
	-- GivenDose REAL,
	-- FOREIGN KEY (OrderNumber) REFERENCES pharma_orders(OrderNumber) ON UPDATE CASCADE ON DELETE SET NULL,
	-- FOREIGN KEY (PharmaID) REFERENCES pharma_names(PharmaID) ON UPDATE CASCADE ON DELETE SET NULL,
	-- FOREIGN KEY (CompID) REFERENCES pharma_compounds(CompID) ON UPDATE CASCADE ON DELETE SET NULL
-- );

/*
CREATE TABLE pharma(
	OrderNumber INTEGER, 
	CompID INTEGER, 
	PharmaID INTEGER, 
	DateTime REAL, 
	RowNumber INTEGER,
	GivenDose REAL
	PRIMARY KEY (OrderNumber, CompID, PharmaID, DateTime, RowNumber),
	FOREIGN KEY (PatientID) REFERENCES patients(PatientID) ON UPDATE CASCADE ON DELETE SET NULL,
	FOREIGN KEY (PharmaID) REFERENCES pharma_names(PharmaID) ON UPDATE CASCADE ON DELETE SET NULL,
	FOREIGN KEY (CompID) REFERENCES pharma_compounds(CompID) ON UPDATE CASCADE ON DELETE SET NULL,
	FOREIGN KEY (OrderNumber) REFERENCES pharma_orders(OrderNumber) ON UPDATE CASCADE ON DELETE SET NULL
) WITHOUT ROWID;
*/


-- CREATE VIEW hv_pharma AS
-- SELECT 
    -- po.PatientID,
    -- p.OrderNumber,
    -- p.PharmaID,
    -- pn.PharmaName,
    -- p.CompID,
    -- pc.CompName,
    -- p.RowNumber,
    -- datetime(p.DateTime, 'unixepoch', 'localtime') AS DateTime,
    -- p.Rate,
    -- p.GivenDose,
    -- doe.Concentration,
    -- doe.ConcentrationUnit,
    -- pp.EstimatedConcentration,
    -- pp.EstimatedUnit
-- FROM 
    -- pharma p
-- LEFT JOIN
	-- pharma_orders po
	-- ON p.OrderNumber = po.OrderNumber
-- LEFT JOIN 
    -- daily_orders_estimated doe 
    -- ON p.OrderNumber = doe.OrderNumber 
    -- AND p.PharmaID = doe.PharmaID
-- LEFT JOIN 
    -- pharma_compounds pc 
    -- ON p.CompID = pc.CompID
-- LEFT JOIN 
    -- pharma_names pn 
    -- ON p.PharmaID = pn.PharmaID
-- LEFT JOIN 
    -- pharma_prep pp 
    -- ON p.OrderNumber = pp.OrderNumber 
    -- AND p.PharmaID = pp.PharmaID 
    -- AND p.CompID = pp.CompID;