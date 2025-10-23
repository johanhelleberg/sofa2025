CREATE TABLE daily_orders( 
				  PatientID INTEGER NOT NULL,
				  OrderNumber INTEGER NOT NULL, 
				  PharmaID INTEGER NOT NULL,
				  Concentration REAL,
				  ConcentrationUnit TEXT,
				  PRIMARY KEY (PatientID, OrderNumber, PharmaID),
					FOREIGN KEY (PatientID) REFERENCES patients(PatientID) ON UPDATE CASCADE ON DELETE SET NULL,
					FOREIGN KEY (PharmaID) REFERENCES pharma_names(PharmaID) ON UPDATE CASCADE ON DELETE SET NULL
                  ) WITHOUT ROWID;

-- CREATE TABLE daily_orders( 
				  -- OrderNumber INTEGER NOT NULL, 
				  -- PharmaID INTEGER,
				  -- Concentration REAL,
				  -- ConcentrationUnit TEXT,
				  -- Dose REAL,
				  -- RowNumber INTEGER,
				  -- PRIMARY KEY (OrderNumber, PharmaID, RowNumber),
					-- FOREIGN KEY (OrderNumber) REFERENCES pharma_orders(OrderNumber) ON UPDATE CASCADE ON DELETE SET NULL,
					-- FOREIGN KEY (PharmaID) REFERENCES pharma_names(PharmaID) ON UPDATE CASCADE ON DELETE SET NULL
                  -- ) WITHOUT ROWID;


-- use this table to only keep the first entry of concentration per ordernumber and pharmaID, for quick joins
-- CREATE TABLE daily_orders_first( 
				  -- OrderNumber INTEGER NOT NULL,
				  -- PharmaID INTEGER NOT NULL,				   
				  -- Concentration REAL,
				  -- ConcentrationUnit TEXT,
				  -- PRIMARY KEY (OrderNumber, PharmaID),
					-- FOREIGN KEY (OrderNumber) REFERENCES pharma_orders(OrderNumber) ON UPDATE CASCADE ON DELETE SET NULL,
					-- FOREIGN KEY (PharmaID) REFERENCES pharma_names(PharmaID) ON UPDATE CASCADE ON DELETE SET NULL
                  -- ) WITHOUT ROWID;					  
				  