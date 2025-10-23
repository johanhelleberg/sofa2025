PRAGMA table_info('daily_orders');				  
SELECT * FROM daily_orders LIMIT 10;
--SELECT COUNT(1) FROM daily_orders;

--thanks chatGPT
CREATE VIEW hv_pharma_infusions_full AS
SELECT
    po.PatientID,
    pi.OrderNumber,
    pi.RowNumber,
    datetime(pi.DateTime, 'unixepoch', 'localtime') AS DateTime,
    pi.Rate,
    pp.PharmaID,
    pn.PharmaName,
    pp.CompID,
    pc.CompName,
    do.Concentration,
    do.ConcentrationUnit
FROM
    pharma_infusions AS pi
JOIN
    pharma_orders AS po ON pi.OrderNumber = po.OrderNumber
JOIN
    pharma_prep AS pp ON pi.OrderNumber = pp.OrderNumber
JOIN
    (
        SELECT
            do.OrderNumber,
            do.PharmaID,
            do.RowNumber,
            do.Concentration,
            do.ConcentrationUnit
        FROM
            daily_orders AS do
        JOIN
            (
                SELECT
                    OrderNumber,
                    PharmaID,
                    MIN(RowNumber) AS MinRowNumber
                FROM
                    daily_orders
                GROUP BY
                    OrderNumber,
                    PharmaID
            ) AS min_do ON do.OrderNumber = min_do.OrderNumber AND do.PharmaID = min_do.PharmaID AND do.RowNumber = min_do.MinRowNumber
    ) AS do ON pp.OrderNumber = do.OrderNumber AND pp.PharmaID = do.PharmaID
JOIN
    pharma_names AS pn ON pp.PharmaID = pn.PharmaID
JOIN
    pharma_compounds AS pc ON pp.CompID = pc.CompID;