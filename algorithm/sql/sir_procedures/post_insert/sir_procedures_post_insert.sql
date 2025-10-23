PRAGMA table_info('sir_procedures');				  
SELECT * FROM sir_procedures LIMIT 10;
select Count(1) from sir_procedures;

CREATE VIEW 'hv_procedures' AS 
SELECT VtfHuvudId, 
KvaCode, 
datetime(StartDateTime, 'unixepoch', 'localtime') AS StartDateTime,
datetime(EndDateTime, 'unixepoch', 'localtime') AS EndDateTime,
DurationMinutes,
HasDuration 
FROM 'sir_procedures';