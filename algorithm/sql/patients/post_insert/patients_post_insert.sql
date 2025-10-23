PRAGMA table_info('patients');				  
SELECT * FROM patients LIMIT 10;
select Count(1) from patients;

CREATE VIEW 'hv_patients' AS 
SELECT 
PatientID, 
StudyID, 
MinAge, 
SexMale, 
datetime(FirstAdmissionDateTime, 'unixepoch', 'localtime') AS FirstAdmissionDateTime,
datetime(LastDischargeDateTime, 'unixepoch', 'localtime') AS LastDischargeDateTime,
DeadInICU
FROM 'patients';