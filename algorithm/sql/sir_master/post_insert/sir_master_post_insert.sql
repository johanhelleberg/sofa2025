
PRAGMA table_info('sir_master');				  
SELECT * FROM sir_master LIMIT 10;
select Count(1) from sir_master;

CREATE VIEW 'hv_master' AS 
SELECT StudyID, 
PatientID, 
VtfHuvudId, 
Unit, 
AdmissionYear, 
Age, 
SexMale, 
CareType, 
MotherClinic, 
EmergencyAdmission, 
UnderwentSurgery, 
AdmittedFrom, 
AdmissionReason,
datetime(AdmissionDateTime, 'unixepoch', 'localtime') AS AdmissionDateTime, 
DischargedTo, 
DischargeReason, 
datetime(DischargeDateTime, 'unixepoch', 'localtime') AS DischargeDateTime,
CareTimeMinutes, 
Outcome, 
datetime(DeathDateTime, 'unixepoch', 'localtime') AS DeathDateTime,
PrimaryICUICD, 
TreatmentStrategy, 
VTS, 
VTS2014
FROM 'sir_master';