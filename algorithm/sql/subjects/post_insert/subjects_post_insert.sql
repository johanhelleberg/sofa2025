PRAGMA table_info('subjects');				  
SELECT * FROM subjects LIMIT 10;
select Count(1) from subjects;

CREATE VIEW 'hv_subjects' AS
SELECT 
StudyID, 
Dead,
datetime(DeathDateTime, 'unixepoch', 'localtime') AS DeathDateTime
FROM 'subjects';