@echo off
:Start 
::Start the program
echo Running runSleepStudy.exe ...
echo Until 8:00 AM, code will restart if it is closed. To prevent it from restarting, press Ctrl-C and answer Y (yes).
start "" "C:\Users\Pragya Sharma\Documents\LabVIEW Data\0_NCS_Code\builds\RunNCS\Fifth Build\RunNCS.exe"

:CheckRunning
::Every 5 seconds, check if runSleepStudy.exe is still running. If it is running and it is past 8 AM but before the afternoon, close the program. 
::If it is running but not past 8 am, let it continue to run. If it is NOT running and past 8 am, do nothing. If it is NOT running and
::past 8 am, restart the program.
choice /d y /t 5 > nul

for /F "tokens=1-2" %%a in ('time /t') do (SET mytime=%%a && SET dayHalf=%%b) 
SET hour=%mytime:~0,2%


SET keepRunning=true
if %dayHalf%==AM 	(
	if %hour% GEQ 08	(
		if %hour% LSS 12 (
		SET keepRunning=false
		)
	)
)

tasklist | find /I "runSleepStudy.exe" 

if %errorlevel%==0 (
    ::running
	if %keepRunning%==true (
		goto CheckRunning
	) else (
		taskkill /f /t /im runSleepStudy.exe
	)
) 
if %errorlevel%==1 (
	::not running
    if %keepRunning%==true (
		goto Start
	)
)

echo a|xcopy C:\Users\Public\SleepCenter\Data\* "C:\Users\SleepCenter2\Box Sync\Data\" /Y /E /Q