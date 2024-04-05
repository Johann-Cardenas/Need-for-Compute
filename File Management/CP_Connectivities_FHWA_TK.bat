@echo off
setlocal enabledelayedexpansion

REM Define source and target directories
set "SOURCE_DIR=E:\FHWA_Thick"
set "TARGET_DIR=C:\Users\johannc2\Box\FAA Data\04_FEM\00_FEM DATA\FHWA_Thick\FHWA_Thick_Connectivity"


REM Iterate through each folder in the source directory
for /d %%D in ("%SOURCE_DIR%\*") do (

  REM Iterate over each subfolder in the current folder
  for /d %%S in ("%%D\*") do (
      REM Capture the subfolder name
      set "SUBFOLDER_NAME=%%~nXS"

      REM Define the specific files to be copied, matching the naming pattern
      set "ELEM_FILE=%%S\!SUBFOLDER_NAME!_Elem.txt"
      set "NODES_FILE=%%S\!SUBFOLDER_NAME!_Nodes.txt"
      
      REM Check if the Elem file exists and copy it
      if exist "!ELEM_FILE!" (
          echo Copying file: !ELEM_FILE!
          copy "!ELEM_FILE!" "%TARGET_DIR%" >nul
      ) else (
          echo Elem file not found in folder: !SUBFOLDER_NAME!
      )
      
      REM Check if the Nodes file exists and copy it
      if exist "!NODES_FILE!" (
          echo Copying file: !NODES_FILE!
          copy "!NODES_FILE!" "%TARGET_DIR%" >nul
      ) else (
          echo Nodes file not found in folder: !SUBFOLDER_NAME!
      )
  )
)

echo All applicable files have been copied.
pause
