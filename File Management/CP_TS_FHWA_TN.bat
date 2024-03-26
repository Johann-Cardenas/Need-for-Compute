@echo off
setlocal enabledelayedexpansion

REM Define source and target base directories
set "SOURCE_DIR=G:\FHWA_Thin"
set "TARGET_BASE_DIR=C:\Users\johannc2\Box\FAA Data\04_FEM\00_FEM DATA\FHWA_Thin\FHWA_Thin_Responses"

REM Iterate through each folder in the source directory
for /d %%D in ("%SOURCE_DIR%\*") do (

  REM Iterate over each subfolder in the current folder
  for /d %%S in ("%%D\*") do (
      REM Capture the subfolder name
      set "SUBFOLDER_NAME=%%~nXS"
      
      REM Construct the target directory path
      set "TARGET_DIR=%TARGET_BASE_DIR%\!SUBFOLDER_NAME!"
        
      REM Check if the corresponding folder exists in the target directory
      if not exist "!TARGET_DIR!" (
          echo Folder !SUBFOLDER_NAME! does not exist in the target directory. Creating...
          mkdir "!TARGET_DIR!"
      )
      
      REM Now that the directory is guaranteed to exist, copy the files
      echo Folder found or created: !SUBFOLDER_NAME!. Copying files...
      
      REM Copy the .inp file
      REM copy "%%S\*.inp" "!TARGET_DIR!" >nul
      
      REM Copy the 3DResponse files
      for /l %%i in (1,1,17) do (
          if exist "%%S\*_3DResponse_Tire%%i.txt" (
              copy "%%S\*_3DResponse_Tire%%i.txt" "!TARGET_DIR!" >nul
            )
        ) 
    )
)

echo Operation completed.
pause