@echo off
setlocal enabledelayedexpansion

REM Define source and target base directories
set "SOURCE_DIR=D:\Single Axle Cases"
set "TARGET_BASE_DIR=C:\Users\johannc2\Box\FAA Data\04_FEM\00_FEM DATA\Single_Thin\Single_Thin_Responses"

REM Iterate through each folder in the source directory
for /d %%D in ("%SOURCE_DIR%\*") do (
    REM Capture the folder name
    set "FOLDER_NAME=%%~nxD"
    
    REM Construct the target directory path
    set "TARGET_DIR=%TARGET_BASE_DIR%\!FOLDER_NAME!"
    
    REM Check if the corresponding folder exists in the target directory
    if exist "!TARGET_DIR!" (
        REM Copy the specific files to the corresponding target directory
        echo Folder found: !FOLDER_NAME!. Copying files...
        
        REM Copy the .inp file
        copy "%%D\!FOLDER_NAME!.inp" "!TARGET_DIR!" >nul
        
        REM Copy the 3DResponse files
        for /l %%i in (1,1,17) do (
            if exist "%%D\!FOLDER_NAME!_3DResponse_tire%%i.txt" (
                copy "%%D\!FOLDER_NAME!_3DResponse_tire%%i.txt" "!TARGET_DIR!" >nul
            )
        )
    ) else (
        echo Folder !FOLDER_NAME! does not exist in the target directory. Skipping...
    )
)

echo Operation completed.
pause