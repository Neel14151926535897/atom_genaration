@echo off
SET RAYLIB_PATH=C:\raylib\raylib
SET COMPILER_PATH=C:\raylib\w64devkit\bin
SET PATH=%COMPILER_PATH%;%PATH%
gcc main.c atom.c atom.res -o atom.exe -Iraylib/src -Lraylib/src -lraylib -lopengl32 -lgdi32 -lwinmm -Oz -s -ffunction-sections -fdata-sections -Wl,--gc-sections,--strip-all -mwindows -fno-ident -fno-asynchronous-unwind-tables -fno-stack-protector -fno-plt
upx --ultra-brute --lzma atom.exe
if %ERRORLEVEL% EQU 0 (
    echo Build successful! Running simulation...
    atom.exe 118 500
) else (
    echo Build failed with error code %ERRORLEVEL%.
)