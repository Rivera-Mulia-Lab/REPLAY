@echo off
setlocal

set "SCRIPT_DIR=%~dp0"

for /f "delims=" %%i in ('wsl wslpath "%SCRIPT_DIR%"') do set "WSL_DIR=%%i"

wsl bash -lc "cd \"$WSL_DIR\" && chmod +x ./replay && ldd ./replay && echo ==== RUNNING ==== && ./replay"
pause

endlocal