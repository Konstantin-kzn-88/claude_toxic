@echo off
cd /d "%~dp0"
python gui.py --secondary
if errorlevel 1 pause
