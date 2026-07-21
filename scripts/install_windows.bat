@echo off
setlocal
cd /d "%~dp0\.."
python -m pip install --upgrade pip
python -m pip install .
echo.
echo XFit installed. Start GUI with: xfit-gui
echo Or run: python main.py
pause
