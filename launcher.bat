@echo off
chcp 65001 >nul
title Advanced Chemistry Engine — Kirill Mikhailovich Special Edition
color 0B

:: УКАЖИТЕ ЗДЕСЬ ВАШУ ССЫЛКУ НА GITHUB RAW (без слэша на конце)
set "REPO_BASE_URL=https://github.com/sv5190228-web/AdvancedChemistryEngine-KM-edition.git"

echo =========================================================================
echo   🎂 Advanced Chemistry Engine: Ляпишев Кирилл Михайлович Edition 🎂
echo =========================================================================
echo.

:: 1. Проверка дистанционных обновлений (OTA)
echo [1/4] Проверка наличия обновлений на сервере разработчика...
powershell -Command "[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; try { Invoke-WebRequest -Uri '%REPO_BASE_URL%/server.py?t=%RANDOM%' -OutFile 'server.py.tmp' -TimeoutSec 5 -UseBasicParsing; Move-Item -Force 'server.py.tmp' 'server.py'; Invoke-WebRequest -Uri '%REPO_BASE_URL%/AdvancedChemistryVisualisator_3.html?t=%RANDOM%' -OutFile 'html.tmp' -TimeoutSec 5 -UseBasicParsing; Move-Item -Force 'html.tmp' 'AdvancedChemistryVisualisator_3.html'; Write-Host '      [+] Обновления успешно загружены и применены!' -ForegroundColor Green } catch { if (Test-Path 'server.py.tmp') { Remove-Item 'server.py.tmp' }; if (Test-Path 'html.tmp') { Remove-Item 'html.tmp' }; Write-Host '      [*] Запуск актуальной локальной версии.' -ForegroundColor Yellow }"

:: 2. Запуск вычислительного сервера Python
echo [2/4] Запуск химического суперкомпьютера (Python + RDKit)...
start "Chemistry Server [Special Edition]" cmd /k "python server.py"

:: 3. Ожидание инициализации Flask
echo [3/4] Инициализация нейронных и квантовых модулей...
timeout /t 2 /nobreak >nul

:: 4. Открытие 3D визуализатора
echo [4/4] Открытие 3D-лаборатории...
start AdvancedChemistryVisualisator_3.html

echo.
echo Лаборатория запущена! С Днем Рождения, Кирилл Михайлович! 🎉
timeout /t 3 >nul
exit
