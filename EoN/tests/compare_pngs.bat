for %%i in (ExpectedPng\*.png) do start "" "C:\Program Files\TortoiseSVN\bin\TortoiseIDiff.exe" ..\%%~nxi ExpectedPng\%%~nxi
