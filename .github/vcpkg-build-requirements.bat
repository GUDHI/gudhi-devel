@ECHO OFF
git clone https://github.com/Microsoft/vcpkg
vcpkg\bootstrap-vcpkg.bat
dir
for /F "delims=" %%a in (vcpkg-requirements.txt) do ( vcpkg\vcpkg install %%a )