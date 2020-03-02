git clone https://github.com/Microsoft/vcpkg
vcpkg\bootstrap-vcpkg.bat
vcpkg\vcpkg integrate install
dir .github
for /F "delims=" %%a in (.github/vcpkg-requirements.txt) do ( vcpkg\vcpkg install %%a )