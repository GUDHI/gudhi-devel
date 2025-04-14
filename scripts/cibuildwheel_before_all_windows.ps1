set VCPKG_BUILD_TYPE=release
vcpkg install eigen3 cgal --triplet x64-windows
vcpkg version
Get-ChildItem "C:\vcpkg\installed\x64-windows\bin\"
