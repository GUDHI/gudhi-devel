import sysconfig

# Get CFLAGS from sysconfig (this is what is done in setuptools for CXXFLAGS)
# and remove '-I...' as this is done by setuptools (even if CFLAGS or CXXFLAGS is set)
lst = [str for str in sysconfig.get_config_var("CFLAGS").split() if not (str.startswith("-I"))]

# Remove also '-isystem xxx'
indices = [i for i, v in enumerate(lst) if v == "-isystem"]
for index in reversed(indices):
    del lst[index : index + 2]

# Get unique compilation flags and print it for CMake
print("'%s'," % "', '".join(set(lst)))
