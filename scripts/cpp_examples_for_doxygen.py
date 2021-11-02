import os
import glob

for gd_mod in glob.glob("src/*/"):
    mod_files = []
    for paths in [gd_mod + 'utilities', gd_mod + 'example']:
        if os.path.isdir(paths):
            for root, dirs, files in os.walk(paths):
                for file in files:
                    if file.endswith(".cpp"):
                        mod_files.append(str(os.path.join(root, file)).split(paths)[1][1:])
    if len(mod_files) > 0:
        mod = str(gd_mod).split('/')[1]
        print(' * \section ' + mod + '_example_section ' + mod)
        for file in mod_files:
            print(' * @example ' + file)
