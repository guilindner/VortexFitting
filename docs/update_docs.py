import os
import sys
import glob, shutil

is_windows = sys.platform.startswith('win')

if is_windows:
    print("Running on Windows ...")
    cmd='make html'
    os.system(cmd)
    files = glob.iglob(os.path.join('./_build/html', "*.html"))
    for file in files:
        if os.path.isfile(file):
            print(file)
            shutil.copy2(file, './')
else:
    print("Running on Linux ...")
    cmd='make html && cp -R _build/html/* .'
    os.system(cmd)



