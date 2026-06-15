#!/usr/bin/env python

import os

source_dir = "/Users/liani/nu_work/plotting/surprise_plotting" # the folder you copied the files from
target_folder = "/Users/liani/nu_work/pelee" # the folder you copied the files to

for root, dirs, files in os.walk(source_dir):
    for name in files:
        try:
            os.remove(target_folder+"/"+name)
        except FileNotFoundError:
            pass

