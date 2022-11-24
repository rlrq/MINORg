#!/usr/bin/python3

import os
import sys

from data_manip import splitlines

def splitlines(fname, ignore_empty_lines = True):
    data = open(fname, 'r').read().split('\n')
    if ignore_empty_lines:
        return [line for line in data if line]
    else:
        return data

whl, whl_ver = sys.argv[1:3]

## read docker file
dockerfiles = ["/mnt/chaelab/rachelle/scripts/minorgpy/Dockerfile-lite",
               "/mnt/chaelab/rachelle/scripts/minorgpy/Dockerfile-full"]
# dockerfile = "/mnt/chaelab/rachelle/scripts/minorgpy/Dockerfile"

for dockerfile in dockerfiles:
    contents = splitlines(dockerfile, ignore_empty_lines = False)
    ## modify wheel version to latest
    for i, c in enumerate(contents):
        ## directory matters
        if c[:3] == "ADD" and c[-9:] == ".tar.gz .":
            contents[i] = f"ADD dist/minorg-{whl_ver}.tar.gz ."
        elif c[:15] == "WORKDIR minorg-":
            contents[i] = f"WORKDIR minorg-{whl_ver}"
        ## wheel
        elif c[:4] == "COPY" and c[-10:] == "minorg.whl":
            contents[i] = f"COPY dist/{os.path.basename(whl)} minorg.whl"
    
    ## write upated Dockerfile
    with open(dockerfile, 'w+') as f:
        f.write('\n'.join(contents))


## exit
