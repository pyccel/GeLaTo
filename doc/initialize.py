# coding: utf-8

# TODO use join for dir and filename

import os

# ...
def convert_nb(filename):
    root = os.path.dirname(os.path.realpath(__file__))

    dirname = "../notebooks/tutorials"
    os.chdir(dirname)

    cmd = "jupyter nbconvert --to=rst --execute " + filename
    os.system(cmd)

    os.chdir(root)

    fname = filename.split(".")[0]

    f_out = dirname + "/" + fname + ".rst"
    cmd = "mv " + f_out + " tutorials/"
    os.system(cmd)

    d_out = dirname + "/" + fname + "_files"
    cmd = "mv " + d_out + " tutorials/"
    os.system(cmd)
# ...

# ...
cmd = "cp ../README.rst ."
os.system(cmd)
# ...

# ...
cmd = "mkdir -p tutorials"
os.system(cmd)
# ...

# ...
convert_nb("getting_started_1.ipynb")
# ...
