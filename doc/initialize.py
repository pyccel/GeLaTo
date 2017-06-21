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

    # ... cleaning
    cmd = "rm -rf tutorials/" + fname + "*"
    os.system(cmd)
    # ...

    # ...
    f_out = dirname + "/" + fname + ".rst"
    cmd = "mv " + f_out + " tutorials/"
    os.system(cmd)
    # ...

    # ...
    d_out = dirname + "/" + fname + "_files"
    cmd = "mv " + d_out + " tutorials/"
    os.system(cmd)
    # ...
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
#convert_nb("tutorial_1.ipynb")
#convert_nb("tutorial_2.ipynb")
#convert_nb("tutorial_3.ipynb")
#convert_nb("tutorial_4.ipynb")
#convert_nb("tutorial_5.ipynb")
convert_nb("tutorial_6.ipynb")
#convert_nb("tutorial_7.ipynb")
# ...
