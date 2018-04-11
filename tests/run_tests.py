
import os
import sys

# Useful for very coarse version differentiation.
PY2 = sys.version_info[0] == 2
PY3 = sys.version_info[0] == 3
PY34 = sys.version_info[0:2] == (3, 4)
PY35 = sys.version_info[0:2] == (3, 5)

if PY34 or PY35:
    cmd = 'py.test'
else:
    cmd = 'pytest'

dirs = ['tests/expressions', 'tests/weak_form', 'tests/glt']
for d in dirs:
    os.system('{cmd} {dir}'.format(cmd=cmd, dir=d))
