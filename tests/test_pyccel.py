# coding: utf-8

# Usage:
#   python tests/test_pyccel.py --execute

import os

from pyccel.commands.console import pyccel
from utils import clean_tests

# ...
def test_1d():
    print('============== testing 1d ================')
    ignored = []

    base_dir = os.getcwd()

    # ...
    def _get_files(path):
        path_dir = os.path.join(base_dir, path)

        files = sorted(os.listdir(path_dir))
        files = [f for f in files if not(f in ignored) and (f.endswith(".py"))]
        return [os.path.join(path_dir, f) for f in files]
    # ...

    # ...
    folders = ['tests/pyccel/1d']
    # ...

    # ...
    files = []
    for r in folders:
        files += _get_files(r)
    # ...

    # ...
    for f_name in files:
        f = os.path.basename(f_name)

        pyccel(files=[f_name])
        print(('> testing {0}: done'.format(str(f))))
    # ...
# ...

# ...
def test_2d():
    print('============== testing 2d ================')
    ignored = ['ex5b.py', 'ex7b.py']

    base_dir = os.getcwd()

    # ...
    def _get_files(path):
        path_dir = os.path.join(base_dir, path)

        files = sorted(os.listdir(path_dir))
        files = [f for f in files if not(f in ignored) and (f.endswith(".py"))]
        return [os.path.join(path_dir, f) for f in files]
    # ...

    # ...
    folders = ['tests/pyccel/2d']
    # ...

    # ...
    files = []
    for r in folders:
        files += _get_files(r)
    # ...

    # ...
    for f_name in files:
        f = os.path.basename(f_name)

        pyccel(files=[f_name])
        print(('> testing {0}: done'.format(str(f))))
    # ...
# ...

# ...
def test_3d():
    print('============== testing 3d ================')
    ignored = []

    base_dir = os.getcwd()

    # ...
    def _get_files(path):
        path_dir = os.path.join(base_dir, path)

        files = sorted(os.listdir(path_dir))
        files = [f for f in files if not(f in ignored) and (f.endswith(".py"))]
        return [os.path.join(path_dir, f) for f in files]
    # ...

    # ...
    folders = ['tests/pyccel/3d']
    # ...

    # ...
    files = []
    for r in folders:
        files += _get_files(r)
    # ...

    # ...
    for f_name in files:
        f = os.path.basename(f_name)

        pyccel(files=[f_name])
        print(('> testing {0}: done'.format(str(f))))
    # ...
# ...

################################
if __name__ == '__main__':
    clean_tests()
    test_1d()
    test_2d()
    test_3d()
    clean_tests()
