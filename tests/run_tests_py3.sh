python3 tests/expressions/test_calculus.py
python3 tests/expressions/test_1d.py
python3 tests/expressions/test_2d.py
python3 tests/expressions/test_3d.py

python3 tests/weak_form/test_1d.py
python3 tests/weak_form/test_2d.py
python3 tests/weak_form/test_3d.py

cd tests/glt/; python3 test_1d.py; cd ../.. 
cd tests/glt/; python3 test_2d.py; cd ../..
cd tests/glt/; python3 test_3d.py; cd ../..

cd tests/fem/; python3 test_1d.py; cd ../.. 
cd tests/fem/; python3 test_2d.py; cd ../..
cd tests/fem/; python3 test_3d.py; cd ../..

