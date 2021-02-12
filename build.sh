# package as source and binary
python3 setup.py sdist bdist_wheel
# load package to test pypi
python3 -m twine upload dist/*
