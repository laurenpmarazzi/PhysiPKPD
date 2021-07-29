# PhysiPKPD
This repository contains the code for PhysiPKPD, a module template for modeling drug treatments in PhysiCell.

### Running instructions
Run the following commands from the PhysiCell root directory (assuming you have this repository as a subfolder within that directory)
```
make clean
make reset
make data-cleanup
make template

cp PhysiPKPD/src/custom.* custom_modules/
cp PhysiPKPD/src/main.cpp ./
cp PhysiPKPD/src/mymodel.xml config/

make
project.exe ./config/mymodel.xml
```
If on Mac, the last command will instead read
```
./project ./config/mymodel.xml
```
