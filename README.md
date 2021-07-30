# PhysiPKPD
This repository contains the code for PhysiPKPD, a module template for modeling drug treatments in PhysiCell.

### Running instructions
Run the following commands from the PhysiCell root directory (assuming you have this repository as a subfolder within that directory)
```
make clean
make reset
make data-cleanup
make template

cp PhysiPKPD/custom_modules/custom.* custom_modules/
cp PhysiPKPD/main.cpp ./
cp PhysiPKPD/config/mymodel.xml config/
cp PhysiPKPD/Makefile ./ # this will overwrite your Makefile and any custom additions you have made!
make
project.exe ./config/mymodel.xml
```
If on Mac, the last command will instead read
```
./project ./config/mymodel.xml
```

### Data visualization 
Run the following commands from the directory in which you ran the project to create a GIF, movie, and/or data plots.
```
make gif
make movie
make data-plots
```
