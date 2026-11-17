#!/bin/bash

set -e

# Delete previously downloaded openms
rm -rf bin/openms/

# Download openms
wget -O openms.deb https://github.com/OpenMS/OpenMS/releases/download/Release2.8.0/OpenMS-2.8.0-Debian-Linux-x86_64.deb

# Extract only the data-part of the deb package
ar x openms.deb data.tar.gz

# Extract the data into the openms_folder
mkdir -p ./bin/openms

# Delete downloaded and extracted deb package
tar -xf data.tar.gz -C ./bin/openms
rm -rf openms.deb data.tar.gz

# Download Pia:
rm -rf bin/pia/
wget -O pia.zip https://github.com/medbioinf/pia/releases/download/1.4.8/pia-1.4.8.zip
mkdir -p ./bin/pia
unzip -d ./bin/pia pia.zip
rm -rf pia.zip


# Download Dinosaur
wget  -O bin/Dinosaur.jar https://github.com/fickludd/dinosaur/releases/download/1.2.0/Dinosaur-1.2.0.free.jar




### Install tdf2mzml
rm -rf bin/tdf2mzml
git clone https://github.com/mafreitas/tdf2mzml bin/tdf2mzml
OLD_PWD=$(pwd)
cd bin/tdf2mzml
git checkout 85c2b258ea418369ee8753294d8b57978d6ec747
pip install -r requirements.txt
cd $OLD_PWD
# You can call tdf2mzml via: LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$pwd/bin/tdf2mzml python bin/tdf2mzml tdf2mzml.py 


### Install python dependencies
python -m pip install numpy==1.23.5
python -m pip install fisher-py==1.0.21
python -m pip install pyopenms
python -m pip install pandas
python -m pip install sqlalchemy
python -m pip install plotly
python -m pip install scikit-learn
python -m pip install alphatims
python -m pip install nextflow
python -m pip install streamlit==1.28.2


# Make the executables visible in nextflow
chmod +x bin/*
