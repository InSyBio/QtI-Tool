#!/bin/bash

#sudo -s

echo "--- find SCRIPTPATH---"
SCRIPTPATH="/opt/insybio-qti/"


echo "--- Do an update - upgrade  ---"
apt update
apt upgrade -y
apt dist-upgrade -y
apt autoremove
apt autoclean

echo "--- Install Git  ---"

apt install git -y

echo "--- Install OpenMS dependencies from binary  ---"
#wget https://cmake.org/files/v3.10/cmake-3.10.0-Linux-x86_64.sh
#chmod a+x cmake-3.10.0-Linux-x86_64.sh
#sudo ./cmake-3.10.0-Linux-x86_64.sh

apt install build-essential cmake autoconf patch libtool automake -y
apt install qt4-default libqtwebkit-dev -y
apt install libeigen3-dev libwildmagic-dev libxerces-c-dev libboost-all-dev libsvn-dev libgsl-dev libbz2-dev -y
apt install libqt4-sql-sqlite
apt install seqan-dev seqan-apps
apt install eigensoft libsvm3 libsvm-dev libsvm-tools
sudo apt install sqlite3 sqlite3-pcre libsqlite3-*

echo "--- Create insybio-qti folder ---"
mkdir $SCRIPTPATH
cd $SCRIPTPATH

echo "--- Get OpenMS ---"
git clone https://github.com/OpenMS/OpenMS.git

echo "--- Get and install OpenMS dependencies from source ---"
git clone  https://github.com/OpenMS/contrib.git
mkdir $SCRIPTPATH/contrib-build
cd $SCRIPTPATH/contrib-build
cmake -DBUILD_TYPE=LIST ../contrib
cmake -DBUILD_TYPE=SEQAN ../contrib
cmake -DBUILD_TYPE=WILDMAGIC ../contrib
cmake -DBUILD_TYPE=EIGEN ../contrib

 
cmake -DBUILD_TYPE=ALL -DNUMBER_OF_JOBS=4 ../contrib

echo "--- Build Open MS ---"

cd ..
mkdir OpenMS-build
cd OpenMS-build
echo "--- HOME:  SCRIPTPATH ---"
cmake -DOPENMS_CONTRIB_LIBS="$SCRIPTPATH/insybio-qti/contrib-build" -DBOOST_USE_STATIC=OFF -DCMAKE_PREFIX_PATH="/usr;/usr/local" ../OpenMS
make


echo "--- Download thirdparty ---"
cd $SCRIPTPATH
git clone  https://github.com/OpenMS/THIRDPARTY.git
mv THIRDPARTY OpenMS-thirdparty


#export LD_LIBRARY_PATH="$SCRIPTPATH/OpenMS-build/lib:$LD_LIBRARY_PATH"
#export PATH=$PATH:$SCRIPTPATH/OpenMS-build/bin
#export PATH=$SCRIPTPATH/OpenMS-thirdparty/Linux/64bit/XTandem:${PATH}
#source ~/.bashrc

echo "--- Install pip python-libsvn libgs-dev---"
sudo apt-get install python-pip python-libsvm libgs-dev numpy

echo "--- Instal pyopenms 2.3.0.3 ---"
pip install pyopenms

