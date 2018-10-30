#!/usr/bin/env bash

sudo apt-get update
sudo apt-get install g++ -y
sudo apt-get install cmake -y
sudo apt-get install libopenblas-dev -y
sudo apt-get install liblapack-dev -y
sudo apt-get install libarpack2-dev -y
sudo apt-get install libarmadillo-dev -y

#optimLib
git clone https://github.com/kthohr/optim.git
cd optim
./configure -i "/usr/local" -p
sudo make
sudo make install
sudo /sbin/ldconfig -v

echo "--------------------------------------------------------------------"
echo "OPTIMLIB DONE"
echo "FREEGLUT STARTING"
echo "--------------------------------------------------------------------"
sudo apt-get install binutils-gold -y
sudo apt-get install freeglut3 freeglut3-dev libgl1-mesa-dev libglu1-mesa-dev -y

echo "--------------------------------------------------------------------"
echo "FREEGLUT DONE"
echo "GUI STARTING"
echo "--------------------------------------------------------------------"

sudo apt-get install -y xfce4 virtualbox-guest-dkms virtualbox-guest-utils virtualbox-guest-x11
sudo VBoxClient --clipboard
sudo VBoxClient --draganddrop
sudo VBoxClient --display
sudo VBoxClient --checkhostversion
sudo VBoxClient --seamless