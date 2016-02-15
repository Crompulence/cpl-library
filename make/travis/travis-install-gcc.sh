#!/bin/sh
set -e
sudo apt-get install build-essential
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get update
case $1 in
  gcc-4.6) set -x;
    sudo apt-get install gcc-4.6 g++-4.6;;
  gcc-4.7) set -x;
    sudo apt-get install gcc-4.7 g++-4.7;;
  gcc-4.8) set -x;
    sudo apt-get install gcc-4.8 g++-4.8;;
  gcc-4.9) set -x;
    sudo apt-get install gcc-4.9 g++-4.9;;
  gcc-5) set -x;
    sudo apt-get install gcc-5 g++-5;;
  *)
    echo "Unknown gcc version:" $1; exit 1;;
esac
