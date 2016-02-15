#!/bin/sh
set -e
sudo apt-get install build-essential
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt-get update -q
case $1 in
  gcc-4.8) set -x;
    sudo apt-get install gcc-4.8 g++-4.8 -y;;
  gcc-4.9) set -x;
    sudo apt-get install gcc-4.9 g++-4.9 -y;;
  gcc-5) set -x;
    sudo apt-get install gcc-5 g++-5 -y;;
  *)
    echo "Unknown gcc version:" $1; exit 1;;
esac
