#!/bin/bash

# Set up ROOT path for python plotting scripts
# Should be a valid ROOT version with PyROOT setup
export ROOTLIB='/home/zhangweibin/bin/root/lib'

# Path for Python Modules
# This is set automatically
export PATH_OnePlot_MODULE1=`pwd`
export PATH_OnePlot_MODULE2=$PATH_OnePlot_MODULE1"/modules/"

# Path for atlas plotting style
export PATH_AtlasStyle=`pwd`"/atlas-style/"
