#!/bin/tcsh
# This file to compile java and mex files in the repository after clones and updates.
# Note that this repository works with a local copy of vision and cell finder jars to avoid multiplying dependencies.

find . -name "*.java" -print | xargs javac -source 1.7 -target 1.7 -cp "./vision/Vision.jar:./cellfinder/Cell-Finder.jar"
find . -name "*.cpp" -print | xargs mex
