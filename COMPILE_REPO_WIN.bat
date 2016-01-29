@ECHO OFF

REM This file to compile java and cpp-mex files in the repository after clones and updates.
REM Note that this repository works with a local copy of vision and cell finder jars to avoid multiplying dependencies.

dir /s /B *.java > .\javacSourcesTemp.txt
javac -source 1.7 -target 1.7 -cp ".\vision\Vision.jar;.\cellfinder\Cell-Finder.jar" @javacSourcesTemp.txt
rm .\javacSourcesTemp.txt

REM If this doesn't work, proceed to "mex -setup" and set compiler appropriately
dir /s /B *.cpp > mexSourcesTemp.txt
start /b /wait mex -silent @mexSourcesTemp.txt
REM mex runs a background thread and it is not actually easy to wait
REM appropriately to remove the txt file.