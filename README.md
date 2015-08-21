# mvision
Matlab implementation of the neuron-finder pipeline of vision

Cloning
New cloning users shall compile all mex and java files within the repository.
In matlab:
mex -setup % If necessary
mex file
...

In UNIX shell
Provides the list of source java files
compile appropriately with either
-cp cellFinder/Cell-Finder.jar
or
-cp vision/Vision.jar
if necessary
find . '-name' '*.java'
