To generate test files:
=======================

Run: "python testgen.py MATLAB ph.test >phtest.m", and then execute phtest.m under MATLAB (of course, after loading the BuTools package by BuToolsInit)

For other environments:

"python testgen.py Python ph.test >phtest.py"
"python testgen.py Mathematica ph.test >phtest.txt" --> then copy the contents of the file to a notebook and evaluate the test function

To generate examples for the documentation:
===========================================

Step 1: create example generator scripts and generate examples

Run: "python testgen.py MATLABDocEx ph.test >phdocex.m", then execute phdocex.m under MATLAB (of course, after loading the BuTools package by BuToolsInit). The output is PH_matlab.docex.



Step 2: add examples to the documentation
