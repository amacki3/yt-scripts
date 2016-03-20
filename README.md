# yt-scripts
scripts for MPhys project

This is a repository of codes used in my 5th year MPhys project to investigate missing baryons in Large scale structure.
The code is in python and depends on yt, numpy and matplotlib.

Convds.py was provided by Britton Smith, filament.py was also provided but has been modified.

There are known issues:
  Currently many of the parallel modules are broken due to some form of bug with yt parallal_analysis_interface,
  a script that should fix this has been included, but the other scripts will need slight reworks.
  
  denstemp.py needs reorganised such that the main function is not as ugly and is decomposed into its seperate functions
  
  A lot of the code needs sanitised of hard coded variables, such as the directory data is stored/loaded from
  
  Some code lacks meaningful commentary
  Code is not uniform in style/syntax or variable naming.
