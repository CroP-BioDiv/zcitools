# zcitools
Tools used in ZCI research

Main idea of the project is to store python code for running calculations and calculation results
in some organized way for simpler reusability, or at least to later find what was calculated in the first place.

Generally, research is done in steps. Each step depends on some previous step(s) and/or human decisions.


## Project structure

Structure of a project results is simple. Project is a directory, and each step is it's subdirectory.
Step subdirectory stores data (mostly calculation result) and information how it was produced.

Step data has type so it is known how to read it.

History information are used for tracking what was done in a research, and to make a script that can
reproduce whole process.


## Code structure

ToDo: ...
