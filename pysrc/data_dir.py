""" 
The purpose of the scipt is to store the directory in which the data will be saved. The data_path holds all of the data for these runs. The
user then specifies what directory in data_path holds the actual data from the command line in my_utils.py with the -o flag. 

For example:
Since, data_dir.data_path -> /Data
from the command line, I would run the code as 
$main.py -o run1
Thus, I would have /Data/run1 as the path that would hold my data

This script also, allows one to get the path of the object directory which holds the .so files needed to build the python modules
that run fortran code.

The user specifies what path should hold the data in this script.
"""

import os

#-----------------------------------------------
# Script
#-----------------------------------------------
"""
# store the python path
temp	= os.getcwd()

# move up the data heirarchy
for i in range(6):
  temp	= os.path.dirname(temp)

# store the path that holds the data
data_path	= "/" + temp [1:]+ "/Data"

# store the obj file path
obj_path	= os.path.dirname( os.getcwd() )
obj_path	= obj_path + "/myobj"
"""

data_path 	= os.environ["MY_DATA_DIR"]
#obj_path	= os.environ["MY_ELECTRO_DIR"] + "/myobj"
#base_path	= os.environ["MY_BASE_DIR"]