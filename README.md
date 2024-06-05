# DavisRead
Code for interpreting Davis Rig .ms8.txt files

Consider all of this a work-in-progress.

Included are two scripts, one for R and one for Python, that each contain a function for interpreting .ms8.txt files.
Using the R file is fairly straightforward; in a separate analysis script, use source("path/to/Davis_Read.r") to run that script, which will add a function to your environment that does the interpretation.

TBD on the python script. I think it's functional, but I don't yet know the best way to use it.

# Getting Started:
To download the code and set it up, open a terminal and run the following:
```
git clone https://github.com/Mraymon5/DavisRead.git #Download the repository
cd <path/to/DavisRead>/requirements #Move your working directory to the repository
conda create --name DavisRead python=3.11 #Make a new environment for the repository
conda activate DavisRead #Activate the environment
pip install -r requirements.txt #Install dependencies
```

Once the code has been downloaded and requirements installed, you can move on to using the script.
It's packaged as a function, so you don't need to open the script itself, unless you want to alter its behavior.
Instead, you can import the function into an analysis script, or in the terminal if you're one of those.
Open a new script window, or open your pre-existing analysis script, and at the top with the rest of your imports, add:
```
from DavisRead.Davis_Read import Davis_Read
```
Then, you can call the function like so:
```
output = Davis_Read(folder="/path/to/your/data/", output_folder="/path/to/save/destination/")
```
