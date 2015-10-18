# EasyPyMOL
Script to facilitate the making of horizontal scripts

-----------------------
roundview.py    
-----------------------
(C) Blaine Mooers,Ph.D.,  Universitry of Oklahoma Health Sciences Center, 2015-2016

##Problem: viewing port matrix from get_view() in PyMOL is too hard to copy and paste quickly onto one line with other commands during horizontal scripting##

	PyMOL>get_view  

	### cut below here and paste into script ###  
	set_view (\
  	 -0.832868993,    0.398498207,    0.383896619,\  
 	 -0.260102808,   -0.894237876,    0.363985002,\  
 	  0.488390923,    0.203309149,    0.848513067,\  
  	  0.000000000,    0.000000000,  -61.396984100,\  
 	 -46.246913910,   -4.663769245,   42.401920319,\  
  	 56.260883331,   66.533096313,  -20.000000000 )  
	### cut above here and paste into script ###

    

##Solution: roundview()##

	PyMOL>roundview  
	set_view(-0.83,0.4,0.38,-0.26,-0.89,0.36,0.49,0.2,0.85,0.0,0.0,-61.4,-46.25,-4.66,42.4,56.26,66.53,-20.0);

    
#### Paste above reformated set_view() command onto the PyMOL command line in the external gui or into a script. ####




--------------------
What this is
------------

This python script includes the short function roundview() to reformat the viewing port setting matrix from six rows to one row for use with other commands separated by a 
command line on a single row.  


--------------------
Requirements
------------

Requires a molecular object loaded into an interactive session of PyMOL. Does not
require any external modules. Tested on Ubuntu 14.04 64 bit with PyMOL 1.7.2.2. 
Tested on Windows 8 32 bit  running PyMOL 1.7.6.2 and PyMOL 1.7.6.6. Tested on 
Mac OSX 10.10.5 64 bit running PyMOL 1.5.0.5 and 1.7.6.6. 

--------------
How to install
--------------

Get the script. Either download the folder from [this link](https://github.com/MooersLab/PymoMooersLab/archive/master.zip) or type in a terminal window in your home directory or the directory where you store your PyMOL python scripts 

	git clone https://github.com/MooersLab/MooersLabPymol/roundview

You need the program <code>git</code> installed on your computer. Git is avialbe via macports on a mac or otherwise [see these instructions for installing git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).


There are at least four ways of loading the script into PyMOL:

1) move roundview.py to the working directory of pymol. In pymol, type 

	run roundview.py


2) load roundview.py with the plugin manager (see plugin pulldown) in PyMOL

3) copy roundview.py to the following hidden folder ~.pymol/startup in your home directory or to the startup folder inside your MacPyMOL.app.

4) create or edit the hidden file <code>.pymolrc</code> text file in the home directory so that it includes the following lines so that roundview.py is always loaded upon startup:


	 import sys
 	 sys.path.append('/Path/To/roundview')
 	 run /Path/To/roundview.py

--------------
How to run
--------------
After loading a pdb file and setting up the molecular scene, type on one of the command lines in PyMOL:


  	PyMOL>roundview



-------------------------
How to get help
---------------
Type the following on ona of the command lines in PyMOL"


  PyMOL> help roundview


Usage: roundview [view, significant digits, outname] # the values in the [] are optional.  
The default view is "0". The default number of significant digits is 2. The "outname" is the name of a plain text file to which the output of roundveiw() is written.

-------------------------------------------------------------------------------------------
Quick test with wildtype bacteriophage T4 lysozyme, 3fa0
----------------------------------------------------

Open up pymol. Copy and paste the entire line below onto the command line in the external gui (the on above the GL viewing port). This is an example of a horizontal script. By hitting the up arrow key, you can recall this command for editing on the command line. This code block is  more agile to edit than opening, editing, saving, and loading an external script file.


  	fetch 3fa0,async=0;orient;turn z,-90;turn y,-5;turn x,10; hide everything; bg_color white; show cartoon;color red, ss h;color yellow, ss s;color green, ss l+'';roundview  
  
To apply the canonical view similar to that in Gassner et al. 2004, copy and paste the following onto the command line:
  
  	set_view (-0.18,-0.69,-0.7,0.98,-0.17,-0.09,-0.06,-0.7,0.71,0.0,0.0,-154.87,34.77,11.27,9.52,121.27,188.47,-20.0); ray 1500,1600; png coanonical.png
  	
You should get back an image that looks like the following:



<img src="https://cloud.githubusercontent.com/assets/15176203/10561743/5f59bb82-74fd-11e5-828c-dbe1dbb2c648.png" width="90%"></img> 






To test some other argument values, copy paste the following command. 

  	PyMOL>roundview 0,1,firstscene.txt

Do a <code>ls *.txt</code> to list the files in the working diretory. The file "firstscene.txt" should be listed. The default filename is roundedview.txt. These file is
appended with each execution of the roundview command.


--------------
License
---------
MIT License

--------------
Reference
---------

Mooers, B. H. M. (submitted) Easier use of PyMOL with horizontal scripts.
