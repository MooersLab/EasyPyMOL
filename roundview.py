from __future__ import division
from pymol import stored,cmd

def roundview(storedView = 0, sigDigits = 2, 
    outputName = "roundedview.txt" ):
    """
DESCRIPTION

Adds the command "roundview" that gets a view (default is 0,
the current view; you can get a stored view assigned to some
other digit with the view command) and rounds to two
significant digits (default) the viewpoint matrix elements
and rewrites the matrix elements on a single line with no
whitespaces and a semicolon at the end. The saved space
eases the making of a single line of PyMOL commands
separated by semicolons. This enables rapid and interactive
editing of chunks of PyMOL commands. The viewpoints are
appended to the bottom of a text file in the present working
directory called "roundedview.txt". The line could be easier
to copy from this file than from the command history window
in the external gui. A semicolon with nothing to the right
of it at the end of a line of grouped commands is harmless.

Usage: roundview [view, significant digits, outname] 

    Note that the values in the [] are optional.

The default argument values are 0,2, roundedview.txt

Examples: "run roundview.py" and then call with "roundview" or
store in .pymol/startup/ and then quit and restart pymol and
you will find the function listed in the list of plugins.

Then invoke by typing roundview

            roundview # typical use 
            
            roundview 0,3  # current view with three significant digits
            
            roundview 0,3,settings.txt  # as above but write to "settings.txt"

By default, the settings are appended to the file "roundedview.txt".

Warning: This program will accept any string as a filename.
Unix will take anything except back and forward slashes in filenames.
Avoid whitespaces in the filename of the output file.


18 elements of the view matrix (0-17)

0 - 8 = column-major 3x3 matrix that rotates the model axes 
to camera axes 

9 - 11 = origin of rotation relative to the camera 
in camera space

12 - 14 = origin of rotation in model space 

15 = front plane distance from the camera 

16 = rear plane distance from the camera 

17 = orthoscopic flag 
(not implemented in older versions of PyMOL)


The MIT License (MIT)

Copyright (c) <2015-2016> <Blaine H. M. Mooers, Ph.D.>
<University of Oklahoma Health Sciences Center>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS
OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 """
    #Check the user input
    try:
        storedView = int(storedView) 
    except ValueError:
        print("Sorry, could not parse storedView. It should be an integer, 0 or larger")

    try:
        sigDigits = int(sigDigits)
    except ValueError:
        print("Sorry, could not parse sigDigits. It should be an integer, 0 or larger")
        

#call the get_view function
    m = cmd.get_view(storedView)

#Make a list of the elements in the orientation matrix. 
    myList = [m[0], m[1], m[2], m[3], m[4], m[5], m[6],
    m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14],
    m[15], m[16], m[17]]

#Round of to two significant digits the matrix elements.
#This rounding approach solved the problem of unwanted
#whitespaces when I tried using a string format statement
    myRoundedList = [ round(elem, sigDigits) for elem in myList ]

#x is the format of the output. The whitespace is required
#between the "set_view" and "(". 
    x = 'set_view ({0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},\
{11},{12},{13},{14},{15},{16},{17});' 
#print to the external gui.
    print x.format(*myRoundedList)

#Write to a text file. 
    with open(outputName, "a") as myFile:
        myFile.write(x.format(*myRoundedList) + "\n" ) 
    return

#The extend command makes roundview into a PyMOL command.
cmd.extend("roundview",roundview)
