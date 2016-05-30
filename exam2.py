from __future__ import print_function
# -*- coding: utf-8 -*-
"""
    DESCRIPTION AND INSTALLATION

        Defines aliases from exam 2. of the OUHSC Macromolecular Systems course.
         
        Create ~/Scripts/PyMOLScripts
        and store the script in this subfolder or store in some other folder of your choosing 
        if you know what you are doing. 
        
        Enter on the command line in PyMOL the following command:
        
        run ~/Scripts/PyMOLScripts/exam2function.py
        
        Now the aliases q1,q2, ..., q8 are active.
        
        Tested on PyMOL versions 1.5.0.5 AND 1.8.0.5. 
                
    Copyright Notice
    ================
    Copyright (C) 2016  Blaine Mooers

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    The source code in this file is copyrighted, but you can
    freely use and copy it as long as you don't change or remove any of
    the copyright notices.

    Blaine Mooers , PhD 
    blaine-mooers@ouhsc.edu
    975 NE 10th St, BRC 466
    University of Oklahoma Health Sciences Center
    Oklahoma City, OK, USA 73104
""" 
from pymol import cmd
__author__ = "Blaine Mooers"
__copyright__ = "Blaine Mooers, University of Oklahoma Health Sciences Center, Oklahoma City, OK, USA 73104"
__license__ = "GPL-3"
__version__ = "1.0.2"
__credits__ = ["William Beasley","Chiedza Kanyumbu"] 
# people who reported bug fixes, made suggestions, etc. 
__date__ = "30 May 2016"
__maintainer__ = "Blaine Mooers"
__email__ = "blaine-mooers@ouhsc.edu"
__status__ = "Production" 

print("2015 OUHSC MSI-Lecture 17 and 18 takehome exam")
print("on nucleic acid structure.")
print(" ")
print("Enter 'q1' for questin 1, 'q2' for question 2, ...,")
print("'q8' for question 8.")
print(" ")
print("Questions (1-4). HIV reverse transcriptase bound to DNA")
print("with AZT incorporated (PDB-ID 3V6D).")
print(" ")
print("Questions (5-8). The full-length hammerhead")
print("ribozyme (PDB-ID 3zp8).")
print(" ")
print("Enter 'help exam2' to learn how to use this script.'")
print("Enter 'help q1' to see quesiton 1 and the corresponding horizontal script.'")
print("All or part of the horizontal script can be copied from the command history window for reuse.'")


def exam2():
    """
    USAGE:
    
        Enter 'run exam2.py' without including the quotation
        marks. Then enter 'q1' for question 1. When you are
        ready for question 2, enter 'q2' for question 2, and
        so on. 

    DESCRIPTION

        Defines the molecular scenes for MSI-Lecture 17 and
        18 takehome exam.    
        
        Each question with a molecular scene has the scene
        defined by a specific function qX(), where X is the
        number of the question. Each function
        defines a series of PyMOL commands. Each
        function also has its own documentation that explains
        the molecular scene. and that is printed to the
        command history window when you type 'help q1' or
        'help q2' and so on. Each row of commands starts
        with 'reinitialize' to remove old molecules and settings. Each
        function prints to the command history window
        an alias to the commands in the function. 
        
    LIMITATIONS
        
        Each alias uses the fetch command, so an active
        internet connection is required. 
        
        Replace the fetch command with the load command to 
        load atomic coordinates from your harddrive. 
        
        Enter "rock" to stop the rocking motion.
        
    REQUIREMENTS
            
        When developing a variant of this program., run the 
        program roundview.py to obtained a shorter version 
        of the viewport settings.

    NOTES
    
        Use gray scroll bar to the right side of the command 
        history window to see all of this documentation. or hit 
        escape after moving the cursor to the viewing port 
        (the GUI with the molecule). 

  Copyright Notice
  ================
  
  The PyMOL function source code in this file is copyrighted, but you can
  freely use and copy it as long as you don't change or remove any of
  the copyright notices.
  
    """ 
    print("The function exam2() stores this docstring.")
cmd.extend( "exam2",exam2)


def q1():
    '''
    DESCRIPTION

    Question 1: Explain how AZT terminates extension of the
    DNA chain. 
    
    This is the definition of the alias "q1": 
    delete all; fetch 3v6d, HIVrt, async=0; rock;preset.ball_and_stick("c. P and i. 822");set_view (-0.99,-0.1,0.06,0.09,-0.39,0.92,-0.07,0.92,0.39,0.0,0.0,-29.2,-10.56,24.72,39.27,23.02,35.38,-20.0); 
    '''
    cmd.reinitialize()
    cmd.fetch('3v6d', type='pdb', name='HIVrt',async='0')
    cmd.rock()
    preset.ball_and_stick('c. P and i. 822')
    cmd.set_view('(-0.99,-0.1,0.06,0.09,-0.39,0.92,-0.07,\
    0.92,0.39,0.0,0.0,-29.2,-10.56,24.72,39.27,\
    23.02,35.38,-20.0)') 
    print('Enter "q1" to set scene for question 1.')
    print('Enter "help q1" for more information.') 
cmd.extend( 'q1',q1)



def q2():
    '''
    DESCRIPTION

    Question 2: Look along the helical axis of the DNA. 
    Is the DNA helix bent?
    
    This is the definition of the alias "q2": 
    delete all; fetch 3v6d, HIVrt, async=0; rock;show cartoon, c. P or c. T;set_view (0.61,0.29,-0.73,0.78,-0.37,0.51,-0.12,-0.88,-0.45,-0.0,-0.0,-192.99,-29.84,8.42,47.76,178.27,207.7,-20.0);  
    '''
    cmd.reinitialize()
    cmd.fetch('3v6d',type='pdb',name='HIVrt',async='0')
    cmd.rock()
    cmd.show_as("cartoon","chain P or chain T")
    cmd.set_view('(0.61,0.29,-0.73,0.78,-0.37,0.51,-0.12,-0.88,-0.45,-0.0, 0.0,-192.99,-29.84,8.42,47.76,178.27,207.7,-20.0);')
    cmd.print('Enter "q2" to set scene for question 2.')    
    cmd.print('Enter "help q2" for more information')
cmd.extend( 'q2',q2)


def q3():
    '''
    DESCRIPTION

    Question 3: In which groove of the DNA 
    is the protein making the most contacts? 
    Is this unusual? 
    
    This is the definition of the alias "q3": 
    delete all; fetch 3v6d, HIVrt, async=0; show cartoon, c. A or c. B;hide lines, c. A or c. B;set_view (0.53,-0.06,-0.84,0.82,-0.21,0.53,-0.21,-0.98,-0.07,-0.0,-0.0,-192.99,-29.84,8.42,47.76,178.27,207.7,-20.0);rock;
    '''
    cmd.reinitialize()
    cmd.fetch('3v6d', type='pdb', name= 'HIVrt', async='0')
    cmd.show_as("cartoon","c. A or c. B") 
    cmd.hide("lines","c. A or c. B") 
    cmd.rock()
    cmd.set_view('(0.53,-0.06,-0.84,0.82,-0.21,0.53,-0.21,-0.98,-0.07,-0.0,-0.0,-192.99,-29.84,8.42,47.76,178.27,207.7,-20.0);')
    print('Enter "q3" to set scene for question 3.')
    print('Enter "help q3" for more information')
cmd.extend( 'q3',q3)   



def q4():
    '''
    DESCRIPTION

    Question 4:  What the dihedral angle about the disulfide
    bond between MRG81 of chain F and Cys258 of chain C?
    This is a cross link between the protein and the DNA. Is
    this a cis or trans conformation of the bonds about the
    S--S bond? Is this conformation energetically favorable
    or unfavorable? 
    
    This is the definition of the alias "q4": 
    delete all; fetch 3v6d, HIVrt, async=0;preset.ball_and_stick("(c. C and i. 258) or (c. F and i. 817 )");set_view (0.21,-0.91,0.34,-0.84,0.01,0.54,-0.5,-0.4,-0.77,0.0,-0.0,-38.82,-39.67,-55.1,10.96,36.31,41.34,-20.0); 
    '''
    cmd.reinitialize()
    cmd.fetch('3v6d',type='pdb',name='HIVrt',async='0')
    cmd.do("preset.ball_and_stick(selection='(chain C and resi 258) or (chain F and resi 817 )')")
    cmd.set_view('(0.21,-0.91,0.34,-0.84,0.01,0.54,-0.5,-0.4,-0.77,0.0,-0.0,-38.82,-39.67,-55.1,10.96,36.31,41.34,-20.0);')
    cmd.rock()
    print('Enter "q4" to set scene for question 4.')
    print('Enter "help q4" for more information')
cmd.extend('q4',q4)



def q5():
    '''
    DESCRIPTION

    Question 5:  How many axial stacks of helices does 
    the ribozyme have?
    
    This is the definition of the alias "q5": 
    delete all;fetch 3zp8, hammer, async=0;show cartoon, hammer;set_view (-0.5,0.18,-0.85,-0.17,-0.98,-0.11,-0.85,0.09,0.52,0.0,0.0,-167.2,-18.45,10.92,-12.11,126.37,208.02,-20.0);rock;
    '''
    cmd.reinitialize()
    cmd.fetch('3zp8', type='pdb', name= 'hammer', async='0')
    cmd.show_as('cartoon','hammer')
    cmd.rock()
    cmd.set_view('(-0.5,0.18,-0.85,-0.17,-0.98,-0.11,-0.85,0.09,0.52,0.0,0.0,-167.2,-18.45,10.92,-12.11,126.37,208.02,-20.0);')
    
    print('Enter "q5" to set scene for question 5.')    
    print('Enter "help q5" for more information')
cmd.extend('q5',q5)


def q6():
    '''
    DESCRIPTION

    Question 6:  What is the average distance of the Na1044
   ligand bonds? Give the residue numbers of the RNA
   nucleotides and the sodium to identify them. How many
   ligands are from RNA?
   
    This is the definition of the alias "q6": 
    delete all;fetch 3zp8, hammer, async=0; rock;preset.ball_and_stick("all");distance ligand1, i. 1044, c. A and i. 22 and n. N7;distance ligand2, i. 1044, c. A and i. 21 and n. OP2;distance ligand3, i. 1044, i. 2121;distance ligand4, i. 1044, i. 2120;distance ligand5, i. 1044, i. 2122;distance ligand6, i. 1044, i. 2130;set_view (-0.87,0.18,-0.46,-0.39,-0.81,0.44,-0.29,0.56,0.78,-0.0,0.0,-20.47,-18.05,14.02,-18.89,17.47,23.47,-20.0);
    '''
    cmd.reinitialize()
    cmd.fetch('3zp8', type='pdb', name= 'hammer', async='0')
    cmd.rock()
    cmd.do('preset.ball_and_stick("all")')
    cmd.distance('Sodiumligand1',' resi 1044', 'chain A and i. 22 and n. N7')
    cmd.distance('Sodiumligand2', 'resi 1044', 'chain A and i. 21 and n. OP2')
    cmd.distance('Sodiumligand3', 'resi 1044', 'resi 2121')
    cmd.distance('Sodiumligand4', 'resi 1044', 'resi 2120')
    cmd.distance('Sodiumligand5', 'resi 1044', 'resi 2122')
    cmd.distance('Sodiumligand6', 'resi 1044', 'resi 2130')
    cmd.do('set label_size, -0.4')
    cmd.set_view('(-0.87,0.18,-0.46,-0.39,-0.81,0.44,-0.29,0.56,0.78,-0.0,0.0,-20.47,-18.05,14.02,-18.89,17.47,23.47,-20.0);')   
    print('Enter "q6" to set scene for question 6.')
    print('Enter "help q6" for more information')
cmd.extend('q6',q6)


def q7():
    '''
    DESCRIPTION

    Question 7:  Measure the longest dimension and the
    shortest dimension of the ribozyme. Enter "rock" to 
    stop the rocking motion. What is the ratio of
    the longest dimension to the shortest dimension? Is it
    globular like a protein?  
    
    This is the definition of the alias "q7":
    delete all;fetch 3zp8, hammer, async=0;show ribbon, hammer;set_view (0.62,0.14,0.78,0.13,-0.99,0.07,0.78,0.05,-0.63,-0.0,-0.0,-169.8,-16.43,9.44,-9.63,143.54,196.05,-20.0); rock;
    '''
    cmd.reinitialize()
    cmd.fetch('3zp8', type='pdb', name= 'hammer', async='0')
    cmd.show_as('ribbon', 'hammer')
    cmd.rock()
    cmd.set_view('(0.62,0.14,0.78,0.13,-0.99,0.07,0.78,0.05,-0.63,-0.0,-0.0,-169.8,-16.43,9.44,-9.63,143.54,196.05,-20.0);')
    print('Enter "q7" to set scene for question 7.')
    print('Enter "help q7" for more information')
cmd.extend('q7',q7)


def q8():
    '''
    DESCRIPTION

   Question 8:  Find the unusual base pair between A21 and
   G36. What is the length of the H-bonds between the bases
   (ignore the H atoms in the distance measurement)? List
   the distance with the residue name, residue number, and
   atom name. What additional H-bond occurs between a base
   in this base pair and a ribose ring of one of the two
   nucleotides in this base pair? 
   
    This is the definition of the alias "q8": 
    fetch 3zp8, hammer, async=0;hide everything;show ribbon;show sticks, resi 21 or resi 36;set_view (-0.9,-0.19,0.39,0.39,-0.74,0.55,0.19,0.65,0.74,0.0,0.0,-37.58,-21.66,15.71,-23.32,35.42,39.74,-20.0);rock; 
    '''
    cmd.reinitialize()
    cmd.fetch('3zp8', type='pdb', name= 'hammer', async='0')
    cmd.show_as('ribbon', 'hammer')
    cmd.rock()
    cmd.hide('everything')
    cmd.show_as( 'ribbon','all') 
    cmd.show_as('sticks', 'resi 21 or resi 36')
    cmd.set_view('(-0.9,-0.19,0.39,0.39,-0.74,0.55,0.19,0.65,0.74,0.0,0.0,-37.58,-21.66,15.71,-23.32,35.42,39.74,-20.0);')
    print('Enter "q8" to set scene for question 8.')
    print('Enter "help q8" for more information')
cmd.extend('q8',q8)

