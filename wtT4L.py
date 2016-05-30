from __future__ import print_function
# -*- coding: utf-8 -*-
""" 

  Copyright Notice
  ================
  
     Copyright (C) 2016  Blaine Mooers

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
    See the GNU General Public License for more details:
    http://www.gnu.org/licenses/.

  The source code in this file is copyrighted, but you can
  freely use and copy it as long as you don't change or 
  remove any of the copyright notices.
  
  Blaine Mooers, PhD 
  blaine-mooers@ouhsc.edu
  975 NE 10th St, BRC 466
  University of Oklahoma Health Sciences Center, 
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

print("Enter the name of the alias to the left of the semicolon.") 
print(" ")
print("site11: the active side glutamate Glu11")
print(" ")
print("site31: the His31-Asp 70 salt bridge in which His31's pKa")
print("     is shifted by 3 orders of magnitude")
print(" ")
print("site96: the site of the Arg96-->His mutation that has a large loss in stability")
print("     that was hard to explain from the structure alone.")
print(" ")
print("site99: the Leu99-->Ala mutant protein has a large cavity")
print("     that can accomodate benzene.")
print(" ")
print("site145: poor H-bond geometery between Arg-145 and ")
print("     Glu-11 suggests that the pKa of Glu11 is shifted upwards ")
print("     by the unfavroable interactions.")
print(" ")
print("Type 'help wtT4L' for more information. ")
print(" ")

def wtT4L():
    """
    DESCRIPTION

    Creates aliases to compound commands. The compound commands
    are executed in PyMOL by typing their name. The command 
    wtT4L.py only needs to be run once in a PyMOL session to 
    create the aliases. The aliases disappear upon ending 
    the PyMOL session. The names of the aliases and a short description 
    of each site is printed to the command history window. The additional
    script 'site31lables.pml' must be present in the present working
    directory for alias site31 to work. 
    
    After changing the molecule's orientation, run the program 
    'roundview.py' to obtained a short version of the viewport settings.

    USAGE

    run wtT4L.py
    
    Then enter the name of an alias from the list above, for example: 
    
    site11
    
    """
cmd.extend('wtT4L',wtT4L)
    
print(wtT4L.__doc__)
cmd.bg_color("white")
cmd.fetch(code="3fa0", name="wtT4L",state = 0,async='0')

cmd.alias('site11', 'zoom resi 11; preset.technical("wtT4L")') 
cmd.alias('site31', 'preset.ball_and_stick("wtT4L");@site31labels.pml;set_view (0.31,-0.93,0.21,0.92,0.24,-0.29,0.22,0.28,0.93,-0.09,-0.05,-9.88,37.55,10.06,30.09,20.0,23.82,-20.0);')
cmd.alias('site96', 'preset.technical("wtT4L");set_view (-0.75,-0.65,0.11,0.62,-0.75,-0.22,0.22,-0.1,0.97,0.0,-0.0,-32.32,29.16,-1.45,6.77,27.32,37.32,-20.0);') 
cmd.alias('site99', 'preset.ball_and_stick("wtT4L");color bluewhite, i. 99; set_view (-0.24,-0.95,-0.21,0.51,0.07,-0.85,0.83,-0.32,0.46,0.03,-0.5,-5.37,22.35,-18.6,18.83,29.89,37.68,-20.0);') 
cmd.alias('site145', 'preset.technical("wtT4L");set_view (0.02,-0.63,-0.78,0.37,0.73,-0.56,0.93,-0.27,0.24,0.19,-0.3,-0.0,24.28,2 .24,13.53,15.23,21.26,-20.0);')



