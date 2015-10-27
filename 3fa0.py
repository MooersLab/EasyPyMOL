import pymol 
from pymol import * 
""" 
DESCRIPTION

    Create aliases to compound commands. The compound commands
    are executed in PyMOL by typing their name. The command 
    3fa0sites.py only needs to be run once in a session to 
    create the aliases. The aliases disappear upon ending 
    a session. The names of the aliases and a short description 
    of each site is printed to the command history window. Run 
    the program roundview.py to obtained a short version of 
    the viewport settings.

USAGE
    run 3Fa0sites.py
    Then enter the name of an alias, for example: site11

26 October 2015

  Copyright Notice
  ================
  
  The PyMOL function source code in this file is copyrighted, but you can
  freely use and copy it as long as you don't change or remove any of
  the copyright notices.
  
  ----------------------------------------------------------------------
  This PyMOL script is Copyright (C) 2015 by 
  Blaine Mooers , PhD <blaine-mooers@ouhsc.edu>
  University of Oklahoma Health Sciences Center, Oklahoma City, OK, USA
  
                         All Rights Reserved
  
  Permission to use, copy, modify, distribute, and distribute modified
  versions of this software and its documentation for any purpose and
  without fee is hereby granted, provided that the above copyright
  notice appear in all copies and that both the copyright notice and
  this permission notice appear in supporting documentation, and that
  the name(s) of the author(s) not be used in advertising or publicity
  pertaining to distribution of the software without specific, written
  prior permission.
  
  THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS
  SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS.  IN NO EVENT SHALL THE
  AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
  CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING
  FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
  CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT
  OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
  SOFTWARE.
  ----------------------------------------------------------------------""" 
print("Enter the name of the alias to the left of the semicolon.") 
print(" ")
print("Alias: description of the molecular scene.")
print(" ")
print("site11: the active side glutamate Glu11")
print("site31: the His31-Asp 70 salt bridge in which His31's pKa")
print("     is shifted by 3 orders of magnitude")
print(" ")
print("site96: the site of the Arg96His mutation that was ")
print("     hard to explain from the structure alone and that first")
print("     published in Nature in 1979.")
print(" ")
print("site99: the leu99ala mutant protein was a large cavity")
print("     that can accomodate benezene.")
print(" ")
print("site145: poor H-bond geometery between Arg-145 and ")
print("     Glu-11 suggests that the pKa of Glu11 is shifted upwards ")
print("     by the unfavroable interaction.")
cmd.bg_color("white")
cmd.fetch(code="3fa0", name="wtT4L", state = 0)
cmd.alias('site11', 'zoom resi 11; preset.technical("wtT4L")') 
cmd.alias('site31', 'preset.ball_and_stick("wtT4L");\
set_view (0.31,-0.93,0.21,0.92,0.24,-0.29,0.22,0.28,0.93,\
-0.09,-0.05,-9.88,37.55,10.06,30.09,20.0,23.82,-20.0);')
cmd.alias('site96', 'preset.technical("wtT4L");set_view(-0.75,\
-0.65,0.11,0.62,-0.75,-0.22,0.22,-0.1,0.97,0.0,-0.0,\
-32.32,29.16, -1.45,6.77,27.32,37.32,-20.0);') 
cmd.alias('site99', 'preset.ball_and_stick("wtT4L");\
color bluewhite, i. 99; set_view (-0.24,-0.95,-0.21,\
0.51,0.07,-0.85,0.83,-0.32,0.46,0.03,-0.5,-5.37,22.35,\
-18.6,18.83,29.89,37.68,-20.0);') 
cmd.alias('site145', 'preset.technical("wtT4L");set_view (0.02,-0.63,-0.78,\
0.37,0.73,-0.56,0.93,-0.27,\
0.24,0.19,-0.3,-0.0,24.28,2 .24,13.53,15.23,21.26,-20.0);')

