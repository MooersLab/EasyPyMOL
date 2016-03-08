from __future__ import print_function
# -*- coding: utf-8 -*-
"""
    DESCRIPTION

        Defines on startup aliases to three categories of commands.
        Aliases are listed here instead of in the pymolrc file
        to avoid clutter of the command history window on start up of 
        PyMOL. Source from  your .pymolrc file on the mac or linux or 
        in your pymolrc.pml file on Windows by adding the command:
        
        run ~/Scripts/PyMOLScripts/StartUpAliases.py    
    
        Tested on PyMOL version 1.8.0.5. 
        
        Requires version 1.6 or higher.  Earlier versions of PyMOL 
        should use StartupAliasesSupplementalOldVersions.py.
        
        Requires quat.py from PyMOL Wiki. Store in ~/Scripts/PyMOLScripts/.

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
  freely use and copy it as long as you don't change or remove any of
  the copyright notices.
  
  Blaine Mooers , PhD 
  blaine-mooers@ouhsc.edu
  975 NWE 10th St, BRC 466
  University of Oklahoma Health Sciences Center, Oklahoma City, OK, USA

""" 
from pymol import cmd
__author__ = "Blaine Mooers"
__copyright__ = "Blaine Mooers, University of \
Oklahoma Health Sciences Center, Oklahoma City, OK, USA"
__license__ = "GPL-3"
__version__ = "1.0.1"
__credits__ = ["William Beasley","Chiedza Kanyumbu"] 
# people who reported bug fixes, made suggestions, etc. 
__date__ = "19 February 2016"
__maintainer__ = "Blaine Mooers"
__email__ = "blaine-mooers@ouhsc.edu"
__status__ = "Production" 

def SA():
    '''
DESCRIPTION
    StartUpAliases.py  Copyright (C) 2016  Blaine Mooers. 
    This script comes with 
    ABSOLUTELY NO WARRANTY; for details, 
    see source file. 

    Format of list below:
    Active alias, description: PDB code, where applicable. 
    
    Molecules in standard orientations: 
    T4L, WT T4 lysozyme (1.09 ang) as a ribbon diagram: 3fa0. 
    U8, 16-mer dsRNA with 8 contiguous Us. U-helix RaNA (1.37 ang): 3nd3.
    WC8, 16-mer RNA with all Watson-Crick base pairs (1.67 ang): 3nd4.
    N9, neuraminidase as cartoon, biological unit (1.55 ang): 4dgr.
    GGT, gamma glutamyl transpeptidase as cartoon (1.67 ang): 4gdx.
    GU, 10-mer RNA with eight GU base pairs (1.32 ang): 4pco.

    Complex figures to serve as templates: 
    BST, Base-stacking figure, (1.32 ang): 4pco. 
    LG, Electron density map of nine sugar glycan,(1.55 ang):, 4dgr. 
    NA, Sodium cation in major groove of 16-mer RNA: 3nd4.


    Complex representations applied to any visible molecular object:
    
    AO, Make ambient occlusion image. Requires global view of protein.
    BS, Make fancy ball and stick representation of visible atoms. 
    BU, Display biological unit. 
    CB, Define color blind compatible coloring scheme. 
    BW, Make black and white ribbon cartoon on white background.
    CSS, Color ribbon and cartoons by secondary structure: red, green and yellow. 
    CBSS, Color ribbon and cartoons with dcolorblind friendly colors. 
    CR, Commands to make colored filled-ring cartoon of nucleic acids..
    FR, Commands to make filled-ring cartoon of nucleic acids.


    Type the alias name to execute the commands. The names are case
    sensitive. Type 'help <AliasName>' (e.g., help T4L) for
    description and the commands. Some aliases require additional
    scripts. The commands can be copied from the command history
    window and pasted onto the command line for code reuse. The first
    set of commands has line breaks for easy selection of code
    fragment selection and copying. The second set of commands is one
    one line for easy copying and pasting of the entire horizontal
    script. 
    
    Type 'SA' to refresh the list of aliases.
    '''
    print(SA.__doc__)
cmd.extend('SA',SA)


def T4L():
    '''
DESCRIPTION
    
    WT T4 lysozyme as ribbon diagram (1.08 Ang):  3Fa0. 
    
USAGE
    Type 'T4L' to activate. Type 'help T4L' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    delete all;fetch 3fa0,type=pdb,async=0;orient;
    turn z,-90;turn y,-5;turn x,10; 
    hide everything; bg_color white; 
    show cartoon;color red, ss H;color yellow, ss S;
    color green, ss L+;set_view (-0.18,-0.69,-0.7,0.98,-0.17,-0.09,
    -0.06,-0.7,0.71,0.0,0.0,-165.67,34.77,11.27,9.52,132.07,
    199.27,-20.0); ray 1500,1600;
    
    Commands without linebreaks:
    delete all;fetch 3fa0,type=pdb,async=0;orient;turn z,-90;turn y,-5;turn x,10; hide everything; bg_color white;show cartoon;color red, ss H;color yellow, ss S;color green, ss L+;set_view (-0.18,-0.69,-0.7,0.98,-0.17,-0.09,-0.06,-0.7,0.71,0.0,0.0,-165.67,34.77,11.27,9.52,132.07,199.27,-20.0); ray 1500,1600; 
    
    '''
    cmd.reinitialize()
    cmd.fetch('3fa0', type='pdb', async='0')
    cmd.orient()
    cmd.turn('z', '-90')
    cmd.turn('y', '-5')
    cmd.turn('x', '10')
    cmd.hide('everything')
    cmd.bg_color('white')
    cmd.show('cartoon')
    cmd.color('red', 'ss H')
    cmd.color('yellow', 'ss S')
    cmd.color('green', 'ss L+')
    cmd.set_view('(-0.18,-0.69,-0.7,0.98,-0.17,-0.09,-0.06,-0.7,0.71,0.0,0.0,-165.67,34.77,11.27,9.52,132.07,199.27,-20.0)')
    cmd.ray('1500', '1600')
cmd.extend('T4L',T4L)
    
    
def U8():
    '''
DESCRIPTION

    16-mer dsRNA with 8 contiguous Us. U-helix RNA (1.37 Ang):  3nd3.
    Has one strand in the asymmetric unit. Uses quat.py to generate
    the second strand. Cartoon with filled rings and bases cartoon.
    
USAGE
    Type 'U8' to activate. Type 'help U8' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    delete all;fetch 3nd3,type=pdb,async=0;
    run /Users/blaine-mooers/Scripts/PyMOLScripts/quat.py;quat 3nd3;
    hide everything;bg_color white; show sticks;set cartoon_ring_mode, 3;
    set cartoon_ring_finder, 1;set cartoon_ladder_mode, 1;
    set cartoon_nucleic_acid_mode, 4;
    set cartoon_ring_transparency, 0.5;as cartoon;
    set_view (-1.0,-0.03,0.06,-0.06,0.01,-1.0,0.04,-1.0,-0.01,-0.09,
    -0.02,-168.02,7.85,15.56,-0.21,137.38,199.33,-20.0);draw; 
    
    Commands without linebreaks:
    delete all;fetch 3nd3,type=pdb,async=0;run /Users/blaine-mooers/Scripts/PyMOLScripts/quat.py;quat 3nd3;hide everything;bg_color white; show sticks;set cartoon_ring_mode, 3;set cartoon_ring_finder, 1;set cartoon_ladder_mode, 1;set cartoon_nucleic_acid_mode, 4;set cartoon_ring_transparency, 0.5;as cartoon;set_view (-1.0,-0.03,0.06,-0.06,0.01,-1.0,0.04,-1.0,-0.01,-0.09,-0.02,-168.02,7.85,15.56,-0.21,137.38,199.33,-20.0);draw; 
    '''
    
    cmd.reinitialize()
    cmd.fetch('3nd3', type='pdb', async='0')
    cmd.do('run /Users/blaine-mooers/Scripts/PyMOLScripts/quat.py')
    cmd.do('quat 3nd3')
    cmd.hide('everything')
    cmd.bg_color('white')
    cmd.show('sticks')
    cmd.set('cartoon_ring_mode', '3')
    cmd.set('cartoon_ring_finder', '1')
    cmd.set('cartoon_ladder_mode', '1')
    cmd.set('cartoon_nucleic_acid_mode', '4')
    cmd.set('cartoon_ring_transparency', '0.5')
    cmd.show_as('cartoon')
    cmd.set_view('(-1.0,-0.03,0.06,-0.06,0.01,-1.0,0.04,-1.0,-0.01,-0.09,-0.02,-168.02,7.85,15.56,-0.21,137.38,199.33,-20.0)')
    cmd.draw()
cmd.extend('U8',U8)
    
    
def WC8():
    '''
DESCRIPTION

    16-mer dsRNA. Watson-Crick helix RNA. 1.55 Angstrom 
    resolution: 3nd4.  Has one strand in the asymmetric unit. 
    Needs quat.py to generate the second strand. Use the 
    BU alias. Cartoon with filled rings and bases cartoon.
    
USAGE
    Type 'WC8' to activate. Type 'help WC8' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    delete all; fetch 3nd4,type=pdb,async=0;hide everything;
    run /Users/blaine-mooers/Scripts/PyMOLScripts/quat.py;
    quat 3nd4;bg_color white; show sticks; set stick_radius, 0.12; 
    set nb_spheres_size, 0.25;show nb_spheres;set stick_ball, on; 
    set stick_ball_ratio, 1.8;
    set_view (-0.99,-0.03,0.17,-0.18,0.02,-0.98,0.03,-1.0,-0.03,
    0.0,0.0,-169.97,8.1,15.62,-1.69,139.24,200.7,-20.0);
    hide everything,name H*;draw 

    Commands without linebreaks:
    delete all; fetch 3nd4,type=pdb,async=0;hide everything; run /Users/blaine-mooers/Scripts/PyMOLScripts/quat.py; quat 3nd4;bg_color white; show sticks; set stick_radius, 0.12; set nb_spheres_size, 0.25; show nb_spheres; set stick_ball, on; set stick_ball_ratio, 1.8;set_view (-0.99,-0.03,0.17,-0.18,0.02,-0.98,0.03,-1.0,-0.03,0.0,0.0,-169.97,8.1,15.62,-1.69,139.24,200.7,-20.0);hide everything, name H*;draw 

    '''
    cmd.reinitialize()
    cmd.fetch('3nd4', type='pdb', async='0')
    cmd.remove('name H*')
    cmd.hide('everything')
    cmd.do('run /Users/blaine-mooers/Scripts/PyMOLScripts/quat.py')
    cmd.do('quat 3nd4')
    cmd.bg_color('white')
    cmd.do("preset.ball_and_stick(selection='all')")
    cmd.set_view('(-0.96,-0.03,0.3,-0.31,0.02,-0.95,0.03,-1.0,-0.03,0.0,0.0,-231.24,8.16,15.68,-1.66,200.47,262.01,-20.0)')
    cmd.rock()
cmd.extend('WC8',WC8)
    
    
def N9():
    '''
DESCRIPTION
    
    Influenza N9 neuraminidase at 1.55 Angstrom resolution, PDB code
    4dgr. The biological unit has four copies of the asymmetric unit.
    View is down the four-fold axis. Requires the quat.py script by
    Thomas Holder and available at the PyMOL Wiki page. Store quat.py
    in ~/Scripts/PyMOLScripts

USAGE
    Type 'N9' to activate. Type 'help N9' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.

    delete all;fetch 4dgr, type=pdb, async=0;
    run /Users/blaine-mooers/Scripts/PyMOLScripts/quat.py;
    quat 4dgr;as cartoon; bg_color white;
    color red, 4dgr_1 and ss H;color yellow,4dgr_1 and ss S;
    color green, 4dgr_1 and ss L+;
    color cyan, (not 4dgr_1 and ss H);
    color magenta, (not 4dgr_1 and ss S);
    color orange, (not 4dgr_1 and ss L+);
    set_view (0.98,-0.22,0.01,0.22,0.98,0.02,-0.01,-0.02,
    1.0,-0.0,0.0,-323.44,1.46,5.33,56.19,274.72,372.15,-20.0);
    draw 
    
    delete all;fetch 4dgr, type=pdb, async=0;run /Users/blaine-mooers/Scripts/PyMOLScripts/quat.py; quat 4dgr;as cartoon; bg_color white;color red, 4dgr_1 and ss H;color yellow,4dgr_1 and ss S;color green, 4dgr_1 and ss L+;color cyan, (not 4dgr_1 and ss H);color magenta, (not 4dgr_1 and ss S);color orange, (not 4dgr_1 and ss L+);set_view (0.98,-0.22,0.01,0.22,0.98,0.02,-0.01,-0.02,1.0,-0.0,0.0,-323.44,1.46,5.33,56.19,274.72,372.15,-20.0); draw 
    '''
    cmd.reinitialize()
    cmd.fetch('4dgr', type='pdb', async='0')
    cmd.do('run /Users/blaine-mooers/Scripts/PyMOLScripts/quat.py')
    cmd.do('quat 3nd4')
    cmd.show_as('cartoon')
    cmd.bg_color('white')
    cmd.color('red', '4dgr_1 and ss H')
    cmd.color('yellow', '4dgr_1 and ss S')
    cmd.color('green', '4dgr_1 and ss L+')
    cmd.color('cyan', '(not 4dgr_1 and ss H)')
    cmd.color('magenta', '(not 4dgr_1 and ss S)')
    cmd.color('orange', '(not 4dgr_1 and ss L+)')
    cmd.set_view('(0.98,-0.22,0.01,0.22,0.98,0.02,-0.01,-0.02,1.0,-0.0,0.0,-323.44,1.46,5.33,56.19,274.72,372.15,-20.0)')
    cmd.draw()
cmd.extend('N9',N9)
    
    
def GGT():
    '''
    DESCRIPTION

    WT human gamma glutamyl transpeptidase at 1.67  Angstrom
    resolution as cartoon. PDB Code 4gdx.
    
USAGE
    Type 'GGT' to activate. Type 'help GGT' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    delete all;fetch 4gdx, type=pdb, async=0;remove name H*;
    as cartoon;bg_color white; hide (name c+o+n);
    set cartoon_side_chain_helper, on;color red, 4gdx and ss H; 
    color yellow,4gdx and ss S;color green,4gdx and ss L+; 
    select ASNNAG, resn NAG or resi 95 or i. 120  or i. 230 or
    i. 266 or i. 344 or i. 511 or i. 381; 
    color red, elem o and ASNNAG; 
    color blue, elem n and ASNNAG;
    color yellow, elem c and ASNNAG;show sticks,ASNNAG;
    disable ASNNAG;
    set_view (0.55,-0.83,0.07,0.5,0.26,-0.82,0.66,0.49,0.56,0.0,0.0,-197.16,-22.42,-22.69,-12.01,155.44,238.88,-20.0); 
    draw 
    
    Commands without linebreaks:
    delete all;fetch 4gdx, type=pdb, async=0;remove  name H*;as cartoon;bg_color white; hide (name c+o+n);set cartoon_side_chain_helper,  on;color red, 4gdx and ss H; color yellow,4gdx and ss S;color green,4gdx and ss L+; select ASNNAG,resn NAG or resi 95 or i. 120  or i. 230 or i. 266 or i. 344 ori. 511 or i. 381; color red, elem o and ASNNAG; color blue, elem n and ASNNAG;color yellow, elem c  and ASNNAG;show sticks,ASNNAG;disable ASNNAG; set_view(0.55,-0.83,0.07,0.5,0.26,-0.82,0.66,0.49,0.56,0.0,0.0,-197.16,-22.42,-22.69,-12.01,155.44,238.88,-20.0); draw 
    
    '''
    cmd.reinitialize()
    cmd.fetch('4gdx', type='pdb', async='0')
    cmd.remove('name H*')
    cmd.show_as('cartoon')
    cmd.bg_color('white')
    cmd.hide('(name c+o+n)')
    cmd.set('cartoon_side_chain_helper', 'on')
    cmd.color('red', '4gdx and ss H')
    cmd.color('yellow', '4gdx and ss S')
    cmd.color('green', '4gdx and ss L+')
    cmd.select('ASNNAG', 'resn NAG or resi 95 or i. 120  or i. 230 or i. 266 or i. 344 or i. 511 or i. 381')
    cmd.color('red', 'elem o and ASNNAG')
    cmd.color('blue', 'elem n and ASNNAG')
    cmd.color('yellow', 'elem c and ASNNAG')
    cmd.show('sticks', 'ASNNAG')
    cmd.disable('ASNNAG')
    cmd.set_view('(0.55,-0.83,0.07,0.5,0.26,-0.82,0.66,0.49,0.56,0.0,0.0,-197.16,-22.42,-22.69,-12.01,155.44,238.88,-20.0)')
    cmd.draw()
cmd.extend('GGT',GGT)
    
    
def GU():
    '''
    DESCRIPTION

    10-mer dsRNA with 8 contiguous Us. U-helix RNA. 
    1.32 Angstrom resolution: 4PCO. Has five strands in 
    the asymmetric unit. Deleted chain E and cobalt 
    hexammine 102. Cartoon with filled rings and
    bases cartoon.
    
USAGE
    Type 'GU' to activate. Type 'help GU' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    delete all;fetch 4PCO,type=pdb,async=0;hide everything; 
    bg_color white; cartoon oval; set cartoon_ring_mode, 3; 
    set cartoon_nucleic_acid_color, blue;select rna_A, resn A; 
    select rna_C, resn C;select rna_G, resn G; 
    select rna_U, resn U;color yellow, rna_A; color red, rna_C;
    color gray40, rna_G; color palecyan, rna_U;as cartoon; 
    disable rna_U; set stick_radius, 0.12;
    set nb_spheres_size, 0.3; show nb_spheres; set stick_ball, on;
    set stick_ball_ratio, 1.8; show sticks, resn NCO; 
    show spheres, name Cl; set_view (0.34,-0.81,0.48,0.89,0.11,
    -0.45,0.31,0.58,0.76,-0.0,0.0,-196.36,-9.82,6.76,15.84,159.01,
    233.71,-20.0);draw') 
    
    Commands without linebreaks: 
    delete all;fetch 4PCO,type=pdb,async=0;hide everything;bg_color white; cartoon oval; set cartoon_ring_mode, 3;set cartoon_nucleic_acid_color, blue;select rna_A, resn A;select rna_C,resn C;select rna_G, resn G;select rna_U, resn U;color yellow, rna_A; color red, rna_C;color gray40, rna_G; color palecyan, rna_U;as cartoon;disable rna_U; set stick_radius, 0.12;set nb_spheres_size, 0.3; show nb_spheres; set stick_ball, on;set stick_ball_ratio, 1.8; show sticks, resn NCO;show spheres, name Cl; set_view (0.34,-0.81,0.48,0.89,0.11,-0.45,0.31,0.58,0.76,-0.0,0.0,-196.36,-9.82,6.76,15.84,159.01,233.71,-20.0);draw') 
    '''
    cmd.reinitialize();
    cmd.fetch('4PCO', type='pdb', async='0')
    cmd.hide('everything')
    cmd.bg_color('white')
    cmd.cartoon('oval')
    cmd.set('cartoon_ring_mode', '3')
    cmd.set('cartoon_nucleic_acid_color', 'blue')
    cmd.select('rna_A', 'resn A')
    cmd.select('rna_C', 'resn C')
    cmd.select('rna_G', 'resn G')
    cmd.select('rna_U', 'resn U')
    cmd.color('yellow', 'rna_A')
    cmd.color('red', 'rna_C')
    cmd.color('gray40', 'rna_G')
    cmd.color('palecyan', 'rna_U')
    cmd.show_as('cartoon')
    cmd.disable('rna_U')
    cmd.set('stick_radius', '0.12')
    cmd.set('nb_spheres_size', '0.3')
    cmd.show('nb_spheres')
    cmd.set('stick_ball', 'on')
    cmd.set('stick_ball_ratio', '1.8')
    cmd.show('sticks', 'resn NCO')
    cmd.show('spheres', 'name Cl')
    cmd.set_view('(0.34,-0.81, 0.48,0.89,0.11,-0.45,0.31,0.58,0.76,-0.0,0.0,-196.36,-9.82,6.76,15.84,159.01,233.71,-20.0)')
    cmd.draw()
cmd.extend('GU',GU)
    
    
#######Commands to display complex scenes. #############
    
def BST():
    '''
    DESCRIPTION
    
    G2G3/U9U8 base step , PDB code 4PCO. 
    From the 1.32 Angstrom resolution structure 
    of the RNA decamer with 8 GU base pairs.
    
USAGE
    Type 'BST' to activate. Type 'help BST' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    

    fetch 4PCO, type=pdb, async=0;
    select G2G3, ( ((resi 2 or resi 3) and chain A) or ((resi 8 or resi 9) and chain B));
    remove not G2G3;
    bg_color white;show sticks;set stick_radius=0.14;
    set stick_ball, on; set stick_ball_ratio,1.9;
    set_view 
    (-0.75,0.09,0.66,-0.2,0.92,-0.35,-0.64,-0.39,-0.67,-0.0,-0.0,-43.7
    ,7. 24,9.55,11.78,29.46,57.91,-20.0);
    hide everything, element H;
    select carbon1, element C and (resi 3 or resi 8) 
    # select lower 
    base pair;select carbon2, element C and (resi 2 or resi 9) 
    #select upper base pair;
    color gray70, carbon1;
    color gray10, carbon2;
    space cmyk;
    distance hbond1, /4PCO//B/U`9/N3,/4PCO//A/G`2/O6;
    distance hbond2, /4PCO//B/U`9/O2,/4PCO//A/G`2/N1;
    distance hbond3, /4PCO//A/U`3/N3,/4PCO//B/G`8/O6;
    distance hbond4, /4PCO//A/U`3/O2,/4PCO//B/G`8/N1;
    color black, hbond1;
    color black, hbond2;
    color gray70, hbond3;
    color gray70, hbond4;
    show nb_spheres;
    set nb_spheres_size, 0.35;
    hide labels;
    ray 1600,1000;png 4PCO.png
    
    Commands without linebreaks: 
    delete all;fetch 4PCO, type=pdb,async=0;select G2G3, ( ((resi 2 or resi 3) and chain A) or ((resi 8 or resi 9) and chain B) );remove not G2G3;bg_color white;show sticks;set stick_radius=0.14;set stick_ball, on; set stick_ball_ratio, 1.9;set_view (-0.75,0.09,0.66,-0.2,0.92,-0.35,-0.64,-0.39,-0.67,-0.0,-0.0,-43.7,7.24,9.55,11.78,29.46,57.91,-20.0);hide everything, element h; select carbon1,element c and (resi 3 or resi 8)# select lower base pair;select carbon2, element c and (resi 2 or resi 9)#select upper base pair;delete all;fetch 4PCO, type=pdb, async=0;\select G2G3, ( ((resi 2 or resi 3) and chain A) or ((resi 8 or resi 9) and chain B));remove not G2G3;bg_color white;show sticks;set stick_radius=0.14;set stick_ball, on; set stick_ball_ratio,1.9;set_view (-0.75,0.09,0.66,-0.2,0.92,-0.35,-0.64,-0.39,-0.67,-0.0,-0.0,-43.7,7. 24,9.55,11.78,29.46,57.91,-20.0);hide everything, element h;select carbon1, element c and (resi 3 or resi 8) # select lower base pair;select carbon2, element c and (resi 2 or resi 9) #select upper base pair;color gray70, carbon1;color gray00,carbon2;space cmyk;distance hbond1,/4PCO//B/U`9/N3,/4PCO//A/G`2/O6;distance hbond2, /4PCO//B/U`9/O2,/4PCO//A/G`2/N1;distance hbond3, /4PCO//A/U`3/N3,/4PCO//B/G`8/O6;distance hbond4, /4PCO//A/U`3/O2,/4PCO//B/G`8/N1;color black, hbond1;color black, hbond2;color gray70, hbond3;color gray70, hbond4;show nb_spheres;set nb_spheres_size, 0.35;hide labels;ray 1600,1000;png 4PCO.png 
    '''
    cmd.reinitialize()
    cmd.fetch('4PCO', type='pdb', async='0')
    cmd.select('G2G3', '( ((resi 2 or resi 3) and chain A)or ((resi 8 or resi 9) and chain B) )')
    cmd.remove('not G2G3')
    cmd.bg_color('white')
    cmd.set('stick_radius', '0.14')
    cmd.set('stick_ball', 'on')
    cmd.set('stick_ball_ratio', '1.9')
    cmd.set_view('(-0.75,0.09,0.66,-0.2,0.92,-0.35,-0.64,-0.39,-0.67,-0.0,-0.0,-43.7,7.24,9.55,11.78,29.46,57.91,-20.0)')
    cmd.remove('name H*')
    cmd.select('carbon1', 'element C and (resi 3 or resi 8)')
    cmd.select('carbon2', 'element C and (resi 2 or resi 9)')
    cmd.color('gray70', 'carbon1')
    cmd.color('gray10', 'carbon2')
    cmd.show('sticks')
    cmd.space('cmyk')
    cmd.distance('hbond1', '/4PCO//B/U`9/N3', '/4PCO//A/G`2/O6')
    cmd.distance('hbond2', '/4PCO//B/U`9/O2', '/4PCO//A/G`2/N1')
    cmd.distance('hbond3', '/4PCO//A/U`3/N3', '/4PCO//B/G`8/O6')
    cmd.distance('hbond4', '/4PCO//A/U`3/O2', '/4PCO//B/G`8/N1')
    cmd.color('black', 'hbond1')
    cmd.color('black', 'hbond2')
    cmd.color('gray70', 'hbond3')
    cmd.color('gray70', 'hbond4')
    cmd.show('nb_spheres')
    cmd.set('nb_spheres_size', '0.35')
    cmd.hide('labels')
    cmd.ray('1600', '1000')
    cmd.png('4PCO.png')
cmd.extend('BST',BST)


def LG():
    '''
    DESCRIPTION
    
    Nine sugar glycan in influenza N9 neuraminidase at 
    1.55 Angstrom  resolution, PDB code 4dgr. 
    The electron density map is contoured at 1.0 sigma. 
    39 commands were used to make this figure.  
    
USAGE
    Type 'LG' to activate. Type 'help LG' to see this documentation
   printed to the command history window. Select from the command
   history individual lines of code to build a new script. Select the
   hortizontal script at the bottom if retaining most of the commands
   in your new script. Copy and paste onto the comand line below.
   Works only with the command line immediately under the command
   history window at the top of the gui.
    
    delete all;fetch 4dgr, async=0;fetch 4dgr, type=2fofc, 
    async=0;select LongGlycan, resi 469:477;
    orient LongGlycan;remove not LongGlycan;
    remove name H*;isomesh 2fofcmap, 4dgr_2fofc, 1, 
    LongGlycan, carve = 1.8;color density, 2fofcmap; 
    show sticks;show spheres;set stick_radius, .07;
    set sphere_scale, .19;set sphere_scale, .13, elem H;
    set bg_rgb=[1, 1, 1];set stick_quality, 50;
    set sphere_quality, 4;color gray85, elem C;
    color red, elem O;color slate, elem N;
    color gray98, elem H;set stick_color, gray50;
    set ray_trace_mode, 1;set ray_texture, 2;
    set antialias, 3;set ambient, 0.5;set spec_count, 5;
    set shininess, 50;set specular, 1;set reflect, .1;
    set dash_gap, 0;set dash_color, black;
    set dash_gap, .15;set dash_length, .05;
    set dash_round_ends, 0;set dash_radius, .05;
    set_view (0.34,-0.72,0.61,0.8,0.56,0.22,-0.51,0.4,
    0.77,0.0,0.0,-81.31,44.64,-9.02,58.62,65.34,97.28,-20.0);
    preset.ball_and_stick("all",mode=1);draw 
    
    Commands without linebreaks:
    delete all;fetch 4dgr, async=0;fetch 4dgr, type=2fofc, async=0;select LongGlycan, resi 469:477;orient LongGlycan;remove not LongGlycan;remove name H*;isomesh 2fofcmap, 4dgr_2fofc, 1, LongGlycan, carve = 1.8;color density, 2fofcmap; show sticks;show spheres;set stick_radius, .07;set sphere_scale, .19;set sphere_scale, .13, elem H;set bg_rgb=[1, 1, 1];set stick_quality, 50;set sphere_quality, 4;color gray85, elem C;color red, elem O;color slate, elem N;color gray98, elem H;set stick_color, gray50;set ray_trace_mode, 1;set ray_texture, 2;set antialias, 3;set ambient, 0.5;set spec_count, 5;set shininess, 50;set specular, 1;set reflect, .1;set dash_gap, 0;set dash_color, black;set dash_gap, .15;set dash_length, .05;set dash_round_ends, 0;set dash_radius, .05;set_view (0.34,-0.72,0.61,0.8,0.56,0.22,-0.51,0.4,0.77,0.0,0.0,-81.31,44.64,-9.02,58.62,65.34,97.28,-20.0);preset.ball_and_stick("all",mode=1);draw 
    '''
    cmd.reinitialize()
    cmd.fetch('4dgr', async='0')
    cmd.fetch('4dgr', type='2fofc', async='0')
    cmd.select('LongGlycan', 'resi 469:477')
    cmd.orient('LongGlycan')
    cmd.remove('not LongGlycan')
    cmd.remove('name H*')
    cmd.isomesh('2fofcmap', '4dgr_2fofc', '1', 'LongGlycan', carve ='1.8')
    cmd.color('density', '2fofcmap')
    cmd.show('sticks')
    cmd.show('spheres')
    cmd.set('stick_radius', '.07')
    cmd.set('sphere_scale', '.19')
    cmd.set('sphere_scale', '.13', 'elem H')
    cmd.set('bg_rgb', '[1, 1, 1]')
    cmd.set('stick_quality', '50')
    cmd.set('sphere_quality', '4')
    cmd.color('gray85', 'elem C')
    cmd.color('red', 'elem O')
    cmd.color('slate', 'elem N')
    cmd.color('gray98', 'elem H')
    cmd.set('stick_color', 'gray50')
    cmd.set('ray_trace_mode', '1')
    cmd.set('ray_texture', '2')
    cmd.set('antialias', '3')
    cmd.set('ambient', '0.5')
    cmd.set('spec_count', '5')
    cmd.set('shininess', '50')
    cmd.set('specular', '1')
    cmd.set('reflect', '.1')
    cmd.set('dash_gap', '0')
    cmd.set('dash_color', 'black')
    cmd.set('dash_gap', '.15')
    cmd.set('dash_length', '.05')
    cmd.set('dash_round_ends', '0')
    cmd.set('dash_radius', '.05')
    cmd.set_view('(0.34,-0.72,0.61,0.8,0.56,0.22,-0.51,0.4,0.77,0.0,0.0,-81.31,44.64,-9.02,58.62,65.34,97.28,-20.0)')
    preset.ball_and_stick("all",mode=1);
    cmd.draw()
cmd.extend('LG',LG)
    
    
def NA():
    '''
DESCRIPTION

    Hydrated sodium cation bound in major groove of a 
    16-mer RNA of Watson-Crick base pairs.
    The sodium is bound to the N7 nitrogen atom of 
    Adenine 3 at 1.55 Angstrom resolution, PDB code 3nd4. 
    57 commands were used to make this figure. 

    More than one label in a horizontal script is not 
    allowed. This one label has to be at the end of the line.
    Labels can be imported from a Labels.pml file.
    Store the label commands one per row in this file.
    Import the file with the @Labels.pml command. 
    Include the path to the file if the labels file is not 
    in the current working directory of PyMOL. 
    
USAGE
    Type 'NA' to activate. Type 'help NA' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    delete all;
    viewport 900,600;
    fetch 3nd4, type=pdb, async=0;
    run ~/Scripts/PyMOLScripts/quat.py; 
    quat 3nd4; 
    show sticks;
    set stick_radius=0.125;
    hide everything, name H*;
    bg_color white;
    create coorCov, (3nd4_1 and (resi 19 or resi 119 or resi 219 or resi 319 or resi 419 or resi 519 or (resi 3 and name N7)));
    bond (coorCov//A/NA`19/NA),(coorCov//A/A`3/N7);
    bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`119/O); 
    bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`219/O); 
    bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`319/O);
    bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`419/O); 
    bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`519/O);
    distance (3nd4_1 and chain A,and resi 19 and name NA), (3nd4_1 and chain A and resi 519);
    distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 419);
    distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 119);
    distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 319);
    distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 219);
    show nb_spheres; 
    set nb_spheres_size, .35;
    distance hbond1,/3nd4_1/1/A/HOH`119/O, /3nd4_1/1/A/A`3/OP2;
    distance hbond2,/3nd4_1/1/A/HOH`319/O, /3nd4_1/1/A/A`3/OP2;
    distance hbond3,/3nd4_1/1/A/HOH`91/O,/3nd4_1/1/A/HOH`119/O;
    distance hbond4,/3nd4_1/1/A/G`4/N7,/3nd4_1/1/A/HOH`91/O;
    distance hbond5,/3nd4_1/1/A/G`4/O6, /3nd4_1/1/A/HOH`419/O;
    distance hbond6,/3nd4_1/1/A/HOH`91/O,/3nd4_1/1/A/G`4/OP2;
    distance hbond7,/3nd4_1/1/A/HOH`319/O,/3nd4_1/1/A/G`2/OP2;
    distance  hbond9,/3nd4_1/1/A/HOH`419/O,/3nd4_2/2/A/HOH`74/O;
    distance hbond10,/3nd4_2/2/A/C`15/O2,/3nd4_1/1/A/G`2/N2;
    distance hbond11, /3nd4_2/2/A/C`15/N3,/3nd4_1/1/A/G`2/N1;
    distance hbond12,/3nd4_2/2/A/C`15/N4,/3nd4_1/1/A/G`2/O6;
    distance hbond13, /3nd4_2/2/A/U`14/N3,/3nd4_1/1/A/A`3/N1;
    distance hbond14,3nd4_2/2/A/U`14/O4,/3nd4_1/1/A/A`3/N6;
    distance hbond15, /3nd4_2/2/A/C`13/N4,/3nd4_1/1/A/G`4/O6;
    distance hbond16,/3nd4_2/2/A/C`13/N3, /3nd4_1/1/A/G`4/N1;
    distance hbond17, /3nd4_1/1/A/G`4/N2,/3nd4_2/2/A/C`13/O2;
    distance hbond18,/3nd4_1/1/A/G`2/N2,/3nd4_2/2/A/C`15/O2;
    distance hbond19,/3nd4_1/1/A/HOH`91/O,/3nd4_1/1/A/G`4/OP2;
    set depth_cue=0;set ray_trace_fog=0;set dash_color, black;
    set label_font_id, 5;set label_size, 36;
    set label_position, (0.5, 1.0, 2.0);
    set label_color, black;set dash_gap, 0.2;
    set dash_width, 2.0;set dash_length, 0.2;
    set label_color, black;set dash_gap, 0.2;
    set dash_width, 2.0;set dash_length, 0.2;
    select carbon, element C; color yellow, carbon;
    disable carbon;set_view
    (-0.9,0.34,-0.26,0.33,0.18,-0.93,-0.27,-0.92,-0.28,-0.07,-0.23,-27.83,8.63,19.85,13.2,16.0,31.63,-20.0)
    
    Commands without linebreaks:
    delete all;viewport 900,600;fetch 3nd4, type=pdb,async=0;run ~/Scripts/PyMOLScripts/quat.py;quat 3nd4; show sticks;set stick_radius=0.125;hide everything, name H*;bg_color white;create coorCov, (3nd4_1 and (resi 19 or resi 119 or resi 219 or resi 319 or resi 419 or resi 519 or (resi 3 and name N7)));bond (coorCov//A/NA`19/NA),(coorCov//A/A`3/N7); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`119/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`219/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`319/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`419/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`519/O);distance (3nd4_1 and chain Aand resi 19 and name NA), (3nd4_1 and chain A and resi 519);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 419);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 119);distance (3nd4_1 and chain A and resi 19 and name NA),(3nd4_1 and chain A and resi 319);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 219);show nb_spheres; set nb_spheres_size, .35;distance hbond1,/3nd4_1/1/A/HOH`119/O, /3nd4_1/1/A/A`3/OP2;distance hbond2,/3nd4_1/1/A/HOH`319/O, /3nd4_1/1/A/A`3/OP2;distance hbond3,/3nd4_1/1/A/HOH`91/O, /3nd4_1/1/A/HOH`119/O;distance hbond4,/3nd4_1/1/A/G`4/N7,/3nd4_1/1/A/HOH`91/O;distance hbond5,/3nd4_1/1/A/G`4/O6, /3nd4_1/1/A/HOH`419/O;distance hbond6,/3nd4_1/1/A/HOH`91/O, /3nd4_1/1/A/G`4/OP2;distance hbond7,/3nd4_1/1/A/HOH`319/O, /3nd4_1/1/A/G`2/OP2;distance  hbond9,/3nd4_1/1/A/HOH`419/O,/3nd4_2/2/A/HOH`74/O;distance hbond10,/3nd4_2/2/A/C`15/O2,/3nd4_1/1/A/G`2/N2;distance hbond11, /3nd4_2/2/A/C`15/N3,/3nd4_1/1/A/G`2/N1;distance hbond12,/3nd4_2/2/A/C`15/N4,/3nd4_1/1/A/G`2/O6;distance hbond13, /3nd4_2/2/A/U`14/N3,/3nd4_1/1/A/A`3/N1;distance hbond14,3nd4_2/2/A/U`14/O4,/3nd4_1/1/A/A`3/N6;distance hbond15, /3nd4_2/2/A/C`13/N4,/3nd4_1/1/A/G`4/O6;distance hbond16,/3nd4_2/2/A/C`13/N3, /3nd4_1/1/A/G`4/N1;distance hbond17, /3nd4_1/1/A/G`4/N2,/3nd4_2/2/A/C`13/O2;distance hbond18,/3nd4_1/1/A/G`2/N2,/3nd4_2/2/A/C`15/O2;distance hbond19,/3nd4_1/1/A/HOH`91/O,/3nd4_1/1/A/G`4/OP2;set depth_cue=0;set ray_trace_fog=0;set dash_color, black;set label_font_id, 5;set label_size, 36;set label_position, (0.5, 1.0, 2.0);set label_color, black;set dash_gap, 0.2;set dash_width, 2.0;set dash_length, 0.2;set label_color, black;set dash_gap, 0.2;set dash_width, 2.0;set dash_length, 0.2;select carbon, element C; color yellow, carbon;disable carbon;set_view (-0.9,0.34,-0.26,0.33,0.18,-0.93,-0.27,-0.92,-0.28,-0.07,-0.23,-27.83,8.63,19.85,13.2,16.0,31.63,-20.0); 
    '''
    cmd.reinitialize();
    cmd.viewport('900','600');
    cmd.fetch('3nd4', type='pdb', async='0');
    cmd.do('run /Users/blaine-mooers/Scripts/PyMOLScripts/quat.py')
    cmd.do('quat 3nd4');
    cmd.show('sticks');
    cmd.set('stick_radius', '0.125');
    cmd.hide('everything', 'name H*');
    cmd.bg_color('white');
    cmd.create('coorCov', '(3nd4_1 and (resi 19 or resi 119 or resi 219 or resi 319 or resi 419 or resi 519 or (resi 3 and name N7)))');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/A`3/N7)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`119/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`219/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`319/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`419/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`519/O)');
    cmd.distance('(3nd4_1 and chain A and resi 19 and name NA)','(3nd4_1 and chain A and resi 519)');
    cmd.distance('(3nd4_1 and chain A and resi 19 and name NA)','(3nd4_1 and chain A and resi 419)');
    cmd.distance('(3nd4_1 and chain A and resi 19 and name NA)','(3nd4_1 and chain A and resi 119)');
    cmd.distance('(3nd4_1 and chain A and resi 19 and name NA)','(3nd4_1 and chain A and resi 319)');
    cmd.distance('(3nd4_1 and chain A and resi 19 and name NA)','(3nd4_1 and chain A and resi 219)');
    cmd.show('nb_spheres');
    cmd.set('nb_spheres_size', '.35');
    cmd.distance('hbond1', '/3nd4_1/1/A/HOH`119/O', '/3nd4_1/1/A/A`3/OP2');
    cmd.distance('hbond2', '/3nd4_1/1/A/HOH`319/O', '/3nd4_1/1/A/A`3/OP2');
    cmd.distance('hbond3', '/3nd4_1/1/A/HOH`91/O', '/3nd4_1/1/A/HOH`119/O');
    cmd.distance('hbond4', '/3nd4_1/1/A/G`4/N7', '/3nd4_1/1/A/HOH`91/O');
    cmd.distance('hbond5', '/3nd4_1/1/A/G`4/O6', '/3nd4_1/1/A/HOH`419/O');
    cmd.distance('hbond6', '/3nd4_1/1/A/HOH`91/O', '/3nd4_1/1/A/G`4/OP2');
    cmd.distance('hbond7', '/3nd4_1/1/A/HOH`319/O', '/3nd4_1/1/A/G`2/OP2');
    cmd.distance('hbond9', '/3nd4_1/1/A/HOH`419/O', '/3nd4_2/2/A/HOH`74/O');
    cmd.distance('hbond10', '/3nd4_2/2/A/C`15/O2', '/3nd4_1/1/A/G`2/N2');
    cmd.distance('hbond11', '/3nd4_2/2/A/C`15/N3', '/3nd4_1/1/A/G`2/N1');
    cmd.distance('hbond12', '/3nd4_2/2/A/C`15/N4', '/3nd4_1/1/A/G`2/O6');
    cmd.distance('hbond13', '/3nd4_2/2/A/U`14/N3', '/3nd4_1/1/A/A`3/N1');
    cmd.distance('hbond14', '/3nd4_2/2/A/U`14/O4', '/3nd4_1/1/A/A`3/N6');
    cmd.distance('hbond15', '/3nd4_2/2/A/C`13/N4', '/3nd4_1/1/A/G`4/O6');
    cmd.distance('hbond16', '/3nd4_2/2/A/C`13/N3', '/3nd4_1/1/A/G`4/N1');
    cmd.distance('hbond17', '/3nd4_1/1/A/G`4/N2', '/3nd4_2/2/A/C`13/O2');
    cmd.distance('hbond18', '/3nd4_1/1/A/G`2/N2', '/3nd4_2/2/A/C`15/O2');
    cmd.distance('hbond19', '/3nd4_1/1/A/HOH`91/O', '/3nd4_1/1/A/G`4/OP2');
    cmd.set('depth_cue', '0');
    cmd.set('ray_trace_fog', '0');
    cmd.set('dash_color', 'black');
    cmd.set('label_font_id', '5');
    cmd.set('label_size', '36')
    cmd.set('label_position', '(0.5, 1.0,2.0)');
    cmd.set('label_color', 'black');
    cmd.set('dash_gap', '0.2');
    cmd.set('dash_width', '2.0');
    cmd.set('dash_length', '0.2');
    cmd.set('label_color', 'black');
    cmd.set('dash_gap', '0.2');
    cmd.set('dash_width', '2.0');
    cmd.set('dash_length', '0.2');
    cmd.select('carbon', 'element C');
    cmd.color('yellow', 'carbon');
    cmd.disable('carbon');
    cmd.set_view('-0.9,0.34,-0.26,0.33,0.18,-0.93,-0.27,-0.92,-0.28,-0.07,-0.23,-27.83,8.63,19.85,13.2,16.0,31.63,-20.0');
cmd.extend('NA',NA)
    
    
    
# ##### Commands applicable to displayed molecules. ##########
def AO():
    '''
DESCRIPTION
    
    Commands to make ambient occlusion image like those in Qutemole. 
    
USAGE
    Type 'AO' to activate. Type 'help AO' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    set_color oxygen, [1.0,0.4,0.4];
    set_color nitrogen, [0.5,0.5,1.0];
    remove solvent;
    as spheres;
    util.cbaw;
    bg white;
    set light_count,10;
    set spec_count,1;
    set shininess, 10;
    set specular,0.25;
    set ambient,0;
    set direct,0;
    set reflect,1.5;
    set ray_shadow_decay_factor, 0.1;
    set ray_shadow_decay_range, 2;
    unset depth_cue;
    ray
    
    Commands without linebreaks:
    set_color oxygen, [1.0,0.4,0.4];set_color nitrogen, [0.5,0.5,1.0];remove solvent;as spheres;util.cbaw;bg white;set light_count,10;set spec_count,1;set shininess, 10;set specular,0.25;set ambient,0;set direct,0;set reflect,1.5;set ray_shadow_decay_factor, 0.1;set ray_shadow_decay_range, 2;unset depth_cue;ray 
    
    '''
    cmd.set_color('oxygen', '[1.0,0.4,0.4]')
    cmd.set_color('nitrogen', '[0.5,0.5,1.0]')
    cmd.remove('solvent')
    cmd.show_as('spheres')
    cmd.util.cbaw()
    cmd.bg_color('white')
    cmd.set('light_count', '10')
    cmd.set('spec_count', '1')
    cmd.set('shininess', '10')
    cmd.set('specular', '0.25')
    cmd.set('ambient', '0')
    cmd.set('direct', '0')
    cmd.set('reflect', '1.5')
    cmd.set('ray_shadow_decay_factor', '0.1')
    cmd.set('ray_shadow_decay_range', '2')
    cmd.unset('depth_cue')
    cmd.ray()
cmd.extend('AO',AO)
    
    
def BS():
    '''
    DESCRIPTION
    
    Commands to make ball-and-stick representation. Applied globally.
    Recycle code to include a selection in the "show sticks" command. 
    
USAGE
    Type 'BS' to activate. Type 'help BS' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    hide everything, name H 
    show sticks;
    set stick_radius, 0.12;
    set stick_ball, on; 
    set stick_ball_ratio, 1.9;
    show nb_spheres;
    set nb_spheres_size=0.33
    
    Commands without linebreaks:
    show sticks;set stick_radius, 0.12;set stick_ball, on; set stick_ball_ratio, 1.9;show nb_spheres; set nb_spheres_size=0.33 
    '''
    cmd.show('sticks')
    cmd.set('stick_radius', '0.12')
    cmd.set('stick_ball', 'on')
    cmd.set('stick_ball_ratio', '1.9')
    cmd.show('nb_spheres')
    cmd.set('nb_spheres_size', '0.33')
    cmd.hide('everything', 'name H*')
cmd.extend('BS',BS)
    
def BW():
    '''
    DESCRIPTION
    
    Commands to make black-and white-ribbon cartoon on a white background.
    Good for avoiding color figure charges. Requires a pdb file. 
    
    USAGE
    Orient struture as desired. Then type 'BW' to execute the function. 
    Type 'help BW' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    show cartoon; 
    hide lines; 
    hide nonbonded; 
    # black and white cartoon;
    #note how the dcomment is on a separate line and not to the right of a command; 
    set ray_trace_mode, 2; 
    bg_color white; note how comment is on a separate line and not include to the right of a command; 
    set antialias, 2; 
    ray 1600,1600; 
    png test.png
    
    Commands without linebreaks:
    show cartoon; hide lines; hide nonbonded; set ray_trace_mode, 2; # black and white cartoon; bg_color white; set antialias, 2; ray 1600,1600; png test.png
    '''
    cmd.show_as("cartoon", "all"); 
    cmd.hide('lines'); 
    cmd.hide('nonbonded'); 
    # black and white cartoon; 
    cmd.set('ray_trace_mode', '2'); 
    cmd.bg_color('white'); 
    cmd.set('antialias', '2'); 
    cmd.ray('600','600'); 
    cmd.png('test.png')
cmd.extend('BW',BW) 

    
def BU():
    '''
    DESCRIPTION
    
    Commands to make biological unit. Requires a pdb file. There are
    other ways of displaying the biological unit in PyMOL. Depends on
    the quat.py script by Thomas Holder.
    
    USAGE
    
    Type 'BU' to activate. Type 'help BU' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.

    Commands without linebreaks:
    run ~/Scripts/PyMOLScripts/quat.py; quat 
    '''
#    cmd.alias('aBU', 'run ~/Scripts/PyMOLScripts/quat.py; quat') 
    cmd.do('run /Users/blaine-mooers/Scripts/PyMOLScripts/quat.py')
#    cmd.run('/Users/blaine-mooers/Scripts/PyMOLScripts/quat.py')
    cmd.do('quat')
cmd.extend('BU',BU)    
    
def CB():
    '''
    DESCRIPTION
    
    Loads Jared Sampson's script "colorblindfriendly.py" from the
    ~/Pymol-script-repo directory. The new colorblind-friendly color
    names are printed to the command history window and are available
    for use like standard colors in PyMOL.

    Read the header of the script file for more information. Although
    this script is in my PyMOL plugin collection, it is out of sight
    and out of mind. The listing of the alias to this script is a
    reminder to use it. 
    
USAGE
    Type 'CB' to activate. Type 'help CB' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    Commands without linebreaks:
    run ~/Pymol-script-repo/colorblindfriendly.py 
    '''
    cmd.run('~/Pymol-script-repo/colorblindfriendly.py')
cmd.extend('CB',CB)
    
    
def CR():
    '''
    DESCRIPTION
    
    Commands to make colored filled-ring cartoon of nucleic acids. May
    need to 'hide everything' first. If asymmetric unit has one strand
    of a dsRNA or dsDNA, remember to apply the BU alias to display the
    second strand. ddd
    
    Adapted from KP Wu's blog post:
    https://kpwu.wordpress.com/2012/05/24/pymol-different-colors-of-
    nucleic-acid-rings/
    
USAGE
    Type 'CR' to activate. Type 'help CR' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    hide everything;
    bg_color white;
    cartoon oval;
    set cartoon_ring_mode,3;
    set cartoon_nucleic_acid_color, blue;
    select rna_A, resn A;
    select rna_C, resn C;
    select rna_G, resn G;
    select rna_U, resn U;
    color yellow, rna_A;
    color red, rna_C;
    color gray40, rna_G;
    color palecyan, rna_U;
    as cartoon 
    
    Commands without linebreaks:
    hide everything;bg_color white;cartoon oval;set cartoon_ring_mode,3;set cartoon_nucleic_acid_color, blue;select rna_A, resn A;select rna_C, resn C;select rna_G, resn G;select rna_U, resn U;color yellow, rna_A;color red, rna_C;color gray40, rna_G;color palecyan, rna_U;as cartoon 
    '''
    cmd.hide('everything')
    cmd.bg_color('white')
    cmd.cartoon('oval')
    cmd.set('cartoon_ring_mode', '3')
    cmd.set('cartoon_nucleic_acid_color', 'blue')
    cmd.select('rna_A', 'resn A')
    cmd.select('rna_C', 'resn C')
    cmd.select('rna_G', 'resn G')
    cmd.select('rna_U', 'resn U')
    cmd.color('yellow', 'rna_A')
    cmd.color('red', 'rna_C')
    cmd.color('gray40', 'rna_G')
    cmd.color('palecyan', 'rna_U')
    cmd.show_as('cartoon')
    cmd.disable('rna_U')    
cmd.extend('CR',CR)
    
    
def CSS():
    '''
    DESCRIPTION
    
    Commands to color ribbon or cartoon representations of proteins by
    secondary structures. 
    
USAGE
    Type 'CSS' to activate. Type 'help CSS' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    Commands without linebreaks:
    as cartoon;
    color red, ss H;
    color yellow,ss S;
    color green, ss L+; 
    '''
    cmd.show_as('cartoon')
    cmd.color('red', 'ss H')
    cmd.color('yellow', 'ss S')
    cmd.color('green', 'ss L+')
cmd.extend('CSS',CSS)
    
    
def CBSS():
    '''
    DESCRIPTION
    
    Commands to do colorblind friendly coloring by secondary
    structures of ribbon or cartoon representations of proteins.
    Depends on colorblindfriendly.py. Script assumed to be stored in
    the home directory in ~/Pymol-script-repo.
    
USAGE
    Type 'CBSS' to activate. Type 'help CBSS' to see this
    documentation printed to the command history window. Select from
    the command history individual lines of code to build a new
    script. Select the hortizontal script at the bottom if retaining
    most of the commands in your new script. Copy and paste onto the
    comand line below. Works only with the command line immediately
    under the command history window at the top of the gui.
    
    run ~/Pymol-script-repo/colorblindfriendly.py;
    as cartoon;
    color cb_red, ss H;
    color cb_yellow,ss S;
    color cb_green, ss L+; 
    
    Commands without linebreaks:
    run ~/Pymol-script-repo/colorblindfriendly.py;as cartoon;color cb_red, ss H;color cb_yellow,ss S;color cb_green, ss L+; 
    '''
    
    cmd.run('~/Pymol-script-repo/colorblindfriendly.py')
    cmd.show_as('cartoon')
    cmd.color('cb_red', 'ss H')
    cmd.color('cb_yellow', 'ss S')
    cmd.color('cb_green', 'ss L+')
cmd.extend('CBSS',CBSS)
    
    
def DU():
    '''
    DESCRIPTION
    
    Commands to make dumbbell (ribbons with rolled edges) cartoon of
    the main chains of nucleic acids and proteins. 50% transparent
    cartoon so it can be combined with lines, sticks, and
    ball-and-sticks (try BS alias).
    
USAGE
    Type 'DU' to activate. Type 'help DU' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    cartoon dumbbell;
    set cartoon_dumbbell_width, 0.2;
    set cartoon_dumbbell_radius, 0.4;
    show cartoon; 
    
    Commands without linebreaks:
    cartoon dumbbell;set cartoon_dumbbell_width, 0.2;set cartoon_dumbbell_radius, 0.4;show cartoon; 
    '''
    cmd.cartoon('dumbbell')
    cmd.set('cartoon_dumbbell_width', '0.2')
    cmd.set('cartoon_dumbbell_radius', '0.4')
    cmd.show('cartoon')
cmd.extend('DU',DU)    
    
    
def FR():
    '''
    DESCRIPTION
    
    Commands to make filled-ring cartoon of nucleic acids. May need to
    'hide everything' first. Adapted from script on
    http://www-cryst.bioc.cam.ac.uk/members/zbyszek/figures_pymol. 
    
USAGE
    Type 'FR' to activate. Type 'help FR' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    show sticks;set cartoon_ring_mode, 3;
    set cartoon_ring_finder, 1;
    set cartoon_ladder_mode, 1;
    set cartoon_nucleic_acid_mode, 4;
    set cartoon_ring_transparency, 0.5;
    as cartoon;
    
    Commands without linebreaks:
    show sticks;set cartoon_ring_mode, 3;set cartoon_ring_finder, 1;set cartoon_ladder_mode, 1;set cartoon_nucleic_acid_mode, 4;set cartoon_ring_transparency, 0.5;as cartoon; 
    '''
    cmd.show('sticks')
    cmd.set('cartoon_ring_mode', '3')
    cmd.set('cartoon_ring_finder', '1')
    cmd.set('cartoon_ladder_mode', '1')
    cmd.set('cartoon_nucleic_acid_mode', '4')
    cmd.set('cartoon_ring_transparency', '0.5')
    cmd.show_as('cartoon')
cmd.extend('FR',FR)    
    
    
def HH():
    '''
    DESCRIPTION
    
    Commands to hide hydrogen atoms. 
    
USAGE
    Type 'HH' to activate. Type 'help HH' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    Commands without linebreaks:
    hide everything, name H* 
    
    '''
    cmd.hide('everything', 'name H*') 
cmd.extend('HH',HH)    
    
    
def PU():
    '''
    DESCRIPTION
    
    Commands to make putty cartoon of main chain of nucleic acids and
    proteins. The radius of the cartoon is inversely proportional to
    the B-factors. 50% transparent cartoon so it can be combined with
    lines, sticks, and ball-and-sticks (try BS alias).
    
USAGE
    Type 'PU' to activate. Type 'help PU' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    cartoon putty;
    set cartoon_ladder_mode, 0;
    set cartoon_transparency,0.5;
    set cartoon_ring_finder, 0;
    show cartoon 
    
    Commands without linebreaks 
    cartoon putty;set cartoon_ladder_mode, 0;set cartoon_transparency,0.5;set cartoon_ring_finder, 0;show cartoon 
    '''
    cmd.cartoon('putty')
    cmd.set('cartoon_ladder_mode', '0')
    cmd.set('cartoon_transparency', '0.5')
    cmd.set('cartoon_ring_finder', '0')
    cmd.show('cartoon')
cmd.extend('PU',PU)
    
""" Print the aliases on startup of PyMOL"""
print(SA.__doc__)
