# PyMOL made easier with horizontal scripting
*Script to facilitate the making of horizontal scripts*

Welcome to the `EasyPyMOL` repository, which implements the PyMOL approaches described in the manuscript:
[Blaine Mooers](http://www.oumedicine.com/department-of-biochemistry-and-molecular-biology/faculty/blaine-mooers-ph-d-) (*under revision*).
"Easier use of PyMOL with horizontal scripts". 

Click on icon to see 4 minute introductory video:

[![IMAGE ALT TEXT HERE](http://img.youtube.com/vi/watch?v=GnRtEhGvPBQ/0.jpg)](http://www.youtube.com/watch?v=XRsAaKq4afs)

#### Ambient occlusion image of phage T4 lysozyme made with the alias "AO" 
17 commands were placed on one line as a horizontal script.


<img src="https://cloud.githubusercontent.com/assets/15176203/13590209/9ad0758c-e4a3-11e5-995a-0ed5fb2cc88f.png" width="75%"/> 

####  Ribbon diagram of the above molecule made with the alias "T4L" 
Another image made with a horizontal script.

<img src="https://cloud.githubusercontent.com/assets/15176203/10561743/5f59bb82-74fd-11e5-828c-dbe1dbb2c648.png" width="60%"/>




### *Problem*: the view  port settings return from [`get_view()`](http://pymolwiki.org/index.php/Get_View) in [PyMOL](https://www.pymol.org/) is too hard to copy and paste quickly onto one line with other commands during horizontal scripting.

```python
PyMOL>get_view  

### cut below here and paste into script ###  
set_view(\
  -0.832868993,    0.398498207,    0.383896619,\
  -0.260102808,   -0.894237876,    0.363985002,\
   0.488390923,    0.203309149,    0.848513067,\
   0.000000000,    0.000000000,  -61.396984100,\
 -46.246913910,   -4.663769245,   42.401920319,\
  56.260883331,   66.533096313,  -20.000000000 )
### cut above here and paste into script ###
```

## *Solution*: The function `roundview()` in the script [`roundview.py`](.\roundview.py) which is available from this website.

```py
PyMOL>roundview  
set_view(-0.83,0.4,0.38,-0.26,-0.89,0.36,0.49,0.2,0.85,0.0,0.0,-61.4,-46.25,-4.66,42.4,56.26,66.53,-20.0);
```
    
Paste the above reformated [`set_view()`](http://pymolwiki.org/index.php/Set_View) command onto the PyMOL [command line](http://pymolwiki.org/index.php/Command_Line_Options) in the external gui or into a [script](http://www.pymolwiki.org/index.php/Running_Scripts) in a plain text editor.

What this is
--------------------------------------------------------------------------------

The python script roundview.py includes the short function `roundview()` that reformats the viewing port settings from seven rows to one row. The more compact format from roundview() is easy to copy and paste onto the command line. Other commands that are separated by semicolons can be added to the command line along with the settings. This defines a horiontal script. The script can include commnents that are isolated by semicolons. The horizontal script can be edited and tested repeatedly within PyMOL for many cycles without using an external text editor. This saves time during the development of a new molecular scene. The cursor can be moved around quickly on the command line with the readline commands:

* **cntrl-a**  moves cursor to the beginning of the line
* **cntrl-e**  moves cursor to the end of the line
* **shift-cntrl-a**  selects everthing to from the cursor to the beginning of the line
* **shift-cntrl-e**  selects everthing to from the cursor to the end of the line
* **command-f**  move forward by one word
* **commend-b**  move backward by one word

Requirements
--------------------------------------------------------------------------------

Requires a molecular object loaded into an interactive session of PyMOL. Does not
require any modules outside of two in PyMOL. Should work on all versions of PyMOL. 
Tested on:
* Ubuntu 14.04 64 bit with PyMOL 1.7.2.2. 
* Windows 8 32 bit  running PyMOL 1.7.6.2 and PyMOL 1.7.6.6. 
* Mac OSX 10.10.5 64 bit running PyMOL 1.5.0.5, 1.7.6.6 (via macports), and 1.8.0.5. 

We tried to make the code backward compatible to PyMOL 1.5. We do not gaurantee that this
code works with earlier versions of PyMOL. Nor do we guarantee that the code will
not fail in future versions of PyMOL. In PyMOL version 1.6, there were
several changes that reduced backward compatibility. We also do not guarantee that
the code will work if you install in a way that is different from the methods that we describe here.

Instructions
--------------------------------------------------------------------------------

#### Quick start instructions for beginning users of Github

Watch this 1 minute video or read on. 

[![IMAGE ALT TEXT HERE](http://img.youtube.com/vi/watch?v=GnRtEhGvPBQ/0.jpg)](http://www.youtube.com/watch?v=GnRtEhGvPBQ)


Copy script from [this link](https://github.com/MooersLab/EasyPyMOL/blob/master/roundview.py) after clicking on "RAW" in the upper right corner and paste into a plain text file (NOT a doc, docx, or rtf file). Name the script [`roundview.py`](./roundview.py). Save to your home directory (e.g., /Users/username or /home/username or C:\Users\username). Start PyMOL. Check that PyMOL's current directory is in the home directory by entering the "pwd" command on the command line in PyMOL. Check for presence of roundview.py with "ls *.py". 

```shell
ls *.py
```
Paste the following horizontal script on the command line just below the command history window in the top or external gui:
```shell
fetch 1lw9, async=0; run roundview.py; roundview 0,1
```
You should see the following in the command history window of the top gui:
```shell
set_view (1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,-155.2,35.1,11.5,9.7,122.3,188.0,-20.0);
```

Type the following to see the default format from get_view.
```shell
get_view
```

Which looks easier to copy from the command history window and paste onto the command line in PyMOL?

#### More advanced ways to install roundview.py

Get the script. Either download the folder from [this link](https://github.com/MooersLab/EasyPyMOL/archive/master.zip) or type the following command in a terminal window in your home directory or the directory where you store your PyMOL python scripts 

```shell
git clone https://github.com/MooersLab/EasyPyMOL.git
```

You need the program [Git](https://git-scm.com/) installed on your computer. Git is available via macports on a mac or otherwise [see these instructions for installing git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).


There are at least four different ways of loading the script into PyMOL:

1. move [`roundview.py`](./roundview.py) to the working directory of pymol. In pymol, type 

   ```py
    run roundview.py
   ```
Please note that the "run" command just loads the script into PyMOL. It does not execute it. Now the roundview command is available by typing "roundview", and the on-line documentation is available by typing "roundview".

2. load [`roundview.py`](./roundview.py) with the plugin manager (see plugin pulldown) in PyMOL

3. copy [`roundview.py`](./roundview.py) to safe folder that will not be deleted when you delete PyMOL. I use `~/Scripts/PyMOLScripts/`. Then load into PyMOL using method 1 or 2. 

4. create or edit the hidden file `.pymolrc` (`pymolrc.pml` and not hidden on Windows) text file in the home directory so that it includes the following lines so that roundview.py is always loaded upon startup. This option also works without the first two lines. 
  
   ```py
    import sys
    sys.path.append('/Path/To/roundview')
    run /Path/To/roundview.py
   ```
For example:
   ```py
    import sys
    sys.path.append('/Users/blaine-mooers/Scripts/PyMOLScripts/')
    run /Users/blaine-mooers/Scripts/Scripts_PyMOL/roundview.py
   ```
Restart pymol. You should see something like the following in the command history window if your path to the script is correct. 
   ```py
    PyMOL>import sys
    PyMOL>sys.path.append('/Users/blaine-mooers/Scripts/PyMOLScripts/')
    PyMOL>run /Users/blaine-mooers/Scripts/Scripts_PyMOL/roundview.py
   ```

Type "roundview" on either command line. You should get back this if no molecule is loaded:
   ```py
    set_view (1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,-50.0,0.0,0.0,0.0,40.0,100.0,-20.0);
   ```
   
## Link to roundview.py installation video for MacPyMOL:

[![IMAGE ALT TEXT HERE](http://img.youtube.com/vi/watch?v=XTwwWgwj4z0/0.jpg)](http://www.youtube.com/watch?v=XTwwWgwj4z0)

## Link to roundview.py installation video for PyMOL om Windows:

[![IMAGE ALT TEXT HERE](http://img.youtube.com/vi/watch?v=XTwwWgwj4z0/0.jpg)](http://www.youtube.com/watch?v=Dh2ihD5mIxY)

## Link to roundview.py installation video for Linux and PyMOLHyridX11 for the mac on YouTube:

[![IMAGE ALT TEXT HERE](http://img.youtube.com/vi/watch?v=XQWQzq48DeA/0.jpg)](http://www.youtube.com/watch?v=XQWQzq48DeA)


#### How to run
--------------------------------------------------------------------------------

After loading a pdb file and setting up the molecular scene, type on one of the command lines in PyMOL:

   ```py
    PyMOL>roundview
   ```

##### How to get help

Type the following on one of the command lines in PyMOL

   ```py
    PyMOL> help roundview
   ```

Usage: `roundview [view, decimal places, outname]` 
* The values in the [ ] are optional.  
* The default view is "0". 
* The default number of `decimal places` is 2. 
* The `outname` is the name of a plain text file to which the output of `roundveiw()` is written.

Quick test with wildtype bacteriophage T4 lysozyme, 3fa0
--------------------------------------------------------------------------------

Open up pymol. Copy and paste the entire line below onto the command line in the external gui (the on above the GL viewing port). This is an example of a horizontal script. By hitting the up arrow key, you can recall this command for editing on the command line. This code block is more agile to edit than opening, editing, saving, and loading an external script file.

```py
fetch 3fa0,async=0;orient;turn z,-90;turn y,-5;turn x,10; hide everything; bg_color white; show cartoon;color red, ss H;color yellow, ss S;color green, ss L+'';roundview  
```
  
To apply the canonical view similar to that in [Gassner et al. 2003](http://www.ncbi.nlm.nih.gov/pubmed/12646375), copy and paste the following onto the command line:
  
```py
set_view (-0.18,-0.69,-0.7,0.98,-0.17,-0.09,-0.06,-0.7,0.71,0.0,0.0,-154.87,34.77,11.27,9.52,121.27,188.47,-20.0); ray 1500,1600; png coanonical.png
```
  	
You should get back an image that looks like the following:

<img src="https://cloud.githubusercontent.com/assets/15176203/10561743/5f59bb82-74fd-11e5-828c-dbe1dbb2c648.png" width="90%"></img> 

To test some other argument values, copy and paste the following command into PyMOL. 

```py
PyMOL>roundview 0,1,firstscene.txt
```

Do a `ls *.txt` to list the files in the working diretory. The file "firstscene.txt" should be listed. The default filename is "roundedview.txt". This file is appended with each execution of the roundview command. You may find it easier to copy the set_view line from this text file than from the command history window in PyMOL.

## Scripts that use aliases to horizontal scripts. Some aliases contain compact scene settings from roundview().

### Exam2function.py

  Defines aliases q1-q8 for questions 1-8 from exam 2 of the OUHSC Macromolecular Systems course. Each alias is 
  is mapped to a number of commands.  
  
  Create ~/Scripts/PyMOLScripts and store the script in this subfolder. 
  
  Enter on the command line in PyMOL the following command:
 
```py  
run ~/Scripts/PyMOLScripts/exam2function.py
```

Now the aliases q1,q2, ..., q8 are active.
  
Type 'q1' to execute the alias assoicated with exam quesitno 1.
  
Type 'help q1' to print the documentation to the PyMOL command history window. 
The bottom of the documentation includes the corresponding horizontal script.
All or parts of the horizontal script can be copied from the command history 
window for reuse of the code in another horizontal script or a traditional 
vertical script.
  
Tested on PyMOL versions 1.5.0.5, 1.7.7.2 (from macports), and 1.8.0.5. 


### StartupAliases.py

Copy to ~/Scripts/PyMOLScripts/.
  
Add this command on one line in your .pymolrc file (pymolrc.pml on Windows):

```py    
run ~/Scripts/PyMOLScripts/StartUpAliases.py
```
Now these aliases will be available whenever you startup PyMOL.

Type the 'alias name' to execute it.

Type 'help alias name' to see the documentation, which includes a vertical list of the commands mapped to the alias to ease the copying of isolated commands from the command history window during code reuse. The corresponding horizontal script without line breaks is also printed. It can be selected in the command history window and pasted onto the command line. 

Format of list below:

Alias name, description: PDB code, where applicable. 
  
#### Molecules in standard orientations: 
  
* T4L, WT T4 lysozyme (1.09 ang) as a ribbon diagram: 3fa0. 
* U8, 16-mer dsRNA with 8 contiguous Us. U-helix RaNA (1.37 ang): 3nd3.
* WC8, 16-mer RNA with all Watson-Crick base pairs (1.67 ang): 3nd4.
* N9, neuraminidase as cartoon, biological unit (1.55 ang): 4dgr.
* GGT, gamma glutamyl transpeptidase as cartoon (1.67 ang): 4gdx.
* GU, 10-mer RNA with eight GU base pairs (1.32 ang): 4pco.


#### Complex figures to serve as templates: 
  
* BST, Base-stacking figure, (1.32 ang): 4pco. 
* LG, Electron density map of nine sugar glycan,(1.55 ang):, 4dgr. 
* NA, Sodium cation in major groove of 16-mer RNA: 3nd4.


#### Complex representations applied to any visible molecular object:
  
* AO, Make ambient occlusion image. Requires global view of protein.
* BU, Display biological unit. 
* CB, Define color blind compatible coloring scheme. 
* BW, Make black and white ribbon cartoon on white background.
* CSS, Color ribbon and cartoons by secondary structure: red, green and yellow. 
* CBSS, Color ribbon and cartoons with dcolorblind friendly colors. 
* CR, Commands to make colored filled-ring cartoon of nucleic acids..
* FR, Commands to make filled-ring cartoon of nucleic acids.

#### Demo 1: 

Type "T4L" on the command line. Now type "AO". You should get an image like the following:

<img src="https://cloud.githubusercontent.com/assets/15176203/13590209/9ad0758c-e4a3-11e5-995a-0ed5fb2cc88f.png" width="90%"></img>


Type 'help AO' on the command line to see the documentation for the AO alias. It is mapped to 17 commands. 

PyMOL>help AO
 
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

#### Demo 2: 

Type "T4L" on the command line. Now type "BW". You should get a black and white  image like the following. These black and white
figures are useful when color figures are not needed:

<img src="https://cloud.githubusercontent.com/assets/15176203/13590390/201de41c-e4a5-11e5-9835-7aced7982306.png" width="90%"></img>


#### Demo 3:

Type "U8" on the command line. Convert this to black and white with "BW".

<img src="https://cloud.githubusercontent.com/assets/15176203/13672122/089167a6-e699-11e5-89bd-19307f6e1181.png" width="90%"></img>

#### Demo 4: 

Type 'help NA' to see a very long script mapped to two a letter command.

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
    set depth_cue=0;
    set ray_trace_fog=0;
    set dash_color, black;
    set label_font_id, 5;
    set label_size, 36;
    set label_position, (0.5, 1.0, 2.0);
    set label_color, black;set dash_gap, 0.2;
    set dash_width, 2.0;set dash_length, 0.2;
    set label_color, black;set dash_gap, 0.2;
    set dash_width, 2.0;set dash_length, 0.2;
    select carbon, element C; 
    color yellow, carbon;
    disable carbon;
    set_view (-0.9,0.34,-0.26,0.33,0.18,-0.93,-0.27,-0.92,-0.28,-0.07,-0.23,-27.83,8.63,19.85,13.2,16.0,31.63,-20.0)
    
    Commands without linebreaks for easy selecting, copying, and pasting onto the command line for code reuse:
    delete all;viewport 900,600;fetch 3nd4, type=pdb,async=0;run ~/Scripts/PyMOLScripts/quat.py;quat 3nd4; show sticks;set stick_radius=0.125;hide everything, name H*;bg_color white;create coorCov, (3nd4_1 and (resi 19 or resi 119 or resi 219 or resi 319 or resi 419 or resi 519 or (resi 3 and name N7)));bond (coorCov//A/NA`19/NA),(coorCov//A/A`3/N7); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`119/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`219/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`319/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`419/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`519/O);distance (3nd4_1 and chain Aand resi 19 and name NA), (3nd4_1 and chain A and resi 519);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 419);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 119);distance (3nd4_1 and chain A and resi 19 and name NA),(3nd4_1 and chain A and resi 319);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 219);show nb_spheres; set nb_spheres_size, .35;distance hbond1,/3nd4_1/1/A/HOH`119/O, /3nd4_1/1/A/A`3/OP2;distance hbond2,/3nd4_1/1/A/HOH`319/O, /3nd4_1/1/A/A`3/OP2;distance hbond3,/3nd4_1/1/A/HOH`91/O, /3nd4_1/1/A/HOH`119/O;distance hbond4,/3nd4_1/1/A/G`4/N7,/3nd4_1/1/A/HOH`91/O;distance hbond5,/3nd4_1/1/A/G`4/O6, /3nd4_1/1/A/HOH`419/O;distance hbond6,/3nd4_1/1/A/HOH`91/O, /3nd4_1/1/A/G`4/OP2;distance hbond7,/3nd4_1/1/A/HOH`319/O, /3nd4_1/1/A/G`2/OP2;distance  hbond9,/3nd4_1/1/A/HOH`419/O,/3nd4_2/2/A/HOH`74/O;distance hbond10,/3nd4_2/2/A/C`15/O2,/3nd4_1/1/A/G`2/N2;distance hbond11, /3nd4_2/2/A/C`15/N3,/3nd4_1/1/A/G`2/N1;distance hbond12,/3nd4_2/2/A/C`15/N4,/3nd4_1/1/A/G`2/O6;distance hbond13, /3nd4_2/2/A/U`14/N3,/3nd4_1/1/A/A`3/N1;distance hbond14,3nd4_2/2/A/U`14/O4,/3nd4_1/1/A/A`3/N6;distance hbond15, /3nd4_2/2/A/C`13/N4,/3nd4_1/1/A/G`4/O6;distance hbond16,/3nd4_2/2/A/C`13/N3, /3nd4_1/1/A/G`4/N1;distance hbond17, /3nd4_1/1/A/G`4/N2,/3nd4_2/2/A/C`13/O2;distance hbond18,/3nd4_1/1/A/G`2/N2,/3nd4_2/2/A/C`15/O2;distance hbond19,/3nd4_1/1/A/HOH`91/O,/3nd4_1/1/A/G`4/OP2;set depth_cue=0;set ray_trace_fog=0;set dash_color, black;set label_font_id, 5;set label_size, 36;set label_position, (0.5, 1.0, 2.0);set label_color, black;set dash_gap, 0.2;set dash_width, 2.0;set dash_length, 0.2;set label_color, black;set dash_gap, 0.2;set dash_width, 2.0;set dash_length, 0.2;select carbon, element C; color yellow, carbon;disable carbon;set_view (-0.9,0.34,-0.26,0.33,0.18,-0.93,-0.27,-0.92,-0.28,-0.07,-0.23,-27.83,8.63,19.85,13.2,16.0,31.63,-20.0); 


Type "NA" to get the resulting image of a sodium cation bound with inner sphere coordination to the N7 nitrogen of an adenine and to five waters. The sodium is in the major groove of a double-strande RNA molecule (PDB-ID 3nd4). The dashed lines represent hydrogen bonds. The numbers are distances in angstroms.

![naadenine](https://cloud.githubusercontent.com/assets/15176203/13609887/853536a4-e520-11e5-9978-9ea1da4f0884.png)


Reference, License, and Copyright
--------------------------------------------------------------------------------

* Mooers, B. H. M. (submitted) Easier use of PyMOL with horizontal scripts.
* GNU General Public License ([GPL-3](http://www.gnu.org/licenses/gpl-3.0.en.html))
* (C) [Blaine Mooers](http://www.oumedicine.com/department-of-biochemistry-and-molecular-biology/faculty/blaine-mooers-ph-d-), Ph.D.,  [University of Oklahoma Health Sciences Center](http://www.ouhsc.edu/), 2015-2016
