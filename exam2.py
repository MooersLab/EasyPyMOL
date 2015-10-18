import pymol
from pymol import cmd

print("2015 OUHSC MSI-Lecture 17 and 18 takehome exam")
print("on nucleic acid structure.")
print( " ")
print("Enter 'q1' for questin 1, 'q2' for question 2, and so on")
print("through 'q8' for question 8.")
print( " ")
print( "Questions (1-4). HIV reverse transcriptase bound to DNA")
print("with AZT incorporated (PDB-ID 3V6D).")
print( " ")
print("Questions (5-8). The full-length hammerhead")
print("ribozyme (PDB-ID 3zp8).")
print(" ")
print( "Enter 'help exam2' to learn how to use this script.'")
print( "Enter 'help q1' to see quesiton 1 and so on.'")

def exam2(question=0):
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
        defined by a specific function quX(), where X is the
        integer number of the question. Each function
        defines an alias to a row of PyMOL commands. Each
        function also has its own docstring that explains
        the molecular scene. and that is printed to the
        command history window when you type 'help qu1' or
        'help qu2' and so on. Each row of commands starts
        with 'reinitialize' to give you a fresh start. Each
        function is followed an alias command that defines
        'q1' as an alias to the combined command 'qu1;aq1'.
        This compound command calls the qu1() function and
        then the alias that is defined by the function qu1().
        
    LIMITATIONS
        
        Each alias uses the fetch command, so an active
        internet connection is required. Replace the fetch
        command with the load command to load coordinates
        from your harddrive. 
        
    REQUIREMENTS
            
        When developing a variant of this program., run the 
        program roundview.py to obtained a shorter version 
        of the viewport settings.

    NOTES
    
        Use gray scroll bar to the right side of the command 
        history window to see all of this documentation. or hit 
        escape after moving the cursor to the viewing port 
        (the GUI with the molecule). 

    Copyright (c) 2015 Blaine Mooers PhD, 
    OUHSC, September 20, 2015
    """ 
    print "The function exam2() stores this docstring." 
cmd.extend( "exam2",exam2)


def qu1(question=1):
    '''
    DESCRIPTION

    Question 1: Explain how AZT terminates extension of the
    DNA chain. 
    '''
#    print question
    cmd.alias('aq1','reinitialize;fetch 3v6d, HIVrt,async=0;\
    rock;preset.ball_and_stick("c. P and i. 822");\
    set_view (-0.99,-0.1,0.06,0.09,-0.39,0.92,-0.07,\
    0.92,0.39,0.0,0.0,-29.2,-10.56,24.72,39.27,\
    23.02,35.38,-20.0)') 
    print('Enter "q1" to set scene for question 1.' )
    print('Enter "help(qu1)" for more information') 
    print('This is the definition of the alias "aq1": \
    fetch 3v6d, HIVrt, async=0; rock;\
    preset.ball_and_stick("c. P and i. 822");\
    set_view (-0.99,-0.1,0.06,0.09,-0.39,0.92,-0.07,\
    0.92,0.39,0.0,0.0,-29.2,-10.56,24.72,39.27,23.02,\
    35.38,-20.0)')
cmd.extend( 'qu1',qu1)
cmd.alias('q1','qu1;aq1')


def qu2(question=2):
    '''
    DESCRIPTION

    Question 2: Look along the helical axis of the DNA. Is
    the DNA helix bent?
    '''
#    print question
    cmd.alias('aq2','reinitialize;\
    fetch 3v6d, HIVrt, async=0;\
    rock;show cartoon, c. P or c. T;\
    set_view (0.61,0.29,-0.73,0.78,-0.37,\
    0.51,-0.12,-0.88,-0.45,-0.0, 0.0,-192.99,\
    -29.84,8.42,47.76,178.27,207.7,-20.0);')
    print('Enter "q2" to set scene for question 2.')    
    print('Enter "help(qu2)" for more information')
    print('This is the definition of the alias "aq2": \
    fetch 3v6d, HIVrt, async=0; rock;show cartoon, c. P or c. T;\
    set_view (0.61,0.29,-0.73,0.78,-0.37,0.51,-0.12,\
    -0.88,-0.45,-0.0,-0.0,-192.99,-29.84,8.42,47.76,\
    178.27,207.7,-20.0);')
cmd.extend( 'qu2',qu2)
cmd.alias('q2','qu2;aq2')
#print 'Enter "reinitialize" before going to next problem"


def qu3(question=3):
    '''
    DESCRIPTION

    Question 3: In which groove is the protein making the
    most contacts? How is this unusual? 
    '''
#    print question
    cmd.alias('aq3','reinitialize;fetch 3v6d, HIVrt,\
    async=0; show cartoon, c. A or c. B;\
    rock;set_view (0.53,-0.06,-0.84,0.82,-0.21,0.53,\
    -0.21,-0.98,-0.07,-0.0,-0.0,-192.99,-29.84,8.42,\
    47.76,178.27,207.7,-20.0);')
    print('Enter "q3" to set scene for question 3.')
    print('Enter "help(qu3)" for more information')
    print('This is the definition of the alias "aq3": \
    fetch 3v6d, HIVrt, async=0; show cartoon, c. A or c. B;rock; \
    set_view (0.53,-0.06,-0.84,0.82,-0.21,0.53,-0.21,\
    -0.98,-0.07,-0.0,-0.0,-192.99,-29.84,8.42,47.76,\
    178.27,207.7,-20.0)')
cmd.extend( 'qu3',qu3)
cmd.alias('q3','qu3;aq3')
#print 'Enter "reinitialize" before going to next problem"


def qu4(question=4):
    '''
    DESCRIPTION

    Question 4:  What the dihedral angle about the disulfide
    bond between MRG81 of chain F and Cys258 of chain C?
    This is a cross link between the protein and the DNA. Is
    this a cis or trans conformation of the bonds about the
    S--S bond? Is this conformation energetically favorable
    or unfavorable? 
    '''
#    print question
    cmd.alias('aq4','reinitialize;fetch 3v6d, HIVrt, async=0;\
    preset.ball_and_stick("(c. C and i. 258) or (c. F and i. 817 )");\
    set_view (0.48,0.84,-0.27,-0.81,0.3,-0.51,-0.35,0.46,0.82,\
    0.0,-0.0,-27.26,-40.77,-53.15,12.15,25.47,29.09,-20.0);')
    print('Enter "q4" to set scene for question 4.')
    print('Enter "help(qu4)" for more information')
    print('This is the definition of the alias "aq4": \
    fetch 3v6d, HIVrt, async=0; \
    preset.ball_and_stick("(c. C and i. 258) or (c. F and i. 817 )"); \
    set_view (0.48,0.84,-0.27,-0.81,0.3,-0.51,-0.35,0.46, \
    0.82,0.0,-0.0,-27.26,-40.77,-53.15,12.15,25.47,29.09,-20.0)')
cmd.extend( 'qu4',qu4)
cmd.alias('q4','qu4;aq4')
#print 'Enter "reinitialize" before going to next problem"


def qu5(question=5):
    '''
    DESCRIPTION

    Question 5:  How many axial stacks of helices does 
    the ribozyme have?
    '''
#    print question
    cmd.alias('aq5', 'reinitialize; fetch 3zp8, hammer,\
    async=0; show cartoon,hammer;rock;set_view \
    (-0.5,0.18,-0.85,-0.17,-0.98,-0.11,-0.85,0.09,0.52,0.0,0\
    .0,-167.2,-18.45,10.92,-12.11,126.37,208.02,-20.0);')
    print('Enter "q5" to set scene for question 5.')    
    print('Enter "help(qu5)" for more information')
    print('This is the definition of the alias "aq5": \
    fetch 3zp8, hammer, async=0; \
    show cartoon, c. A or c. B;rock;\
    set_view (0.53,-0.06,-0.84,0.82,-0.21,0.53,-0.21,\
    -0.98,-0.07,-0.0,-0.0,-192.99,-29.84,8.42,47.76,\
    178.27,207.7,-20.0)')
cmd.extend( 'qu5',qu5)
cmd.alias('q5','qu5;aq5')
#print 'Enter "reinitialize" before  before going to next
#problem"
#

def qu6(question=6):
    '''
    DESCRIPTION

    Question 6:  What is the average distance of the Na1044
   ligand bonds? Give the residue numbers of the RNA
   nucleotides and the sodium to identify them. How many
   ligands are from RNA?
    '''
#    print question
    cmd.alias('aq6','reinitialize;fetch 3zp8, hammer,async=0;\
    rock; hide cartoon;\
    distance ligand1, i. 1044, c. A and i. 22 and n. N7;\
    distance ligand2, i. 1044, c. A and i. 21 and n. OP2;\
    distance ligand3, i. 1044, i. 2121;\
    distance ligand4, i. 1044, i. 2120; \
    distance ligand5, i. 1044, i. 2122; \
    distance ligand6, i. 1044, i. 2130;\
    set_view (-0.87,0.18,-0.46,-0.39,-0.81,0.44,-0.29,\
    0.56,0.78,-0.0,0.0,-20.47,-18.05,14.02,-18.89,\
    17.47,23.47,-20.0);')   
    print 'Enter "q6" to set scene for question 6.'
    print 'Enter "help(qu6)" for more information' 
    print 'This is the definition of the alias "aq6": \
    fetch 3zp8, hammer, async=0; rock; hide cartoon;\
    distance ligand1, i. 1044, c. A and i. 22 and n. N7;\
    distance ligand2, i. 1044, c. A and i. 21 and n. OP2;\
    distance ligand3, i. 1044, i. 2121;\
    distance ligand4, i. 1044, i. 2120;\
    distance ligand5, i. 1044, i. 2122;\
    distance ligand6, i. 1044, i. 2130;\
    set_view (-0.87,0.18,-0.46,-0.39,-0.81,0.44,\
    -0.29,0.56,0.78,-0.0,0.0,-20.47,-18.05,14.02,\
    -18.89,17.47,23.47,-20.0)'
cmd.extend('qu6',qu6)

cmd.alias('q6','qu6;aq6')
#print 'Enter "reinitialize" before  before going to next
#problem"


def qu7(question=7):
    '''
    DESCRIPTION

    Question 7:  Measure the longest dimension and the
    shortest dimension of the ribozyme. What is the ratio of
    the longest dimension to the shortest dimension? Is it
    globular like a protein?
    '''
    
#    print question
    cmd.alias('aq7','reinitialize;fetch 3zp8, hammer,\
    async=0; show ribbon, hammer;rock;\
    set_view (0.62,0.14,0.78,0.13,-0.99,0.07,0.78,\
    0.05,-0.63,-0.0,-0.0,-169.8,-16.43,9.44,-9.63,\
    143.54,196.05,-20.0);')
    print 'Enter "q7" to set scene for question 7.'
    print 'Enter "help(qu7)" for more information' 
    print 'This is the definition of the alias "aq7":\
    reinitialize;fetch 3zp8, hammer, async=0;\
    show ribbon,3zp8;\
    set_view (0.62,0.14,0.78,0.13,-0.99,0.07,\
    0.78,0.05,-0.63,-0.0,-0.0,-169.8,-16.43,\
    9.44,-9.63,143.54,196.05,-20.0)'
cmd.extend('qu7',qu7)

cmd.alias('q7','qu7;aq7')
#print 'Enter "reinitialize" before  before going to next
#problem"


def qu8(question=8):
    '''
    DESCRIPTION

   Question 8:  Find the unusual base pair between A21 and
   G36. What is the length of the H-bonds between the bases
   (ignore the H atoms in the distance measurement)? List
   the distance with the residue name, residue number, and
   atom name. What additional H-bond occurs between a base
   in this base pair and a ribose ring of one of the two
   nucleotides in this base pair?   
    '''
#    print question
    cmd.alias('aq8','reinitialize;\
    fetch 3zp8, hammer,async=0;  rock; hide everything;\
    show ribbon; show sticks, resi 21 or resi 36;\
    set_view (-0.9,-0.19,0.39,0.39,-0.74,0.55,0.19,0.65,0.74,\
    0.0,0.0,-37.58,-21.66,15.71,-23.32,35.42,39.74,-20.0);')
    print('Enter "q8" to set scene for question 8.')
    print('Enter "help(qu8)" for more information')
    print('This is the definition of the alias "aq8": \
    fetch 3zp8, hammer, async=0; rock; hide everything; show ribbon;\
    show sticks, resi 21 or resi 36;\
    set_view (-0.9,-0.19,0.39,0.39,-0.74,0.55,0.19,0.65,\
    0.74,0.0,0.0,-37.58,-21.66,15.71,-23.32,35.42,39.74,-20.0)')
cmd.extend('qu8',qu8)

cmd.alias('q8','qu8;aq8')
#print 'Enter "reinitialize" before  before going to 
# next problem"
