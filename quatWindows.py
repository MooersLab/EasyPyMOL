'''
(c) 2010-2011 Thomas Holder, MPI for Developmental Biology
 
Module for reading REMARK records from PDB files and in particular
generate quaterny structure from REMARK 350.
'''
 
import sys, os
from pymol import cmd, stored
 
local_mirror_divided = '/mnt/bio/db/pdb.divided'
 
def pdbremarks(filename):
    '''
    Read REMARK lines from PDB file. Return dictionary with remarkNum as key
    and list of lines as value.
    '''
    remarks = dict()
    if not isinstance(filename, basestring):
        f = filename
    elif filename[-3:] == '.gz':
        import gzip
        f = gzip.open(filename)
    else:
        f = open(filename)
    for line in f:
        recname = line[0:6]
        if recname == 'REMARK':
            num = int(line[7:10])
            lstring = line[11:]
            remarks.setdefault(num, []).append(lstring)
    return remarks
 
def quat350(rem350):
    '''
    Get transformation matrices for biomolecule 1 from REMARK 350.
    '''
    biomt = dict()
    chains = tuple()
    seenbiomolecule = False
    for line in rem350:
        if line.startswith('BIOMOLECULE:'):
            if seenbiomolecule:
                break
            seenbiomolecule = True
        elif line.startswith('APPLY THE FOLLOWING TO CHAINS:'):
            chains = tuple(chain.strip() for chain in line[30:].split(','))
        elif line.startswith('                   AND CHAINS:'):
            chains += tuple(chain.strip() for chain in line[30:].split(','))
        elif line.startswith('  BIOMT'):
            row = int(line[7])
            num = int(line[8:12])
            vec = line[12:].split()
            vec = map(float, vec)
            biomt.setdefault(chains, dict()).setdefault(num, []).extend(vec)
    return biomt
 
def quat(name=None, filename=None, prefix=None, quiet=0):
    '''
DESCRIPTION
 
    Read REMARK 350 from `filename` and create biological unit
    (quaternary structure)
 
USAGE
 
    quat [name [, filename [, prefix]]]
 
ARGUMENTS
 
    name = string: name of object and basename of PDB file, if
    filename is not given {default: first loaded object}
 
    filename = string: file path {default: <name>.pdb}
 
    prefix = string: prefix for new objects {default: <name>}
 
EXAMPLE
 
    fetch 1rmv
    quat 1rmv
    '''
    quiet = int(quiet)
    if name is None:
        name = cmd.get_object_list()[0]
    if prefix is None:
        prefix = name
    if filename is None:
        candidates = [
            '%s.pdb' % (name),
            '%s/%s.pdb' % (cmd.get('fetch_path'), name),
            '%s/%s/pdb%s.ent.gz' % (local_mirror_divided, name[1:3], name),
        ]
        for filename in candidates:
            if os.path.exists(filename):
                break
        else:
            print 'please provide filename'
            return
        if not quiet:
            print 'loading from %s' % (filename)
    remarks = pdbremarks(filename)
    if 350 not in remarks:
        print 'There is no REMARK 350 in', filename
        return
    quat = quat350(remarks[350])
    for chains in quat:
        matrices = quat[chains]
        for num in matrices:
            mat = matrices[num][0:12]
            mat.extend([0,0,0,1])
            copy = '%s_%d' % (prefix, num)
            if not quiet:
                print 'creating %s' % (copy)
            cmd.create(copy, '/%s//%s' % (name, '+'.join(chains)))
            cmd.alter(copy, 'segi="%d"' % (num))
            cmd.transform_object(copy, mat)
    cmd.disable(name)
    cmd.group('%s_quat' % (prefix), '%s_*' % (prefix))
 
cmd.extend('quat', quat)
 
# vi:expandtab:smarttab
