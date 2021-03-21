#!/usr/bin/env python

import sys,os
import argparse
import numpy

atom_dtype = numpy.dtype( [
    ("atomname", numpy.unicode_, 4),
    ("resname", numpy.unicode_, 3),
    ("chnid", numpy.unicode_, 1),
    ("resid", numpy.int32),
    ("X", numpy.float32, 3),
] )

def parse_pdb(pdbfile):
    allatoms = []
    for line in pdbfile:
        if line[:4] == 'ATOM' or line[:6] == "HETATM":
            split_line = (
                line[12:16], line[17:20], line[21], line[22:26], 
                (line[30:38], line[38:46], line[46:54])
            )
            allatoms.append(split_line)
    return (numpy.array(allatoms, dtype=atom_dtype))


# main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Symmetric alignment from pointgroups.')
    parser.add_argument('--pdb', help='input PDB file containing ASU only', required=True)
    parser.add_argument('--residues',  help='residues to constrain', type=str, default='')
    parser.add_argument('--width',  help='Constraint width', type=float, default=1.0)
    parser.add_argument('--weight',  help='Constraint width', type=float, default=1.0)
    parser.add_argument('--atoms',  help='Atoms to constrain', default="N,CA,C")
    parser.add_argument('--vrt',  help='vrtid', type=int, default="0")

    args = parser.parse_args()

    # residues
    constrain = dict()
    resranges = args.residues.split(',')
    for resrange in resranges:
        res = resrange.split('-')
        if (len(res)==1):
            constrain[int(res[0])] = 1
        elif (len(res)==2):
            for i in range(int(res[0]),int(res[1])+1):
                constrain[i] = 1
        else:
            assert False, "Error parsing "+resrange

    # pdb
    with open(args.pdb) as pdbfile:
        pdb = parse_pdb( pdbfile )

    atoms = args.atoms.split(',')

    stripnames = [a.strip() for a in pdb['atomname']]
    mask = numpy.array([i in atoms for i in stripnames])
    pdb = pdb[mask]

    # atom pair constrain
    asu = pdb[pdb['chnid']=='A']

    for i in range(len(asu)):
        a_i = asu[i]
        if (a_i['resid'] in constrain and a_i['atomname'].strip() in args.atoms):
            x_i = a_i['X']
            weight = args.weight

            if (weight != 0):
                print ("CoordinateConstraint",a_i['atomname'],a_i['resid'],"ORIG",args.vrt,x_i[0],x_i[1],x_i[2], 
                    "SCALARWEIGHTEDFUNC",weight,"HARMONIC",0.0,args.width)
