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
    parser.add_argument('--seqsep',  help='min sequence separation', type=int, default=0)
    parser.add_argument('--residues',  help='residues to constrain', type=str, default='')
    parser.add_argument('--toall', dest='toall', action='store_true')
    parser.add_argument('--maxdis',  help='maximum distance', type=float, default=9.0)
    parser.add_argument('--width',  help='Constraint width', type=float, default=1.0)
    parser.add_argument('--weight',  help='Constraint width', type=float, default=1.0)
    parser.add_argument('--wt_helix',  help='Constraint width', type=float, default=2.0)
    parser.add_argument('--atoms',  help='Atoms to constrain', default="N,CA,C")

    args = parser.parse_args()

    # pdb
    with open(args.pdb) as pdbfile:
        pdb = parse_pdb( pdbfile )

    # residues
    constrain = dict()
    if (args.residues != ''):
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

    atoms = args.atoms.split(',')

    stripnames = [a.strip() for a in pdb['atomname']]
    mask = numpy.array([i in atoms for i in stripnames])
    pdb = pdb[mask]

    asu = pdb[pdb['chnid']=='A']
    if (len(constrain) > 0):
        asu = numpy.array(
            [x for x in asu if x['resid'] in constrain], dtype=atom_dtype)
        if (not args.toall):
            pdb = numpy.array(
                [x for x in pdb if x['resid'] in constrain], dtype=atom_dtype)

    atomdists = numpy.linalg.norm( asu['X'][:,None,:] - pdb['X'][None,:,:], axis=-1 )
    neighbors = (atomdists < args.maxdis).nonzero()
    mask = asu[neighbors[0]]['resid'] < pdb[neighbors[1]]['resid']-args.seqsep
    neighbors = (neighbors[0][mask],neighbors[1][mask])

    ncsts = len(neighbors[0])
    for i in range(ncsts):
        a_i = asu[neighbors[0][i]]
        a_j = pdb[neighbors[1][i]]
        d_ij = atomdists[neighbors[0][i], neighbors[1][i]]
        weight = args.weight

        if (weight != 0):
            print ("AtomPair",a_i['atomname'],a_i['resid'],a_j['atomname'],a_j['resid'], 
                "SCALARWEIGHTEDFUNC",weight,"HARMONIC",d_ij,args.width)
                #"BOUNDED",d_ij-0.1,d_ij+0.1,args.width,'X')
