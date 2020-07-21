'''
Converts a ligand .pdb file to a .pdbqt file for AutoDock Vina

--------------------------------------------------------------------------------
MIT License

Copyright (c) 2020 Brent Larzalere

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

'''

import sys
import argparse
from pymol.cgo import *
import os


def process_pdb(pdb_input, name):
    # loads the .pdb file in pymol then adds the polar hydrogens and temporarily
    # saves the pdb file; it then converts the .pdb file to a .pdbqt file while
    # preserving the hydrogens
    cmd.load(pdb_input)
    cmd.h_add(selection='acceptors or donors')
    nm = '{0}.pdb'.format(name)
    cmd.save(nm)
    one = 'obabel {0}.pdb -O {1}.pdbqt -xh'.format(name, name)
    os.system(one)
    os.remove(nm)


parser = argparse.ArgumentParser(description='Prep a .pdb ligand for AutoDock Vina')
parser.add_argument('-p', '--process_pdb', nargs='+', help='pdb_file ligand_name')
args = parser.parse_args()


def main():
    if args.process_pdb:
        process_pdb(sys.argv[2], sys.argv[3])


if __name__ == '__main__':
    main()
