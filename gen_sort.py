'''
Sorts smiles and .pdbqt files against a list.

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

import argparse
import os
import sys
import shutil
import pandas as pd
from tqdm import tqdm


def pdbqt_sort(stuff, dr):
    dest = './next_gen'
    os.mkdir(dest)
    sort_dir = '{0}'.format(dr)

    inlist = '{0}'.format(stuff)
    with open(inlist, 'r') as f:
        data = [r.rstrip() for r in f]

    for i in tqdm(range(len(data))):
        temp = '{0}.pdbqt'.format(data[i])
        for f_name in os.listdir(sort_dir):
            name = '{0}'.format(f_name)
            if name == temp:
                src = os.path.join(sort_dir, name)
                shutil.copy(src, dest)
    print('Sort Completed!')


def smile_sort(stuff, smi_map):
    inlist = '{0}'.format(stuff)
    with open(inlist, 'r') as f:
        cklist = [r.rstrip() for r in f]

    smi_list = './{0}'.format(smi_map)
    dataset = pd.read_csv(smi_list, engine='python')

    smiles =[]
    for i in tqdm(range(len(cklist))):
        temp = '{0}'.format(cklist[i])
        for j in range(dataset.shape[0]):
            name = '{0}'.format(dataset.Name[j])
            if name == temp:
                smi = dataset.SMILE[j]
                smiles.append(smi)
                continue

    outfile = './nextgen_smiles.csv'
    with open(outfile, 'w') as f:
        for smi in smiles:
            f.write(smi + '\n')


parser = argparse.ArgumentParser(description='Sorts smiles and .pdbqt files against a list')
parser.add_argument('-p', '--pdbqt_sort', nargs='+', help='-p list_name source_directory')
parser.add_argument('-s', '--smile_sort', nargs='+', help='-s list_name smiles_map_file')
args = parser.parse_args()


def main():
    if args.pdbqt_sort:
        pdbqt_sort(sys.argv[2], sys.argv[3])
    elif args.smile_sort:
        smile_sort(sys.argv[2], sys.argv[3])


if __name__ == '__main__':
    main()
