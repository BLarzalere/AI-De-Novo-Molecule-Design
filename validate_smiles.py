'''
Program that validates and outputs only valid RDKIT compatible molecules.

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
from tqdm import tqdm
from rdkit import Chem, RDLogger


RDLogger.DisableLog('rdApp.*')


def main(input_file, output_file):

    # read in the input file
    with open(input_file, 'r') as f:
        smiles = [r.rstrip() for r in f]

    print(f'Input SMILES num: {len(smiles)}')
    print('Performing SMILES validation ...')

    smiles_val, boom = [], []
    for i in tqdm(range(len(smiles))):
        smi = smiles[i]
        mol = Chem.MolFromSmiles(smi)
        if mol:  # check if the smile is a valid molecule
            boom.append(smi)
        else:
            continue
    smiles_val.extend(boom)

    # save smiles output file
    print("Creating SMILES output file ...")
    with open(output_file, 'w') as f:
        for smi in smiles_val:
            f.write(smi + '\n')

    print(f'Output SMILES num: {len(smiles_val)}')
    print("Smiles Validation Complete")

    return


if __name__ == '__main__':
    print('Validate that the SMILES are RDKIT valid molecules')
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='Input file name')
    parser.add_argument('output_file', help='Output file name')

    args = parser.parse_args()
    main(args.input_file, args.output_file)
