'''
Program to clean the input SMILES data

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
from rdkit.Chem import MolStandardize


RDLogger.DisableLog('rdApp.*')


class MolClean(object):
    def __init__(self):
        self.normizer = MolStandardize.normalize.Normalizer()
        self.lfc = MolStandardize.fragment.LargestFragmentChooser()
        self.uc = MolStandardize.charge.Uncharger()

    def clean(self, smi):
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mol = self.normizer.normalize(mol)
            mol = self.lfc.choose(mol)
            mol = self.uc.uncharge(mol)
            smi = Chem.MolToSmiles(mol,  isomericSmiles=False, canonical=True)
            return smi
        else:
            return None


def main(input_file, output_file, **kwargs):
    mc = MolClean()

    with open(input_file, 'r') as f:
        smiles = [r.rstrip() for r in f]

    smiles = [x[:-6:] for x in smiles]

    print(f'input SMILES num: {len(smiles)}')
    print('Clean Up Started ...')

    mc_smiles = [mc.clean(smi) for smi in tqdm(smiles)]

    print('Cleanup Completed')
    print(f'output SMILES num: {len(mc_smiles)}')

    with open(output_file, 'w') as f:
        for smi in mc_smiles:
            f.write(smi + '\n')

    return


if __name__ == '__main__':
    print('Standardize the SMILES - remove salts, charges and stereochemical information')
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='input file')
    parser.add_argument('output', help='output file')

    args = parser.parse_args()
    main(args.input, args.output)
