'''
Performs pre-processing of SMILES data

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
import json
from sklearn.model_selection import train_test_split
from tqdm import tqdm
from rdkit import Chem, RDLogger
from SmilesEnumerator import SmilesEnumerator

RDLogger.DisableLog('rdApp.*')


class PreProcessSmiles(object):

    def enumerate_smile(self, input_smiles):
        # extract the character set and add beginning / ending characters '!' and 'E'
        charset = set("".join(list(input_smiles)) + "\!E@/")
        n_vocab = len(charset)
        max_mol = max([len(x) for x in input_smiles])
        print("Number of unique smiles characters:", n_vocab)
        print("Max molecule length:", max_mol)
        print("Smiles character set:\n", str(charset))

        # create Python dictionaries for mapping characters-to-integers and integers-to-characters
        char_to_int = dict((c, i) for i, c in enumerate(charset))
        int_to_char = dict((i, c) for i, c in enumerate(charset))

        # save the Python dictionaries
        with open('char_to_int.json', 'w') as d:
            json.dump(char_to_int, d)

        with open('int_to_char.json', 'w') as d2:
            json.dump(int_to_char, d2)

        return

    # perform data augmentation of the smiles sequences
    def smile_augment(self, input_smiles, times):
        print('Begin SMILES augmentation ...')
        sme = SmilesEnumerator(pad=0, isomericSmiles=False)

        smiles_aug = []
        times = int(times)
        for i in tqdm(range(len(input_smiles))):
            smi = input_smiles[i]
            boom = []
            for j in range(times):
                temp = sme.randomize_smiles(smi)
                mol = Chem.MolFromSmiles(temp)
                if mol:  # check if the randomized smile is a valid molecule
                    boom.append(temp)
                else:
                    continue
            smiles_aug.extend(boom)

        print("Smiles augmentation completed")
        print("Length of augmented smiles set:", len(smiles_aug))
        print("Sample augmented SMILES:\n", smiles_aug[2:10])

        return smiles_aug


def main(input_file, output_file, **kwargs):

    test_file = args.test_file
    times = args.times

    # read in the input file
    with open(input_file, 'r') as f:
        smiles = [r.rstrip() for r in f]

    print(f'Input SMILES num: {len(smiles)}')

    pp_smi = PreProcessSmiles()

    # extract the smile character set and create char-to-int / int-to-char mapping
    pp_smi.enumerate_smile(smiles)

    # split into train/test sets
    if test_file:
        print("Splitting into Train & Test sets ...")
        smiles_train, smiles_test = train_test_split(smiles, test_size=0.2)
        print("Training set size:", len(smiles_train))
        print("Test set size:", len(smiles_test))

        # perform data augmentation on the smiles sequences
        if times:
            print("Performing data augmentation on training dataset ....")
            smiles_train = pp_smi.smile_augment(smiles_train, times)

        # save training and test sets
        print("Creating SMILES output files ...")
        with open(output_file, 'w') as f:
            for smi in smiles_train:
                f.write(smi + '\n')

        with open(test_file, 'w') as f:
            for smi in smiles_test:
                f.write(smi + '\n')

    else:
        # perform data augmentation on the smiles sequences
        if times:
            print("Performing data augmentation ....")
            smiles = pp_smi.smile_augment(smiles, times)

        # save smiles output file
        print("Creating SMILES output file ...")
        with open(output_file, 'w') as f:
            for smi in smiles:
                f.write(smi + '\n')

        print(f'Output SMILES num: {len(smiles)}')

    print("Pre-Processing Complete")

    return


if __name__ == '__main__':
    print('Preprocess the SMILES and optionally split into train & test sets')
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='Input file name')
    parser.add_argument('output_file', help='Training set or output file name')
    parser.add_argument('--test_file', help='Test set output file name if splitting into train/test sets (optional)')
    parser.add_argument('--times', help='Number times the training smiles are augmented (optional)')

    args = parser.parse_args()
    main(args.input_file, args.output_file)
