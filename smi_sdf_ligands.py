'''
Converts SMILES to SDF files with Excel mapping output

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
import pandas as pd
from pymol.cgo import *
from rdkit import Chem, RDLogger
import rdkit.Chem.PropertyMol

RDLogger.DisableLog('rdApp.*')


def convert_smi_to_mol(smi_list):
    # convert our smiles to mols
    smiles = pd.read_csv(smi_list, sep='\n', header=None)
    mol_list = []
    for i in range(smiles.shape[0]):
        smi = smiles[0][i]
        mol = Chem.MolFromSmiles(smi)
        mol_list.append(mol)
    return mol_list


def mol_props(mols, ext):
    # add the name property to each mol and create a pandas dataframe associating the names to the smiles
    i = 1
    mols_to_export, smile_export = [], []
    for mol in mols:
        name = str(ext) + '-' + str(i)
        pm = Chem.PropertyMol.PropertyMol(mol)
        pm.SetProp('Name', name)
        mols_to_export.append(pm)

        mol_dict = {}
        mol_dict['Name'] = name
        smile = Chem.MolToSmiles(mol)
        mol_dict['SMILE'] = smile
        smile_export.append(mol_dict)
        i += 1

    df = pd.DataFrame(smile_export)
    return df, mols_to_export


def write_sdf(mols, output):
    # convert our mols to .sdf files
    sdf_dir = '{0}/'.format(output)
    os.mkdir(sdf_dir)
    for m in mols:
        pm = Chem.PropertyMol.PropertyMol(m)
        name = pm.GetProp('Name')
        w = Chem.SDWriter(str(sdf_dir) + str(name) + '.sdf')
        w.write(m)


def smi_to_sdf(smi_list, ext, df_name, output):
    # convert smiles to mols
    mols = convert_smi_to_mol(smi_list)

    # add the name property to mols
    smiles_df, mol_list = mol_props(mols, str(ext))

    # save the name-smile dataframe to excel
    smiles_df.to_excel(str(df_name) + '.xlsx')

    # create the .sdf files form our mols
    write_sdf(mol_list, output)
    print('Done creating SDF files')


parser = argparse.ArgumentParser(description='Prep ligands for AutoDock Vina')
parser.add_argument('-s', '--smi_to_sdf', nargs='+', help='smile_list ligand_name exel_file_name output_directory_name')
args = parser.parse_args()


def main():
    if args.smi_to_sdf:
        smi_to_sdf(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])


if __name__ == '__main__':
    main()
