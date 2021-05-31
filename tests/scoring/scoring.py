#Written by Nitin Murthy
#Demo file to show how to use the scoring functions that I've encapsulated. Currently supported are SAScore and SCScore.

import pandas as pd
from rdkit import Chem
import os
from chemlg.sascore import sascore
from chemlg.scscore import scscore

#Demo section
direct = os.getcwd()
# print(direct)
#Tip for the read_csv function, use r'directory\csvfilename.csv' as the argument. Just putting in the csv file name doesn't work for some reason.
df1 = pd.read_csv(r'C:\Users\nitin\Documents\UB\Hachmann Group\Thesis\Retro_synth\integration_lg\tests\scoring\final_smiles.csv')
# print(df1.head(5))

smi = df1['SMILES']
mols = [Chem.MolFromSmiles(s) for s in smi] #SAScore needs a mol input
# print(mols[:5])
#
df_sa = sascore(mols)
#
df_sc = scscore(smi)
