# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 15:17:31 2020

@author: uni21
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from mordred import Calculator,descriptors
import numpy as np
import pandas as pd
import math
# データ可視化ライブラリ
import matplotlib.pyplot as plt

df = pd.read_csv(r'C:\Users\uni21\OneDrive\デスクトップ\Data science\Takeda\sample-data-master\sample-data-master\Materials Project 2.csv')
df = df[df.columns[1:5]]
smiles = df["Smiles"]
df.head()

mols = [Chem.MolFromSmiles(smile) for smile in smiles]
mols

# 記述子を計算
calc = Calculator(descriptors, ignore_3D=False)
df_descriptors_mordred = calc.pandas(mols)

df_descriptors = df_descriptors_mordred.astype(str)
masks = df_descriptors.apply(lambda d: d.str.contains('[a-zA-Z]' ,na=True))
df_descriptors = df_descriptors[~masks]
df_descriptors = df_descriptors.astype(float)
df_descriptors = df_descriptors.fillna("NA")
df_descriptors.head()

df_summary = pd.concat([df["Smiles"], df["EA (eV)"], df_descriptors],axis=1)

df_summary.head(1)

df_summary.to_csv(r'C:\Users\uni21\OneDrive\デスクトップ\Data science\Takeda\sample-data-master\sample-data-master\Materials Project 3.csv', index = False)
