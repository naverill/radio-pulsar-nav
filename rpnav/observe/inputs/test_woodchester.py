from pathlib import Path
import csv
import os

import astropy.units as u

import pandas as pd

FILE_DIR = Path(__file__).parent
VELA_P = 89.33 * u.ms

def test_vela():
    fpath = f"{FILE_DIR}/"

    strong_df = pd.read_csv(fpath + "Vela_2p3_Strong_Pol_171024.csv", index_col=0)

    weak_df = pd.read_csv(fpath + "Vela_2p3_Weak_Pol_171024.csv", index_col=0)

    print(strong_df["S/N"].mean())    
    print((strong_df["Delta f0/f0"].multiply(strong_df["Best BCP"]) * u.ms).mean())    

    print(weak_df["SNR"].mean())
    print()
    print((weak_df["delta F0/F0"].multiply(weak_df["PP (B)"]) * u.ms).mean())    
    # print((weak_df["delta F0/F0"].multiply(weak_df["PP (B)"]) * u.ms).mean())    
    # print((weak_df["delta F0/F0"].multiply(weak_df["PP (B)"]).mean() * u.ms).to(u.second))    
