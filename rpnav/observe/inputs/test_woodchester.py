from pathlib import Path
import csv
import os
import numpy as np

import astropy.units as u

import pandas as pd

FILE_DIR = Path(__file__).parent
VELA_P = 89.33 * u.ms

def test_vela():
    fpath = f"{FILE_DIR}/"

    strong_df = pd.read_csv(fpath + "Vela_2p3_Strong_Pol_171024.csv", index_col=0)

    weak_df = pd.read_csv(fpath + "Vela_2p3_Weak_Pol_171024.csv", index_col=0)

    print(strong_df["S/N"].mean())
    # dt_strong = (abs(strong_df["BCP Error"]).multiply(strong_df["Best BCP"]) * u.ms) 
    dt_strong = (abs(strong_df["BCP Error"]) * u.ms) 
    print(dt_strong.max())    
    print(dt_strong.min())    

    dt_weak = (abs(weak_df["PP (B) Err."]) * u.ms)
    print(weak_df["SNR"].mean())
    print(dt_weak.max())   
    print(dt_weak.min())   
