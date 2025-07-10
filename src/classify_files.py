import compute_r
import classify_r
import pandas as pd

def classify_files(df):
    ran_files = glob.glob('data/coord/*_rand.ecsv')
    dat_files = glob.glob('data/coord/*_data.ecsv')
    df = compute_r.compute_r(df)
    df = classify_r.classify_r(df)
    return df