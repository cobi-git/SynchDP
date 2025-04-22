import glob
import pickle
import pandas as pd
import numpy as np

# ==================================================
# Generating the query_set.pkl input file to SynchDP
# ==================================================

def generate_covid_dict(files):
    df_dict = {}
    for file in files:
        columns = pd.read_csv(file, nrows=1).columns
        df = pd.read_csv(file, usecols=columns)

        df['DATE'] = pd.to_datetime(df['DATE'])
        first_date = df['DATE'].iloc[0]
        df['DATE_DIFF'] = (df['DATE'] - first_date).dt.days

        df = df.loc[:, ['NEWS', 'DATE_DIFF']]
        df.name = file.split('/')[-2]
        df.set_index('DATE_DIFF', inplace=True)

        merge_df = pd.DataFrame(df).T
        df_dict[df.name] = merge_df
    return df_dict

# file_path = './patient_data/*/SeverityScore.csv'
file_path = '/data1/home/inukj/mysoftware/synchdp/tutorial/pair_data/*/SeverityScore.csv'
file_list = glob.glob(file_path)
df_dict = generate_covid_dict(file_list)
selected_sample = list(pd.read_csv('./input/omics_info.csv', index_col=0).index)
selected_sample.sort()

new_data = {key: df_dict[key] for key in selected_sample if key in df_dict}

with open('./query_set.pkl','wb') as fw:
    pickle.dump(new_data, fw)


# ================================
# Generate reference_seq_set.pkl
# ================================

def continuous_norms(params_list, name):
    x_total = []
    y_total = []

    for params in params_list:
        mu, sigma, mul, cnt_per_dist = params
        x = np.linspace(-4, 4, cnt_per_dist)
        y = (1 / np.sqrt(2 * np.pi * sigma**2)) * np.exp(-(x - mu)**2 / (2 * sigma**2)) * mul

        if x_total:
            x = x[1:]
            y = y[1:]

        x_total.extend(x)
        y_total.extend(y)

    y_total = [num * 10 for num in y_total]

    seq_df = pd.DataFrame([y_total], index=['score'])
    seq_df.columns.name = 'DATE_DIFF'
    seq_df.name = name

    return seq_df

params_list_sev = [(0, 1, 1.75, 31)]
params_list_nor = [(0, 1, 0.75, 31)]

sev_seq_df = continuous_norms(params_list_sev, 'Rs')
nor_seq_df = continuous_norms(params_list_nor, 'Rm')

ref_seq = {
    'Rs': sev_seq_df,
    'Rm': nor_seq_df}

with open('reference_seq_set.pkl', 'wb') as f:
    pickle.dump(ref_seq, f)