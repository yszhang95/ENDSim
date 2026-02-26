import ROOT
import pandas as pd
import numpy as np

t_signal = ROOT.RDataFrame("doms", "../c++/pmt_bdt_signal_numu_100000_rootracker_randomVtx_aframe_spacing5m_hex.root")

df_signal = pd.DataFrame(t_signal.AsNumpy())

colagg = {}
info_cols = ['rock_wgt', 'vtxX', 'vtxY', 'vtxZ', 'momX', 'momY', 'momZ', 'muE']
for col in df_signal.columns:
    if col == "event_id":
        continue
    if col not in info_cols:
        colagg[col] = list
    else:
        colagg[col] = np.mean

df_signal = df_signal.groupby("event_id").agg(colagg).reset_index()

max_len = 2790  # hard coded

def pad_trunc(l, length, pad_value=0):
    if len(l) >= length:
        return l[:length]
    return l + [pad_value] * (length - len(l))

# Apply to all columns except event_id and various truth info columns
for col in df_signal.columns:
    if col == "event_id" or col in info_cols:
        continue
    df_signal[col] = df_signal[col].apply(lambda l: pad_trunc(l, max_len, pad_value=-1))


def df_to_numpy(df):
    arrays = []
    for c in df.columns:
        if c == 'event_id' or c in info_cols:
            continue
        if type(df[c][0]) is not list:
            arrays.append(np.array(df[c].to_list()).reshape(-1, 1))
        else:
            arrays.append(np.array(df[c].to_list()))
    return np.hstack(arrays)

X_signal = df_to_numpy(df_signal)  # shape (n_signal,  L * n_cols)
df_signal_info = df_signal[info_cols]

print(X_signal.shape)

np.save('signal_aframe_spacing5m_hex_DU_v2.npy', X_signal)
df_signal_info.to_csv('signal_aframe_spacing5m_hex_DU_v2_info.csv', index=False)

with open('info_aframe_spacing5m_hex_DU_v2.txt', 'w') as info:
    info.write('Columns : %s'%str(list(df_signal.columns)))
    info.write('\nArray Length: %d'%max_len)
