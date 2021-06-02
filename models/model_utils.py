import numpy as np
import pandas as pd
import random

#encode sequence into input tensor
base_dict = {
  'A': np.array([1,0,0,0]),
  'C': np.array([0,1,0,0]),
  'G': np.array([0,0,1,0]),
  'T': np.array([0,0,0,1]),
  'N': np.array([0,0,0,0])
}

def encode_seq(sequences, max_len):
  n = len(sequences)
  a = np.empty((n,max_len,4))
  for i in range(n):
    seq = sequences.iloc[i]
    pad = max_len - len(seq)
    seq += pad*'N'
    seq_array = list(seq)
    for j in range(max_len):
      base = seq_array[j]
      a[i][j] = base_dict[base]
  return a

#one hot encode annotation types
annot_dict = {
  'hg19_genes_introns': np.array([1,0,0,0,0,0]),
  'hg19_genes_exons': np.array([0,1,0,0,0,0]),
  'hg19_genes_1to5kb': np.array([0,0,1,0,0,0]),
  'hg19_genes_promoters': np.array([0,0,0,1,0,0]),
  'hg19_genes_3UTRs': np.array([0,0,0,0,1,0]),
  'hg19_genes_5UTRs': np.array([0,0,0,0,0,1]),
  'no_annot': np.array([0,0,0,0,0,0])
}

def encode_annot(annotations):
  n = len(annotations)
  a = np.empty((n,6))
  for i in range(n):
    a[i] = annot_dict[annotations.iloc[i]]
  return a

#format data correctly to input into model
def get_formatted_data(X_df, Y_df, max_len):
  seq_array = encode_seq(X_df['sequence'], max_len)
  gene_array = X_df['lnTPM'].to_numpy()
  annot_array = encode_annot(X_df['annot'])
  label_array = Y_df['label'].to_numpy()
  return seq_array, gene_array, annot_array, label_array

#balance traning dataset across RBPs for complex tests
def balance(df_all):
  #find RBPs in complex, get list of RBP_dfs, and pick out longest one
  RBP_list = list(set(df_all['RBP']))
  RBP_df_list = []
  max_size = 0
  for RBP in RBP_list:
    df_RBP = df_all[df_all['RBP'] == RBP]
    if len(df_RBP) > max_size:
      max_size = len(df_RBP)
    RBP_df_list.append(df_RBP)

  #balance data for RBP
  for RBP in RBP_list:
    df_RBP = df_all[df_all['RBP'] == RBP]
    RBP_len = len(df_RBP)
    if RBP_len == max_size:
      continue

  #first get integer multiple of smaller df
  integer_mult = max_size // RBP_len
  df_extra = pd.concat([df_RBP for i in range(integer_mult)], ignore_index=True)
  
  #randomly sample the rest to get exact data size match
  remainder = max_size % RBP_len
  df_more = df_RBP.iloc[random.sample(range(0,RBP_len), remainder)]
  
  #combine to get full new dataset for RBP
  df_new = pd.concat([df_extra, df_more])
  RBP_ix = RBP_list.index(RBP)
  RBP_df_list[RBP_ix] = df_new
    
  return pd.concat(RBP_df_list)
