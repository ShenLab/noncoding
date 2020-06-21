import numpy as np
import pandas as pd

#encode sequence into tensor
base_dict = {
    'A': np.array([1,0,0,0]),
    'C': np.array([0,1,0,0]),
    'G': np.array([0,0,1,0]),
    'T': np.array([0,0,0,1]),
    'N': np.array([0,0,0,0])
}

def encode_seq(sequences):
  n = len(sequences)
  a = np.empty((n,1000,4))
  for i in range(n):
    seq = sequences.iloc[i]
    seq_array = list(seq)
    for j in range(1000):
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
def get_formatted_data(X_df, Y_df):
  seq_array = encode_seq(X_df['sequence'])
  gene_array = X_df['lnTPM'].to_numpy()
  annot_array = encode_annot(X_df['annot'])
  label_array = Y_df['label'].to_numpy()
  return seq_array, gene_array, annot_array, label_array