import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import roc_curve, auc, plot_roc_curve
import statistics
import os.path
import csv

import keras
from keras import optimizers
from keras import backend as tf

# custom modules
import cnn_models
import model_utils

cell_lines = ['HepG2','K562']
max_len = 1000
num_channels = 4 # A G C T

dilated_out_file = '5cv_complex_dilated_auc.txt'
regular_out_file = '5cv_complex_auc.txt'

for cell_line in cell_lines:
  data_file = '../complexes/{}_complex_data.csv'.format(cell_line)
  complex_file = '../complexes/{}_complexes.txt'.format(cell_line)

  data = pd.read_csv(data_file)
  data = data.fillna('no_annot')
  data = data.drop(['seqnames','start','end','strand'], axis=1)

  with open(complex_file,'r') as f:
    complex_list = f.readlines()
  complex_list = [x.split('\t') for x in complex_list]
  for cpx in complex_list:
    if cpx[-1][-1] == '\n':
      cpx[-1] = cpx[-1][:-1]

  for test_complex in complex_list:
    print(test_complex)
    df_complex = data[data['RBP'].isin(test_complex)]

    std_cols = ['RBP','sequence','lnTPM','annot']
    complex_cols = [RBP+'_label' for RBP in test_complex]
    keep_cols = std_cols + complex_cols
    df_complex = df_complex[keep_cols]

    complex_len = len(test_complex)
    kf = KFold(n_splits=5, shuffle=True, random_state=13)
    dilated_auc_list = []
    regular_auc_list = []
    for i in range(complex_len):
      dilated_auc_list.append([])
      regular_auc_list.append([])

    for i, (train_index, test_index) in enumerate(kf.split(df_complex)):
      train_df = df_complex.iloc[train_index]
      test_df = df_complex.iloc[test_index]
      train_df, valid_df = train_test_split(train_df, test_size=0.1, shuffle=True, random_state=21)

      #balance training and validation data
      train_df = model_utils.balance(train_df)
      valid_df = model_utils.balance(valid_df)

      #get data in correct format
      trainX, trainY = train_df[['sequence','lnTPM','annot']], train_df[complex_cols]
      validX, validY = valid_df[['sequence','lnTPM','annot']], valid_df[complex_cols]
      testX, testY = test_df[['sequence','lnTPM','annot']], test_df[complex_cols]

      trainX, trainG, trainA, trainY = model_utils.get_formatted_data(trainX, trainY, max_len)
      validX, validG, validA, validY = model_utils.get_formatted_data(validX, validY, max_len)
      testX, testG, testA, testY = model_utils.get_formatted_data(testX, testY, max_len)

      #train dilated model
      dilated_model = cnn_models.get_dilated_model((max_len, num_channels), len(complex_cols))
      dilated_model.compile(optimizer='sgd', loss='binary_crossentropy')
      dilated_model.fit([trainX, trainG], trainY,
        batch_size=16,
        epochs=16,
        verbose=0,
        validation_data=([validX, validG], validY))
      dilated_predictY = dilated_model.predict([testX, testG])

      for i in range(complex_len):
        fpr, tpr, thres = roc_curve(testY[:,i], dilated_predictY[:,i])
        roc_auc = auc(fpr, tpr)
        dilated_auc_list[i].append(roc_auc)

      #train regular model
      regular_model = cnn_models.get_gene_annot_model((max_len, num_channels), len(complex_cols))
      regular_model.compile(optimizer='sgd', loss='binary_crossentropy')
      regular_model.fit([trainX, trainG, trainA], trainY,
        batch_size=16,
        epochs=16,
        verbose=0,
        validation_data=([validX, validG, validA], validY))
      regular_predictY = regular_model.predict([testX, testG, testA])

      for i in range(complex_len):
        fpr, tpr, thres = roc_curve(testY[:,i], regular_predictY[:,i])
        roc_auc = auc(fpr, tpr)
        regular_auc_list[i].append(roc_auc)

    dilated_auc_means = [str(np.mean(l)) for l in dilated_auc_list]
    regular_auc_means = [str(np.mean(l)) for l in regular_auc_list]

    with open(dilated_out_file, 'a') as f:
      f.write('\t'.join(test_complex))
      f.write('\n')
      f.write('\t'.join(dilated_auc_means))
      f.write('\n')

    with open(regular_out_file, 'a') as f:
      f.write('\t'.join(test_complex))
      f.write('\n')
      f.write('\t'.join(regular_auc_means))
      f.write('\n')
