import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import roc_curve, auc, plot_roc_curve
from statistics import stdev
import os.path
import csv

import keras
from keras import optimizers
from keras import backend as tf

# custom modules
import cnn_models
import model_utils
import rnn_models

with open('../RBP_list.txt', 'r') as f:
    RBP_list = f.readlines()

if RBP_list[-1][-1] != '\n':
    RBP_list[-1] = RBP_list[-1] + '\n'
RBP_list = [x[:-1] for x in RBP_list]

print('reading data...')
pos = pd.read_csv('../pos_rbp_seqs_1kb_SEQ.csv')
neg = pd.read_csv('../neg_rbp_seqs_1kb_SEQ.csv')

print('done reading data!')
pos = pos[pos['RBP'].isin(RBP_list)]
neg = neg[neg['RBP'].isin(RBP_list)]

#fill in annotation if empty
pos = pos.fillna('no_annot')
neg = neg.fillna('no_annot')

#add pos/neg labels to sequence data
pos['label'] = 1
neg['label'] = 0

#keep only relevant columns of df
pos = pos[['RBP','sequence','lnTPM','annot','label']]
neg = neg[['RBP','sequence','lnTPM','annot','label']]

auc_out_file = '5cv_cnn_gene_only_auc.csv'
max_len = 1000
num_channel = 4
model_out_dir = 'saved_gene_models'
if not os.path.exists(model_out_dir):
  os.makedirs(model_out_dir)

for RBP in RBP_list:
  print(RBP)
  df_pos = pos[pos['RBP'] == RBP]
  df_neg = neg[neg['RBP'] == RBP]
  df_all = pd.concat([df_pos, df_neg], ignore_index=True)

  X = df_all[['sequence','lnTPM','annot']]
  y = df_all[['label']]
  kf = KFold(n_splits=5, shuffle=True, random_state=13)

  auc_list = []
  max_auc = 0
  model_file_name = '{}/{}_binding_model.h5'.format(model_out_dir, RBP)
  for i, (train_index, test_index) in enumerate(kf.split(X, y)):
    #create train and test sets
    trainX, testX = X.iloc[train_index], X.iloc[test_index]
    trainY, testY = y.iloc[train_index], y.iloc[test_index]

    #create validation set
    (trainX, validX, trainY, validY) = train_test_split(trainX, trainY, test_size=0.1, random_state=21, shuffle=True)

    #get data in correct format
    trainX, trainG, trainA, trainY = model_utils.get_formatted_data(trainX, trainY, max_len)
    validX, validG, validA, validY = model_utils.get_formatted_data(validX, validY, max_len)
    testX, testG, testA, testY = model_utils.get_formatted_data(testX, testY, max_len)

    model = cnn_models.get_gene_model((max_len, num_channel), 1)
    model.compile(optimizer='sgd', loss='binary_crossentropy')
    model.fit([trainX, trainG], trainY,
      batch_size=16,
      epochs=16,
      verbose=0,
      validation_data=([validX, validG], validY))
    predictY = model.predict([testX, testG]).ravel()

    try:
      fpr, tpr, thres = roc_curve(testY, predictY)
      roc_auc = auc(fpr, tpr)
      print(roc_auc)
      auc_list.append(roc_auc)

      if roc_auc > max_auc:
        model.save(model_file_name)
    except:
      continue

  out_data = [RBP, str(np.mean(auc_list)), str(stdev(auc_list))]
  if os.path.exists(auc_out_file):
    with open(auc_out_file, 'a') as result_file:
      writer = csv.writer(result_file)
      writer.writerow(out_data)
    continue

  with open(auc_out_file, 'w') as result_file:
    writer = csv.writer(result_file)
    writer.writerow(out_data)

