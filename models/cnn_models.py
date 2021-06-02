import keras
from keras.models import Sequential, Model
from keras.layers import Dense, Flatten, MaxPooling1D, Conv1D, Dropout, Input, Concatenate
from keras.regularizers import l2
from keras import backend as tf

# base model without gene expression or annotation input
def get_base_model(input_shape, out_dim):
  model_in = Input(shape=input_shape)

  conv1 = Conv1D(100, kernel_size=8, activation='relu', input_shape=input_shape, kernel_regularizer=l2(.01))(model_in)
  pooling = MaxPooling1D(pool_size=4, strides=None, padding='valid')(conv1)
  drop1 = Dropout(.2)(pooling)

  conv2 = Conv1D(100, kernel_size=8, activation='relu', kernel_regularizer=l2(.01))(drop1)
  pooling2 = MaxPooling1D(pool_size=4, strides=None, padding='valid')(conv2)
  drop2 = Dropout(.2)(pooling2)

  flat = Flatten()(drop2)

  dense1 = Dense(600, activation='relu')(flat)
  dense2 = Dense(300, activation='relu')(dense1)
  dense3 = Dense(64, activation='relu')(dense2)

  model_out = Dense(out_dim, activation="sigmoid")(dense3)
  model = Model(inputs=model_in, outputs=model_out)

  return model

# model with gene expression input
def get_gene_model(input_shape, out_dim):
  model_in = Input(shape=input_shape)
  rank_in = Input(shape=(1,))

  conv1 = Conv1D(100, kernel_size=8, activation='relu', input_shape=input_shape, kernel_regularizer=l2(.01))(model_in)
  pooling = MaxPooling1D(pool_size=4, strides=None, padding='valid')(conv1)
  drop1 = Dropout(.2)(pooling)

  conv2 = Conv1D(100, kernel_size=8, activation='relu', kernel_regularizer=l2(.01))(drop1)
  pooling2 = MaxPooling1D(pool_size=4, strides=None, padding='valid')(conv2)
  drop2 = Dropout(.2)(pooling2)

  flat = Flatten()(drop2)

  dense1 = Dense(600, activation='relu')(flat)
  dense2 = Dense(300, activation='relu')(dense1)
  dense3 = Dense(64, activation='relu')(dense2)

  dense_rank = Dense(1,)(rank_in)

  mult_layer_gene = Concatenate()([dense3, dense_rank])
  model_gene_out = Dense(out_dim, activation="sigmoid")(mult_layer_gene)
  model_gene = Model(inputs=[model_in, rank_in], outputs=model_gene_out)
  return model_gene


# model with annotation data
def get_annot_model(input_shape, out_dim):
  model_in = Input(shape=input_shape)
  annot_in = Input(shape=(6,))

  conv1 = Conv1D(100, kernel_size=8, activation='relu', input_shape=input_shape, kernel_regularizer=l2(.01))(model_in)
  pooling = MaxPooling1D(pool_size=4, strides=None, padding='valid')(conv1)
  drop1 = Dropout(.2)(pooling)

  conv2 = Conv1D(100, kernel_size=8, activation='relu', kernel_regularizer=l2(.01))(drop1)
  pooling2 = MaxPooling1D(pool_size=4, strides=None, padding='valid')(conv2)
  drop2 = Dropout(.2)(pooling2)

  flat = Flatten()(drop2)

  dense1 = Dense(600, activation='relu')(flat)
  dense2 = Dense(300, activation='relu')(dense1)
  dense3 = Dense(64, activation='relu')(dense2)

  dense_annot = Dense(6,)(annot_in)

  mult_layer_annot = Concatenate()([dense3, dense_annot])
  model_annot_out = Dense(out_dim, activation="sigmoid")(mult_layer_annot)
  model_annot = Model(inputs=[model_in, annot_in], outputs=model_annot_out)

  return model_annot


# model with gene expression and annotation data
def get_gene_annot_model(input_shape, out_dim):
  model_in = Input(shape=input_shape)
  rank_in = Input(shape=(1,))
  annot_in = Input(shape=(6,))

  conv1 = Conv1D(100, kernel_size=8, activation='relu', input_shape=input_shape, kernel_regularizer=l2(.01))(model_in)
  pooling = MaxPooling1D(pool_size=4, strides=None, padding='valid')(conv1)
  drop1 = Dropout(.2)(pooling)

  conv2 = Conv1D(100, kernel_size=8, activation='relu', kernel_regularizer=l2(.01))(drop1)
  pooling2 = MaxPooling1D(pool_size=4, strides=None, padding='valid')(conv2)
  drop2 = Dropout(.2)(pooling2)

  flat = Flatten()(drop2)

  dense1 = Dense(600, activation='relu')(flat)
  dense2 = Dense(300, activation='relu')(dense1)
  dense3 = Dense(64, activation='relu')(dense2)

  dense_rank = Dense(1,)(rank_in)
  dense_annot = Dense(6,)(annot_in)

  mult_layer_gene_annot = Concatenate()([dense3, dense_rank, dense_annot])
  model_gene_annot_out = Dense(out_dim, activation="sigmoid")(mult_layer_gene_annot)
  model_gene_annot = Model(inputs=[model_in, rank_in, annot_in], outputs=model_gene_annot_out)

  return model_gene_annot

#model with dilation step in second convolution layer
def get_dilated_model(input_shape, out_dim):
  model_in = Input(shape=input_shape)
  rank_in = Input(shape=(1,))

  conv1 = Conv1D(100, kernel_size=8, activation='relu', input_shape=input_shape, kernel_regularizer=l2(.01))(model_in)
  pooling = MaxPooling1D(pool_size=4, strides=None, padding='valid')(conv1)
  drop1 = Dropout(.2)(pooling)

  conv2 = Conv1D(100, kernel_size=8, activation='relu', kernel_regularizer=l2(.01), dilation_rate=2)(drop1)
  pooling2 = MaxPooling1D(pool_size=4, strides=None, padding='valid')(conv2)
  drop2 = Dropout(.2)(pooling2)

  flat = Flatten()(drop2)

  dense1 = Dense(600, activation='relu')(flat)
  dense2 = Dense(300, activation='relu')(dense1)
  dense3 = Dense(64, activation='relu')(dense2)

  dense_rank = Dense(1,)(rank_in)

  mult_layer_gene_annot = Concatenate()([dense3, dense_rank])
  model_gene_annot_out = Dense(out_dim, activation="sigmoid")(mult_layer_gene_annot)
  model_gene_annot = Model(inputs=[model_in, rank_in], outputs=model_gene_annot_out)

  return model_gene_annot