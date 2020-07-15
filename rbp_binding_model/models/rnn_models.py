import keras
from keras.models import Sequential, Model
from keras.layers import Dense, Flatten, MaxPooling1D, Conv1D, Dropout, LSTM, Input, Concatenate, Layer
from keras.regularizers import l2
from keras import backend as tf

# CNN-LSTM combined model
def get_cnn_rnn_model(input_shape, out_dim):
  model_in = Input(shape=input_shape)
  rank_in = Input(shape=(1,))
  annot_in = Input(shape=(6,))

  conv = Conv1D(100, kernel_size=8, activation='relu', input_shape=input_shape, 
                    kernel_regularizer=l2(.01))(model_in)
  pooling = MaxPooling1D(pool_size=4, strides=None, padding='valid')(conv)
  drop1 = Dropout(.2)(pooling)

  lstm = LSTM(50, return_sequences=True, activation='relu',
                    kernel_regularizer=l2(.01))(drop1)
  pooling2 = MaxPooling1D(pool_size=4, strides=None, padding='valid')(lstm)
  drop2 = Dropout(.2)(pooling2)

  flat = Flatten()(drop2)

  dense1 = Dense(300, activation='relu')(flat)
  dense2 = Dense(64, activation='relu')(dense1)

  dense_rank = Dense(1,)(rank_in)
  dense_annot = Dense(6,)(annot_in)

  mult_layer_gene_annot = Concatenate()([dense2, dense_rank, dense_annot])
  model_gene_annot_out = Dense(out_dim, activation="sigmoid")(mult_layer_gene_annot)
  model_gene_annot = Model(inputs=[model_in, rank_in, annot_in], outputs=model_gene_annot_out)

  return model_gene_annot


# Simple RNN model
def get_rnn_model(input_shape, out_dim):
  model_in = Input(shape=input_shape)
  rank_in = Input(shape=(1,))
  annot_in = Input(shape=(6,))

  lstm = LSTM(50, return_sequences=True, activation='relu', input_shape=input_shape, kernel_regularizer=l2(.01))(model_in)
  drop2 = Dropout(.5)(lstm)

  flat = Flatten()(drop2)

  dense1 = Dense(300, activation='relu')(flat)
  dense2 = Dense(64, activation='relu')(dense1)

  dense_rank = Dense(1,)(rank_in)
  dense_annot = Dense(6,)(annot_in)

  mult_layer_gene_annot = Concatenate()([dense2, dense_rank, dense_annot])
  model_gene_annot_out = Dense(out_dim, activation="sigmoid")(mult_layer_gene_annot)
  model_gene_annot = Model(inputs=[model_in, rank_in, annot_in], outputs=model_gene_annot_out)

  return model_gene_annot


class attention(Layer):
    def __init__(self,**kwargs):
        super(attention,self).__init__(**kwargs)

    def build(self,input_shape):
        self.W=self.add_weight(name="att_weight",shape=(input_shape[-1],1),initializer="normal")
        self.b=self.add_weight(name="att_bias",shape=(input_shape[1],1),initializer="zeros")        
        super(attention, self).build(input_shape)

    def call(self,x):
        et=tf.squeeze(tf.tanh(tf.dot(x,self.W)+self.b),axis=-1)
        at=tf.softmax(et)
        at=tf.expand_dims(at,axis=-1)
        output=x*at
        return tf.sum(output,axis=1)

    def compute_output_shape(self,input_shape):
        return (input_shape[0],input_shape[-1])

    def get_config(self):
        return super(attention,self).get_config()


def get_rnn_attention_model(input_shape, out_dim):
  model_in = Input(shape=input_shape)
  rank_in = Input(shape=(1,))
  annot_in = Input(shape=(6,))

  conv = Conv1D(100, kernel_size=8, activation='relu', input_shape=input_shape, kernel_regularizer=l2(.01))(model_in)
  pooling = MaxPooling1D(pool_size=4, strides=None, padding='valid')(conv)
  drop1 = Dropout(.2)(pooling)

  att_in = LSTM(50, return_sequences=True, activation='relu', kernel_regularizer=l2(.01))(drop1)
  att_out = attention()(att_in)
  drop2 = Dropout(.2)(att_out)

  dense1 = Dense(300, activation='relu')(drop2)
  dense2 = Dense(64, activation='relu')(dense1)

  dense_rank = Dense(1,)(rank_in)
  dense_annot = Dense(6,)(annot_in)

  mult_layer_gene_annot = Concatenate()([dense2, dense_rank, dense_annot])
  model_gene_annot_out = Dense(out_dim, activation="sigmoid")(mult_layer_gene_annot)
  model_gene_annot = Model(inputs=[model_in, rank_in, annot_in], outputs=model_gene_annot_out)

  return model_gene_annot
