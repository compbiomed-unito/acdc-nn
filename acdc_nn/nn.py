from silence_tensorflow import silence_tensorflow
silence_tensorflow()

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.layers import Lambda, Dense, Input, Dropout, Activation, Flatten, BatchNormalization, Dropout
from tensorflow.keras import backend
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dot, Conv1D, concatenate, GlobalAveragePooling1D
from tensorflow.keras import backend as K


# Build the network with global average.
def ACDC(num_H, d, num_3d):
    Coding = 20 # size of the coding, 20 residues, 
    #
    num_Env3D = num_3d   # number of 3D neighbours
    num_Env3D_i = num_3d   # number of 3D neighbours
    num_Env1D = 5   # number of sequeunce neighbours
    num_Env1D_i = 5 # number of sequence neighbours
    num_IMut = Coding
    num_IMut_i = Coding
    #
    # Create an input for the 3D neighbours
    input_env3D = Input(shape=(num_Env3D,Coding), name="input_env3D") # direct 
    input_env3D_i = Input(shape=(num_Env3D_i,Coding), name="input_env_i3D") # inverse
    #
    conv_3D = Conv1D(filters=Coding, kernel_size=1, strides=1, use_bias=False, name="net_envConv3D")
    #
    net_env3D = conv_3D(input_env3D) # direct 
    net_env3D_i = conv_3D(input_env3D_i) # inverse
    # create the input for the mutation for 3D
    # direct
    inputMut = Input(shape=(num_IMut,), name="inputMut")
    pool = GlobalAveragePooling1D()(net_env3D)###
    dot_prod = Dot(axes=(1,1))([pool,inputMut]) 
    # inverse
    inputMut_i = Input(shape=(num_IMut_i,), name="inputMut_i")
    pool_i = GlobalAveragePooling1D()(net_env3D_i)###
    dot_prod_i = Dot(axes=(1,1))([pool_i,inputMut_i])
    # Create an input for the Sequence neighbours
    input_env1D = Input(shape=(num_Env1D,Coding), name="input_env1D") # direct
    input_env1D_i = Input(shape=(num_Env1D_i,Coding), name="input_env_i1D") # inverse
    #
    conv_1D = Conv1D(filters=Coding, kernel_size=1, strides=1, use_bias=False, name="net_envConv1D")
    #
    net_env1D = conv_1D(input_env1D) # direct
    net_env1D_i = conv_1D(input_env1D_i) # inverse
    #
    # create the input for the mutation for 1D
    dot_prod1D = Dot(axes=(2,1))([net_env1D,inputMut])
    dot_prod1D_i = Dot(axes=(2,1))([net_env1D_i,inputMut_i])
    #
    # flat dense layer after the convolutional
    net = Sequential(name="TopLayer")
    #num_H=[32,16] # hidden nodes for each layer starting from the bottom
    i = 1
    for hidden in num_H:
        net.add(Dense(hidden, activation="relu", name="Dense_"+str(i)))
        net.add(Dropout(d)) #0.2
        i+=1
    net.add(Dense(1, activation='linear', name="Dense_"+str(i)))
    # Concatenate 3D and 1D inputs
    # direct
    conc_dir = concatenate([dot_prod,dot_prod1D, inputMut]) ###
    net_dir = net(conc_dir)
    # inverse
    conc_inv = concatenate([dot_prod_i,dot_prod1D_i, inputMut_i])###
    net_inv = net(conc_inv)
    # difference and average layers
    out_d = Lambda(diff, output_shape=diff_shape, name="out_diff")([net_dir, net_inv])
    out_av = Lambda(avout, output_shape=avout_shape, name="out_av")([net_dir, net_inv])
    # final training model
    net_dir_inv = Model([input_env3D, input_env1D, inputMut, input_env3D_i, input_env1D_i, inputMut_i],[out_d, out_av])
    # partial model, to use only one flow either direct or inverse
    net_DDG = Model([input_env3D, input_env1D, inputMut],[net_dir])
    net_dir_inv.compile(optimizer='adam',
              loss=['logcosh', abs_loss],
              loss_weights=[1., 1.] )
    return net_dir_inv, net_DDG



def diff(vects):
    x, y = vects
    return (x-y)/2.0


def diff_shape(shapes):
    shape1, shape2 = shapes
    return (shape1[0], 1)

def avout(vects):
    x, y = vects
    return (x+y)/2.0


def avout_shape(shapes):
    shape1, shape2 = shapes
    return (shape1[0], 1)
    

def abs_loss(y_true, y_pred):
    return K.mean(K.abs(y_pred), axis=-1)



def mkInp(X, col_3d):
    Xm = X[:, 0:20]
    X1D = X[:, 20:120].reshape(X.shape[0], 100//20, 20)
    X3D = X[:, 120:].reshape(X.shape[0], col_3d//20, 20)
    return Xm, X1D, X3D    
