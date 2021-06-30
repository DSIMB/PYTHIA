#!/usr/bin/env python3

"""
Copyright (c) 2019 Cretin Gabriel, Gelly Jean-Christophe
"""

import argparse
import datetime
import math
import os
import sys
import numpy as np
import pandas as pd
import tensorflow_addons as tfa
from tensorflow.keras.models import model_from_json
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras.layers import Activation
from tensorflow.keras.optimizers import Adam

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f','--fasta', help="Path to the fasta file containing the input sequence", type=str, required=True)
parser.add_argument('-mjson','--model-json', help="Path to model.json containing the model's architecture", type=str, required=True)
parser.add_argument('-mh5','--model-h5', help="Path to model.h5 containing the model's weights", type=str, required=True)
parser.add_argument('-m', '--merge', help="Path to the merge feature vector file", type=str, required=True)
parser.add_argument('-o', '--output', help="Path to output prediction directory", type=str, required=True)

args = parser.parse_args()

# Set input variables
CLASSES = 16
CLASS_TO_PB = np.array(['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p'])
PATH_MODEL_JSON = args.model_json
PATH_MODEL_H5 = args.model_h5
PATH_FASTA = args.fasta
PATH_MERGE = args.merge
PATH_OUTPUT = args.output

class Mish(Activation):

    def __init__(self, activation, **kwargs):

        super(Mish, self).__init__(activation, **kwargs)

        self.__name__ = 'Mish'

get_custom_objects().update({'Mish': Mish(tfa.activations.mish)})

def lr_schedule(epoch):
    """Learning Rate Schedule

    Learning rate is scheduled to be reduced after 80, 120, 160, 180 epochs.
    Called automatically every epoch as part of callbacks during training.

    # Arguments
        epoch (int): The number of epochs

    # Returns
        lr (float32): learning rate
    """
    if epoch > 65:
        lr = 1e-5
    elif epoch > 60:
        lr = 2.5e-5
    elif epoch > 50:
        lr = 5e-5
    else:
        lr = 5e-4
    print("Learning rate (epoch %d): %f" %(epoch, lr))
    return lr


############################

# Check args

if not os.path.isfile(PATH_FASTA):
    sys.exit(f"{PATH_FASTA} does not exist.\nPlease enter a valid fasta file.")

if not os.path.isfile(PATH_MODEL_JSON):
    sys.exit(f"{PATH_MODEL_JSON} does not exist.\nPlease enter a valid json model file.")

if not os.path.isfile(PATH_MODEL_H5):
    sys.exit(f"{PATH_MODEL_H5} does not exist.\nPlease enter a valid h5 model file.")

if not os.path.isfile(PATH_MERGE):
    sys.exit(f"{PATH_MERGE} does not exist.\nPlease enter a valid input features merge file.")

if not os.path.exists(PATH_OUTPUT):
    os.makedirs(PATH_OUTPUT, exist_ok=True)

### LOAD MODEL

# windows, nb features
ROWS, COLS = 61, 100

# Load JSON arch
with open(PATH_MODEL_JSON, 'r') as json_file:
    loaded_model_json = json_file.read()
model = model_from_json(loaded_model_json)
print("Loaded model architecture from: " + PATH_MODEL_JSON)

# Load weights from H5
model.load_weights(PATH_MODEL_H5)
print("Loaded weights from: " + PATH_MODEL_JSON)

# Compile model (required to make predictions)
model.compile(loss='categorical_crossentropy', optimizer=Adam(lr=lr_schedule(0)), metrics=['accuracy'])

# Load full fasta
input_features = np.loadtxt(PATH_MERGE)
input_features = np.reshape(input_features, (len(input_features), ROWS, COLS))

# Load fasta file
seq = np.loadtxt(PATH_FASTA, dtype="str", skiprows=1)
seq_AA = seq.tolist()

y_pred = model.predict(input_features, verbose=0)
y_pred_classes = np.argmax(y_pred, axis=1)
y_pred_pb = CLASS_TO_PB[y_pred_classes]
prob_y_pred = np.amax(y_pred, axis=1)

# First 2 and last 2 PBs cannot be predicted
# because of definition of Protein Blocks
# So we set 'ZZ...ZZ' extremities
np.put(y_pred_pb, [0, 1, -2, -1], ["Z", "Z", "Z", "Z"])
np.put(prob_y_pred, [0, 1, -2, -1], [0, 0, 0, 0])
y_pred[0, :] = 0
y_pred[1, :] = 0
y_pred[-2, :] = 0
y_pred[-1, :] = 0

# Create a dataframe from zipped list
zippedList =  list(zip(seq_AA, 
                       y_pred_pb, 
                       prob_y_pred,
                       y_pred[:, 0],
                       y_pred[:, 1], 
                       y_pred[:, 2], 
                       y_pred[:, 3], 
                       y_pred[:, 4], 
                       y_pred[:, 5], 
                       y_pred[:, 6], 
                       y_pred[:, 7], 
                       y_pred[:, 8], 
                       y_pred[:, 9], 
                       y_pred[:, 10], 
                       y_pred[:, 11], 
                       y_pred[:, 12], 
                       y_pred[:, 13], 
                       y_pred[:, 14], 
                       y_pred[:, 15]))
df = pd.DataFrame(zippedList, columns = ["Residue", "Predicted_PB", "P_max", "P_a", "P_b", "P_c", "P_d", "P_e", "P_f", "P_g", "P_h", "P_i", "P_j", "P_k", "P_l", "P_m", "P_n", "P_o", "P_p"])
df.to_csv(os.path.join(PATH_OUTPUT, "PB_prediction.csv"), index=False, float_format="%.2f", sep="\t")

# Get seq id from input fasta
with open(PATH_FASTA) as f:
    seq_id = f.readline().strip()
# Concatenate the predicted Protein Blocks
seq = "".join([pb for pb in y_pred_pb])
# Write the PB sequence as fasta output
with open(os.path.join(PATH_OUTPUT, "predicted_PB.fasta"), "w") as f:
    f.write(f"{seq_id}\n{seq}")



















