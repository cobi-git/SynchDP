import sys
import os
import shutil
import glob
import pickle
import numpy as np
import pandas as pd

# with open('./reference_seq_set.pkl','rb') as f:
#     dat = pickle.load(f)

with open('./query_set.pkl','rb') as f:
    dat = pickle.load(f)


# COV-CCO-035
# COV-CCO-079