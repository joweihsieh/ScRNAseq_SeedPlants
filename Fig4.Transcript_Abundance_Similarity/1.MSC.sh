################ 
# setting and testing code by toy example
################

#launch Jupyter Notebook
#pip3 install matplotlib
#pip3 install pandas
#pip install tk (for vscode)

################
# Testing - Import packages
################

#import pandas
#import tkinter as tk
#from tkinter import filedialog
#import MSClustering as MSC
#import numpy as np
#import matplotlib.pyplot as plt 
#import time

################
# Generate distance matrix of UMI counts for MSC
################

nohup Rscript DistMatrix.R &

################
# open terminal
################

cd /Users/joweihsieh/Dropbox/YCL/Single_Cell_Cla/other_clustering/MSC_python_Hu
python MSC_demo.py

################
# select the excel table of similarity matrix on GUI
################