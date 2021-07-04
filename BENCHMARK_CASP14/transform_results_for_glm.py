#!/usr/bin/env python3

"""
FILES
Generate 2 files to create GLM figures:
    1 - Output probas + true (1) or false (0): 2 columns: balanced_for_glm.txt + global_for_glm.txt
    2 - Neq for a residue position + true (1) or (false): 2 columns: balanced_neq.txt + global_neq_txt

FIGURES
    1 - glm of PYTHIA models output probas vs. TPR
    2 - glm of Neq vs. TPR
"""


import os
import sys
import numpy as np
from sklearn.metrics import accuracy_score

res_casp_balanced = ["results/balanced/6uf2A/PYTHIA.o3P1-20210630150307/prediction/PB_prediction.csv",
            "results/balanced/6xc0C/PYTHIA.3Os3-20210630170825/prediction/PB_prediction.csv",
            "results/balanced/6y4fA/PYTHIA.wr5s-20210630150605/prediction/PB_prediction.csv",
            "results/balanced/6ya2A/PYTHIA.KQjS-20210630150723/prediction/PB_prediction.csv",
            "results/balanced/6ya2B/PYTHIA.MdY5-20210630150847/prediction/PB_prediction.csv",
            "results/balanced/6ya2C/PYTHIA.dl2Z-20210630151013/prediction/PB_prediction.csv",
            "results/balanced/6zycA/PYTHIA.os28-20210630151137/prediction/PB_prediction.csv",
            "results/balanced/7d2oA/PYTHIA.NoOj-20210630151312/prediction/PB_prediction.csv",
            "results/balanced/7jtlA/PYTHIA.6O2b-20210630151433/prediction/PB_prediction.csv",
            "results/balanced/7m7aA/PYTHIA.Pv0z-20210630213955/prediction/PB_prediction.csv",
            "results/balanced/7k7wA/PYTHIA.uoT9-20210630214250/prediction/PB_prediction.csv",
            "results/balanced/7m7aB/PYTHIA.Zjwp-20210630151741/prediction/PB_prediction.csv",
            "results/balanced/7m7aC/PYTHIA.TghN-20210630152014/prediction/PB_prediction.csv"]
res_casp_global = ["results/global/6uf2A/PYTHIA.VqX2-20210630133531/prediction/PB_prediction.csv",
            "results/global/6xc0C/PYTHIA.TuZ8-20210630133701/prediction/PB_prediction.csv",
            "results/global/6y4fA/PYTHIA.Dlkj-20210630133829/prediction/PB_prediction.csv",
            "results/global/6ya2A/PYTHIA.f75G-20210630133946/prediction/PB_prediction.csv",
            "results/global/6ya2B/PYTHIA.F8fB-20210630134112/prediction/PB_prediction.csv",
            "results/global/6ya2C/PYTHIA.KrAQ-20210630134237/prediction/PB_prediction.csv",
            "results/global/6zycA/PYTHIA.yM0B-20210630134403/prediction/PB_prediction.csv",
            "results/global/7d2oA/PYTHIA.PVyX-20210630134522/prediction/PB_prediction.csv",
            "results/global/7jtlA/PYTHIA.ABvP-20210630134645/prediction/PB_prediction.csv",
            "results/global/7m7aA/PYTHIA.VyQ7-20210630214852/prediction/PB_prediction.csv",
            "results/global/7k7wA/PYTHIA.Q6hg-20210630215149/prediction/PB_prediction.csv",
            "results/global/7m7aB/PYTHIA.LHV0-20210630134954/prediction/PB_prediction.csv",
            "results/global/7m7aC/PYTHIA.gaFD-20210630135229/prediction/PB_prediction.csv"]

pb_dict = {'a':0,'b':1,'c':2,'d':3,'e':4,'f':5,'g':6,'h':7,'i':8,'j':9,'k':10,'l':11,'m':12,'n':13,'o':14,'p':15}
c = 0
a = 0
b = 0
with open("balanced_for_glm.txt", "w") as f, open("balanced_neq.txt", "w") as f2:
    for res in res_casp_balanced:
        target = res.split("/")[2]
        probs = np.loadtxt(res, skiprows=1, usecols = tuple(range(3,19)))
        with open("data/pdb_pbfasta/"+target+".pb_adjust") as f1:
            f1.readline()
            true_pb_seq = np.array(list(f1.readline().strip()))
        ind = np.where(true_pb_seq == "Z")
        probs = np.delete(probs, ind, axis=0)
        true_pb_seq = np.delete(true_pb_seq, ind)
        for i, pb in enumerate(true_pb_seq):
            c += 1
            tmp = np.zeros(16, dtype=int)
            tmp[pb_dict[pb]] = 1
            neq = 0
            for j in range(16):
                if probs[i][j] == 0:
                    neq += 0
                else:
                    neq += probs[i][j] * np.log(probs[i][j])
                f.write(str(probs[i][j])+" "+str(tmp[j])+"\n")
            if np.argmax(probs[i]) == pb_dict[pb]:
                a += 1
                f2.write(f" {np.exp(-neq):.2f} 1\n")
            else:
                b += 1
                f2.write(f" {np.exp(-neq):.2f} 0\n")
print(c, a, b)

with open("global_for_glm.txt", "w") as f, open("global_neq.txt", "w") as f2:
    for res in res_casp_global:
        target = res.split("/")[2]
        probs = np.loadtxt(res, skiprows=1, usecols = tuple(range(3,19)))
        with open("data/pdb_pbfasta/"+target+".pb_adjust") as f1:
            f1.readline()
            true_pb_seq = np.array(list(f1.readline().strip()))
        ind = np.where(true_pb_seq == "Z")
        probs = np.delete(probs, ind, axis=0)
        true_pb_seq = np.delete(true_pb_seq, ind)
        for i, pb in enumerate(true_pb_seq):
            tmp = np.zeros(16, dtype=int)
            tmp[pb_dict[pb]] = 1
            neq = 0
            for j in range(16):
                if probs[i][j] == 0:
                    neq += 0
                else:
                    neq += probs[i][j] * np.log(probs[i][j])
                f.write(str(probs[i][j])+" "+str(tmp[j])+"\n")
            if np.argmax(probs[i]) == pb_dict[pb]:
                f2.write(f" {np.exp(-neq):.2f} 1\n")
            else:
                f2.write(f" {np.exp(-neq):.2f} 0\n")
