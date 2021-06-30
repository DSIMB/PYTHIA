#!/usr/bin/env python3

import os
import numpy as np
from sklearn.metrics import accuracy_score

res_casp_balanced = ["results/balanced/6uf2A/PYTHIA.o3P1-20210630150307/prediction/predicted_PB.fasta",
            "results/balanced/6xc0C/PYTHIA.3Os3-20210630170825/prediction/predicted_PB.fasta",
            "results/balanced/6y4fA/PYTHIA.wr5s-20210630150605/prediction/predicted_PB.fasta",
            "results/balanced/6ya2A/PYTHIA.KQjS-20210630150723/prediction/predicted_PB.fasta",
            "results/balanced/6ya2B/PYTHIA.MdY5-20210630150847/prediction/predicted_PB.fasta",
            "results/balanced/6ya2C/PYTHIA.dl2Z-20210630151013/prediction/predicted_PB.fasta",
            "results/balanced/6zycA/PYTHIA.os28-20210630151137/prediction/predicted_PB.fasta",
            "results/balanced/7d2oA/PYTHIA.NoOj-20210630151312/prediction/predicted_PB.fasta",
            "results/balanced/7jtlA/PYTHIA.6O2b-20210630151433/prediction/predicted_PB.fasta",
            "results/balanced/7m7aA/PYTHIA.Pv0z-20210630213955/prediction/predicted_PB.fasta",
            "results/balanced/7k7wA/PYTHIA.uoT9-20210630214250/prediction/predicted_PB.fasta",
            "results/balanced/7m7aB/PYTHIA.Zjwp-20210630151741/prediction/predicted_PB.fasta",
            "results/balanced/7m7aC/PYTHIA.TghN-20210630152014/prediction/predicted_PB.fasta"]
res_casp_global = ["results/global/6uf2A/PYTHIA.VqX2-20210630133531/prediction/predicted_PB.fasta",
            "results/global/6xc0C/PYTHIA.TuZ8-20210630133701/prediction/predicted_PB.fasta",
            "results/global/6y4fA/PYTHIA.Dlkj-20210630133829/prediction/predicted_PB.fasta",
            "results/global/6ya2A/PYTHIA.f75G-20210630133946/prediction/predicted_PB.fasta",
            "results/global/6ya2B/PYTHIA.F8fB-20210630134112/prediction/predicted_PB.fasta",
            "results/global/6ya2C/PYTHIA.KrAQ-20210630134237/prediction/predicted_PB.fasta",
            "results/global/6zycA/PYTHIA.yM0B-20210630134403/prediction/predicted_PB.fasta",
            "results/global/7d2oA/PYTHIA.PVyX-20210630134522/prediction/predicted_PB.fasta",
            "results/global/7jtlA/PYTHIA.ABvP-20210630134645/prediction/predicted_PB.fasta",
            "results/global/7m7aA/PYTHIA.VyQ7-20210630214852/prediction/predicted_PB.fasta",
            "results/global/7k7wA/PYTHIA.Q6hg-20210630215149/prediction/predicted_PB.fasta",
            "results/global/7m7aB/PYTHIA.LHV0-20210630134954/prediction/predicted_PB.fasta",
            "results/global/7m7aC/PYTHIA.gaFD-20210630135229/prediction/predicted_PB.fasta"]


print("Accuracy PYTHIA balanced model:")
for res in res_casp_balanced:
    target = res.split("/")[2]
    with open(res) as f:
        f.readline()
        pred_pb_seq = np.array(list(f.readline().strip()))[2:-2]
    with open("data/pdb_pbfasta/"+target+".pb_adjust") as f:
        f.readline()
        true_pb_seq = np.array(list(f.readline().strip()))[2:-2]
    ind = np.where(true_pb_seq == "Z")
    pred_pb_seq = np.delete(pred_pb_seq, ind)
    true_pb_seq = np.delete(true_pb_seq, ind)
    acc = accuracy_score(true_pb_seq, pred_pb_seq)
    print("{} {:.3f}".format(target, acc))


print("Accuracy PYTHIA global model:")
for res in res_casp_global:
    target = res.split("/")[2]
    with open(res) as f:
        f.readline()
        pred_pb_seq = np.array(list(f.readline().strip()))[2:-2]
    with open("data/pdb_pbfasta/"+target+".pb_adjust") as f:
        f.readline()
        true_pb_seq = np.array(list(f.readline().strip()))[2:-2]
    ind = np.where(true_pb_seq == "Z")
    pred_pb_seq = np.delete(pred_pb_seq, ind)
    true_pb_seq = np.delete(true_pb_seq, ind)
    acc = accuracy_score(true_pb_seq, pred_pb_seq)
    print("{} {:.3f}".format(target, acc))












