#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 17:06:40 2017

@author: qiaochug
"""

import time
import argparse
import numpy as np
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
from sklearn.metrics import precision_score, recall_score, accuracy_score, roc_curve
import pickle
import matplotlib.pyplot as plt

def normalize(X):
    mu = np.zeros(np.size(X,1))
    sigma = np.zeros(np.size(X,1))
    count = 0
    trans = X.T
    for col in trans:
        mu[count] = np.mean(col)
        sigma[count] = np.std(col)
        innercount = 0
        for ele in col:
            trans[count][innercount] = (ele - mu[count])/sigma[count]
            innercount = innercount + 1
        count = count + 1
    
    normalized_X = trans.T
    return mu, sigma, normalized_X

def normalizeUsing(X, mu, sigma):
    trans= X.T
    count = 0
    for col in trans:
        innercount = 0
        for ele in col:
            trans[count][innercount] = (ele - mu[count])/sigma[count]
            innercount = innercount + 1
        count = count + 1
    
    normalized_X = trans.T
    return normalized_X

def trainpredictsave(testdata, loadclf, musigma):

    testdata = pd.read_csv(testdata)
    
    testy = testdata['y'].values
    testX = testdata.drop('y', axis = 1).values
    musigma = pickle.load(open(musigma, "rb"))
    mu = musigma[0]
    sigma = musigma[1]
    #testX = normalizeUsing(testX,mu,sigma)
    
    #pickle.dump([mu,sigma], open(musigma, "wb"))
    clf = pickle.load(open(loadclf, "rb"))
    
    t0 = time.time()
    predicty = clf.predict_proba(testX)
    t1 = time.time()
    
    print(clf.classes_)
    fpr,tpr,thresholds = roc_curve(testy, predicty[:,1])
    plt.plot(fpr, tpr)
    plt.show()
    plt.title("ROC_CURVE RanForest AllData Balanced")
    plt.close()
    print(fpr, tpr, thresholds, '\n')
    
    predicty = clf.predict(testX)
    print('Time: ', (t1-t0), ' sec\n')
    accuracy = accuracy_score(testy, predicty)*100
    print('Accuracy: ', accuracy, '%\n')
    recall = recall_score(testy,predicty)*100
    print('Recall: ', recall, '%\n')
    precision = precision_score(testy, predicty)*100
    print('Precision: ', precision, '%\n')
    
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('testdata', nargs = 1, default = None,
                         help = 'the testing set')
    parser.add_argument('loadclf', nargs = 1, default = None, 
                        help = 'the pickle file that saves the trained parameters')
    parser.add_argument('loadmusigma', nargs = 1, default = None, 
                        help = 'the pickle file that saves the mu and sigma')
    
    args = parser.parse_args()
    
    trainpredictsave(args.testdata[0], args.loadclf[0], args.loadmusigma[0])

if __name__ == '__main__':
    main()