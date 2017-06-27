#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 17:06:40 2017

@author: qiaochug
"""

import argparse
import numpy as np
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
from sklearn.metrics import precision_score, recall_score, accuracy_score
import pickle

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

def trainpredictsave(traindata, testdata, clfOut, musigmaOut):
    traindata = pd.read_csv(traindata)
    testdata = pd.read_csv(testdata)
    trainy = traindata['y'].values
    trainX = traindata.drop('y', axis = 1).values
    mu,sigma,trainX = normalize(trainX)
    
    testy = testdata['y'].values
    testX = testdata.drop('y', axis = 1).values
    testX = normalizeUsing(testX,mu,sigma)
    
    clf = RandomForestClassifier(n_estimators=10)
    clf = clf.fit(trainX, trainy)
    
    pickle.dump(clf, open(clfOut,"wb"))
    pickle.dump([mu, sigma], open(musigmaOut,"wb"))
    
    predicty = clf.predict(testX)
    accuracy = accuracy_score(testy, predicty)*100
    print('Accuracy: ', accuracy, '%\n')
    recall = recall_score(testy,predicty)*100
    print('Recall: ', recall, '%\n')
    precision = precision_score(testy, predicty)*100
    print('Precision: ', precision, '%\n')
    
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('traindata', nargs = 1, default = None,
                         help = 'the training set')
    parser.add_argument('testdata', nargs = 1, default = None,
                         help = 'the testing set')
    parser.add_argument('saveClf', nargs = 1, default = None, 
                        help = 'the pickle file that saves the trained parameters')
    parser.add_argument('savemusigma', nargs = 1, default = None, 
                        help = 'the pickle file that saves the trained parameters')
    
    args = parser.parse_args()
    
    trainpredictsave(args.traindata[0], args.testdata[0], args.saveClf[0], args.savemusigma[0])

if __name__ == '__main__':
    main()