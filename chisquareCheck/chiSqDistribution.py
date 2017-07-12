import sys
sys.path.insert(0, './linkerCode/')

import numpy as np
import argparse
import pickle
from LinkerLib import Triplet
import matplotlib.pyplot as plt   

def main():
    args = argparse.ArgumentParser()
    args.add_argument('goodtripletFile',  nargs=1, help='path to pickle file of triplets')
    args.add_argument('missingtripletFile',  nargs=1, help='path to pickle file of triplets')
    args.add_argument('badtripletFile', nargs=1, help='path to pickle file of triplets')
    args = args.parse_args()
    goodtriplets = list(pickle.load(open(args.goodtripletFile[0], 'rb')))
    missingtriplets = list(pickle.load(open(args.missingtripletFile[0], 'rb')))
    badtriplets = list(pickle.load(open(args.badtripletFile[0], 'rb')))
    
    savename = args.goodtripletFile[0].split('+')[-1].split('.')[0]
    #the x_axis stores the semi-major axis
    x_axisg = np.zeros(len(goodtriplets))
    x_axism = np.zeros(len(missingtriplets))
    x_axisb = np.zeros(len(badtriplets))
    #the y_axis is the chi-sqaures
    y_axisg = np.zeros(len(goodtriplets))
    y_axism = np.zeros(len(missingtriplets))
    y_axisb = np.zeros(len(badtriplets))

    for y in range(len(goodtriplets)):
        #fill the x_axis with the semi-major axis
        goodtriplets[y].calcOrbit()
        x_axisg[y] = goodtriplets[y].elements['a']
        #fill the y_axis with the chi sqaures
        y_axisg[y] = goodtriplets[y].getChiSq()
    print("good", len(goodtriplets))

    for y in range(len(missingtriplets)):
        missingtriplets[y].calcOrbit()
        x_axism[y] = missingtriplets[y].elements['a']
        y_axism[y] = missingtriplets[y].getChiSq()
    print("missing", len(missingtriplets))

    for y in range(len(badtriplets)):
        badtriplets[y].calcOrbit()
        x_axisb[y] = badtriplets[y].elements['a']
        y_axisb[y] = badtriplets[y].getChiSq()
    print("bad", len(badtriplets))

    #setup for scatter plot
    colors = ["black", "blue", "red"]
    x_axes = [x_axisg, x_axism, x_axisb]
    y_axes = [y_axisg, y_axism, y_axisb]
    for i in np.arange(0,3):
        y_axes[i] = np.log10(y_axes[i])
    labels = ["good", "missing", "bad"]
    area = np.pi*2
    for i in np.arange(0,3):
        plt.scatter(x_axes[i], y_axes[i],s=area,c=colors[i],alpha=0.5,label=labels[i])
    plt.title('ChiSq_' + savename)
    plt.xlabel('semi-major axis')
    plt.ylabel('ChiSq')
    plt.savefig('ChiSq_' + savename + 'log'+'.png')
    plt.close()

if __name__=='__main__':
    main()
