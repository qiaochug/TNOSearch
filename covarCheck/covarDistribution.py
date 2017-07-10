import sys
sys.path.insert(0, './linkerCode/')

import numpy as np
import argparse
import pickle
from LinkerLib import Triplet
import matplotlib.pyplot as plt

#checks the r coefficients of a 6x6 matrix
def checkCovar(covar):
    rcoefs = np.zeros((6,6))
    for x in range(6):
        for y in range(x+1):
            diag1 = covar[x][x]
            diag2 = covar[y][y]
            offdiag = covar[x][y]
            rcoefs[x][y] = round(((offdiag)/(np.sqrt(abs(diag1*diag2))))**2, 8)
    np.set_printoptions(precision=3)
    print(rcoefs)
    return rcoefs    

def main():
    args = argparse.ArgumentParser()
    args.add_argument('tripletFile', nargs=1, help='path to pickle file of triplets')
    args = args.parse_args()
    triplets = list(pickle.load(open(args.tripletFile[0], 'rb')))
    savename = args.tripletFile[0].split('+')[-1].split('.')[0]
    #the x_axis stores the semi-major axis
    x_axis = np.zeros(len(triplets))
    #the y_axes are the Rs
    #there is an array of arrays
    #[y10,y20,y21,y30,y31,y32,y40,y41,y42,y43,y50,y51,y52,y53,y54]
    y_axes = np.zeros((15, len(triplets)))
    
    for y in range(len(triplets)):
        #fill the x_axis with the semi-major axis
        triplets[y].calcOrbit()
        x_axis[y] = triplets[y].elements['a']
        
        covar = triplets[y].getCovar()
        print(triplets[y].toStr())
        rceofs = checkCovar(covar)
        #make list of indexes of covariance values of interest
        ind1 = np.array([1,2,2,3,3,3,4,4,4,4,5,5,5,5,5])
        ind2 = np.array([0,0,1,0,1,2,0,1,2,3,0,1,2,3,4])
        #fill the y_axes with covariance values, from the 0 y_axis to the 14
        for i in np.arange(0,15):
            y_axes[i][y] = rceofs[ind1[i]][ind2[i]]

    #setup for scatter plot
    colors = (0,0,0)
    area = np.pi*2
    for i in np.arange(0,15):
        plt.scatter(x_axis, y_axes[i],s=area,c=colors,alpha=0.5)
        plt.title(savename + '_R' + str(ind1[i]) + str(ind2[i]))
        plt.xlabel('semi-major axis')
        plt.ylabel('r')
        plt.savefig(savename + '_R' + str(ind1[i]) + str(ind2[i]) + '.png')
            

if __name__=='__main__':
    main()
