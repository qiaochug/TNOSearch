
import time
from LinkerLib import Triplet
from LinkerLib import Detection
from LinkerLib import printPercentage
from LinkerLib import writeTriplets
from LinkerLib import pickleTriplets
from GraphNodes import connectGraph
from GraphNodes import graphLinks
import argparse
import pickle
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
import itertools
   
#removes detections one by one to test if orbits fit
def removeDets(triplets):
    goodTriplets = []
    badTriplets = []
    time0 = time.time()
    counter = 0
    for trip in triplets:
        printPercentage(counter, len(triplets), time.time()-time0)
        elements, errs= trip.calcOrbit()
        if(elements['a'] > 2 and elements['e'] < 1):
            goodTriplets.append(trip)
        elif(len(trip.dets) > 3):
            found = False
            for x in range(len(trip.dets)):
                tripTest = Triplet(trip.dets[:x] + trip.dets[(x+1):])
                elements, errs = tripTest.calcOrbit()
                if(elements['a'] > 2 and elements['e'] < 1):
                    goodTriplets.append(tripTest)
                    found = True
                    continue

            if(not found):
                badTriplets.append(trip)
            
        else:
            badTriplets.append(trip)
        counter+=1
    return goodTriplets, badTriplets            

#given a list of triplets, return a list of triplets where no two detections have the same exposure
def filterSameExp(triplets):
    print('\n')
    result = []
    for trip in triplets:
        dets = trip.dets
        #place each detection in its respective night
        mjds = set([x.mjd for x in dets]) 
        #dictionary where mjd is key, and value is list of detections
        mjdCount = dict([(x,[]) for x in mjds])
        for det in dets:
            mjdCount[det.mjd].append(det)
        #convert dictionary to list
        countList = mjdCount.values()
        #get every combination of those triplets
        result += (list(itertools.product(*countList)))
    return result

def worthadding(newblob, listofblob):
    for oldblob in listofblob:
        inthisold = True
        for nodeb in newblob.dets:
            if nodeb not in oldblob.dets:
                inthisold = False
                continue
        if inthisold:
            return False
    return True
        

def growtrees(graph, d):
    good = []
    for CC in d:
        CCgood = []
        for trip in d[CC]:
            blobsnow = [trip]
            while blobsnow:
                for blob in blobsnow:
                    successors = graph.successors(blob.dets[-1])
                    if not successors:
                        blobsnow.remove(blob)
                        if worthadding(blob, CCgood):
                            CCgood.append(blob)
                        continue
                    for next in successors:
                        trydets = blob.dets[:]
                        trydets.append(next)
                        tryblob = Triplet(trydets)
                        elements, errs = tryblob.calcOrbit()
                        if(elements['a'] > 2 and elements['e'] < 1):
                            if blob in blobsnow:
                                blobsnow.remove(blob)
                            blobsnow.append(tryblob)
                    if worthadding(blob,CCgood):
                        CCgood.append(blob)
                    if blob in blobsnow:
                        blobsnow.remove(blob)
    good = good + CCgood
    return good

#given a graph, finds the connected components and attempts to identify objects from them
def organizetriplets(graph, triplets):
    components = list(nx.connected_component_subgraphs(graph.to_undirected(),copy=False))
    # This dictionary associates triplets with the connect components they belong to
    d = dict([(com,[]) for com in components])
    for trip in triplets:
        for com in components:
            if com.has_node(trip.dets[0]):
                #print(com, ' has ', trip, '\n')
                d[com].append(trip)
                break
    return d        
   
def main():
    args = argparse.ArgumentParser()
    args.add_argument('triplets', nargs=1, help='path to pickle file')
    args = args.parse_args()
    triplets = pickle.load(open(args.triplets[0],'rb'))
    saveName = args.triplets[0].split('+')[-1].split('.')[0]
    # merges all triplets onto a graph
    G, dict = connectGraph(triplets)
    
    # know the tripelts inside each connected components
    d  = organizetriplets(G, triplets)

    good = growtrees(G,d)

    #saves every merged triplets 
    #saveGoodName = 'goodTriplets+' + saveName
    #saveBadName = 'badTriplets+' + saveName 
    #writeTriplets(good, saveGoodName + '.txt')
    #writeTriplets(bad, saveBadName + '.txt')
    #pickleTriplets(good, saveGoodName + '.pickle')
    #pickleTriplets(bad, saveBadName + '.pickle')
    
    #graphs original merged graph
    #graphLinks(G, dict, saveName)
    #G, dict= connectGraph(good)
    #graphLinks(G, dict, saveName)
    #, [3.4,5.2,-1,0])
    '''
    merged = findObjects(G)
    
    savename = 'mergedTrips+' + saveName + '.txt'
    writeTriplets(merged, savename)
    
    print('writing pickle file')
    savename = 'mergedTrips+' + saveName + '.pickle'
    with open(savename, 'wb') as f:
        pickle.dump(merged,f)
    '''
if __name__ == '__main__':
    main()
