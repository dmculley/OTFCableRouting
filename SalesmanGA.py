# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 11:34:03 2013

@author: dmc13
"""

# Fixed Start Open Multiple Traveling Salesmen Problem (M-TSP) Genetic Algorithm (GA)

import numpy as np
import math


#### RAW INPUTS ####

turbine_locations = [[1,2],[1,3],[1,3],[1,4],[2,2],[2,3],[2,3],[2,4],[3,2],[3,3],[3,3],[3,4],[4,2],[4,3],[4,3],[4,4]]
substation_location = [[0,9]]

Capacity = 6

PopSize = 80

NumIter = 5000

ShowProg = False


#### PROCESS RAW INPUTS ####

vertices = substation_location + turbine_locations

NRoutes = int(math.ceil(float(len(vertices)) / Capacity))


# Check this - it needs to be as flexible as possible while forcing the minimum number or routes.
MinRoute = len(vertices) / NRoutes

def construct_cost_matrix(vertices):
    '''Constructs a matrix of costs for every potential edge - connections of vertices to themselves are set to infinity'''
    distances = []
    grouped_distances = []
    for i in range(len(vertices)):
        for j in range(len(vertices)):
            if i == j:
                dist = np.inf
            else:
                dist = (np.sqrt((vertices[i][0] - vertices[j][0])**2 + (vertices[i][1] - vertices[j][1])**2))
            distances.append(dist)
    for i in range(0, len(distances), len(vertices)):
        grouped_distances.append(tuple(distances[i:i+len(vertices)]))
    C = np.array(grouped_distances)
    return C
    
DMat = construct_cost_matrix(vertices)


#### Initialisations ####

# Initialise Route Break Selection
NBreaks = NRoutes - 1
dof = len(turbine_locations) - MinRoute * NRoutes
addto = np.ones((1, dof+1))
for k in range(2, NBreaks):
    addto = np.cumsum(addto)
cumProb = np.cumsum(addto)/sum(addto)


def rand_breaks():
    if MinRoute == 1:      # No constraints on breaks
        TempBreaks = np.random.permutation(range(len(turbine_locations)-1))
        Breaks = TempBreaks[0:NBreaks].sort
    else:                  # Force breaks to be at least the minimum tour length
        NAdjust = find(rand < cumprob, 1) - 1
        Spaces = ceil(NBreaks * np.random.rand(1,NAdjust))
        Adjust = np.zeros((1, NBreaks))
        for j in range(NBreaks):
#            Adjust[j] = sum(Spaces == j)
        Breaks = MinTour * range(NBreaks) + np.cumsum(Adjust)


# Initialise the Populations

PopRoute = np.zeros((PopSize, len(turbine_locations)))
PopBreak = np.zeros((PopSize, NBreaks))
PopRoute[0] = range(1, len(turbine_locations)+1)
PopBreak[0] = rand_breaks()
for k in range(1, PopSize):
    PopRoute[k] = np.random.rand(len(turbine_locations)) + 1
    PopBreak[k] = rand_breaks()


#### Run the Genetic Algorithm ####

GlobalMin = np.inf

TotalDist = np.zeros((1, PopSize))

DistHistory = np.zeros((1, NumIter))

TempPopRoute = np.zeros((8, len(turbine_locations)))

TempPopBreaks = np.zeros((8, NBreaks))

NewPopRoute = np.zeros((PopSize, len(turbine_locations)))

NewPopBreaks = np.zeros((PopSize, NBreaks))


def rand_breaks():
    if MinRoute == 1:      # No constraints on breaks
        TempBreaks = np.random.permutation(range(len(turbine_locations)-1))
        Breaks = TempBreaks[0:NBreaks].sort
    else:                  # Force breaks to be at least the minimum tour length
        NAdjust = find(rand < cumprob, 1) - 1
        Spaces = ceil(NBreaks * np.random.rand(1,NAdjust))
        Adjust = np.zeros((1, NBreaks))
        for j in range(NBreaks):
#            Adjust[j] = sum(Spaces == j)
        Breaks = MinTour * range(NBreaks) + np.cumsum(Adjust)




if ShowProg:
    print 'Enter script displaying progress of the algorithm'
    
for iteration in range(NumIter):
    for p in range(PopSize):
        d = 0
        pRoute = PopRoute[p]
        pBreak = PopBreak[p]
#        rng = [[1 pBreak + 1];[PBreak n]]
        for r in range(NRoutes):
            d = DMat[0][pRoute[rng[r][1]]]
            for k in range(rng[r][1], rng[r][2]-1):
                d = d + DMat[pRoute[k]][pRoute[k+1]]
        TotalDist[p] = d
        
    ## Find the Best Route in the Population ##

    MinDist = min(TotalDist)
    MinDistIndex = TotalDist.index(MinDist)
    DistHistory[iteration] = MinDist
    if MinDist < GlobalMin:
        GlobalMin = MinDist
        OptRoute = PopRoute[MinDistIndex]
        OptBreak = PopBreak[MinDistIndex]
#        rng = [[1 OptBreak + 1]; [OptBreak n]]
        if ShowProg:
            print 'Enter script displaying progress of the algorithm'
        
    ## Genetic Algorithm Operators ##
    
    RandomOrder = np.random.permutation(range(PopSize))
    for p in range(8, PopSize, 8):
        rtes = PopRoute[RandomOrder[p-7:p]]
        brks = PopBreak[RandomOrder[p-7:p]]
        dists = TotalDist[RandomOrder[p-7:p]]
        Ignore = min(dists)
        IgnoreIndex = dists.index(dists)
        BestOf8Route = rtes[IgnoreIndex]
        BestOf8Break = brks[IgnoreIndex]
        RouteInsertionPoints = (math.ceil(len(turbine_locations) * np.random.rand(1,2))).sort
        I = RouteInsertionPoints[0]
        J = RouteInsertionPoints[1]
        for k in range(8):
            TempPopRoute[k] = BestOf8Route
            TempPopBreaks[k] = BestOf8Break


            #### GENERATE NEW SOLUTIONS ??? ####
            
        NewPopRoute[p-7][p] = TempPopRoute
        NewPopBreak[p-7][p] = TempPopBreak
        
    PopRoute = NewPopRoute
    PopBreak = NewPopBreak

if ShowResult:
    print 'Enter script to display result'
    
    

        














