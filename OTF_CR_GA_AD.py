# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 10:11:39 2013

@author: dmc13
"""

#import scipy as sp
from matplotlib import pyplot as plt
import math
import copy
import numpy as np
from ad import adnumber
import ad
from ad.admath import *


## 0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0 USEFUL FUNCTIONS 0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0 ##

def construct_cost_matrix(vertices):
    '''Constructs a matrix of costs for every potential edge - connections of vertices to themselves are set to infinity'''
    distances = []
    grouped_distances = []
    for i in range(len(vertices)):
        for j in range(len(vertices)):
            if i == j:
                dist = 99
            else:
                dist = (sqrt((vertices[i][0] - vertices[j][0])**2 + (vertices[i][1] - vertices[j][1])**2))
            distances.append(dist)
    for i in range(0, len(distances), len(vertices)):
        grouped_distances.append(tuple(distances[i:i+len(vertices)]))
    C = np.array(grouped_distances)
    return C
    
    
    
def rand_breaks(n, MinRoute, Nroutes, NBreaks):
    '''Produces a list of random, but valid, route-breaks'''
    RB = [np.random.random_integers(MinRoute, Capacity)]
    for i in range(1, NBreaks):
        RB.append(np.random.random_integers(MinRoute, Capacity) + RB[i-1])
    if RB[-1] < (n - Capacity):
        short = (n - Capacity) - RB[-1]
        add_each = int(np.ceil(0.5 + short / len(RB)))
        for i in range(len(RB)):
            RB[i] = RB[i] + add_each * (i+1)
    RB = np.array((RB))
    return RB
    
    
    
def produce_plot(R, V, dist):
    '''Display with matplotlib'''
    plt.title('Total distance='+str(dist))
    plt.plot(V[:,0],V[:,1],'o')
    for i in range(len(R)):
        plt.plot(R[i][:,0], R[i][:,1], '-')
    for i in range(len(V)):
        plt.text(V[i][0], V[i][1], '%s' % (str(i)))
    plt.axis('equal')
    plt.show()
    
    
    
def differentiate(turbine_locations, C, Route, Break, n):
    '''Differentiate the length of the routing w.r.t. the position of the turbines, produces a n x 2 array, ((dC/dx1, dC/dy1), ...)'''
    print 'Determining dC/dX: Performing automatic differentiation...'
    rting = routing(Route, Break)
    TotalDist = routing_distance(rting, C)
    dC_dX = np.zeros((n,2)) 
    for i in range(n):
        dC_dX[i][0] = TotalDist.d(turbine_locations[i][0])
        dC_dX[i][1] = TotalDist.d(turbine_locations[i][1])
    return dC_dX
    
    
    
def routing(Route, Break):
    '''Combine route and break into an array of the routes described as vertices in the order in which they are toured'''
    Route = Route.tolist()
    Break = Break.tolist()    
    rting = [[0] + Route[0:Break[0]]]
    if len(Break) > 1:
        for f in range(1, len(Break)):
            rting.append([0] + Route[Break[f-1]:Break[f]])
    rting.append([0] + Route[Break[-1]:])
    return rting
    
    
    
def routing_coordinates(vertices, rting):
    '''Convert a routing expressed in indexes to a routing expressed in coordinates'''
#    rting = routing(Route, Break)    
    for i in range(len(rting)):
        for j in range(len(rting[i])):
            rting[i][j] = vertices[rting[i][j]]
        rting[i] = np.array(rting[i])
    return rting 

def routing_distance(rting, C):
    '''Return the geometric length of the routing'''
    d = 0        
    for r in range(len(rting)):
        for v in range(len(rting[r]) - 1):
            d += C[rting[r][v]][rting[r][v+1]]
    return d



def convert_to_adnumber(coordinate_list):
    '''Convert the location vectors from floats into adnumbers to enable differentiation'''
    adnumber_coordinate_list = []
    for i in range(len(coordinate_list)):
        adnumber_coordinate_list.append([adnumber(coordinate_list[i][0]), adnumber(coordinate_list[i][1])])
    coordinate_list = adnumber_coordinate_list
    return coordinate_list



def GA(turbine_locations, substation_location, Capacity):
    
    ## ALGORITHMIC PARAMETERS    
    
    PopSize = 16
    NumIter = 12
    ConvergenceDefinition = 2 # number of iterations after which if there is no change, the model is determined to have coverged
    ConvergenceCounter = 0
    Converged = False
    ShowProg = False
    ShowResult = True
    
    ## PROCESS INPUTS
    
    vertices = substation_location + turbine_locations
    NRoutes = int(math.ceil(float(len(vertices)) / Capacity))
    MinRoute = len(vertices) / NRoutes
    n = len(turbine_locations)
    C = construct_cost_matrix(vertices)
    NBreaks = NRoutes - 1
    
    print 'Initialising Genetic Algorithm...'    
    
    ## INITIALISE POPULATION
    
    PopRoute = np.zeros((PopSize, n), dtype = int)
    PopBreaks = np.zeros((PopSize, NBreaks), dtype = int)
    
    ## GENERATE STARTING POPULATION    
    
    PopRoute[0] = range(1, n+1)
    PopBreaks[0] = rand_breaks(n, MinRoute, NRoutes, NBreaks)
    for i in range(1, PopSize):
        PopRoute[i] = np.random.permutation(n) + 1
        PopBreaks[i] = rand_breaks(n, MinRoute, NRoutes, NBreaks)
    
    ## INITIALISE GENETIC ALGORITHM
    
    GlobalMin = np.inf
    TotalDist = np.zeros((1, PopSize), dtype = ad.ADV)
    DistHistory = np.zeros((1, NumIter), dtype = ad.ADV)
    TempPopRoute = np.zeros((8, len(turbine_locations)), dtype = int)
    TempPopBreaks = np.zeros((8, NBreaks), dtype = int)
    NewPopRoute = np.zeros((PopSize, len(turbine_locations)), dtype = int)
    NewPopBreaks = np.zeros((PopSize, NBreaks), dtype = int)
    
    ## RUN THE GENETIC ALGORITHM    
    
    print 'Running Genetic Algorithm...'     

    for i in range(NumIter):
        while not Converged:
            for p in range(PopSize):
                d = 0
                pRoute = PopRoute[p]
                pBreak = PopBreaks[p]
                rting = routing(pRoute, pBreak)
                d = routing_distance(rting, C)
                TotalDist[0][p] = d
            MDidx = np.argmin(TotalDist)
            MinDist = TotalDist[0][MDidx]
            DistHistory[0][i] = MinDist
            if MinDist < GlobalMin: 
                GlobalMin = copy.deepcopy(MinDist)
                OptRoute = copy.deepcopy(PopRoute[MDidx])
                OptBreak = copy.deepcopy(PopBreaks[MDidx])
                if ShowProg:
                    print 'Current Routing Length', GlobalMin
                
                ConvergenceCounter = 0
            else:
                ConvergenceCounter += 1    
            RandomOrder = np.random.permutation(PopSize)
            for p in range(8, PopSize+1, 8):
                rtes = PopRoute[RandomOrder[p-8:p]]
                brks = PopBreaks[RandomOrder[p-8:p]]
                dists = TotalDist[0][RandomOrder[p-8:p]]
                idx = np.argmin(dists)    
                bestof8Route = rtes[idx]
                bestof8Break = brks[idx]
                selector = np.random.permutation(n)
                randlist = [selector[n/3], selector[2*n/3]]
                I = min(randlist)
                J = max(randlist)
                for k in range(8):
                    TempPopRoute[k] = bestof8Route
                    TempPopBreaks[k] = bestof8Break
                # Transformation 1        
                Temp = TempPopRoute[1][I]; TempPopRoute[1][I:J] = TempPopRoute[1][J:I:-1]; TempPopRoute[1][J] = Temp
                # Transformation 2        
                TempPopRoute[2][I], TempPopRoute[2][J] = TempPopRoute[2][J], TempPopRoute[2][I]
                # Transformation 3        
                TempPopRoute[3] = np.array(TempPopRoute[3].tolist()[0:I] + TempPopRoute[3].tolist()[I+1:J] + [TempPopRoute[3].tolist()[I]] + TempPopRoute[3].tolist()[J:])
                # Transformation 4
                TempPopBreaks[4] = rand_breaks(n, MinRoute, NRoutes, NBreaks)
                # Transformation 5        
                Temp = TempPopRoute[5][I]; TempPopRoute[5][I:J] = TempPopRoute[5][J:I:-1]; TempPopRoute[5][J] = Temp
                TempPopBreaks[5] = rand_breaks(n, MinRoute, NRoutes, NBreaks)
                # Transformation 6        
                TempPopRoute[6][I], TempPopRoute[6][J] = TempPopRoute[6][J], TempPopRoute[6][I]
                TempPopBreaks[6] = rand_breaks(n, MinRoute, NRoutes, NBreaks)
                # Transformation 7        
                TempPopRoute[7] = np.array(TempPopRoute[7].tolist()[0:I] + TempPopRoute[7].tolist()[I+1:J] + [TempPopRoute[7].tolist()[I]] + TempPopRoute[7].tolist()[J:])
                TempPopBreaks[7] = rand_breaks(n, MinRoute, NRoutes, NBreaks)        
                
                NewPopRoute[p-8:p] = TempPopRoute
                NewPopBreaks[p-8:p] = TempPopBreaks   
            PopRoute = NewPopRoute
            PopBreaks = NewPopBreaks
            
            if ConvergenceCounter >= ConvergenceDefinition:
                Converged = True
    
    if ShowResult:
        V = np.array((vertices))
        produce_plot(routing_coordinates(vertices, routing(OptRoute, OptBreak)), V, GlobalMin)
    
    print 'Alogorithm completed'
    
    return OptRoute, OptBreak, GlobalMin
    




## 0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0 RUN PROGRAM 0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0 ##


## INPUTS

turbine_locations = [[1,1],[1.3,3.4],[1.5,4.6],[2.2,2.8],[2.1,3.1],[2.2,4.3],[3.4,2.5],[3.6,3.7],[3.8,4.9],[4.1,2.2],[4.3,3.4],[4.5,4.6],[5.7,2.6],[5.9,3.1],[5.2,4.3],[6.4,2.5],[6.6,3.7],[6.8,4.9]]

substation_location = [[5,9]]

Capacity = 5

turbine_locations = convert_to_adnumber(turbine_locations)

substation_location = convert_to_adnumber(substation_location)


## RUN
    
Optimised_Routing = GA(turbine_locations, substation_location, Capacity)

dC_dX = differentiate(turbine_locations, construct_cost_matrix(substation_location + turbine_locations), Optimised_Routing[0], Optimised_Routing[1], len(turbine_locations))

print 'Completed'

print 'Minimimum Routing Length Discovered: ', int(Optimised_Routing[2])























