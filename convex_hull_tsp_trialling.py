# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 12:09:03 2013

@author: dmc13
"""

'''Hueristic for the capcitated VRP problem using convex hulls'''

import numpy as np
import scipy as sp
import scipy.spatial 
import itertools
from matplotlib import pyplot as plt

### INPUTS ###

Substation = np.array([[0,0]])
Turbines = np.array([[1,2],[1,3],[1,4],[2,2],[2,3],[2,4],[3,2],[3,3],[3,4],[4,2],[4,3],[4,4],[9,9]])
C = 5

V = np.concatenate((Substation, Turbines))


def edge_length_V(V, A, B):
    '''measure the length of an edge joining two vertices, indexed A and B in V'''
    dist = (np.sqrt((V[A][0] - V[B][0])**2 + (V[A][1] - V[B][1])**2))
    return dist
    
def edge_length_R(R, i, A, B):
    '''measure the length of an edge joining two vertices, indexed A and B in R'''
    dist = (np.sqrt((R[i][A][0] - R[i][B][0])**2 + (R[i][A][1] - R[i][B][1])**2))
    return dist


def produce_hull(V, C):
    '''produce the convex hull of the vertices, sort using (sparse) adjacency matrix and convert into a route'''
    # at least 4 points are required to determine the convex hull    
    if len(V) < 4:
        return range(len(V))
    # if there are sufficient vertices, calculate the hull, starting and ending at the substation
    hull = sp.spatial.qhull.Delaunay(V).convex_hull
    g = sp.sparse.csr_matrix((np.ones(hull.shape[0]),hull.T), shape=(hull.max()+1,)*2)
    sorted_hull = sp.sparse.csgraph.depth_first_order(g,hull[0,0],directed=False)[0]
    idx = np.where(sorted_hull == 0)[0]
    A = sorted_hull[idx:len(sorted_hull)]
    B = sorted_hull[0:idx]
    ordered_hull = np.concatenate((A, B))
    # break the longest return edge to the substation, break the route at max capacity
    if len(ordered_hull) > C:
        if edge_length_V(V, ordered_hull[0], ordered_hull[1]) > edge_length_V(V, ordered_hull[0], ordered_hull[-1]):
            inverted_hull = ordered_hull[::-1]
            new_turbine_order = inverted_hull[0: -1]        
            ordered_hull = np.concatenate((np.array([inverted_hull[-1]]), new_turbine_order))
            ordered_hull = ordered_hull[0:C]
            return ordered_hull
        else:
            ordered_hull = ordered_hull[0:C]
            return ordered_hull
    else:
        if edge_length_V(V, ordered_hull[0], ordered_hull[1]) > edge_length_V(V, ordered_hull[0], ordered_hull[-1]):
            inverted_hull = ordered_hull[::-1]
            new_turbine_order = inverted_hull[0: -1]        
            ordered_hull = np.concatenate((np.array([inverted_hull[-1]]), new_turbine_order))
            return ordered_hull
        else: 
            return ordered_hull


def populate_routing(V, Substation, C):
    '''produce and array or routes, each itself an array of coordinates in order of connection'''
    R = []
    while len(V) > 0:
        sorted_hull = produce_hull(V, C)
        r = []    
        for i in sorted_hull:
            v = V.tolist()
            r.append(v[i])
        R.append(np.array(r))
        V = np.array(list(itertools.compress(V, [i not in sorted_hull for i in range(len(V))])))
        if len(V) > 0:
            V = np.concatenate((Substation, V))
        else: break
    print R
    R = np.array(R)
    return R


def produce_plot(R, V, dist):
    '''display with matplotlib'''
    plt.title('Total distance='+str(dist))
    plt.plot(V[:,0],V[:,1],'o')
    for i in range(len(R)):
        plt.plot(R[i][:,0], R[i][:,1], '-')
    for i in range(len(V)):
        plt.text(V[i][0], V[i][1], '%s' % (str(i)))
    plt.axis('equal')
    plt.show()
    
def find_routing_length(R):
    routing_length = []
    for i in range(len(R)):
        for j in range(0, len(R[i])-1):
            routing_length.append(edge_length_R(R, i, j, j+1))
    routing_length = sum(routing_length)
    return routing_length



R = populate_routing(V, Substation, C)
print R
RL = find_routing_length(R)

produce_plot(R, V, RL)

