# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 12:09:03 2013

@author: dmc13
"""

import numpy as np
import scipy as sp
import scipy.spatial 
import itertools

### INPUTS ###

Substation = np.array([[0,0]])
Turbines = np.array([[1,2],[1,3],[1,4],[2,2],[2,3],[2,4],[3,2],[3,3],[3,4],[4,2],[4,3],[4,4],[9,9]])
C = 8

V = np.concatenate((Substation, Turbines))


def edge_length(V, A, B):
    '''measure the length of an edge joining two vertices, indexed A and B in V'''
    dist = (np.sqrt((V[A][0] - V[B][0])**2 + (V[A][1] - V[B][1])**2))
    return dist


## sort hull indices using (sparse) adjacency matrix graph stuff
#
def produce_hull(V):
    '''produce the convex hull of the vertices and convert into a route'''
    print V    
    hull = sp.spatial.qhull.Delaunay(V).convex_hull
    g = sp.sparse.csr_matrix((np.ones(hull.shape[0]),hull.T), shape=(hull.max()+1,)*2)
    sorted_hull = sp.sparse.csgraph.depth_first_order(g,hull[0,0],directed=False)[0]
    idx = np.where(sorted_hull == 0)[0]
    A = sorted_hull[idx:len(sorted_hull)]
    B = sorted_hull[0:idx]
    ordered_hull = np.concatenate((A, B))
    if edge_length(V, ordered_hull[0], ordered_hull[1]) > edge_length(V, ordered_hull[0], ordered_hull[-1]):
        inverted_hull = ordered_hull[::-1]
        new_turbine_order = inverted_hull[0: -1]        
        ordered_hull = np.concatenate((np.array([inverted_hull[-1]]), new_turbine_order))
        return ordered_hull
    else: return ordered_hull

R = []
    
while len(V) > 2:
#    print 'len', len(V)
    sorted_hull = produce_hull(V)
    r = []    
    for i in sorted_hull:
        v = V.tolist()
        r.append(v[i])
    R.append(r)
    V = np.array(list(itertools.compress(V, [i not in sorted_hull for i in range(len(V))])))
#    print 'len', len(V)
#    print 'S', Substation
#    print 'V', V
    if len(V) > 0:
        V = np.concatenate((Substation, V))
    else: break

print 'R', np.array(R)


#dist = (np.sqrt((coord_list[i][0] - coord_list[j][0])**2 + (coord_list[i][1] - coord_list[j][1])**2))


#if len(sorted_hull) < C:
#    Y = np.array(list(itertools.compress(V, [i not in sorted_hull for i in range(len(V))])))
#    Y = np.concatenate((Substation, Y))
#
#hull = sp.spatial.qhull.Delaunay(Y).convex_hull
#g = sp.sparse.csr_matrix((np.ones(hull.shape[0]),hull.T), shape=(hull.max()+1,)*2)
#sorted_hullY = sp.sparse.csgraph.depth_first_order(g,hull[0,0],directed=False)[0]
#sorted_hullY = np.append(sorted_hullY, sorted_hullY[0])
V = np.concatenate((Substation, Turbines))

### display with matplotlib
from matplotlib import pyplot as plt
plt.plot(V[:,0],V[:,1],'o')

#for i in range(len(R)):
#    for j in range(len(R[i])):
#        plt.plot(R[i][j][::0], R[i][j][::1], '-o')
    
for i in range(len(V)):
    plt.text(V[i][0], V[i][1], '%s' % (str(i)))



#plt.plot(Y[sorted_hullY,0],Y[sorted_hullY,1])

plt.show()