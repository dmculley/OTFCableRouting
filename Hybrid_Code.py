# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 14:01:35 2013

@author: dmc13
"""

from matplotlib import pyplot as plt
import math
import copy
import numpy as np
from ad import adnumber
import ad
from ad.admath import *
import random


class CableCostGA(object):
    
    def __init__(self, turbine_locations, substation_location = [[0,0]], capacity = 7, pop_size = 80, num_iter = 200, convergence_definition = 20, converged = False, show_prog = False, show_result = False): 
        self.turbine_locations = turbine_locations
        self.substation_location = substation_location
        self.capacity = capacity
        self.pop_size = pop_size
        self.num_iter = num_iter
        self.convergence_definition = convergence_definition
        self.show_prog = show_prog
        self.show_result = show_result


    def convert_to_adnumber(self, coordinate_list):
        '''Convert the location vectors from floats into adnumbers to enable differentiation'''
        adnumber_coordinate_list = []
        for i in range(len(coordinate_list)):
            adnumber_coordinate_list.append([adnumber(coordinate_list[i][0]), adnumber(coordinate_list[i][1])])
        coordinate_list = adnumber_coordinate_list
        return coordinate_list
        
    
    def differentiate(self, best_chromosome):
        '''Differentiate the length of the routing w.r.t. the position of the turbines, produces a n x 2 array, ((dC/dx1, dC/dy1), ...)'''
        print 'Determining dC/dX: Performing automatic differentiation...'
        turbine_locations = self.convert_to_adnumber(self.turbine_locations)
        substation_location = self.convert_to_adnumber(self.substation_location)
        vertices = substation_location + turbine_locations
        C = self.construct_cost_matrix(vertices)
        route = best_chromosome[0]
        breaks = best_chromosome[1]    
        rting = self.routing(route, breaks)
        total_dist = self.routing_distance(rting, C)
        n = len(turbine_locations)
        dC_dX = np.zeros((n,2)) 
        for i in range(n):
            dC_dX[i][0] = total_dist.d(turbine_locations[i][0])
            dC_dX[i][1] = total_dist.d(turbine_locations[i][1])
        print 'Automatic differentiation complete...'
        return dC_dX
    
    
    def construct_cost_matrix(self, vertices):
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
    
    
    def produce_plot(self, vertices, best_chromosome, global_min):
        '''Display with matplotlib'''
        opt_route = best_chromosome[0]
        opt_breaks = best_chromosome[1]
        R = self.routing_coordinates(vertices, self.routing(opt_route, opt_breaks))
        print R
        print self.routing(opt_route, opt_breaks)
        plt.title('Total distance='+str(global_min))
        V = np.array((vertices))
        plt.plot(V[:,0],V[:,1],'o')
        for i in range(len(R)):
            plt.plot(R[i][:,0], R[i][:,1], '-')
        for i in range(len(V)):
            plt.text(V[i][0], V[i][1], '%s' % (str(i)))
        plt.axis('equal')
        plt.show()
        
    
    def rand_breaks(self, n, min_route, n_routes, n_breaks):
        '''Produces a list of random, but valid, route-breaks'''
        RB = [np.random.random_integers(min_route, self.capacity)]
        for i in range(1, n_breaks):
            RB.append(np.random.random_integers(min_route, self.capacity) + RB[i-1])
        if RB[-1] < (n - self.capacity):
            short = (n - self.capacity) - RB[-1]
            add_each = int(np.ceil(0.5 + short / len(RB)))
            for i in range(len(RB)):
                RB[i] = RB[i] + add_each * (i+1)
        return RB    
    
        
    def initialise_population(self, n, min_route, n_routes, n_breaks):
        '''Randomly produces an initial population'''
        popbreaks = []
        poproutes = []
        for i in range(self.pop_size):
            popbreaks.append(self.rand_breaks(n, min_route, n_routes, n_breaks))
            poproutes.append(random.sample(range(1, n+1), n))
        return poproutes, popbreaks
    
    
    def routing(self, route, breaks):
        '''Combine route and breaks into an array of the routes described as vertices in the order in which they are toured'''
        rting = [[0] + route[0:breaks[0]]]
        if len(breaks) > 1:
            for f in range(1, len(breaks)):
                rting.append([0] + route[breaks[f-1]:breaks[f]])
        rting.append([0] + route[breaks[-1]:])
        return rting
    
    
    def routing_distance(self, routing, C):
        '''Return the geometric length of the routing'''
        d = 0        
        for r in range(len(routing)):
            for v in range(len(routing[r]) - 1):
                d += C[routing[r][v]][routing[r][v+1]]
        return d
    
    
    def routing_coordinates(self, vertices, rting):
        '''Convert a routing expressed in indexes to a routing expressed in coordinates''' 
        for i in range(len(rting)):
            for j in range(len(rting[i])):
                rting[i][j] = vertices[rting[i][j]]
            rting[i] = np.array(rting[i])
        return rting
    
    
    def evaluate_population(self, pop, C, n):
        '''produce a list of the functional evaluations of each member of the population'''
        D = []
        for i in range(self.pop_size):
            D.append(self.routing_distance(self.routing(pop[0][i], pop[1][i]), C))
        return D
        
        
    def clarke_wright(self, route, breaks, vertices):
        '''apply clarke and wright algorithm to chromosome'''
        rting = self.routing(route, breaks)
        rting_coords = self.routing_coordinates(vertices, rting)      
        rting = []
        for i in range(len(rting_coords)):
            temp_r = rting_coords[i].tolist()
            del temp_r[0]
            print temp_r
            CW = Clarke_Wright(temp_r, self.substation_location)
            print temp_r
            temp_r = CW.run()
            print temp_r
            del temp_r[0]
            for j in range(len(temp_r)):
                temp_r[j] = self.turbine_locations.index(temp_r[j]) + 1
            rting += temp_r
        return rting
    
    
    def transformations(self, bestof8route, bestof8breaks, n, min_route, n_routes, n_breaks, vertices, iteration):
        '''Returns a copy of the original chromosome and 7 mutations'''
        selector = random.sample(range(1, n), n-1)
        randlist = [selector[n/3], selector[2*n/3]]
        I = min(randlist)
        J = max(randlist)
        temp_pop_route = []
        temp_pop_breaks = []
        temp_pop_route.append(bestof8route)
        temp_pop_breaks.append(bestof8breaks)
        # transformation 1  
        trans_1 = copy.deepcopy(bestof8route)
        trans_1[I], trans_1[J] = trans_1[J], trans_1[I]
        temp_pop_route.append(trans_1)
        temp_pop_breaks.append(bestof8breaks)
        ## transformation 2
        trans_2 = copy.deepcopy(bestof8route)
        Temp = trans_2[I:J]; Temp = Temp[::-1]; trans_2 = trans_2[:I] + Temp + trans_2[J:]
        temp_pop_route.append(trans_2)
        temp_pop_breaks.append(bestof8breaks)
        ## transformation 3
        trans_3 = copy.deepcopy(bestof8route)
        trans_3.remove(I); trans_3.insert(J, I)
        temp_pop_route.append(trans_3)
        temp_pop_breaks.append(bestof8breaks)
        ## transformation 4
        if iteration < 15:
            trans_4 = copy.deepcopy(bestof8route)
            print 't4', trans_4
            trans_4 = self.clarke_wright(trans_4, bestof8breaks, vertices)
            temp_pop_route.append(trans_4)
            temp_pop_breaks.append(bestof8breaks)
        else:
            temp_pop_route.append(bestof8route)
            temp_pop_breaks.append(self.rand_breaks(n, min_route, n_routes, n_breaks))
        ## transformation 5
        if iteration < 15:
            trans_5 = random.sample(range(1, n+1), n)
            print trans_5
            trans_5 = self.clarke_wright(trans_5, bestof8breaks, vertices)
            temp_pop_route.append(trans_5)
            temp_pop_breaks.append(bestof8breaks)
        else:
            temp_pop_route.append(trans_1)
            temp_pop_breaks.append(self.rand_breaks(n, min_route, n_routes, n_breaks))
        ## transformation 6        
        temp_pop_route.append(trans_2)
        temp_pop_breaks.append(self.rand_breaks(n, min_route, n_routes, n_breaks))
        ## transformation 7        
        temp_pop_route.append(trans_3)
        temp_pop_breaks.append(self.rand_breaks(n, min_route, n_routes, n_breaks))        
        
        return temp_pop_route, temp_pop_breaks
        
     
    def breed_population(self, pop, D, n, min_route, n_routes, n_breaks, vertices, iteration):
        '''Create a new population based upon the best performing specimens of the old population'''
        n = len(self.turbine_locations)        
        chrome_pop = []
        for i in range(self.pop_size):
            chrome_pop.append(pop[0][i] + pop[1][i] + [D[i]])
        random.shuffle(chrome_pop)
        new_pop_route = []
        new_pop_breaks = []
        for i in range(0, self.pop_size, 8):
            Dists = [] 
            for j in range(8):
                Dists.append(chrome_pop[i+j][-1])
            MIndx = min(xrange(len(Dists)),key=Dists.__getitem__)
            bestof8Chrome = chrome_pop[MIndx+i]
            bestof8route = bestof8Chrome[:n]
            bestof8breaks = bestof8Chrome[n:-1]
            new_8 = self.transformations(bestof8route, bestof8breaks, n, min_route, n_routes, n_breaks, vertices, iteration)
            new_pop_route += new_8[0]
            new_pop_breaks += new_8[1]
        return new_pop_route, new_pop_breaks
        
            
    def run_GA(self):
        '''Run the genetic algorithm'''
        np.random.seed(100)
        random.seed(100)        
        vertices = self.substation_location + self.turbine_locations
        n_routes = int(math.ceil(float(len(vertices)) / self.capacity))
        min_route = len(vertices) / n_routes
        n = len(self.turbine_locations)
        n_breaks = n_routes - 1
        converged = False
        convergence_counter = 0
        print 'Running Genetic Algorithm...'
        pop = self.initialise_population(n, min_route, n_routes, n_breaks)
        C = self.construct_cost_matrix(vertices)
        D = self.evaluate_population(pop, C, n)    
        MIndx = min(xrange(len(D)),key=D.__getitem__)    
        global_min = copy.deepcopy(D[MIndx])
        best_chromosome = [copy.deepcopy(pop[0][MIndx]), copy.deepcopy(pop[1][MIndx])]
        dist_history = [global_min]
    
        for iteration in range(self.num_iter):
            if not converged:
                if self.show_prog:
                    print 'Current Routing Length', global_min
                pop = self.breed_population(pop, D, n, min_route, n_routes, n_breaks, vertices, iteration)
                D = self.evaluate_population(pop, C, n)
                MIndx = min(xrange(len(D)),key=D.__getitem__)
                if D[MIndx] < global_min:
                    global_min = copy.deepcopy(D[MIndx])
                    best_chromosome = [copy.deepcopy(pop[0][MIndx]), copy.deepcopy(pop[1][MIndx])]
                    dist_history.append(global_min)
                    convergence_counter = 0
                else: convergence_counter += 1
                if convergence_counter > self.convergence_definition:
                    converged = True
        if self.show_result:
            self.produce_plot(vertices, best_chromosome, global_min)
        print 'GA Complete, Routing Length is: ', global_min
        return best_chromosome, global_min
        
    def compute_cable_cost(self):
        out = self.run_GA()
        self.best_chromosome = out[0]
        return out[1]
        
    def compute_cable_cost_derivative(self):
        dC_dX = self.differentiate(self.best_chromosome)
        return dC_dX
        
        
        
        
## 0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0o0 ##




class Clarke_Wright(object):
    
    def __init__(self, turbine_locations, substation_location, plot_output = False):
        self.turbine_locations = turbine_locations
        self.substation_location = substation_location
        self.plot_output = plot_output
      
    def routing_coordinates(self, vertices, rting):
        '''Convert a routing expressed in indexes to a routing expressed in coordinates'''        
        for i in range(len(rting)):
            rting[i] = vertices[rting[i]]
        return rting
    

    def routing_indexes(self, R, vertices):
        '''Convert a routing expressed in tuples to a routing expressed in indexes'''
        chromosome = [0]
        while len(chromosome) <= len(R):
            for i in range(len(R)):
                if R[i][0] == chromosome[-1]:
                    chromosome.append(R[i][1])
                if R[i][1] == chromosome[-1]:
                    chromosome.append(R[i][0])
        return chromosome
    
    
    def produce_plot(self, R, vertices):
        '''Display with matplotlib'''
        R = self.routing_indexes(R, vertices)
        R = self.routing_coordinates(vertices, R)
        R = np.array(([R])) 
        V = np.array((vertices))
        plt.plot(V[:,0],V[:,1],'o')
        plt.plot(R[0][:,0], R[0][:,1], '-')
        plt.text(V[0][0], V[0][1], '%s' % (str(0)))
        plt.axis('equal')
        plt.show() 
        

    def construct_cost_matrix(self, vertices):
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
        

    def make_graph(self, Rgraph, no_depot = True):
        '''Produces a dictionary representing the graph of the form (vertex: [connected vertices], ...) can remove the depot so as to seperate the routes from one another'''
        if no_depot:
            remove_list=[]
            for edge in Rgraph:
                if edge[0] == 0 or edge[1] == 0:
                    remove_list.append(edge)
            for edge in remove_list: 
                Rgraph.remove(edge)
        vertices = []
        for a_tuple in Rgraph:
            vertices.extend(list(a_tuple))
        vertices = list(set(vertices))
        nVertices = len(vertices)
        G = {}
        for i in range(0,nVertices):
            G[vertices[i]]=[]
        for edge in Rgraph:
            G[edge[1]].append(edge[0])
            G[edge[0]].append(edge[1])
        for vertex in G:
            G[vertex] = list(set(G[vertex]))
        return G
        

    def initialise_R(self, n):
        '''produce the initial edge set - all clients connected to their nearest depot'''
        R = []
        # For each client vertex, identify the nearest depot and add (v,d) to initial edge list R
        for index in range(1, n):
            R.append((index, 0))
        return R
        

    def initialise_G(self, R, n):
        '''produce the graph of the initial edge set'''
        G = []
        Rgraph = copy.deepcopy(R)
        G = self.make_graph(Rgraph)
        return G
        

    def construct_savings_list(self, C, n):
        '''Constructs a list of the savings in making every connection and orders them'''
        S = []    
        # For each possible pair of clients in the graph, calculate the cost saving from joining with edge, and add to list S
        R = self.initialise_R(n)        
        for vertex1 in range(1, n):
            for vertex2 in range(1, n):
                if not vertex1 == vertex2:
                    depot = [j[1] for i, j in enumerate(R) if j[0] == vertex1]
                    saving = C[vertex1, depot] - C[vertex1, vertex2]
                    S.append((vertex1, vertex2, saving))
        # Sort the list S of cost savings in decreasing order
        S.sort(key = lambda s:s[2], reverse = True)
        return S


    def find_path_length(self, graph, start, end, path=[]):
        '''Returns a list of vertices on a path in order that they are connected'''        
        path = path + [start]
        if start == end:
            return path
        if not graph.has_key(start):
            return None
        for node in graph[start]:
            if node not in path:
                newpath = self.find_path_length(graph, node, end, path)
                if newpath: return newpath
        return None


    def neighbour_depot(self, edge, R, Vd):
        '''Ensures that k is neighbouring a depot - ensures no short-circuiting'''      
        nhbrDepot = False
        if (edge[0], 0) or (0, edge[0]) in R:
            nhbrDepot = True
        return nhbrDepot

        
    def one_neighbour(self, edge, R):
        '''Prevent route branching by checking that u has only one neighbour'''
        oneNhbr = True
        for edge_temp in R:
            if edge_temp[1] != 0 and edge_temp[0] == edge[0]:
                oneNhbr = False
            if edge_temp[1] == edge [1]:
                oneNhbr = False
        return oneNhbr

    
    def perform_checks(self, S, G, R, n, edge):
        '''Performs viability checks on proposed saving'''
        if not (self.find_path_length(G, edge[0], edge[1])==None):           
            return False
        if not self.neighbour_depot(edge, R, range(1, n)):
            return False
        if not self.one_neighbour(edge, R):
            return False
        return True

        
    def k_d(self, edge, R):
        '''Defines the edge whose removal is proposed'''
        if (edge[0],0) in R:
            k_d = (edge[0],0)
        elif (0,edge[0]) in R:
            k_d = (0,edge[0])
        else: print 'ERROR', edge, R
        return k_d

        
    def perform_merge(self, edge, R, k_d):
        '''Removes the superfluous old edge and inserts the new edge & updates the routing graph'''
        R.remove(k_d)
        R.append((edge[0], edge[1]))
        Rgraph = copy.deepcopy(R)
        G = self.make_graph(Rgraph)
        return G

        
    def run(self):
        vertices = self.substation_location + self.turbine_locations
        n = len(self.turbine_locations)+1
        C = self.construct_cost_matrix(vertices)
        S = self.construct_savings_list(C, n)
        R = self.initialise_R(n)
        G = self.initialise_G(R, n)
        while (not S == []) and (S[0][2][0] > 0):
            edge = (S[0][0],S[0][1])
            allowable = self.perform_checks(S, G, R, n, edge)
            if allowable:
                G = self.perform_merge(edge, R, self.k_d(edge, R))
            else:                 
                del S[0]
        if self.plot_output:
            self.produce_plot(R, vertices)
        R = self.routing_indexes(R, vertices)
        R = self.routing_coordinates(vertices, R)
        return R

        
        
if __name__ == '__main__':
    turbine_locations = [[1,1],[1.3,3.4],[1.5,4.6],[2.2,2.8],[2.1,3.1],[2.2,4.3],[3.4,2.5],[3.6,3.7],[3.8,4.9],[4.1,2.2],[4.3,3.4],[4.5,4.6],[5.7,2.6],[5.9,3.1],[5.2,4.3],[6.4,2.5],[6.6,3.7],[6.8,4.9]]
    show_prog = True
    
    CC = CableCostGA(turbine_locations, show_prog = True, show_result = True)
    
    print CC.compute_cable_cost()
    print CC.compute_cable_cost_derivative()








