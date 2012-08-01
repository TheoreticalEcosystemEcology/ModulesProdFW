#!/usr/bin/env python
# encoding: utf-8

import networkx as nx
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

class species:
    def __init__(self,id,n):
        self.id = id
        self.n = n
        self.dn = 0
        self.pos = (0,0)

def DiGraphFromList(fname):
    G = nx.DiGraph()
    Dat = np.loadtxt(fname)
    # We loop through the file and record the species names 
    SpNames = []
    for Line in Dat:
        SpNames.append(str(int(Line[0])))
        SpNames.append(str(int(Line[1])))
    # We keep the unique species
    UniqueSp = list(set(SpNames))
    # Each unique species becomes a node in the graph
    SpObj = []
    for Sp in UniqueSp:
        SpObj.append(species(Sp,np.random.random()))
    G.add_nodes_from(SpObj)
    # We loop through the link list
    for Line in Dat:
        Pred = str(int(Line[0]))
        Prey = str(int(Line[1]))
        for sp in G:
            if sp.id == Pred:
                PredObj = sp
            if sp.id == Prey:
                PreyObj = sp
        G.add_edge(PredObj, PreyObj)            
    return G
    
G = DiGraphFromList('webs/w12.txt')

## Get position
for n in G:
        if G.out_degree(n) == 0:
            n.pos = (np.random.random(),0)
        else:
            Paths= []
            for alt in G:
                if G.out_degree(alt) == 0:
                    if nx.has_path(G, n, alt) == True:
                        Paths.append(nx.shortest_path_length(G, n, alt))
            n.pos = (np.random.random(),min(Paths))

pos = {}
for n in G.nodes():
        pos[n] = n.pos

PopSize = [n.n for n in G]

nx.draw_shell(G,node_color=PopSize,with_labels=False,cmap=plt.cm.Greys,vmin=min(PopSize),vmax=max(PopSize))
plt.show()