#!/usr/bin/env python
# encoding: utf-8

import networkx as nx
import numpy as np
import scipy as sp
import scipy.stats as spst

import getopt
import sys
import os

import progressbar as pb

def GetProd(graph):
    bm = 0
    nprod = 0
    bmprod = 0
    ncons = 0
    bmcons = 0
    for sp in graph:
        bm += sp.n
        if graph.out_degree(sp) == 0:
            bmprod += sp.n
            nprod += 1
        else:
            bmcons += sp.n
            ncons += 1
    return [len(graph),bm,bm/float(len(graph)),bmprod/float(nprod),bmcons/float(ncons)]

class species:
    def __init__(self,uname,n):
        self.id = uname
        self.n = n
        self.dn = 0
        self.tl = 1
        self.pos = (0,0)
        self.s1 = 0
        self.s2 = 0
        self.s3 = 0
        self.s4 = 0
        self.s5 = 0

def UpdateTrophicLevel(graph):
    for n in graph:
        if graph.out_degree(n) == 0:
            ## Primary producer
            for alt in graph:
                if graph.out_degree(alt) > 0:
                    if nx.has_path(graph, alt, n):
                        TL = nx.shortest_path_length(graph, alt, n) + 1
                        if (alt.tl == 1) or (alt.tl > TL):
                            alt.tl = TL
    return 0

def DiGraphFromList(fname):
    G = nx.DiGraph()
    Dat = np.loadtxt(fname+'.txt')
    Bio = np.loadtxt(fname+'.biomass')
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
        SpObj.append(species(Sp,Bio[int(Sp)-1]))
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
    UpdateTrophicLevel(G)            
    return G

def MotifCount(graph):
    # Number of each motif, third order, simple and double link
    s1 = 0
    s2 = 0
    s3 = 0
    s4 = 0
    s5 = 0
    d1 = 0
    d2 = 0
    d3 = 0
    d4 = 0
    d5 = 0
    d6 = 0
    d7 = 0
    d8 = 0
    for n1 in graph:
        for n2 in graph:
            for n3 in graph:
                #### Simple linkage motifs
                # Motif S1 : Linear food chain,
                if (graph.has_edge(n1,n2) and (not graph.has_edge(n2,n1))) and (graph.has_edge(n2,n3) and (not graph.has_edge(n3,n2))) and (not graph.has_edge(n3,n1)) and (not graph.has_edge(n1,n3)):
                    s1 += 1
                    n1.s1 += 1
                    n2.s1 += 1
                    n3.s1 += 1
                # Motif S2 : Omnivory
                if (graph.has_edge(n1,n2) and (not graph.has_edge(n2,n1))) and (graph.has_edge(n2,n3) and (not graph.has_edge(n3,n2))) and (not graph.has_edge(n3,n1)) and (graph.has_edge(n1,n3)):
                    s2 += 1
                    n1.s2 += 1
                    n2.s2 += 1
                    n3.s2 += 1
                # Motif S3 : Directed loop
                if (graph.has_edge(n1,n2) and (not graph.has_edge(n2,n1))) and (graph.has_edge(n2,n3) and (not graph.has_edge(n3,n2))) and (graph.has_edge(n3,n1)) and (not graph.has_edge(n1,n3)):
                    s3 += 1 
                    n1.s3 += 1
                    n2.s3 += 1
                    n3.s3 += 1
                # Motif S4 : Exploitative competition
                if (graph.has_edge(n1,n2) and (not graph.has_edge(n2,n1))) and (graph.has_edge(n3,n2) and (not graph.has_edge(n2,n3))) and (not graph.has_edge(n3,n1)) and (not graph.has_edge(n1,n3)):
                    s4 += 1
                    n1.s4 += 1
                    n2.s4 += 1
                    n3.s4 += 1
                # Motif S5 : Apparent competition
                if (graph.has_edge(n1,n2) and (not graph.has_edge(n2,n1))) and (graph.has_edge(n1,n3) and (not graph.has_edge(n3,n1))) and (not graph.has_edge(n3,n2)) and (not graph.has_edge(n2,n3)):
                    s5 += 1
                    n1.s5 += 1
                    n2.s5 += 1
                    n3.s5 += 1
                #### Double linkage motifs
                # Motif D1
                if (graph.has_edge(n1,n2) and graph.has_edge(n2,n1)) and (graph.has_edge(n1,n3) and (not graph.has_edge(n3,n1))) and (graph.has_edge(n2,n3) and (not graph.has_edge(n3,n2))):
                    d1 += 1
                # Motif D2
                if (graph.has_edge(n1,n2) and graph.has_edge(n1,n3)) and (graph.has_edge(n2,n3) and graph.has_edge(n3,n2)) and (not graph.has_edge(n2,n1) and (not graph.has_edge(n3,n1))):
                    d2 += 1
                # Motif D3
                if (graph.has_edge(n1,n2) and graph.has_edge(n2,n1)) and (graph.has_edge(n1,n3) and (not graph.has_edge(n3,n1))) and (not graph.has_edge(n2,n3) and (not graph.has_edge(n3,n2))):
                    d3 += 1
                # Motif D4
                if (graph.has_edge(n1,n2) and (not graph.has_edge(n1,n3))) and (graph.has_edge(n2,n3) and graph.has_edge(n3,n2)) and (not graph.has_edge(n2,n1) and (not graph.has_edge(n3,n1))):
                    d4 += 1
                # Motif D5
                if (graph.has_edge(n1,n2) and (graph.has_edge(n2,n1))) and (graph.has_edge(n2,n3) and (not graph.has_edge(n3,n2))) and (graph.has_edge(n3,n1)) and (not graph.has_edge(n1,n3)):
                    d5 += 1
                # Motif D6
                if (graph.has_edge(n1,n2) and graph.has_edge(n2,n1)) and (graph.has_edge(n1,n3) and graph.has_edge(n3,n1)) and (graph.has_edge(n3,n2) and graph.has_edge(n2,n3)):
                    d6 += 1
                # Motif D7
                if (graph.has_edge(n1,n2) and (not graph.has_edge(n2,n1))) and (graph.has_edge(n1,n3) and graph.has_edge(n3,n1)) and (graph.has_edge(n3,n2) and graph.has_edge(n2,n3)):
                    d7 += 1
                # Motif D8
                if ((not graph.has_edge(n1,n2)) and (not graph.has_edge(n2,n1))) and (graph.has_edge(n1,n3) and graph.has_edge(n3,n1)) and (graph.has_edge(n3,n2) and graph.has_edge(n2,n3)):
                    d8 += 1
    return [s1,s2,s3,s4,s5,d1,d2,d3,d4,d5,d6,d7,d8]

## Analysis

widgets = ['Analyses: ', pb.Percentage(), ' ', pb.Bar(),' ', pb.ETA()]
progress = pb.ProgressBar(widgets=widgets)
for webID in progress(xrange(1,114)):
    fname= 'webs/EP'+str(webID)
    Gr = DiGraphFromList(fname)
    UpdateTrophicLevel(Gr)
    Motifs = MotifCount(Gr)
    Motifs = map(str,Motifs)
    Out = GetProd(Gr)
    Fname = 'output/EPweb'+str(webID)+'-WEB.dat'
    f = open(Fname, 'w')
    Out = map(str,Out)
    f.write(str(webID)+' ')
    f.write(str(0)+' ')
    f.write(str(0)+' ')
    f.write(str(len(Gr.edges()))+' ')
    for motif in Motifs:
        f.write('{0} '.format(motif))
    for record in Out:
        f.write('{0} '.format(record))
    TLs = [n.tl for n in Gr]
    Deg = [Gr.degree(n) for n in Gr]
    OutDeg = [Gr.out_degree(n) for n in Gr]
    InDeg = [Gr.in_degree(n) for n in Gr]
    nPP = np.sum([1 for od in OutDeg if od == 0])/float(len(OutDeg))
    nTP = np.sum([1 for ind in InDeg if ind == 0])/float(len(InDeg))
    TLinfo = [np.mean(TLs),np.var(TLs),np.median(TLs),np.max(TLs),spst.skew(Deg),spst.skew(OutDeg),spst.skew(InDeg),nPP,nTP]
    for record in TLinfo:
        f.write('{0} '.format(record))
    f.write('\n')
    f.close()