#!/usr/bin/env python
# encoding: utf-8

import networkx as nx
import numpy as np
import scipy as sp

import getopt
import sys
import os

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
    UpdateTrophicLevel(G)            
    return G

def SimulWeb(graph,p):
    Record = [0,0,0,0,0]
    nRecord = 0
    for t in xrange(p['timesteps']):
        for sp in graph:
            if graph.out_degree(sp) == 0:
                ## Primary producer
                sp.dn += sp.n * (1 - sp.n / float(p['K']))
                ## Find all of its predators
                for pred in graph.pred[sp]:
                    degPred = graph.out_degree(pred)
                    Omega = 1 / float(degPred)
                    preyBiomass = 0
                    ## Get the prey relative consumption for all preys of this predator
                    for prey_of_pred in graph.succ[pred]:
                        preyBiomass += Omega * np.power(prey_of_pred.n, p['h'])
                    Fij = (Omega * np.power(sp.n,p['h'])) / float( (np.power(p['Bo'],p['h']) + p['c'] * pred.n * np.power(p['Bo'],p['h'])) + preyBiomass )
                    Mc = np.power(p['Z'], pred.tl)
                    Mp = np.power(p['Z'], sp.tl)
                    if Mp == 0:
                        MRatio = 1
                    else:
                        MRatio = Mc/float(Mp)
                    sp.dn -= (p['xi'] * np.power( MRatio , p['scexpo']) * p['y'] * pred.n * Fij) / float(p['eij'])
            else:
                ## If not a primary producer
                Mc = np.power(p['Z'], sp.tl)
                Mp = np.power(p['Z'], 1)
                if Mp == 0:
                    MRatio = 1
                else:
                    MRatio = Mc/float(Mp)
                sp.dn -= p['xi'] * np.power( MRatio , p['scexpo']) * sp.n
                if graph.in_degree(sp) > 0:
                    ## Find all its predators
                    for pred in graph.pred[sp]:
                        degPred = graph.out_degree(pred)
                        Omega = 1 / float(degPred)
                        preyBiomass = 0
                        ## Get the prey relative consumption for all preys of this predator
                        for prey_of_pred in graph.succ[pred]:
                            preyBiomass += Omega * np.power(prey_of_pred.n, p['h'])
                        Fij = (Omega * np.power(sp.n,p['h'])) / float( (np.power(p['Bo'],p['h']) + p['c'] * pred.n * np.power(p['Bo'],p['h'])) + preyBiomass )
                        Mc = np.power(p['Z'], pred.tl)
                        Mp = np.power(p['Z'], sp.tl)
                        if Mp == 0:
                            MRatio = 1
                        else:
                            MRatio = Mc/float(Mp)
                        sp.dn -= (p['xi'] * np.power( MRatio, p['scexpo']) * p['y'] * pred.n * Fij) / float(p['eij'])
                ## Find all of its preys
                Omega = 1 / float(graph.out_degree(sp))
                ## Get the prey relative consumption for all preys of this predator
                preyBiomass = 0
                for prey_of_pred in graph.succ[sp]:
                    preyBiomass += Omega * np.power(prey_of_pred.n, p['h'])
                for prey_of_pred in graph.succ[sp]:
                    Fij = (Omega * np.power(prey_of_pred.n,p['h'])) / float( (np.power(p['Bo'],p['h']) + p['c'] * sp.n * np.power(p['Bo'],p['h'])) + preyBiomass )
                    Mc = np.power(p['Z'], prey_of_pred.tl)
                    Mp = np.power(p['Z'], sp.tl)
                    if Mp == 0:
                        MRatio = 1
                    else:
                        MRatio = Mc/float(Mp)
                    sp.dn += (p['ax'] / float(p['ar'])) * np.power( MRatio, p['scexpo']) * p['y'] * sp.n * Fij
        for n in graph:
            n.n += n.dn * p['SCALAR']
            n.dn = 0
        if (p['timesteps']-t) <= p['record']:
            tRecord = GetProd(graph)
            nRecord += 1
            for i in xrange(len(tRecord)):
                Record[i] += tRecord[i]
    for i in xrange(len(tRecord)):
                Record[i] = Record[i] / float(nRecord) 
    return Record

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

def main(argv=None):
    ## PARAMS
    p = {'y':4.0, 'eij':0.85, 'K':1.0, 'r':1.0, 'ar':1.0, 'ax':0.88, 'xi':0.88, 'scexpo':-0.25, 'Z':2.0, 'c':1.0, 'h':1.0, 'Bo':0.5, 'SCALAR':0.02, 'timesteps':10000, 'record':1000}
    webID = 1
    repltodo = 5
    ## FIX OPTIONS
    opts, args = getopt.getopt(sys.argv[1:], "h", ["w=","t=","rec=","repl="])
    for o in opts:
        if o[0] in ('--w'):
            webID = int(o[1])
        elif o[0] in ('--t'):
            p['timesteps'] = int(o[1])
        elif o[0] in ('--rec'):
            p['record'] = int(o[1])
        elif o[0] in ('--repl'):
            repltodo = int(o[1])
    ##
    fname= 'webs/w'+str(webID)+'.txt'
    Gr = DiGraphFromList(fname)
    Motifs = MotifCount(Gr)
    Motifs = map(str,Motifs)
    for repl in xrange(repltodo):
        for Z in [0.0, 2.0]:
            p['Z'] = Z
            for n in Gr:
                n.n = np.random.uniform(low=0.05, high=1.0, size=1)[0]
            #print "Web "+str(webID)+", Z = "+str(p['Z'])+", replicate "+str(repl+1)
            Out = SimulWeb(Gr, p)
            ## Write the web-wide results
            Fname = 'output/web'+str(webID)+'-WEB.dat'
            f = open(Fname, 'a')
            Out = map(str,Out)
            f.write(str(webID)+' ')
            f.write(str(repl+1)+' ')
            f.write(str(p['Z'])+' ')
            f.write(str(len(Gr.edges()))+' ')
            for motif in Motifs:
                f.write('{0} '.format(motif))
            for record in Out:
                f.write('{0} '.format(record))
            f.write('\n')
            f.close()

if __name__ == "__main__":
    main()
