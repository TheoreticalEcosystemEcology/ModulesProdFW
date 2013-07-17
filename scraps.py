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

nx.draw(G,node_color=PopSize,with_labels=False,cmap=plt.cm.Greys,pos=pos,vmin=min(PopSize),vmax=max(PopSize))
plt.show()


plt.plot(Time,BioMass)
plt.show()

for n in G:
    print str(n.id)+' -- '+str(n.n)

print GetProd(G)