from Bio import Phylo
from io import StringIO
import pandas as pd
import re
import numpy as np
import argparse

def parse_args():
    """Define the Arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('--RelateSamples',type=str, default=None)
    parser.add_argument('--DerivedFile',type=str, default=None)
    parser.add_argument('--out',type=str,default=None)
    return parser.parse_args()

def newick_to_list(newick_string):
    name_regex = re.compile(r"[^(),:]+")  # Matches node names

    def traverse(node, parent=None):
        if node is None:
            return
        for child in node:
            if isinstance(child, str):
                df.loc[len(df)] = [parent, child]
            else:
                name_match = name_regex.search(child.name)
                node_name = name_match.group() if name_match else None
                df.loc[len(df)] = [parent, node_name]
                traverse(child.clades, parent=node_name)

    df = pd.DataFrame(columns=["parent", "child"])

    tree = Phylo.read(StringIO(newick_string), "newick")

    traverse(tree.clade, parent=None)

    Parents = []
    Childs = []
    p = (df.iloc[:,0])
    c = (df.iloc[:,1])
    for i in p: Parents.append(i)
    for i in c: Childs.append(i)
    for k in range(len(Parents)):
        if Parents[k] == None:
            Parents[k] = "root"
    return([Childs, Parents])

def CalculateDerivedBelow(List, nodee, nodeeindex):
    if not (nodee in List[1]):
        List[3][nodeeindex] = List[2][nodeeindex]
    else:
        child1index = (np.where(np.array(List[1]) == nodee))[0][0]
        child2index = (np.where(np.array(List[1]) == nodee))[0][1]
        CalculateDerivedBelow(List, List[0][child1index], child1index)
        CalculateDerivedBelow(List, List[0][child2index], child2index)
        List[3][nodeeindex] = List[3][child1index] + List[3][child2index]

def CalculateAncestralBelow(List, nodee, nodeeindex):
    if not (nodee in List[1]):
        List[4][nodeeindex] = 1 - List[2][nodeeindex]
    else:
        child1index = (np.where(np.array(List[1]) == nodee))[0][0]
        child2index = (np.where(np.array(List[1]) == nodee))[0][1]
        CalculateAncestralBelow(List, List[0][child1index], child1index)
        CalculateAncestralBelow(List, List[0][child2index], child2index)
        List[4][nodeeindex] = List[4][child1index] + List[4][child2index]

def OneTreeToList2(Newick, IsDerived):
    from Bio import Phylo
    from io import StringIO
    import numpy as np

    # 1) Leggi l’albero e calcola le profondità normalizzate
    treee = Phylo.read(StringIO(Newick), 'newick')
    depths_dict = treee.depths()
    max_depth = max(depths_dict.values())
    depths = [max_depth - treee.distance(clade) for clade in treee.get_nonterminals()]

    # 2) Individua quali cladi sono “derived” vs “ancestral”
    derived_count = sum(IsDerived)
    IsDerivedIndices = [i for i, v in enumerate(IsDerived) if v == 1]
    AncestralI = []
    DerivedI  = []
    DERIVEDTMRCA = None

    for idx, clade in enumerate(treee.get_nonterminals()):
        descendants = [int(leaf.name) for leaf in clade.get_terminals()]
        if any(d not in IsDerivedIndices for d in descendants):
            AncestralI.append(idx)
        else:
            DerivedI.append(idx)
            if len(descendants) == derived_count:
                DERIVEDTMRCA = clade

    # 3) Cerca l’indice del clade “pieno di derived” in modo sicuro
    MIXEDLIENAGEindex = None
    if DERIVEDTMRCA is not None:
        for idx, clade in enumerate(treee.get_nonterminals()):
            if clade is DERIVEDTMRCA:
                MIXEDLIENAGEindex = idx
                break

    # 4) Se trovato, rimuovi dalla lista ancestrali e aggiungi ai derived
    if MIXEDLIENAGEindex is not None:
        AncestralI.pop(MIXEDLIENAGEindex)
        DerivedI.append(MIXEDLIENAGEindex)
    else:
        print("[WARN] no full-derived clade found; skipping mixed-lineage adjustment")

    # 5) Costruisci le stringhe di output (comma-delimited + newline)
    Ancestralbranches = np.sort(np.array(depths)[AncestralI])
    Derivedbranches  = np.sort(np.array(depths)[DerivedI])

    AncString = ",".join(map(str, Ancestralbranches)) + "\n"
    DerString = ",".join(map(str, Derivedbranches))  + "\n"

    return (DerString, AncString)

def calculateageofsamples(newickstring, NewIsDerived):
    # First, read in the tree
    tree = Phylo.read(StringIO(newickstring), 'newick')
    # Get the leaves
    leaves = tree.get_terminals()
    L = len(leaves)
    # Create a L by L matrix to store the pairwise TMRCA
    # For each pair of leaves, get the TMRCA
    MaxDepth = -1
    for i in tree.depths().keys():
        if tree.depths()[i] > MaxDepth:
            MaxDepth = tree.depths()[i]

    LeafNames = []
    for i in range(L):
        LeafNames.append(int(leaves[i].name))
    # Return the matrix

    LeafDepths = []
    for i in range(L):
        LeafDepths.append((MaxDepth - ((tree.depths())[tree.common_ancestor(leaves[i].name, leaves[i].name)])))
    IsDerivedReOrder = []
    for i in range(len(LeafNames)):
        IsDerivedReOrder.append(NewIsDerived[LeafNames[i]])

    DERIVEDTIMES = []
    ANCESTRALTIMES = []
    for i in range(len(LeafDepths)):
        if LeafDepths[i] > 1.0:
            if IsDerivedReOrder[i] == 1:
                DERIVEDTIMES.append(float(LeafDepths[i]))
            else:
                ANCESTRALTIMES.append(float(LeafDepths[i]))
    return DERIVEDTIMES,ANCESTRALTIMES

if __name__ == "__main__":
    args = parse_args()
    isder = np.loadtxt(args.DerivedFile)
    IsDerived = isder
    f = open(args.RelateSamples, "r")

    TotalStrings = []
    Lines = (f.readlines())[1:]
    newicktree = (Lines[0].split("\t"))
    newicktree = newicktree[len(newicktree)-1]
    newicktree = re.sub(r':\d+(\.\d+)?,', ',', newicktree)
    newicktree = re.sub(r':\d+(\.\d+)?\)', ')', newicktree)

    letter = "a"
    alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
    def nextletter(let):
        lastletter = let[len(let) - 1]
        if lastletter != "z":
            ii = alphabet.index(lastletter)
            return(let[0:(len(let)-1)] + alphabet[ii+1])
        else:
            return(let + "a")

    while True:
        try:
            aaaa = (newicktree.index("),"))
            newicktree = newicktree[0:aaaa]  + ")" + letter + ","  + newicktree[(aaaa+2):len(newicktree)]
            letter = nextletter(letter)
        except ValueError:
            try:
                bbbb = (newicktree.index("))"))
                newicktree = newicktree[0:bbbb] + ")" + letter + ")" + newicktree[(bbbb+2):len(newicktree)]
                letter = nextletter(letter)
            except ValueError:
                break

    List = newick_to_list(newicktree)

    TotalDerived = np.sum(IsDerived)
    TotalAncestral = len(IsDerived) - TotalDerived
    #remove branch length by substituting :____, by , nad replacing :____) by )
    #Then, add a letter after each closed parenthesis.
    List.append([-1] * len(List[0]))

    nnnnn = np.array(List[0])
    for i in range(len(IsDerived)):
        List[2][(np.where(nnnnn == str(i)))[0][0] ] = IsDerived[i]

    List.append ( [0] * len(List[0])) # number of derived below
    List.append ( [0] * len(List[0])) # number of ancestral below

    child1index = (np.where(np.array(List[1]) == "root"))[0][0]
    child2index = (np.where(np.array(List[1]) == "root"))[0][1]
    CalculateDerivedBelow(List,  List[0][child1index], child1index)
    CalculateDerivedBelow(List,  List[0][child2index], child2index)

    child1index = (np.where(np.array(List[1]) == "root"))[0][0]
    child2index = (np.where(np.array(List[1]) == "root"))[0][1]
    CalculateAncestralBelow(List,  List[0][child1index], child1index)
    CalculateAncestralBelow(List,  List[0][child2index], child2index)

    TotalError = [0] * len(List[0])
    for i in range(len(TotalError)):
        TotalError[i] = TotalDerived - List[3][i] + List[4][i]
    minflips = min(TotalError)
    if minflips == 0:
        print("Infinite sites assumption satisfied. No allele flips necessary.")
        NewIsDerived = IsDerived
    else:
        minnn = min(TotalError)
        if minnn == 1:
            print("Infinite sites assumption not satisfied. Flipping " + str(int(minnn)) + " leaf out of "  + str(len(IsDerived)) + " leaves.")
        else:
            print("Infinite sites assumption not satisfied. Flipping " + str(int(minnn)) + " leaves out of " + str(len(IsDerived)) + " leaves.")
        FlippingIndex = TotalError.index(minnn)

        NewIsDerived = [0] * len(IsDerived)
        Descendants = []

        def FindDescendants(List, nodee, nodeeindex):
            if not (nodee in List[1]):
                Descendants.append(int(List[0][nodeeindex]))
            else:
                child1index = (np.where(np.array(List[1]) == nodee))[0][0]
                child2index = (np.where(np.array(List[1]) == nodee))[0][1]
                FindDescendants(List, List[0][child1index], child1index)
                FindDescendants(List, List[0][child2index], child2index)

        FindDescendants(List,  List[0][FlippingIndex], FlippingIndex)

        for i in Descendants:
            NewIsDerived[i] = 1
        for i in range(len(NewIsDerived)):
            if NewIsDerived[i] > IsDerived[i]:
                print("Flipped leaf " + str(i) + " from ancestral (0) to derived (1).")
            if NewIsDerived[i] < IsDerived[i]:
                print("Flipped leaf " + str(i) + " from derived (1) to ancestral (0).")


    newicktreetemp = (Lines[0].split("\t"))
    newicktreetemp = newicktreetemp[len(newicktreetemp)-1]
    DERIVEDTIMES,ANCESTRALTIMES = calculateageofsamples(newicktreetemp[:-1], NewIsDerived)

    if len(DERIVEDTIMES) + len(ANCESTRALTIMES) > 0:
        if len(DERIVEDTIMES) == 1:
            print(str(len(DERIVEDTIMES)) + " ancient haplotype found with derived allele." )
        else:
            print(str(len(DERIVEDTIMES)) + " ancient haplotypes found with derived allele." )
        if len(ANCESTRALTIMES) == 1:
            print(str(len(ANCESTRALTIMES)) + " ancient haplotype found with ancestral allele." )
        else:
            print(str(len(ANCESTRALTIMES)) + " ancient haplotypes found with ancestral allele." )
        DERIVEDTIMES.sort()
        ANCESTRALTIMES.sort()
        derivedtimestring=""
        for iiivv in range(len(DERIVEDTIMES)):
            derivedtimestring = derivedtimestring + str(DERIVEDTIMES[iiivv]) + ";"
        derivedtimestring = derivedtimestring[:-1] + "\n"
        ancestraltimestring=""
        for iiivv in range(len(ANCESTRALTIMES)):
            ancestraltimestring = ancestraltimestring + str(ANCESTRALTIMES[iiivv]) + ";"
        ancestraltimestring = ancestraltimestring[:-1] + "\n"
        TotalStrings.append(derivedtimestring)
        TotalStrings.append(ancestraltimestring)
    for i in range(len(Lines)):
        newicktree = (Lines[i].split("\t"))
        newicktree = newicktree[len(newicktree)-1]
        TotalStrings.extend(OneTreeToList2((newicktree[:-1]), NewIsDerived))
    
    if minflips > 0:
        tempderived = []
        for i in NewIsDerived: tempderived.append(str(i) + "\n")
        f = open(args.out + "_derived.txt", "w+")
        f.writelines(tempderived)
        f.close()
    f = open(args.out + "_times.txt", "w+")
    f.writelines(TotalStrings)
    f.close()
