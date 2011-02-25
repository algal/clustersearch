
# global variables defining defaults.
alphabet = "CGAU"
colors = ["white", "black", "red"]

def mut(point):
    """Returns set of all neighbors of POINT"""
    neighbors = list()
    for pos in xrange(len(point)): # for every position  ...
        ch = point[pos]            # ... find its alternative values
        others = alphabet.replace(ch,'')
        for newch in others:       # ... and build resulting mutants
            variant = newch.join((point[:pos],point[pos+1:]))
            neighbors.append(variant)
    return neighbors

def colorOf(point):
    """Given a point, computes its color"""
    import random
    return random.choice(colors)

def search(root):
    """Searches breadth-first from root node, recording color of every
    node observed, but traversing onward only from "white" nodes.

    Nodes x and y are connected if (x in mut(y)). Each node x has
    color colorOf(x).

    Returns a dictionary of the observed nodes and their colors."""
    SPECIAL_COLOR = colors[0]
    observed_nodes = dict()
    observed_nodes[root]=SPECIAL_COLOR
    traversed = list()
    to_traverse = list()
    to_traverse.append(root)
    while to_traverse:
        # process head of queue of nodes still to be traversed
        cursor = to_traverse.pop(0)
        print("Traversing from node %s" % cursor)
        # compute the neighbors not previously observed
        neighbors = set(mut(cursor)) - set(observed_nodes.keys())
        print("\tFound %s had neighbors %s of which the previously unobserved were %s" % (cursor, mut(cursor), neighbors))
        if neighbors:
            # else, assign them colors
            newnodes = dict([(s,colorOf(s)) for s in neighbors])
            print("\tObserved these nodes to be colored: %s" % newnodes)
            # plan to visit only the special ones later
            to_traverse.extend([i for i in neighbors if newnodes[i]==SPECIAL_COLOR])
            # add them to the db of observed nodes
            observed_nodes.update(newnodes)
            # mark the current node as traversed
            traversed.append(cursor)
    return (observed_nodes,traversed)


def dorun(length,alphabetsize,numOfColors):
    """Does a search run.

    NOTE:
    - modifies global variables alphabet and colors
    - alphabetsize and numOfColors must be < 26"""
    global alphabet
    global colors
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[0:alphabetsize]
    colors = "abcdefghijklmnopqrstuvwxyz"[0:numOfColors]
    origin = ''.join([alphabet[0] for i in range(length)])
    return search(origin)[0]

def calcStats(nodes):
    results = dict()
    SPECIAL_COLOR = colors[0]
    results["cluster_size"] = len([i for i in nodes.values() if i==SPECIAL_COLOR])
    results["perimeter_of_the_cluster"] = len([i for i in nodes.values() if i != SPECIAL_COLOR])
    results["number_of_unique_colors"] = len(set([i for i in nodes.values() if i != SPECIAL_COLOR]))
    return results


