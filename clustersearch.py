
# global variables defining defaults.
alphabet = "CGAU"
colors = ["white", "black", "red"]

debug = False
def log(str):
    if debug:
        print(str)

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
    node observed, but traversing onward only from 'white' nodes.

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
        log("Traversing from node %s" % cursor)
        # compute the neighbors not previously observed
        neighbors = set(mut(cursor)) - set(observed_nodes.keys())
        log("\tFound %s had neighbors %s of which the previously unobserved were %s" % (cursor, mut(cursor), neighbors))
        if neighbors:
            # "discover" their colors
            newnodes = dict([(s,colorOf(s)) for s in neighbors])
            log("\tObserved these nodes to be colored: %s" % newnodes)
            # plan to visit only the special ones later
            to_traverse.extend([i for i in neighbors if newnodes[i]==SPECIAL_COLOR])
            # add them to the db of observed nodes
            observed_nodes.update(newnodes)
            # mark the current node as traversed
            traversed.append(cursor)
    #we ignore traversal path for now.
    return observed_nodes


def doRun(length,alphabetsize,numOfColors):
    """Do a search run.

    Returns dictionary over all points in a cluster (i.e., neutral
    network) and its mutational neighborhood, 'picked' from random
    string graph with string LENGTH, built from an alphabet of
    ALPHABETSIZE, where every node is painted one color out of
    NUMOFCOLORS.

    'Picks' this cluster by doing a breadth-first search, assigning
    colors randomly as it progresses.

    NOTE:
    - modifies global variables alphabet and colors
    - alphabetsize and numOfColors must be < 26"""
    global alphabet
    global colors
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[0:alphabetsize]
    colors = list(range(numOfColors))
    origin = ''.join([alphabet[0] for i in range(length)])
    return search(origin)

def calcStats(nodes):
    "Calculates various descriptive statistics for one search"
    SPECIAL_COLOR = colors[0]
    results = dict()
    results["cluster_size"] = len([i for i in nodes.values() if i==SPECIAL_COLOR])
    results["perimeter_size"] = len([i for i in nodes.values() if i != SPECIAL_COLOR])
    results["perimeter_color_count"] = len(set([i for i in nodes.values() if i != SPECIAL_COLOR]))
    results["exits_size"] = len([(x,y) for x in nodes.keys() for y in nodes.keys() if (x < y) and (x in mut(y)) and (nodes[x] != nodes[y])]);
    return results


def makeRunRecord(length,alphabetsize,numOfColors):
    "Returns inputs and resulting statistics for a single search -- i.e., a row."
    results = calcStats(doRun(length,alphabetsize,numOfColors))
    row = list()
    row = [length,alphabetsize,numOfColors]
    row.append(results["cluster_size"])
    row.append(results["perimeter_size"])
    row.append(results["perimeter_color_count"])
    row.append(results["exits_size"])
    return row

def do_many_runs(length,alphabetsize,numOfColors,runcount):
    """Does multiple identical runs with a set of input params"""
    results = list()
    for runid in range(1,runcount+1):
        row = list()
        row.append(runid)
        row.extend( makeRunRecord(length,alphabetsize,numOfColors) )
        results.append(row)
    return results

def mean(seq):
    return float(sum(seq)) / float(len(seq))

def calcAverages(length,alphabetsize,numOfColors,runcount=200):
    """Calculates average stats over multiple runs.
    [length,alphabetsize,numOfColors,
    avg(cluster_size),avg(perimiter_size),avg(perimiter_color_count),avg(exits_size)]
    """
    rows = do_many_runs(length,alphabetsize,numOfColors,runcount)
    result = list()
    result.extend(rows[0][1:4])    # drop runcount, cp length,alphabetsize,numOfColors
    result.append(mean( [row[4] for row in rows] ) )
    result.append(mean( [row[5] for row in rows] ) )
    result.append(mean( [row[6] for row in rows] ) )
    result.append(mean( [row[7] for row in rows] ) )
    return result

def calcAveragesOverRangedInputs(length,alphabetsize,numOfColors,runcount):
    """Calculates average stats separately for specified inputs or ranges of inputs.

    For example, calcAveragesOverRangedInputs(3,4,5,100)
    does 100 runs with length=3,alphabetsize=4,numOfColors=5.

    and, calcAveragesOverRangedInputs(range(3,11),4,5,100)
    does 100 runs each for every length in range(3,11) with alphabetsize=4,numOfColors=5.    
    """
    def rangify(x):
        if type(x) == type(1):
            return range(x,x+1)
        else:
            return x
    length = rangify(length)
    alphabetsize = rangify(alphabetsize)
    numOfColors = rangify(numOfColors)

    return [calcAverages(mylength,myalphabetsize,mynumOfColors,runcount)
            for mylength in length
            for myalphabetsize in alphabetsize
            for mynumOfColors in numOfColors]
    
def saveAsCSV(table, filename):
    """Saves a list of lists as csv.

    For instance, can be used with output of
    calcAveragesOverNumOfColors or of do_many_runs
    """
    import csv
    fileobj = open(filename, "w")
    my_writer  = csv.writer(fileobj)
    my_writer.writerows(table)
