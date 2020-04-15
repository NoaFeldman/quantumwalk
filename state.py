# from tree import Graph, getGraphFromBaseGraph
import math
import numpy as np
import cmath
import random as rand


def getTreeNeighbors(site: int, lenTree: int):
    result = []
    if 2 * site < lenTree: # site is not a leaf
        result.append(2 * site)
        result.append(2 * site + 1)
    if site > 1:
        result.append(int(site / 2))
    return result


def getNumberOfLeaves(treeLen: int):
    return 2**int(math.log(treeLen, 2))

class Node:
    def __init__(self, region: str, name: int):
        self.region = region
        self.name = name

#                      1000 = 2 * 100 + 0
#                  100
#                      1001
#                10
#                      1010
#                  101
#                      1011
#  ... -2 -1 0 1
#                      1100
#                  110
#                      1101
#                11
#                      1110
#                  111
#                      1111
class Graph:
    def __init__(self, branchLen, tree1HopWeights, tree1Potentials, \
                 tree2HopWeights, tree2Potentials, permutation, connectingHopWeights):
        self.branchLen = branchLen
        self.tree1HopWeights = tree1HopWeights
        self.tree1Potentials = tree1Potentials
        self.tree2HopWeights = tree2HopWeights
        self.tree2Potentials = tree2Potentials
        self.permutation = permutation
        self.connectingHopWeights = connectingHopWeights


    def basicGraph(branchLen: int, treesDepth: int):
        minHopWeight = -5
        maxHopWeight = 5
        potentialResolution = 0.01
        potentials = [i * potentialResolution - 1 for i in range(int(2 / potentialResolution) + 1)]

        tree1HopWeights = []
        tree1Potentials = []
        for i in range(treesDepth):
            for j in range(2**i):
                tree1HopWeights.append(1 + rand.randint(minHopWeight, maxHopWeight) / (maxHopWeight - minHopWeight))
                tree1Potentials.append(potentials[rand.randint(0, len(potentials) - 1)])
        tree2HopWeights = []
        tree2Potentials = []
        for i in range(treesDepth):
            for j in range(2**i):
                tree2HopWeights.append(1 + rand.randint(minHopWeight, maxHopWeight) / (maxHopWeight - minHopWeight))
                tree2Potentials.append(potentials[rand.randint(0, len(potentials) - 1)])
        permutation = [i for i in np.random.permutation(getNumberOfLeaves(len(tree1HopWeights)) + 1) if i != 0]
        connectingHopWeights = []
        for i in range(2**(treesDepth - 1)):
            connectingHopWeights.append(1 + rand.randint(minHopWeight, maxHopWeight) / (maxHopWeight - minHopWeight))
        return Graph(branchLen, tree1HopWeights, tree1Potentials, \
                     tree2HopWeights, tree2Potentials, permutation, connectingHopWeights)

    def nearestNeighborsPairs(self):
        result = []
        for j in range(self.branchLen - 1):
            result.append(['inBranch', j, 'inBranch', j+1, 1])
        if self.branchLen > 0:
            result.append(['inBranch', -1, 'tree1', 0, 1])
        for j in range(len(self.tree1HopWeights)):
            # Sites in tree are numbered from 1 so I can use the rule in the drawing above
            neighbors = [n for n in getTreeNeighbors(site=j + 1, lenTree=len(self.tree1HopWeights)) if (n - 1) > j]
            for neighbor in neighbors:
                result.append(['tree1', j, 'tree1', neighbor - 1, self.tree1HopWeights[neighbor-1]])
        numberOfLeaves = getNumberOfLeaves(treeLen=len(self.tree1HopWeights))
        for l in range(1, numberOfLeaves + 1):
            result.append(['tree1', -l, 'tree2', -self.permutation[l - 1], self.connectingHopWeights[-l]])
        for j in [len(self.tree1HopWeights) - 1 - i for i in range(len(self.tree1HopWeights))]:
            neighbors = [n for n in getTreeNeighbors(site=j + 1, lenTree=len(self.tree2HopWeights)) if (n - 1) < j]
            for neighbor in neighbors:
                result.append(['tree2', j, 'tree2', neighbor - 1, self.tree2HopWeights[j]])
        if self.branchLen > 0:
            result.append(['tree2',  0, 'outBranch', 0, 1])
        for j in range(self.branchLen - 1):
            result.append(['outBranch', j, 'outBranch', j+1, 1])
        return result


def getGraphFromBaseGraph(graph: Graph, weightsDisorder: float, potentialDisorder: float):
    return Graph(graph.branchLen, [(1 + (w - 1) * weightsDisorder) for w in graph.tree1HopWeights], \
                 [v * potentialDisorder for v in graph.tree1Potentials],
                 [(1 + (w - 1) * weightsDisorder) for w in graph.tree2HopWeights], \
                 [v * potentialDisorder for v in  graph.tree2Potentials], \
                 graph.permutation, \
                 [(1 + (w - 1) * weightsDisorder) for w in graph.connectingHopWeights])


class State:
    def __init__(self, graph: Graph, k = 0.0, inBranchCoeffs = [], tree1Coeffs = [], tree2Coeffs = [], outBranchCoeffs = []):
        self.graph = graph
        if inBranchCoeffs != []:
            self.inBranchCoeffs = inBranchCoeffs
        else:
            self.inBranchCoeffs = [cmath.exp(1j * k * j) / math.sqrt(graph.branchLen) for j in range(graph.branchLen)]
        if tree1Coeffs != []:
            self.tree1Coeffs = tree1Coeffs
        else:
            self.tree1Coeffs = [0] * len(graph.tree1HopWeights)
        if tree2Coeffs != []:
            self.tree2Coeffs = tree2Coeffs
        else:
            self.tree2Coeffs = [0] * len(graph.tree2HopWeights)
        if outBranchCoeffs !=[]:
            self.outBranchCoeffs = outBranchCoeffs
        else:
            self.outBranchCoeffs = [0] * graph.branchLen

    def getPassingAmplitude(self):
        return self.tree2Coeffs[0]


def trotterStep(state: State, tau: float):
    localTerms = state.graph.nearestNeighborsPairs()
    # Start from originals and not zeros
    # use [i, j] = M[i, j]
    # try ilist = newTree1Coeffs, jList = ... and see if it works since list is mutable
    newInBranchCoeffs = state.inBranchCoeffs
    newTree1Coeffs = state.tree1Coeffs
    newTree2Coeffs = state.tree2Coeffs
    newOutBranchCoeffs = state.outBranchCoeffs

    newTree1Coeffs[0] *= cmath.exp(1j * tau / 2 * state.graph.tree1Potentials[0])

    for term in localTerms + [localTerms[- 1 - i] for i in range(len(localTerms))]:
        diag = math.cos(term[4] * tau / 2)
        offDiag = 1j * math.sin(term[4] * tau / 2)
        if term[0] == 'inBranch':
            regionI = newInBranchCoeffs
        elif term[0] == 'tree1':
            regionI = newTree1Coeffs
        elif term[0] == 'tree2':
            regionI = newTree2Coeffs
        else:
            regionI = newOutBranchCoeffs
        if term[2] == 'inBranch':
            regionJ = newInBranchCoeffs
        elif term[2] == 'tree1':
            regionJ = newTree1Coeffs
        elif term[2] == 'tree2':
            regionJ = newTree2Coeffs
        else:
            regionJ = newOutBranchCoeffs
        newI = diag * regionI[term[1]] + offDiag * regionJ[term[3]]
        newJ = diag * regionJ[term[3]] + offDiag * regionI[term[1]]
        [regionI[term[1]], regionJ[term[3]]] = [newI, newJ]

        if term[0] == 'tree1' and term[2] == 'tree1':
            newTree1Coeffs[term[3]] *= cmath.exp(1j * tau / 2 * state.graph.tree1Potentials[term[3]])
        if term[0] == 'tree2' and term[2] == 'tree2':
            newTree2Coeffs[term[3]] *= cmath.exp(1j * tau / 2 * state.graph.tree2Potentials[term[3]])

    return State(state.graph, \
                 inBranchCoeffs=newInBranchCoeffs, tree1Coeffs=newTree1Coeffs, \
                 tree2Coeffs=newTree2Coeffs, outBranchCoeffs=newOutBranchCoeffs)

def multiplyByConstant(state: State, c: float):
    result =  state
    for i in range(len(result.inBranchCoeffs)):
        result.inBranchCoeffs[i] *= c
    for i in range(len(result.outBranchCoeffs)):
        result.outBranchCoeffs[i] *= c
    for i in range(len(state.tree1Coeffs)):
        result.tree1Coeffs[i] *= c
    for i in range(len(state.tree2Coeffs)):
        result.tree2Coeffs[i]
    return result


def getOverlap(state1: State, state2: State):
    result =  0
    for i in range(len(state1.inBranchCoeffs)):
        result += np.conj(state1.inBranchCoeffs[i]) * state2.inBranchCoeffs[i]
    for i in range(len(state1.outBranchCoeffs)):
        result += np.conj(state1.outBranchCoeffs[i]) * state2.outBranchCoeffs[i]
    for i in range(len(state1.tree1Coeffs)):
        result += np.conj(state1.tree1Coeffs[i]) * state2.tree1Coeffs[i]
    for i in range(len(state1.tree2Coeffs)):
        result += np.conj(state1.tree2Coeffs[i]) * state2.tree2Coeffs[i]
    return result

import sys
import pickle

d = 6 #int(sys.argv[1])
graph = Graph.basicGraph(branchLen=0, treesDepth=d)
with open('basicRandomTree0_' + str(d) + '.data', 'wb') as filehandle:
    pickle.dump(graph, filehandle)
#
# with open('basicRandomTree_' + str(d) + '.data', 'rb') as filehandle:
#     graph = pickle.load(filehandle)
hoppingDisorder = 0 #float(sys.argv[2])
potentialDisorder = 0.2 #float(sys.argv[3])
graph = getGraphFromBaseGraph(graph, hoppingDisorder, potentialDisorder)
k = math.acos(0.5) #math.acos(float(sys.argv[4] / 2))
state = State(graph = graph, k=k)
state.tree1Coeffs[0] = 1
state = multiplyByConstant(state,  1/math.sqrt(np.real(getOverlap(state, state))))


if potentialDisorder != 0:
    tau = 1e-2 * potentialDisorder
else:
    tau = 1e-2
amps = []
for i in range(int(1 * d * math.pi/tau)):
    state = trotterStep(state, tau)
    amps.append(state.getPassingAmplitude())
with open('randomTrees_' + str(d) +
          '_' + str(hoppingDisorder) + '_' + str(potentialDisorder) + '_' + \
          str(k) + '.data', 'wb') as filehandle:
    pickle.dump(max([abs(amp) for amp in amps]), filehandle)
