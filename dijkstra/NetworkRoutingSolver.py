#!/usr/bin/python3


from CS312Graph import *
import time


class Array:
    nodes = []

    # Time complexity: O(1)
    def insert(self, node):
        self.nodes.append(node)

    # Time complexity: O(1)
    def decreaseKey(self, i, newDist):
        self.nodes[i].dist = newDist

    # Time complexity: O(n) where n is the number of nodes
    def deleteMin(self):
        minNode = self.nodes[0]
        minDist = 100000

        for i in self.nodes:
            if i.dist < minDist:
                minDist = i.dist
                minNode = i

        self.nodes.remove(minNode)
        return minNode


class Heap:
    def __init__(self):
        self.nodes = [0]
        self.currentSize = 0

    def bubbleUp(self, i):
        while i // 2 > 0:
            j = i // 2
            if self.nodes[i].dist < self.nodes[j].dist:
                temp = self.nodes[j]
                self.nodes[j] = self.nodes[i]
                self.nodes[i] = temp
            i = j


    def insert(self, node):
        self.nodes.append(node)
        self.currentSize += 1
        self.bubbleUp(self.currentSize)

    def bubbleDown(self, i):
        while (i * 2) <= self.currentSize:
            mc = self.minChild(i)
            if self.nodes[i].dist > self.nodes[mc].dist:
                temp = self.nodes[i]
                self.nodes[i] = self.nodes[mc]
                self.nodes[mc] = temp
            i = mc

    def minChild(self, i):
        i *= 2
        if i + 1 > self.currentSize:
            return i
        else:
            if self.nodes[i].dist < self.nodes[i + 1].dist:
                return i
            else:
                return i + 1


    def decreaseKey(self, i, newDist):
        self.nodes[i].dist = newDist

    def deleteMin(self):
        toReturn = self.nodes[1]
        self.nodes[1] = self.nodes[self.currentSize]
        self.currentSize = self.currentSize - 1
        self.nodes.pop()
        self.bubbleDown(1)
        return toReturn


class Node(object):
    id = None
    prev = None
    dist = 10000


def makeNode(id, prev, dist):
    node = Node()
    node.id = id
    node.prev = prev
    node.dist = dist
    return node


class NetworkRoutingSolver:
    def __init__(self):
        pass

    def initializeNetwork(self, network):
        assert (type(network) == CS312Graph)
        self.network = network

    def getShortestPath(self, destIndex):
        self.dest = destIndex
        # TODO: RETURN THE SHORTEST PATH FOR destIndex
        #       INSTEAD OF THE DUMMY SET OF EDGES BELOW
        #       IT'S JUST AN EXAMPLE OF THE FORMAT YOU'LL
        #       NEED TO USE
        shortestPaths = self.toBeAccessed
        path_edges = []
        total_length = 0

        startNode = self.network.nodes[self.source]
        currentNode = None
        for i in self.network.nodes:
            if i.node_id == destIndex:
                currentNode = i
                break
        if currentNode is None:
            return {'cost' : float('inf'), 'path' : path_edges}

        reachedSource = False

        while not reachedSource:
            for i in shortestPaths:
                if currentNode.node_id == i.id:

                    previousNode = None
                    for j in self.network.nodes:
                        if j.node_id == i.prev:
                            previousNode = j
                            break

                    if previousNode is None:
                        return {'cost': float('inf'), 'path': path_edges}

                    addEdge = None
                    for e in previousNode.neighbors:
                        if e.src == previousNode and e.dest == currentNode:
                            addEdge = e

                    path_edges.append((previousNode.loc, currentNode.loc, '{:.0f}'.format(addEdge.length)))
                    total_length += addEdge.length

                    if previousNode == startNode:
                        reachedSource = True
                    else:
                        currentNode = previousNode

        #edges_left = 3
        #while edges_left > 0:
        #    edge = node.neighbors[2]
        #    path_edges.append((edge.src.loc, edge.dest.loc, '{:.0f}'.format(edge.length)))
        #    total_length += edge.length
        #    node = edge.dest
        #    edges_left -= 1
        return {'cost': total_length, 'path': path_edges}

    def computeShortestPaths(self, srcIndex, use_heap=False):
        self.source = srcIndex
        t1 = time.time()
        # TODO: RUN DIJKSTRA'S TO DETERMINE SHORTEST PATHS.
        #       ALSO, STORE THE RESULTS FOR THE SUBSEQUENT
        #       CALL TO getShortestPath(dest_index)

        shortestPaths = []

        if use_heap:
            myQueue = Heap()

            # makeQueue functionality; all are set to maximum distance except the source
            for i in range(len(self.network.nodes)):
                tempNode = self.network.nodes[i]

                if tempNode.node_id == srcIndex:
                    addNode = makeNode(tempNode.node_id, tempNode.node_id, 0)
                else:
                    addNode = makeNode(tempNode.node_id, None, 100000)

                myQueue.insert(addNode)
                shortestPaths.append(addNode)

            # run dijkstra on the entire queue (priority BFS)
            while len(myQueue.nodes) > 1:
                minNode = myQueue.deleteMin()
                graphCounterpart = self.network.nodes[minNode.id]

                for i in graphCounterpart.neighbors:
                    totalDist = minNode.dist + i.length

                    # find node in distance array
                    for j in range(len(shortestPaths)):
                        if i.dest.node_id == shortestPaths[j].id:
                            indexToCompare = j
                            break

                    # compare distances
                    if totalDist < shortestPaths[indexToCompare].dist:
                        shortestPaths[indexToCompare].dist = totalDist
                        shortestPaths[indexToCompare].prev = i.src.node_id
                        for j in range(1, len(myQueue.nodes)):
                            if myQueue.nodes[j].id == i.dest.node_id:
                                myQueue.decreaseKey(j, totalDist)
                                break


        else:
            myQueue = Array()

            # makeQueue functionality; all are set to maximum distance except the source
            for i in range(len(self.network.nodes)):
                tempNode = self.network.nodes[i]

                if tempNode.node_id == srcIndex:
                    addNode = makeNode(tempNode.node_id, tempNode.node_id, 0)
                else:
                    addNode = makeNode(tempNode.node_id, None, 100000)

                myQueue.insert(addNode)
                shortestPaths.append(addNode)

            #run dijkstra on the entire queue (priority BFS)
            while len(myQueue.nodes) > 0:
                minNode = myQueue.deleteMin()
                graphCounterpart = self.network.nodes[minNode.id]

                for i in graphCounterpart.neighbors:
                    totalDist = minNode.dist + i.length

                    #find node in distance array
                    for j in range(len(shortestPaths)):
                        if i.dest.node_id == shortestPaths[j].id:
                            indexToCompare = j
                            break

                    #compare distances
                    if totalDist < shortestPaths[indexToCompare].dist:
                        shortestPaths[indexToCompare].dist = totalDist
                        shortestPaths[indexToCompare].prev = i.src.node_id
                        for j in range(len(myQueue.nodes)):
                            if myQueue.nodes[j].id == i.dest.node_id:
                                myQueue.decreaseKey(j, totalDist)
                                break

        #save for getShortestPath()
        self.toBeAccessed = shortestPaths

        t2 = time.time()
        return (t2 - t1)
