#!/usr/bin/python3
import copy

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time
import numpy as np
from TSPClasses import *
import heapq
import itertools


# Reduces a cost matrix to ensure that each row and column either contains a 0 or is all math.inf
# Time complexity: O(n^2) because each for loop iterates n times, and each calls min() which has O(n) complexity
# (https://stackoverflow.com/questions/35386546/big-o-of-min-and-max-in-python). I'm not sure what the time complexity
# of rcm[i, :] or rcm[:, i] is, but I think it would at most be O(n) as well, so the whole function is O(n^2) overall.
# Space complexity: Doesn't use any new storage other than ints
def reduceMatrix(rcm, bound, n):
    for i in range(n):
        rowMin = min(rcm[i, :])
        if rowMin == 0 or rowMin == math.inf:
            continue
        rcm[i, :] -= rowMin
        bound += rowMin
    for i in range(n):
        colMin = min(rcm[:, i])
        if colMin == 0 or colMin == math.inf:
            continue
        rcm[:, i] -= colMin
        bound += colMin
    return {'bound': bound, 'rcm': rcm}


class TSPSolver:
    def __init__(self, gui_view):
        self._scenario = None

    def setupWithScenario(self, scenario):
        self._scenario = scenario

    ''' <summary>
        This is the entry point for the default solver
        which just finds a valid random tour.  Note this could be used to find your
        initial BSSF.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of solution, 
        time spent to find solution, number of permutations tried during search, the 
        solution found, and three null values for fields not used for this 
        algorithm</returns> 
    '''

    def defaultRandomTour(self, time_allowance=60.0):
        results = {}
        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        bssf = None
        start_time = time.time()
        while not foundTour and time.time() - start_time < time_allowance:
            # create a random permutation
            perm = np.random.permutation(ncities)
            route = []
            # Now build the route using the random permutation
            for i in range(ncities):
                route.append(cities[perm[i]])
            bssf = TSPSolution(route)
            count += 1
            if bssf.cost < np.inf:
                # Found a valid route
                foundTour = True
        end_time = time.time()
        results['cost'] = bssf.cost if foundTour else math.inf
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        return results

    ''' <summary>
        This is the entry point for the greedy solver, which you must implement for 
        the group project (but it is probably a good idea to just do it for the branch-and
        bound project as a way to get your feet wet).  Note this could be used to find your
        initial BSSF.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution, 
        time spent to find best solution, total number of solutions found, the best
        solution found, and three null values for fields not used for this 
        algorithm</returns> 
    '''

    def greedy(self, time_allowance=60.0):
        pass

    ''' <summary>
        This is the entry point for the branch-and-bound algorithm that you will implement
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution, 
        time spent to find best solution, total number solutions found during search (does
        not include the initial BSSF), the best solution found, and three more ints: 
        max queue size, total number of states created, and number of pruned states.</returns>
    '''

    def branchAndBound(self, time_allowance=60.0):
        # create all necessary variables
        cities = self._scenario.getCities()
        ncities = len(cities)
        startTime = time.time()
        bssf = self.defaultRandomTour()['soln']
        lastUpdate = time.time()

        solutionCount = 0
        totalStateCount = 1
        prunedStateCount = 0
        maxSize = 1
        priorityq = []

        # Create cost matrix and get initial reduced cost matrix (rcm)
        # Time: O(n^2) where n = num of cities. The first for loop will iterate n times, and the inner loop will also
        # iterate n times, meaning n^2 total.
        # Space: numpy.zeros() returns an array, and while I couldn't find any definitive source about the space
        # complexity of python arrays, I did find a StackOverflow answer which says python arrays are complexity O(n)
        # where n is the size of the array.
        # (https://stackoverflow.com/questions/45330006/what-is-the-time-complexity-and-space-complexity-of-array-1)
        # So if that is correct, the space complexity here would be O(n^2) where n = num of cities.
        rcm = np.zeros((ncities, ncities))
        for i in range(ncities):
            city = cities[i]
            for j in range(ncities):
                if i == j:
                    rcm[i, j] = math.inf
                else:
                    rcm[i, j] = city.costTo(cities[j])
        # Time complexity: O(n^2) (detailed in function definition)
        reducedMatrix = reduceMatrix(rcm, 0, ncities)
        bound = reducedMatrix['bound']
        rcm = reducedMatrix['rcm']

        # Push first city on the queue
        # Time complexity: O(n*logn)
        # https://python.plainenglish.io/python-for-interviewing-an-overview-of-the-core-data-structures-666abdf8b698
        heapq.heappush(priorityq, self.State([0], bound, rcm))

        # While queue is not empty and there is still time, try another state
        # I'm not sure how to analyze the time complexity of this loop. Without the time limit, the worst case
        # scenario would be O((n-1)!). With the time limit, the worst case scenario is that it runs for 60 seconds,
        # but I don't know how many states it can analyze in a second so I'm not sure how to narrow it down from there.
        # I guess a really nonspecific complexity would be O(60n) where n is the number of states it can analyze in a
        # second. I can't think of how to define it with anything more specific than that.
        # Within each loop, the time complexity is O(logn + n + n + m * n^2), which simplifies to O(m * n^2). This makes
        # sense intuitively because the rest of the algorithm is basically just comparisons, pruning, and new variable
        # assignments, while the section that creates the child states is a complex loop with reduced cost matrices.
        # Space complexity: O(n^2) because of the rcm deepcopy in the final child states loop
        while len(priorityq) != 0 and time.time() - startTime < time_allowance:
            # Pop the next state off the queue
            # Time complexity: O(logn)
            # https://stackoverflow.com/questions/48978388/why-heappop-time-complexity-is-ologn-not-on-in-python
            currentState = heapq.heappop(priorityq)

            # Prune states whose bounds were below the BSSF when queued but are now above the updated BSSF
            if bssf.cost <= currentState.bound:
                prunedStateCount += 1
                continue

            # Create a list of the cities which have not yet been visited (in the current state)
            # Time: O(n) because it loops n times
            # Space: O(n) because remainingStates is n items long before subtraction, and at most n-1 after subtraction
            remainingStates = []
            for i in range(ncities):
                remainingStates.append(i)
            remainingStates = list(set(remainingStates) - set(currentState.sequence))

            # If all cities have been reached, potentially update BSSF. Does not need to compare currentState.bound and
            # BSSF.cost because that was already done above on line 169
            # Time: O(n) because the for loop will iterate n times
            # Space: O(n) because the path is n items long
            if len(remainingStates) == 0 and currentState.rcm[currentState.sequence[-1], currentState.sequence[0]] \
                    is not math.inf:
                solutionCount += 1
                path = []
                for i in currentState.sequence:
                    path.append(cities[i])
                bssf = TSPSolution(path)

            # Check all children of currentState and either prune or queue
            # Time: combining the time complexities listed below, each loop iteration is O(n + n^2 + n*logn), which
            # simplifies to O(n^2) where n = num of cities. The loop iterations m times where m is the number of items
            # in remainingStates, so the time complexity is O(n^2 * m)
            # Space: O(n^2) because of the RCM deepcopy
            for state in remainingStates:
                totalStateCount += 1

                newBound = currentState.bound
                newRCM = copy.deepcopy(currentState.rcm)
                cityFrom = currentState.sequence[-1]
                cityTo = state

                # Add newest edge to bound and set relevant edge values to infinity
                # Time complexity: O(n) due to the for loop
                # Space complexity: O(n^2) because the RCM is O(n^2) and deepcopy creates a new instance
                newBound += newRCM[cityFrom, cityTo]
                newRCM[cityTo, cityFrom] = math.inf
                for i in range(ncities):
                    newRCM[cityFrom, i] = math.inf
                    newRCM[i, cityTo] = math.inf

                # Time complexity: O(n^2) (detailed in function definition)
                reducedMatrix = reduceMatrix(newRCM, newBound, ncities)
                newBound = reducedMatrix['bound']
                newRCM = reducedMatrix['rcm']

                newSequence = currentState.sequence[:]
                newSequence.append(state)

                # Queue child state if bond is below the BSSF cost, prune otherwise
                # Time: heappush() is O(n*logn); source referenced above
                if newBound < bssf.cost:
                    childState = self.State(newSequence, newBound, newRCM)
                    # if (tempFlag):
                    #     childState.sortFlag = True
                    #     tempFlag = False
                    heapq.heappush(priorityq, childState)
                else:
                    prunedStateCount += 1
            maxSize = max(maxSize, len(priorityq))
        prunedStateCount += len(priorityq)
        end_time = time.time()

        results = {'cost': bssf.cost, 'time': end_time - startTime, 'count': solutionCount, 'soln': bssf,
                   'max': maxSize, 'total': totalStateCount, 'pruned': prunedStateCount}

        return results

    class State:
        def __init__(self, sequence, bound, rcm):
            self.sequence = sequence
            self.bound = bound
            self.rcm = rcm
            self.sortFlag = False

        def __lt__(self, compare):
            if self.sortFlag:
                return True
            if self.bound != compare.bound:
                return self.bound < compare.bound
            return len(self.sequence) > len(compare.sequence)

    ''' <summary>
        This is the entry point for the algorithm you'll write for your group project.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution,
        best solution found.  You may use the other three field however you like.
        algorithm</returns> 
    '''

    def fancy(self, time_allowance=60.0):
        pass
