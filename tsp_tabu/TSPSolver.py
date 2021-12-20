#!/usr/bin/python3

from typing import List, Tuple, Optional
from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
# elif PYQT_VER == 'PYQT4'
# 	from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time
import numpy as np
from TSPClasses import *
import heapq
import itertools

import TSPClasses as classes


class Branch:
    def __init__(self, matrix: np.ndarray, cost: float, route: list):
        self.matrix = matrix
        self.cost = cost
        self.route = route

    def __lt__(self, other):
        # Prefer depth over breadth
        if len(self.route) != len(other.route):
            return len(self.route) > len(other.route)
        return self.cost < other.cost
        # return self.cost / len(self.route) < other.cost / len(other.route)


class TSPSolver:
    def __init__(self, gui_view):
        self._scenario: classes.Scenario = None

    def setupWithScenario(self, scenario):
        self._scenario = scenario

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

    def greedy(self, time_allowance=60.0, costs=None, lower_bound=None):
        """
        Perform a greedy search for a solution.

        Note, this algorithm will sometimes fail. To guarantee some solution, run defaultRandomTour if greedy fails.

        :param time_allowance: allowed time to run
        :param costs: optional costs matrix
        :param lower_bound: optional lower bound

        :time complexity: O(n^3)
        :space complexity: O(n^3)
        """
        t_start = time.time()
        results = {
            'total': 0
        }
        if costs is None or lower_bound is None:
            costs = self.makeCostMatrix()
            costs, lower_bound = self.reduceMatrix(costs)
        h, w = np.shape(costs)

        def find_zeros(arr: list):
            # O(n)
            for i, v in enumerate(arr):
                if v == 0:
                    yield i

        # Find a starting row (Doesn't have multiple 0s in the row)
        start = 0
        # O(n)
        for start in range(h):
            if len(list(find_zeros(costs[start]))) == 1:
                break

        choices = [{
            'route': [start],
            'costs': costs,
            'cost': lower_bound,
        }]

        # O(n^2)
        def run_path(route: list, costs: np.ndarray, cost: int):
            if len(route) == h:
                # print(f"Possible solution: {route}")
                return {
                    'route': route,
                    'cost': cost
                }
            r = route[-1]
            # (n)
            zeros = list(find_zeros(costs[r]))
            if not zeros:
                return False
            # O(n^2) On average there are between 1 and 2 zeros
            for z in zeros:
                tmp = self.remove_backtracking(costs, r, z)
                tmp[z][start] = np.inf
                tmp, cst = self.reduceMatrix(tmp)
                results['total'] += 1
                choices.append({
                    'route': route + [z],
                    'costs': tmp,
                    'cost': cost + cst
                })
            return None

        soln = None
        iters = 0
        max_queue = len(choices)
        # O(n^3) Not O(n^2 + b^n) because the loop breaks as soon as a valid solution is made
        while choices and time.time() - t_start < time_allowance:
            max_queue = max(max_queue, len(choices))
            choice = choices.pop()
            # O(n^2)
            result = run_path(**choice)
            if result is not None:
                iters += 1
                if not result:
                    continue
                cities = self._scenario.getCities()
                soln = TSPSolution([cities[i] for i in result['route']])
                # Stop when we find a valid solution
                if result['cost'] == soln.cost:
                    break
                else:
                    soln = None

        t_end = time.time()

        results['cost'] = soln.cost if soln else np.inf
        results['time'] = t_end - t_start
        results['count'] = iters
        results['soln'] = soln
        results['max'] = max_queue
        results['pruned'] = None
        return results

    def branchAndBound(self, time_allowance=60.0):
        """
        Perform a branch and bound algorithm to find the best TSP solution.

        This internally will run the greedy algorithm to find a starting point

        :param time_allowance: allowed time to run
        :return:

        :time complexity: Omega(n^4)
        :space complexity: Omega(n^3)
        """
        t_start = time.time()
        costs = self.makeCostMatrix()
        costs, lower_bound = self.reduceMatrix(costs)
        h, w = np.shape(costs)
        results = {
            'max': 0,
            'total': 1,
            'pruned': 0,
            'count': 0,
        }

        upper_bound = None

        greedy_results = self.greedy(time_allowance, costs, lower_bound)
        t_greedy = time.time()
        if greedy_results['soln']:
            soln = greedy_results['soln']
            upper_bound = greedy_results['cost']
        else:
            print("Couldn't find a solution using the greedy algorithm, falling back to random tour")
            default_results = self.defaultRandomTour(time_allowance - (t_greedy - t_start))
            t_default = time.time()
            if default_results['soln']:
                soln = default_results['soln']
                upper_bound = default_results['cost']
            else:
                results['cost'] = np.inf
                results['time'] = t_default - t_start
                results['count'] = 0
                results['soln'] = None
                results['max'] = None
                results['total'] = None
                results['pruned'] = None
                return results

        start = 0
        branches = [Branch(costs, lower_bound, [start])]

        print(f"Starting solution of cost {soln.cost} with lower bound {lower_bound}")

        def branch(matrix: np.ndarray, cost: float, route: list):
            if cost > upper_bound:
                results['pruned'] += 1
                return None
            if len(route) == h:
                cities = self._scenario.getCities()
                tmp = TSPSolution([cities[i] for i in route])
                if tmp.cost == cost:
                    return tmp
                return None
            r = route[-1]
            # O(n^3)
            for c, cst in enumerate(matrix[r]):
                if cst == np.inf:
                    continue
                tmp = self.remove_backtracking(matrix, r, c)
                tmp[c][start] = np.inf
                # O(n^2)
                tmp, reduce_cst = self.reduceMatrix(tmp)
                results['total'] += 1
                if cost + reduce_cst > upper_bound:
                    # Skip this branch as it is not worth our time
                    results['pruned'] += 1
                    continue
                # O(log n)
                heapq.heappush(branches, Branch(tmp, cost + reduce_cst, route + [c]))
            results['max'] = max(results['max'], len(branches))
            return None

        loop = 0
        new_line = True

        # Omega(n^4), O(n!) It will never reach n! because of the bounding function, but we have no way to assert that
        # it is better than O(n!)
        while branches and time.time() - t_start < time_allowance:
            loop += 1
            curr = heapq.heappop(branches)
            # O(n^3)
            new_sol = branch(curr.matrix, curr.cost, curr.route)
            if new_sol and new_sol.cost < soln.cost:
                results['count'] += 1
                soln = new_sol
                upper_bound = soln.cost
                print(f"New solution of cost {soln.cost} found!")
                new_line = True
            # Add an info message for every 1000 loops
            if loop == 1000:
                loop = 0
                begin = '\033[F' if not new_line else ''
                print(
                    f"{begin}Discovered {results['total']} nodes and pruned {results['pruned']} nodes so far ({round(results['pruned'] / results['total'] * 100)}% pruned)")
                new_line = False

        t_stop = time.time()

        print(f"Final solution of cost {soln.cost}")

        results['cost'] = soln.cost if soln else np.inf
        results['time'] = t_stop - t_start
        results['soln'] = soln
        return results

    def fancy(self, time_allowance=60.0):
        """
        Run the tabu search algorithm

        :time complexity: O(n^3)
        """
        t_start = time.time()
        stop_time = t_start + time_allowance

        costs = self.makeCostMatrix()
        costs, lower_bound = self.reduceMatrix(costs)
        results = {
            'max': 0,
            'total': 1,
            'pruned': 0,
            'count': 0,
        }

        # O(n^3)
        greedy_results = self.greedy(time_allowance, costs, lower_bound)
        t_greedy = time.time()
        if greedy_results['soln']:
            soln = greedy_results['soln']
        else:
            print("Couldn't find a solution using the greedy algorithm, falling back to random tour")
            default_results = self.defaultRandomTour(time_allowance - (t_greedy - t_start))
            t_default = time.time()
            if default_results['soln']:
                soln = default_results['soln']
            else:
                results['cost'] = np.inf
                results['time'] = t_default - t_start
                results['count'] = 0
                results['soln'] = None
                results['max'] = None
                results['total'] = None
                results['pruned'] = None
                return results

        tabu_list = []
        bssf = soln
        cbssf = soln
        curr = soln
        solution_count = 0

        longterm = []

        # Because runningTimes is constant no matter what n, the main loop is O(1)
        runningTimes = 200
        last_i = 0
        first_iteration = True
        # O(200n^3) = O(n^3)
        for i in range(runningTimes):
            if time.time() - t_start > time_allowance:
                break
            # O(n^3)
            neighbors = self.find_neighbors(curr, tabu_list + longterm, stop_time)
            curr = neighbors.pop(0)

            # When we are on the first iteration, add our selected route to the long-term memory
            if first_iteration:
                longterm.append(curr.route)
                first_iteration = False
            # O(n)
            self.add_shortterm_item(tabu_list, curr.route, 50)

            # If the current solution is better than our bssf, then log it
            if curr.cost < cbssf.cost:
                solution_count += 1
                print(f"Found a better solution: {curr.cost} in {i - last_i} iterations")
                cbssf = curr
                last_i = i

            # If we are taking too long to improve our solution, start over
            if i - last_i > 15:
                # We have gotten stuck
                print("We have gotten stuck")
                if cbssf.cost < bssf.cost:
                    bssf = cbssf
                    print("We have found an overall better solution than our current best one")
                curr = soln
                cbssf = soln
                last_i = i
                first_iteration = True

        if cbssf.cost < bssf.cost:
            bssf = cbssf

        # Perhaps a use of intermediate would be:
        # when we get the list of neighbors, we only follow the best one as a potential solution
        # BUT we store the next best 5-10 (or whatever number we choose) neighbors in intermediate
        # memory so that, in the future, if we reach a point where there are no neighbors that aren't
        # taboo, we can backtrack and find the next best neighbor and follow that

        t_stop = time.time()

        print(f"Final solution is {bssf.cost}")

        results['time'] = t_stop - t_start
        results['soln'] = bssf
        results['cost'] = bssf.cost
        results['count'] = solution_count
        return results

    def find_neighbors(self, solution: classes.TSPSolution, tabu_list: List[TSPSolution], stop_time: float) -> List[TSPSolution]:
        """
        Create a list of neighbors that would be made of only one change.

        :param solution: The current solution to find neighbors for
        :param tabu_list: The list of 'taboo' solutions already searched
        :param stop_time: Time when the calculation should be stopped

        :return: a list of neighbors

        :time complexity: O(n^3)
        :space complexity: O(n^2)
        """

        route = solution.route
        neighborhood = []

        # Time complexity: after some trial and error, found the complexity function: n^2 - (n^2 + n)/2
        # This simplifies to O(n^2), unfortunately

        # I also calculated the time complexity for this (It ends up being the same n^2) but creating TSPSolution
        # calculates the cost on creation which is O(n), so the overall complexity should be n^3 - Ben
            # ah makes sense - Haile

        # \sum(\sum 1 j=i+1 to n) i=1 to n-1 = 0.5(n - 1)n
        # O((0.5(n - 1)n)n) = O(0.5(n^3 - n^2)) = O(n^3)
        # Space complexity O(n^2)
        for i in range(len(route) - 1):
            if time.time() > stop_time:
                break
            for j in range(i + 1, len(route)):
                if self.swappable(route, i, j):
                    newRoute = route[:]
                    newRoute[i] = route[j]
                    newRoute[j] = route[i]
                    # O(n)
                    neighbor = TSPSolution(newRoute)
                    neighborhood.append(neighbor)

        # O(n^2)
        neighborhood = self.prune_by_tabu(tabu_list, neighborhood)

        # Sort the neighborhood in place (nlog n)
        neighborhood.sort(key=lambda x: x.cost)

        # Instead of pruning by cost, keeping only the top n might be better for the intermediate memory
        neighborhood = neighborhood[:5]

        return neighborhood

    def swappable(self, route: List[classes.City], i: int, j: int) -> bool:
        """
        Check if route[i] and route[j] can be swapped

        :param route: route
        :param i: node index
        :param j: node index
        :return: whether the two nodes can be swapped

        :time complexity: O(1)
        :space complexity: O(1)
        """
        firstCity = (i == 0)
        lastCity = (j == len(route) - 1)
        adjacent = (j - 1 == 1) or (firstCity and lastCity)
        city1 = route[i]
        city2 = route[j]

        totalCost = 0

        # If first city, i-1 == -1 which python automatically handles
        # If last city, j+1 % (len(route) -1) == 0
        totalCost += route[i-1].costTo(city2)                       # to city2
        totalCost += city1.costTo(route[(j+1) % (len(route) - 1)])  # from city1

        # This is just a weird edge case because the first and last city are adjacent. The index for the last city
        #  comes after the index for the first city (in the array), but technically in visitation order, because they
        #  are adjacent and share an edge, the first city should come after the last city for the calculation in the
        #  adjacent branch.
        # ie for the path 1 > 3 > 2 > 4, to switch 1 and 4 we need to check edges 2 > 1, 4 > 3, and 1 > 4. The above
        #  lines of code will check 2 > 1 and 4 > 3, but because 4 comes after 1 in the array, 1 is city1 and 4 is
        #  city2. So the code for the adjacent branch below would calculate edge 4 > 1 which already exists in the
        #  current solution (looping back around). It needs to calculate 1 > 4 instead.
        # So, we switch the two just for the purposes of the adjacent calculation:
        if firstCity and lastCity:
            city1 = route[j]
            city2 = route[i]

        # If the two cities are not right next to each other, we have to check four edges: to and from both
        # nodes. But if they're right next to each other, we only have to check three edges, because the edge
        # from the first and to the second is the same edge.
        if adjacent:
            totalCost += city2.costTo(city1)        # between cities
        else:
            totalCost += city2.costTo(route[i + 1]) # from city2
            totalCost += route[j - 1].costTo(city1) # to city1

        return totalCost != math.inf                # returns true if all new edges are valid, false otherwise

    def prune_by_tabu(self, tabu_list: List[classes.TSPSolution], neighborhood: List[classes.TSPSolution]) -> List[classes.TSPSolution]:
        """
        Prune any paths that are contained in the tabu list

        :param tabu_list: tabu list
        :param neighborhood: neighborhood
        :return: the pruned route list

        :time complexity: O(n^2)
        :space complexity: O(n)
        """
        # space O(n)
        newNeighborhood = []
        # O(n^2)
        for neighbor in neighborhood:
            # O(n)
            if neighbor.route not in tabu_list:
                newNeighborhood.append(neighbor)
        return newNeighborhood

    def add_shortterm_item(self, short_memory: list, item, max_size: int):
        """
        Add an item to the short term memory in place

        :param short_memory: short term memory object
        :param item: item to add
        :param max_size: maximum size of the memory

        :time complexity: O(n)
        :space complexity: O(1)
        """

        short_memory.append(item)
        # O(n)
        while len(short_memory) > max_size:
            short_memory.pop(0)

    def makeCostMatrix(self):
        """
        Create a matrix of the cost for each city to get to the other

        It costs infinity for a city to go to itself.

        :return: np matrix of the cost

        :time complexity: O(n^2)
        :space complexity: O(n^2)
        """
        cities = self._scenario.getCities()
        rows = []
        for i, cityA in enumerate(cities):
            row = []
            rows.append(row)
            for j, cityB in enumerate(cities):
                if i == j:
                    row.append(float("inf"))
                    continue
                row.append(cityA.costTo(cityB))
        return np.array(rows)

    def reduceMatrix(self, matrix: np.ndarray) -> Tuple[np.ndarray, float]:
        """
        Reduce a cost matrix as much as possible and return the cost to reduce
        that matrix

        :param matrix: matrix

        :return: (reduced cost matrix, cost to reduce the matrix)

        :time complexity: O(n^2)
        :space complexity: O(n^2)
        """
        matrix = np.copy(matrix)
        cost = 0
        # Reduce rows - O(n^2)
        for i, row in enumerate(matrix):
            # O(n)
            lowest = min(row)
            if lowest == float('inf'):
                continue
            # O(n)
            matrix[i] = row - lowest
            # O(1)
            cost += lowest
        # Reduce the cols - O(n^2)
        matrix = np.transpose(matrix)
        for j, col in enumerate(matrix):
            lowest = min(col)
            if lowest == float('inf'):
                continue
            matrix[j] = col - lowest
            cost += lowest
        return np.transpose(matrix), cost

    def remove_backtracking(self, matrix: np.ndarray, row: int, col: int) -> np.ndarray:
        """
        Remove the backtracking of a matrix. This will invalidate any paths that should be inaccessible in the future.

        :param matrix: matrix to remove backtracking from
        :param row: row to remove
        :param col: column to remove
        :return: the matrix with backtracking removed

        :time complexity: O(n^2)
        :space complexity: O(n^2)
        """
        val = matrix[row][col]
        mat = np.copy(matrix)
        h, w = np.shape(matrix)
        mat[row] = [np.inf] * w
        mat = np.transpose(mat)
        mat[col] = [np.inf] * h
        mat = np.transpose(mat)
        mat[row][col] = val

        return mat
