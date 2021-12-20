from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF, QObject
# elif PYQT_VER == 'PYQT4':
#    from PyQt4.QtCore import QLineF, QPointF, QObject
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time

# Some global color constants that might be useful
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)

# Global variable that controls the speed of the recursion automation, in seconds
#
PAUSE = 0.25


#
# This is the class you have to complete.
#
class ConvexHullSolver(QObject):

    # Class constructor
    def __init__(self):
        super().__init__()
        self.pause = False

    # Some helper methods that make calls to the GUI, allowing us to send updates
    # to be displayed.

    def showTangent(self, line, color):
        self.view.addLines(line, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseTangent(self, line):
        self.view.clearLines(line)

    def blinkTangent(self, line, color):
        self.showTangent(line, color)
        self.eraseTangent(line)

    def showHull(self, polygon, color):
        self.view.addLines(polygon, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseHull(self, polygon):
        self.view.clearLines(polygon)

    def showText(self, text):
        self.view.displayStatusText(text)

    # This is the method that gets called by the GUI and actually executes
    # the finding of the hull
    def compute_hull(self, points, pause, view):
        self.pause = pause
        self.view = view
        assert (type(points) == list and type(points[0]) == QPointF)

        t1 = time.time()
        points.sort(key=lambda point: point.x())
        t2 = time.time()

        t3 = time.time()
        convexHull = self.dac_convex_hull_solver(points)
        t4 = time.time()

        # when passing lines to the display, pass a list of QLineF objects.  Each QLineF
        # object can be created with two QPointF objects corresponding to the endpoints
        polygon = [QLineF(convexHull[i], convexHull[(i + 1) % len(convexHull)]) for i in range(len(convexHull))]
        self.showHull(polygon, RED)
        self.showText('Time Elapsed (Convex Hull): {:3.3f} sec'.format(t4 - t3))

    # divide and conquer algorithm splits the points in half until reaching the base case, then recurses back up
    # Time complexity: O(log n) because it will run log(n) times where n is the size of the list of points
    # Space complexity: O(n) because it has two arrays that are a total length of n
    def dac_convex_hull_solver(self, points):
        if (len(points) == 1) or (len(points) == 2): return points
        middleIndex = len(points) // 2
        left = points[:middleIndex]
        right = points[middleIndex:]
        return self.merge(self.dac_convex_hull_solver(left), self.dac_convex_hull_solver(right))

    # merge algorithm to combine the two sets of points, keeping only those that form the combined set's convex hull
    # Time complexity: O(n) in the worst case, where n is the total number of points, because that is when every single
    # point in both point sets is part of the convex hull, which means the entire set gets traversed and appended once
    # Space complexity: O(n) as well, because mergedHull will never be larger than the total size of both distinct sets
    # even in the worst case
    def merge(self, leftPoints, rightPoints):
        upperLeft, upperRight = self.find_tangent(leftPoints, rightPoints, True)
        lowerLeft, lowerRight = self.find_tangent(leftPoints, rightPoints, False)

        mergedHull = []

        # keep left side of leftPoints (ie delete right side))
        i = lowerLeft
        while i != upperLeft:
            mergedHull.append(leftPoints[i])
            i = (i + 1) % len(leftPoints)
        mergedHull.append(leftPoints[upperLeft])

        # keep right side of rightPoints (ie delete left side)
        i = upperRight
        while i != lowerRight:
            mergedHull.append(rightPoints[i])
            i = (i + 1) % len(rightPoints)
        mergedHull.append(rightPoints[lowerRight])

        return mergedHull

    # find the tangent; circle through each side clockwise and counterclockwise until the line with the best slope is
    # found, then return that line
    # Time complexity: Worst case time complexity would be n^2, where n is the number of points in each set (they may be
    # slightly different sizes by one point but this difference is negligible). The inner while loops will iterate a
    # maximum of n times, though will rarely (if ever) actually iterate through all points in the set. The outer while
    # loop could run a maximum of n times, if each time only saw one adjustment from each set. Therefore, the worst case
    # time complexity is n^2.
    # Space complexity: O(1) because the only memory used are ints and booleans
    def find_tangent(self, leftPoints, rightPoints, upper):
        leftIndex = self.get_innermost_index(leftPoints, True)
        rightIndex = self.get_innermost_index(rightPoints, False)

        tangentFound = False
        while not tangentFound:
            tangentFound = True

            # traverse the left hull
            currSlope = self.calc_slope(leftPoints[leftIndex], rightPoints[rightIndex])
            nextSlope = None
            while currSlope != nextSlope:
                nextLeft = ((leftIndex - 1) if upper else (leftIndex + 1)) % len(leftPoints)
                nextSlope = self.calc_slope(leftPoints[nextLeft], rightPoints[rightIndex])
                if (nextSlope < currSlope) if upper else (nextSlope > currSlope):
                    tangentFound = False

                    leftIndex = nextLeft
                    currSlope = nextSlope
                    nextSlope = None
                else:
                    nextSlope = currSlope

            # traverse the right hull
            currSlope = self.calc_slope(leftPoints[leftIndex], rightPoints[rightIndex])
            nextSlope = None
            while currSlope != nextSlope:
                nextRight = ((rightIndex + 1) if upper else (rightIndex - 1)) % len(rightPoints)
                nextSlope = self.calc_slope(leftPoints[leftIndex], rightPoints[nextRight])
                if (nextSlope > currSlope) if upper else (nextSlope < currSlope):
                    tangentFound = False

                    rightIndex = nextRight
                    currSlope = nextSlope
                    nextSlope = None
                else:
                    nextSlope = currSlope
        return leftIndex, rightIndex

    # Helper function to find the rightmost or leftmost index
    # Time complexity: O(n) because it traverses the array once
    # Space complexity: O(1), no memory used except an int
    def get_innermost_index(self, points, left):
        index = -1
        for i in range(len(points)):
            if index < 0:
                index = i
                continue
            elif (points[i].x() > points[index].x()) if left else (points[i].x() < points[index].x()):
                index = i
        return index

    # Helper function to calculate the slope
    # Time complexity: O(1)
    # Space complexity: Constant, no memory used
    def calc_slope(self, pointA, pointB):
        return (pointB.y() - pointA.y()) / (pointB.x() - pointA.x())
