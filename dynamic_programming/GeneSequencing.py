#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3
BANDWIDTH = (MAXINDELS * 2) + 1
BANDED_NUM_OF_COLS = BANDWIDTH + 1

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


class GeneSequencing:

    def __init__(self):
        pass

    # This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean
    # that tells you whether you should compute a banded alignment or full alignment, and _align_length_ tells you how
    # many base pairs to use in computing the alignment
    def align(self, seq1, seq2, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length
        self.seq1 = seq1
        self.seq2 = seq2

        self.numOfRows = min(len(self.seq1), align_length) + 1
        self.numOfCols = min(len(self.seq2), align_length) + 1

        if banded:
            r = self.bandedAlignment()
        else:
            r = self.unbandedAlignment()

        score = r['score']
        alignment1 = r['align1']
        alignment2 = r['align2']

        ###################################################################################################
        # your code should replace these three statements and populate the three variables: score, alignment1
        # and alignment2
        # 		score = random.random()*100;
        # 		alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
        # 			len(seq1), align_length, ',BANDED' if banded else '')
        # 		alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
        # 			len(seq2), align_length, ',BANDED' if banded else '')
        ###################################################################################################

        return {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}

    # Time complexity: Adding up the separate time complexities commented below results in O(mn + mn + (m + n) +
    # (m + n)). This simplifies to 2mn + 2m + 2n, then mn + m + n after dropping the constant coefficients. m and n are
    # negligible compared to mn in Big-O notation, so it simplifies to O(mn).
    # Space complexity: The largest users of space are the two arrays initialized in the first two lines, editTable
    # and typeTracker. These are both of size mxn, where m = numOfRows and n = numOfCols. The only other remotely
    # comparable use of space is with the whichIsUsed array, which is mx1 so not really comparable at all.
    def unbandedAlignment(self):
        # Initialize 2D lists
        # Time complexity: O(mn) where m = numOfRows and n = numOfCols
        editTable = [[0 for _ in range(self.numOfCols)] for _ in range(self.numOfRows)]
        typeTracker = [[0 for _ in range(self.numOfCols)] for _ in range(self.numOfRows)]
        # Initialize first row and column to i * INDEL (or 1/2 for editType)
        for i in range(1, self.numOfCols):
            editTable[0][i] = i * INDEL
            typeTracker[0][i] = 1
        for i in range(1, self.numOfRows):
            editTable[i][0] = i * INDEL
            typeTracker[i][0] = 2

        # Calculate edit distance and type for all possible paths
        # Time complexity: O(mn) where m = numOfRows and n = numOfCols
        for i in range(1, self.numOfRows):
            for j in range(1, self.numOfCols):
                r = self.minEditDistance(editTable, i, j, False)
                editTable[i][j] = r['distance']
                typeTracker[i][j] = r['type']

        # Backtrace to determine what edits were made between each step
        # Time complexity: difficult to know exactly what the upper bound is. After the loop finishes, it's clear that
        # the loop runs O(n) where n = len(whichWasUsed). But before, it's a bit harder. Theoretically, the longest an
        # alignment string would be is m + n when m = len(seq1) and n = len(seq2), if the alignment was entirely inserts
        # and deletes, but the optimal alignment will never be that. I'm not quite sure how to calculate the actual
        # upper bound, but I know that it's somewhere between max(m, n) and m + n, so for simplicity's sake I'll say
        # this is O(m + n).
        whichWasUsed = []
        i = self.numOfRows - 1
        j = self.numOfCols - 1
        while j != 0 and i != 0:
            editType = typeTracker[i][j]
            whichWasUsed.insert(0, editType)
            if editType == 1:
                j -= 1
            if editType == 2:
                i -= 1
            if editType == 3:
                i -= 1
                j -= 1

        # Create alignment strings based on the nx1 list of edits
        # Time complexity: O(n) where n = len(whichWasUsed) (calculated at the function definition). This could also be
        # written as O(m + n) where m = len(seq1) and n = len(seq2), because as determined in the previous comment, the
        # alignments will definitely be shorter than m + n.
        r = self.buildAlignStrings(whichWasUsed)
        return {'score': editTable[self.numOfRows - 1][self.numOfCols - 1], 'align1': r['align1'],
                'align2': r['align2']}

    # Time complexity: Adding up the separate time complexities commented below results in O(n + n + n + n). This
    # equals 4n, which simplifies to O(n) in Big-O notation.
    # Space complexity: The largest users of space are the two arrays initialized in the first two lines, editTable
    # and typeTracker. These are both of size nx7, where n = numOfRows. The only other remotely comparable use of space
    # is with the whichIsUsed array, which is nx1. Both of these (7n and 1n) simplify to O(n) in Big-O notation.
    def bandedAlignment(self):
        # Automatically reject strings with a massive size discrepancy
        if abs(len(self.seq1) - len(self.seq2)) > (.8 * self.MaxCharactersToAlign):
            return {'score': math.inf, 'align1': "No Alignment Possible", 'align2': "No Alignment Possible"}

        # Initialize 2D lists
        # Time complexity: O(7n) where n = numOfRows – drops the coefficient and simplifies to O(n)
        editTable = [[0 for _ in range(BANDWIDTH)] for _ in range(self.numOfRows)]
        # Initialize first row and column to i * INDEL (or 1/2 for editType) with some special index arithmetic to
        # make sure that the correct indexes are initialized at the start (because the first j in the row will decrease
        # as i increases)
        typeTracker = [[0 for _ in range(BANDWIDTH)] for _ in range(self.numOfRows)]
        for i in range(1, MAXINDELS + 1):
            editTable[0][i + MAXINDELS] = i * INDEL
            typeTracker[0][i + MAXINDELS] = 1
        for i in range(1, MAXINDELS + 1):
            editTable[i][MAXINDELS - i] = i * INDEL
            typeTracker[i][MAXINDELS - i] = 2

        # Calculate edit distance and type for all possible paths. Special arithmetic at the beginning to ensure
        # that the top left corner isn't incorrectly included.
        # Time complexity: O(7n) where n = numOfRows – drops the coefficient and simplifies to O(n)
        for i in range(1, self.numOfRows):
            for j in range(0, BANDWIDTH):
                if i > MAXINDELS:
                    r = self.minEditDistance(editTable, i, j, True)
                    editTable[i][j] = r['distance']
                    typeTracker[i][j] = r['type']
                else:
                    if j > MAXINDELS - i:
                        r = self.minEditDistance(editTable, i, j, True)
                        editTable[i][j] = r['distance']
                        typeTracker[i][j] = r['type']

        # Depending on the size difference between seq1 and seq2, the final cost may not be in the center of the row.
        # For example, polynomial -> exponential finds the cost at j = 4, not j = 3. This calculates what the index of
        # the final cost will be.
        # Time complexity: Should all be O(1), no complex function calls or arithmetic used
        if self.numOfRows == self.MaxCharactersToAlign + 1:
            endingIndex = MAXINDELS
        elif len(self.seq1) > len(self.seq2):
            endingIndex = MAXINDELS - min(abs(len(self.seq1) - len(self.seq2)), MAXINDELS)
        else:
            endingIndex = MAXINDELS + min(abs(len(self.seq1) - len(self.seq2)), MAXINDELS)

        # Backtrace to determine what edits were made between each step
        # Time complexity: easier to calculate than in the unbanded, because it's a given that the final alignment will
        # be nowhere near as long as mxn. The max length of the alignment strings will be n + 7, which simplifies to
        # O(n) in Big-O notation.
        whichWasUsed = []
        i = self.numOfRows - 1
        j = endingIndex
        while i != 0:
            editType = typeTracker[i][j]
            whichWasUsed.insert(0, editType)
            if editType == 1:
                j -= 1
            if editType == 2:
                i -= 1
                j += 1
            if editType == 3:
                i -= 1

        # Create alignment strings based on the nx1 list of edits
        # Time complexity: O(n) where n = len(whichWasUsed) (calculated at the function definition below). This could
        # also be written as O(n + 7) where n = len(seq1), because as determined in the previous comment, the
        # alignments will definitely be no longer then n + 7. This is of course simplified to O(n).
        r = self.buildAlignStrings(whichWasUsed)

        return {'score': editTable[self.numOfRows - 1][endingIndex], 'align1': r['align1'], 'align2': r['align2']}

    # Time complexity: O(1), because aside from the function calls, it only uses simply arithmetic, boolean comparisons,
    # and accessing a specific index which is also O(1) regardless of list size. Taking the function calls into account,
    # time complexity is still only O(1), because self.calcValsUnbanded() and self.calcValsBanded() are both O(1) too.
    # Space complexity: O(1) because the local variables are small ints – should always be less than -10000, which is
    # 14 bits + a sign bit
    def minEditDistance(self, E, i, j, banded):
        if not banded:
            r = self.calcValsUnbanded(E, i, j)
        else:
            r = self.calcValsBanded(E, i, j)

        subValue = r['sub']
        insertValue = r['insert']
        deleteValue = r['delete']

        minValue = insertValue
        editType = 1
        if deleteValue < minValue:
            minValue = deleteValue
            editType = 2
        if subValue < minValue:
            minValue = subValue
            editType = 3
        return {'distance': minValue, 'type': editType}

    # Time complexity: O(1), only uses simple arithmetic and accessing a specific index which is also O(1) regardless
    # of list size. The call to self.diff() is also only O(1).
    # Space complexity: O(1) because the local variables are small ints – should always be less than 5000, which is
    # 13 bits
    def calcValsUnbanded(self, E, i, j):
        subValue = E[i - 1][j - 1] + self.diff(i - 1, j - 1)
        insertValue = E[i][j - 1] + INDEL
        deleteValue = E[i - 1][j] + INDEL
        return {'sub': subValue, 'insert': insertValue, 'delete': deleteValue}

    # Time complexity: O(1), only uses simple arithmetic, boolean comparisons, and accessing a specific index which
    # is also O(1) regardless of list size. The call to self.diff() is also only O(1).
    # Space complexity: O(1) because the local variables are small ints – should always be less than -10000, which is
    # 14 bits + a sign bit
    def calcValsBanded(self, E, i, j):
        if j == 0:
            insertValue = math.inf
        else:
            insertValue = E[i][j - 1] + INDEL
        if j == BANDWIDTH - 1:
            deleteValue = math.inf
        else:
            deleteValue = E[i - 1][j + 1] + INDEL
        subValue = E[i - 1][j] + self.diff(i - 1, j - (MAXINDELS - i) - 1)
        return {'sub': subValue, 'insert': insertValue, 'delete': deleteValue}

    # Time complexity: O(1), only uses simple arithmetic, boolean comparisons, and accessing a specific index which
    # is also O(1) regardless of list size
    # Space complexity: O(1) because it uses no new space (would this be O(0)? I assume not.)
    def diff(self, i, j):
        if i >= len(self.seq1) or j >= len(self.seq2):
            return SUB
        if self.seq1[i] == self.seq2[j]:
            return MATCH
        else:
            return SUB

    # Time complexity: O(n) where n = len(typeList). The major contributor of time complexity is the for loop iterating
    # through the list of edits and building the strings, which is O(n) where n = len(typeList). Other than that, the
    # function only uses simple arithmetic, boolean comparisons, and accessing a specific index in a list.
    # Space complexity: O(n) where n = len(typeList). I wasn't quite sure how to calculate this, because I'm not sure
    # whether the space used by typeList counts because it was already allotted before the function call. If the space
    # of typeList is included, then the space complexity is O(n) where n = len(typeList), because typeList is an nx1
    # array of very small ints (2 bits at most). Excluding the space of typeList, the space complexity is actually
    # probably still O(n), where n = len(align1) or len(align2) (the len should be the same at the end), which should
    # also be equal to n = len(typeList). If there are 36 elements in typeList, there should also be 36 elements in
    # each string. Technically, this would be 2n, but constant coefficients are negligible in Big-O notation.
    def buildAlignStrings(self, typeList):
        align1, align2 = "", ""
        seq1index, seq2index = 0, 0
        for i in range(len(typeList)):
            if typeList[i] == 1:
                if seq2index < len(self.seq2):
                    align2 += self.seq2[seq2index]
                    seq2index += 1
                align1 += "-"
            if typeList[i] == 2:
                if seq1index < len(self.seq1):
                    align1 += self.seq1[seq1index]
                    seq1index += 1
                align2 += "-"
            if typeList[i] == 3:
                if seq1index < len(self.seq1):
                    align1 += self.seq1[seq1index]
                    seq1index += 1
                if seq2index < len(self.seq2):
                    align2 += self.seq2[seq2index]
                    seq2index += 1
        return {'align1': align1, 'align2': align2}
