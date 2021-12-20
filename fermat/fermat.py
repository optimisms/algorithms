import random
import math

def prime_test(N, k):
    # This is main function, that is connected to the Test button. You don't need to touch it.
    return fermat(N, k), miller_rabin(N, k)

# Functionality: simple modular exponentiation recursive algorithm. Recurse through, halving y until it reaches 0,
# then return one by one, increasing the exponentiation with each return. On return, at each level, if the exponent y
# is even, square the base and modulo by N. If the exponent y is odd, square the base like normal, but also multiply
# by x, then modulo by N. The result is x^y % N.

# Time complexity is determined as follows: Mult/divide/modulo are each n^2, where n = # of digits in x and N,
# which adds up to 5n^2; depth of the tree is m = # of digits in y. Therefore, the total complexity is approximately
# 5n^2*m. However, it can be assumed that y is of a similar length to x and M, therefore m = n and it becomes 5n^3.
# Finally, as n gets larger, all but the largest order become negligible. Therefore, this simplifies to O(n^3).

# This has a space complexity of n. The last return statement would take the most space. If x and z are both n-digits
# long, multiplying them together will create a number of 3n-digits length. However, as n grows, the 3 becomes
# negligible, so it's just n.
def mod_exp(x, y, N):
    if y == 0: return 1
    z = mod_exp(x, math.floor(y / 2), N)
    if y % 2 == 0:
        return (z * z) % N
    else:
        return (x * z * z) % N


# Functionality: very simple. The larger k gets, the less likely the test is to be wrong. This probability decreases
# exponentially by a power of two. Therefore, the probability that the test is correct increases, and to find that,
# we simply subtract the probability of failure from 1.
# If I remember right, exponentiation has a time complexity of n^3, where n = # of digits in k. There's nothing else
# here to dominate that, so I believe the complexity is O(n^3). I believe the space complexity should just be n or
# n+1, which simplifies to just n. This is a very small algorithm.
def fprobability(k):
    return 1 - (1 / (2 ** k))


def mprobability(k):
    return 1 - (1 / (4 ** k))
    # You will need to implement this function and change the return value.   
    # return 0.0


# Functionality: Fermat's little theorem dictates that "If p is prime, then for every 1 ≤ a < p, a^(p-1) % p = 1."
# This theorem is implemented below. The more a's are tested, the more certain one can be of the accuracy of the
# result.
# First, a random number between 2 and N-1 is generated, and then checked against previous random numbers to ensure
# there are k unique tests. Then, using mod_exp(), modular exponentiation is completed. If the result is not 1,
# the number definitely isn't prime, and the function returns 'composite'. If the program reaches the end of the
# k-length loop without returning, Fermat's theorem has determined that N is most likely prime and returns 'prime'.

# Excluding mod_exp, I believe the time complexity is k, because it will iterate k times. There are no operations
# larger than O(n) (where n = # of digits in N) that I can see. However, I don't see why the complexity of mod_exp
# wouldn't factor in, so the time complexity is really n^3.

# The space complexity is not very large because the largest space complexity, which would come from the modular
# exponentiation, is kept within the modular range. Therefore, I think the space complexity would be n where n = # of
# digits in N.
def fermat(N, k):
    testVals = []
    for i in range(k):
        rand = random.randint(2, N - 1)
        while rand in testVals:
            rand = random.randint(2, N - 1)
        testVals.append(rand)
        result = mod_exp(rand, N-1, N)
        if result != 1:
            return 'composite'
    return 'prime'


# def miller_rabin(N, k):
    # You will need to implement this function and change the return value, which should be
    # either 'prime' or 'composite'.
    #
    # To generate random values for a, you will most likley want to use
    # random.randint(low,hi) which gives a random integer between low and
    #  hi, inclusive.
    # return 'composite'

# def miller_rabin(N, k):
#     testVals = []
#     d = N - 1
#     s = 0
#     while d % 2 == 0:
#         d /= 2
#         s += 1
#     for i in range(k):
#         rand = random.randint(2, N - 1)
#         while rand in testVals:
#             rand = random.randint(2, N - 1)
#         testVals.append(rand)
#
#         result = mod_exp(rand, d, N)
#
#         if result != 1:
#             for j in range(s):
#                 result = mod_exp(result, 2, N)
#                 if result == (N - 1):
#                     break
#                 return 'composite'
#     return 'prime'