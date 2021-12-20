# Project 6 -- Tabu Search Optimization

For this project, we will be implementing the tabu search optimization algorithm
to find an optimal solution to the Traveling Salesperson Problem (TSP). From the
wikipedia page:

> The word tabu comes from the Tongan word to indicate things that cannot be
> touched because they are sacred.

The idea behind the tabu optimization search is that we perform a local search
to find a local minima, then we can remember those possible local minima and
avoid them as we know that we cannot get any better going down that path.

This allows us to back track a small amount to possibly find a better algorithm
that the greediness of local search would prevent us from searching.

## Memory

There are 3 categories of memory that we use to store previous solutions. 

1. Short-term
2. Intermediate-term
3. Long-term

We only need _short-term_ memory in order to improve on the local search
algorithm, but intermediate and long term memory is necessary for solving harder
problems.

If it is possible, we could implement both intermediate, and long-term memory,
but this isn't required if we don't have enough time.

### Short-term Memory

A list of the most recent solutions we have found. We are not allowed to revisit
any nodes in this list--They are considered to be _tabu banned_.

The solutions in short-term memory have a short life-span.

### Intermediate-term Memory

Memory that we can use to encourage or discourage certain paths.

> a memory structure that prevents or induces certain moves (e.g. based on
> frequency memory applied to solutions sharing features in common with
> unattractive or attractive solutions found in the past).

### Long-term Memory

Rules that help move the search into new areas that we haven't visited before.
Usually used when the algorithm gets stuck in a dead end.

## Algorithm

There are a few things that will need to be implemented.

1. findNeighbors
2. short-term memory
3. intermediate-term memory
4. long-term memory
5. Main algorithm to tie it all together.


Below is an algorithm that excludes intermediate/long term memory provided by
the wikipedia article. It may be useful as a starting point to implement our
algorithm.

```
sBest ← s0
bestCandidate ← s0
tabuList ← []
tabuList.push(s0)
while (not stoppingCondition())
    sNeighborhood ← getNeighbors(bestCandidate)
    bestCandidate ← sNeighborhood[0]
    for (sCandidate in sNeighborhood)
        if ( (not tabuList.contains(sCandidate)) and (fitness(sCandidate) > fitness(bestCandidate)) )
            bestCandidate ← sCandidate
        end
    end
    if (fitness(bestCandidate) > fitness(sBest))
        sBest ← bestCandidate
    end
    tabuList.push(bestCandidate)
    if (tabuList.size > maxTabuSize)
        tabuList.removeFirst()
    end
end
return sBest
```


All of the information about this algorithm has come from the wikipedia page
[tabu search][tabu].

[tabu]: https://en.wikipedia.org/wiki/Tabu_search

## Useful links

[Paper 1](http://www.ijarse.com/images/fullpdf/1519813763_NMCOE4098IJARSE.pdf)

[Paper 2](https://www.researchgate.net/publication/331585233_Tabu_Search_Method_for_Solving_the_Traveling_salesman_Problem)
