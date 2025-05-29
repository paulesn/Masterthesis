
**Idee**: for the construction, whenever a edge is added between two well seperated pairs, it does not add the direct edge, it adds the shortest path in the visibility graph.

## Problem
 - The visibility graph is a subgraph of the full graph. the WSPD construction assumes that every edge can be used, which here is not true. it would be possible that the WSPD construction uses (and relies on) one of the edges that are filtered in the visibility graph.
 - The visibility graph could have removed so many edges that two nodes that have a short euklidian distance are far apart in the graph. Then the nodes of one set in a well seperared pair are no longer guranteed to be close enough for the spanner property to hold
