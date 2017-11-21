#!/usr/bin/python2
"""Run the DP on a decomposition with large treedepth."""

import cProfile
import argparse
from collections import deque
from lib.pattern_counting.dp import MemoizedBVKPattern
from lib.pattern_counting.double_count import InclusionExclusion
from lib.decomposition import DFSSweep
from lib.pattern_counting.pattern_counter import PatternCounter
from lib.graph.graphformats import load_graph
from lib.run_pipeline import (coloring_from_file, get_pattern_from_generator,
                              printProfileStats)


def count_large_treedepth(G, H, coloring, p, td_lower):
    table_hints = {"forward": False, "reuse": False}
    counter = PatternCounter(G,
                             [H],
                             [td_lower],
                             coloring,
                             pattern_class=MemoizedBVKPattern,
                             combiner_class=InclusionExclusion,
                             table_hints=table_hints)
    gen = DFSSweep(G,
                   coloring,
                   p,
                   td_lower,
                   len(H),
                   [counter.combiners[0].before_color_set],
                   [counter.combiners[0].after_color_set])
    counter.decomp_generator = gen

    return counter.count_patterns()[0]


def null_function(*args):
    pass

    
def count_depths(G, H, coloring, p, td_lower):
    gen = DFSSweep(G,
                   coloring,
                   p,
                   td_lower,
                   len(H),
                   [null_function],
                   [null_function])
    total_depths = 0
    total_v = 0
    for decomp in gen:
        print decomp.root, len(decomp.vertexRecords)
        # create a post order traversal ordering with a DFS to use in the DP
        ordering = []
        q = deque([decomp.root])
        # print decomp.root, len(decomp),
        # print [(i+1,self.coloring[i]) for i in decomp]
        while q:
            curr = q.pop()
            ordering.append(curr)
            if not decomp.hasLeaf(curr):
                q.extend(reversed(decomp.children(curr)))
        ordering.reverse()
        num_colors = set(coloring[v] for v in ordering)

        if len(num_colors) == p:
            depths = [v.depth for v in ordering]
            total_depths += sum(depths)
            total_v += len(depths)

    print "Average depth in TDDs", float(total_depths) / total_v
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("graph", type=str, help="Graph file")
    parser.add_argument("pattern", type=str,
                        help="Name of pattern to search for")
    parser.add_argument("coloring", type=str, help="Coloring file")
    parser.add_argument("p", type=int,
                        help="Maximum numbers to consider in decomposition")
    parser.add_argument("-d", "--depth_only", action='store_true',
                        help="average the depths of all vertices")
    args = parser.parse_args()
    G = load_graph(args.graph)
    coloring = coloring_from_file(args.coloring, None, None, None, None)
    H, td_lower = get_pattern_from_generator(args.pattern)
    p = args.p

    count_depths(G, H, coloring, p, td_lower)
    if not args.depth_only:
        prof = cProfile.Profile()
        prof.enable()
        count = count_large_treedepth(G, H, coloring, p, td_lower)
        prof.disable()

        print "Number of occurrences of H in G: {}".format(count)

        name = "{} {}".format(args.graph, args.coloring)
        printProfileStats(name, prof)
