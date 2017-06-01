#!/usr/bin/python2
"""Run the DP on a decomposition with large treedepth."""

import cProfile
import argparse
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("graph", type=str, help="Graph file")
    parser.add_argument("pattern", type=str,
                        help="Name of pattern to search for")
    parser.add_argument("coloring", type=str, help="Coloring file")
    parser.add_argument("p", type=int,
                        help="Maximum numbers to consider in decomposition")
    args = parser.parse_args()
    G = load_graph(args.graph)
    coloring = coloring_from_file(args.coloring, None, None, None, None)
    H, td_lower = get_pattern_from_generator(args.pattern)
    p = args.p

    prof = cProfile.Profile()
    prof.enable()
    count = count_large_treedepth(G, H, coloring, p, td_lower)
    prof.disable()

    print "Number of occurrences of H in G: {}".format(count)

    name = "{} {}".format(args.graph, args.coloring)
    printProfileStats(name, prof)
