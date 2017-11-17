#!/usr/bin/env python2.7
#
# This file is part of CONCUSS, https://github.com/theoryinpractice/concuss/,
# and is Copyright (C) North Carolina State University, 2015. It is licensed
# under the three-clause BSD license; see LICENSE.
#


import sys
import os
import argparse
import cProfile
import pstats

from lib.coloring.generate_coloring import start_coloring


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("graph", help="filename of the graph", type=str)
    parser.add_argument("treeDepth", help="tree depth", type=int)
    parser.add_argument("config", help="filename of the config file", type=str,
                        nargs='?', default='config/default.cfg')
    parser.add_argument("-o", "--output", help="filename of the result",
                        type=str, nargs='?', default=None)
    parser.add_argument("-i", "--intermediate",
                        help="filename of intermediate colorings",
                        type=str, nargs='?', default=None)
    args = parser.parse_args()

    start_coloring(args.graph, args.treeDepth, args.config, args.output,
                   args.intermediate)
