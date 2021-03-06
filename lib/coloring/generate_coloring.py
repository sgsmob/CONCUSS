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

from lib.graph.graphformats import load_graph as load_graph
from lib.util.parse_config_safe import parse_config_safe
from lib.graph.graph import Coloring
from lib.coloring.basic.merge_colors import merge_colors

#Config = ConfigParser.ConfigParser()


def save_file(col, filename, override=False, verbose=False):
    num = len(col)
    if not override:
        if os.path.isfile(filename):
            with open(filename, 'r') as f:
                before = int(f.readline())
                if (before > num):
                    override = True
                    if verbose:
                        print "coloring is better, override result (before:",
                        print before, ", now:", num, ")"
                else:
                    if verbose:
                        print "don't override (before:", before, ", now:", num,
                        print ")"
        else:
            override = True
            if verbose: 
                print "there is still no coloring, create result"

    if override:
        with open(filename, 'w') as f:
            f.write(str(num) + '\n')
            for v in col:
                f.write("{0}: {1} \n".format(v, col[v]))


def import_colmodules(name):
    if not name:
        return None

    funcname = name.split('.')[-1]
    modname = "lib.coloring."+name
    module = __import__(modname, fromlist=[funcname])
    return getattr(module, funcname)


def ccalgorithm_factory(cfgfile, silent, execdata):
    Config = parse_config_safe(cfgfile)

    func_ldo = import_colmodules(Config.get('color',
                                            'low_degree_orientation'))
    func_step = import_colmodules(Config.get('color', 'step'))
    func_col = import_colmodules(Config.get('color', 'coloring'))
    func_ctd = import_colmodules(Config.get('color', 'check_tree_depth'))
    if Config.has_option('color', 'optimization'):
        func_opt = import_colmodules(Config.get('color', 'optimization'))
    else:
        func_opt = None

    if Config.has_option('color', 'preprocess'):
        func_preprocess = import_colmodules(Config.get('color',
                                                       'preprocess'))
    else:
        func_preprocess = None

    return CCAlgorithm(
            preprocess=func_preprocess,
            ldo=func_ldo,
            step=func_step,
            col=func_col,
            ctd=func_ctd,
            opt=func_opt,
            silent=silent,
            execdata=execdata)


class CCAlgorithm(object):

    def __init__(self, preprocess=None, ldo=None, step=None, col=None,
                 ctd=None, opt=None, silent=False, execdata=False,
                 profile=False):
        self.preprocess = preprocess
        self.ldo = ldo
        self.step = step
        self.col = col
        self.ctd = ctd
        self.opt = opt
        self.td = None
        self.silent = silent
        self.execdata=execdata
        self.profile = profile

    def echo(self, *msg):
        if not self.silent:
            print " ".join(map(str, msg))

    def start(self, rawgraph, treeDepth):
        self.td = treeDepth

        col = Coloring()

        trans = {}
        frat = {}

        if self.preprocess:
            if self.profile:
                preProfile = cProfile.Profile()
                preProfile.enable()
            self.echo("Preprocess coloring optimizations")
            pp_graph, postprocess = self.preprocess(rawgraph)
            if self.profile:
                preProfile.disable()
                printProfileStats( "preprocessing", preProfile)

        # Normalize graph so that its vertices are named 0, ..., n-1
        pp_graph.remove_loops()
        orig, mapping = pp_graph.normalize()
        for i in xrange(0, len(orig)):
            assert i in orig, "Vertex {0} not contained in " \
                    "norm. graph of size {1}".format(i, len(orig))

        g = self.ldo(orig)
        col = self.col(orig, g, trans, frat, col)
        # print col.color

        # User wants to output execution data
        if self.execdata:
            col_path = 'execdata/color/'
            if not os.path.exists(col_path):
                os.makedirs(col_path)

            os.makedirs(col_path + 'colorings')

            removed_vertices = rawgraph.nodes - pp_graph.nodes
            with open(col_path + 'colorings/0', 'w') as coloring_zero:
                # Write colors for vertices in the preprocessed graph to the
                # coloring file
                for vertex, clr in col.color.iteritems():
                    coloring_zero.write(
                        str(mapping[vertex]) + ": " + str(clr) + '\n')

        correct, nodes = self.ctd(orig, g, col, treeDepth)

        i = 0
        while (not correct):
            if self.profile:
                stepProfile = cProfile.Profile()
                stepProfile.enable()

            i += 1
            self.echo("step", i)

            g, trans, frat = self.step(orig, g, trans, frat, col, nodes, i,
                                       treeDepth, self.ldo)

            col = self.col(orig, g, trans, frat, col)
            # print col.color
            # for key in col:
            #     print str(mapping[key]) + ": " + str(col[key]),
            # User wants to output execution data
            if self.execdata:
                with open(col_path + 'colorings/' + str(i), 'w') as coloring_i:
                    # Write colors for vertices in the preprocessed graph to
                    # the coloring file
                    for vertex, clr in col.color.iteritems():
                        coloring_i.write(
                            str(mapping[vertex]) + ": " + str(clr) + '\n')

            correct, nodes = self.ctd(orig, g, col, treeDepth)

            if self.profile:
                stepProfile.disable()
                printProfileStats("step {0}".format(i), stepProfile)

            if correct:
                self.echo("  step", i, "is correct")
                break

        # end while
        self.echo("number of colors:", len(col))

        if self.opt:
            if self.profile:
                optProfile = cProfile.Profile()
                optProfile.enable()
            self.echo("Optimizing...")
            col = self.opt(orig, g, trans, frat, col, i, treeDepth, self)
            self.echo("number of colors:", len(col))
            if self.profile:
                optProfile.disable()
                printProfileStats( "optimizing", optProfile)

        # Map coloring back to original vertex labels
        colrenamed = Coloring()
        for v in col:
            colrenamed[mapping[v]] = col[v]


        if self.preprocess:
            if self.profile:
                postProfile = cProfile.Profile()
                postProfile.enable()
            self.echo("Postprocessing")
            col_restored = postprocess(colrenamed)

            self.echo("number of colors:", len(col_restored))
            if self.profile:
                postProfile.disable()
                printProfileStats("optimizing", postProfile)

        if self.profile:
            mergeProfile = cProfile.Profile()
            mergeProfile.enable()
        self.echo("Merging color classes")
        col_merged = merge_colors(rawgraph, col_restored, treeDepth)
        self.echo("number of colors:", len(col_merged))

        if self.execdata:
            for idx in range(i + 1):
                with open(col_path + 'colorings/' + str(idx), 'a') \
                        as coloring_i:
                    next_color = len(col)
                    for vertex in removed_vertices:
                        # If degree is 0, assign color 0
                        if rawgraph.degree(vertex) == 0:
                            coloring_i.write(str(vertex) + ": " + '0' + '\n')
                        # Degree not zero, hence assign color equal to length
                        # of coloring
                        else:
                            coloring_i.write(
                                str(vertex) + ": " + str(next_color) + '\n')
                            next_color += 1

            with open(col_path + 'colorings/' + str(i + 1), 'w') as coloring_i:
                # Write colors for vertices in the preprocessed graph to the
                # coloring file
                for vertex, clr in col_merged.color.iteritems():
                    coloring_i.write(str(vertex) + ": " + str(clr) + '\n')

        if self.profile:
            mergeProfile.disable()
            printProfileStats( "merging", mergeProfile)

        return col_merged


def start_coloring(filename, td, cfgfile, output):
    m = ccalgorithm_factory(cfgfile, False, None)

    p, fn = os.path.split(filename)
    graphname, e = os.path.splitext(fn)

    col = m.start(load_graph(filename), td)

    # Store results in common folder and, if supplied, in output file
    save_file(col, 'colorings/' + graphname + str(td), False)
    if output:
        save_file(col, output, True)


def printProfileStats(name, profile, percent=1.0):
    """
    Prints out the function call statistics using the cProfile and
    pstats libraries

    Arguments:
        name:  string labelling the purpose of the statistics
        profile:  cProfile to print
        percent:  decimal proportion of list to print.  Default prints
                all (1.0)
    """
    sortby = 'time'
    restrictions = ""
    ps = pstats.Stats(profile).strip_dirs().sort_stats(sortby)
    print "Stats from {0}".format(name)
    ps.print_stats(restrictions, percent)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("graph", help="filename of the graph", type=str)
    parser.add_argument("treeDepth", help="tree depth", type=int)
    parser.add_argument("config", help="filename of the config file", type=str,
                        nargs='?', default='config/default.cfg')
    parser.add_argument("-o", "--output", help="filename of the result",
                        type=str, nargs='?', default=None)
    args = parser.parse_args()

    start_coloring(args.graph, args.treeDepth, args.config, args.output)
