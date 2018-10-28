#
# This file is part of CONCUSS, https://github.com/theoryinpractice/concuss/,
# and is Copyright (C) North Carolina State University, 2015. It is licensed
# under the three-clause BSD license; see LICENSE.
#


import sys
from copy import copy, deepcopy

from dfs_sweep import DFSSweep
from lib.graph.color_set_td_decomposition import ColorSetTDDecomposition


class SaveTDDDFSSweep(DFSSweep):

    def add(self, data, color):
        """
        Add one color to the color combination and find all new
        connected components using a union find data structure

        :param G: The parent graph
        :param data: The WCSData object
        :param color: The color to add

        """
        # Add the color to the combination
        data.colors_in_combi.add(color)
        # Increment current depth since we have added a color
        data.current_depth += 1
        # Find how many vertices we have
        n = data.max + 1
        # If we are at looking at the first color
        if data.current_depth == 0:
            # Make a ufs structure and initialize it with 0s
            ufs = [0] * n

            # Initialize a component store to store our components
            # Components are stored as sets in a dictionary
            comps = {}

            # For all the vertices of that color
            for v in data.wscolorlist[color].nodes:
                # Initialize the union find entry
                # Make each vertex a root
                ufs[v] = (1 << 2) | self.UFS_TYPE_ROOT
                # Since each vertex is a component by itself
                # Make a set containing only v and add it to the dictionary
                comps[v] = ColorSetTDDecomposition(self.coloring, v)

        # We are looking at at least two colors
        else:
            for x in data.component_store:
                print x
            print "Lengths", data.current_depth, len(data.component_store)
            # Copy the previous dictionary of components
            comps = copy(data.component_store[data.current_depth - 1])
            # Copy the previous union find structure
            ufs = copy(data.union_find[data.current_depth - 1])
            # Initialize new colors in the UFS and comps dictionary
            for v in data.wscolorlist[color].nodes:
                ufs[v] = (1 << 2) | self.UFS_TYPE_ROOT
                comps[v] = ColorSetTDDecomposition(self.coloring, v)

            for v in data.wscolorlist[color].nodes:
                # For each newly added vertex, find its neighbors
                # and check to see if it already belongs to some component
                for u in self.G.neighbours(v):
                    if self.coloring[u] in data.colors_in_combi:
                        root1 = self.ufs_find(ufs, v)
                        root2 = self.ufs_find(ufs, u)
                        if root1 != root2:
                            # Get the UFS entries of these roots
                            ufs_root1 = ufs[root1]
                            ufs_root2 = ufs[root2]
                            # Find which root to append to and which root gets
                            # appended.  'a' is the one that 'd' will be
                            # appended to.
                            a, d = (root1, root2) if ((ufs_root1 >> 2) >
                                                      (ufs_root2 >> 2)) \
                                else (root2, root1)

                            # Increment the count of vertices in 'a'
                            total_count = (ufs_root1 >> 2) + (ufs_root2 >> 2)
                            ufs[a] = (total_count << 2) | self.UFS_TYPE_ROOT
                            # Set 'd' as a child of a
                            # union operation (set parent)
                            ufs[d] = self.UFS_TYPE_CHILD | (a << 2)

                            print "Before"
                            print "\ta ({}) nodes:".format(a), comps[a].nodes
                            print "\td ({}) nodes:".format(d), comps[d].nodes
                            # Join the treedepth decompositions
                            new_comps = comps[root2].join(comps[root1], v,
                                                          self.G.neighbours(v))
                            print "After"
                            print "\ta ({}) nodes:".format(a), comps[a].nodes
                            print "\td ({}) nodes:".format(d), comps[d].nodes
                            comps[a].verify()
                            comps[d].verify()
                            new_comps.verify()
                            comps[a] = new_comps
                            # Delete key since we don't need it anymore
                            # del comps[d]

        # Append the new union find to our stack of union finds
        data.union_find.append(ufs)
        # Append the new component dictionary to our
        # stack of component dictionaries
        data.component_store.append(comps)

    def induced_color_components(self):
        for comp in self.walk_color_space():
            yield comp
