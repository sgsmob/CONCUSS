import copy
from collections import deque, defaultdict

from lib.graph.graph import Graph
from lib.graph.td_decomposition import TDDecomposition
from lib.util.recordtype import recordtype

# Define a record type of information about the vertices
# Attributes:
#     parent:  vertex that is the parent in the decomposition.  None iff vertex
#         is the root
#     children:  list of vertices that are children of this vertex
#     depth:  integer of the depth of the vertex.  Root has depth 0.
VertexInfo = recordtype('VertexInfo', 'parent children depth centers_beneath '
                        'colors_beneath ancestor_N')


class ColorSetTDDecomposition(TDDecomposition):
    """
    A subclass of Graph with additional data structures to describe a
    treedepth decomposition.

    Attributes:
        maxDepth:  The depth of the decomposition
        vertexRecords:  A dictionary mapping a vertex to its VertexInfo
            record (see above)
        root:  the root of the decomposition
    """

    def __init__(self, coloring, v=None):
        super(ColorSetTDDecomposition, self).__init__()
        self.coloring = coloring
        if v is not None:
            self.add_node(v)
            self.update_parent_child(None, v)

    def __default_vertex_record(self, v):
        # print "The right one"
        return VertexInfo(parent=None,
                          children=[],
                          depth=None,
                          centers_beneath=set([v]),
                          colors_beneath=set([self.coloring[v]]),
                          ancestor_N=set())

    def add_node(self, u):
        Graph.add_node(self, u)
        self.vertexRecords.extend(
            [None for x in range(len(self.vertexRecords), u + 1)])
        self.vertexRecords[u] = self.__default_vertex_record(u)

    def add_nodes_from(self, U):
        """Adds an iterable of nodes at once to save time"""
        Graph.add_nodes_from(self, U)
        # Make a local reference to self.vertexRecords for use in the loop
        vertexRecords = self.vertexRecords
        for u in U:
            self.vertexRecords.extend(
                [None for x in range(len(vertexRecords), u + 1)])
            vertexRecords[u] = self.__default_vertex_record(u)

    # def update_parent_child(self, parent, v):
    #     if self.vertexRecords[v].parent is not None:
    #         self.vertexRecords[v].parent
    #     super(ColorSetTDDecomposition, self).update_parent_child(parent, v)

    def update_centers(self, parent, v):
        vr = self.vertexRecords
        centers_to_add = set()
        for center in vr[v].centers_beneath:
            if self.coloring[center] not in vr[parent].colors_beneath:
                centers_to_add.add(center)
        centers_to_remove = set()
        existing_centers = vr[parent].centers_beneath - vr[v].centers_beneath
        for center in existing_centers:
            if self.coloring[center] in vr[v].colors_beneath:
                centers_to_remove.add(center)
        vr[parent].centers_beneath -= centers_to_remove
        vr[parent].centers_beneath |= centers_to_add
        vr[parent].colors_beneath |= vr[v].colors_beneath

    def update_ancestor_N(self, parent, v):
        vr = self.vertexRecords
        vr[parent].ancestor_N |= vr[v].ancestor_N
        vr[parent].ancestor_N.discard(parent)

    def level_order_iter(self, start=None):
        """
        Iterate level-wise (from lowest depth (nearest to the root) to the
        greatest depth) through the subtree rooted at start.  If start is None,
        begin at the root.
        """
        if start is None:
            start = self.root
        Q = deque([start])
        while len(Q) > 0:
            curr = Q.popleft()
            Q.extend(self.vertexRecords[curr].children)
            yield curr

    def subtree_depths(self):
        """Compute the depth of the subdecomposition rooted at each vertex."""
        subtree_depth = dict()
        for u in reversed(self.level_order_iter()):
            children_depths = [subtree_depth[w] for w in
                               self.vertexRecords[u].children]
            children_depths.append(0)
            subtree_depth[u] = max(children_depths) + 1
        return subtree_depth

    def swap(self, v, u=None):
        """
        Rearrange the tree structure so that v is the parent of its current
        parent.  This function does *NOT* verify that the selected swap
        preserves the treedepth decomposition invariant that all edges join
        ancestor-descendant pairs.

        :param u: An ancestor of v.  If provided, repeat the swapping until u
                  is a child of v.
        """
        curr = v
        # alias for shorter lines below
        vr = self.vertexRecords
        finished = False
        while not finished:
            old_parent = vr[curr].parent
            if u is None or u == old_parent:
                finished = True
            elif old_parent is None:
                raise ValueError("Cannot swap root of treedepth decomposition "
                                 "{}".format(curr))
            # update the actual parent child relation
            grandparent = vr[old_parent].parent
            print "swapping {} with {} "\
                "whose parent is {}".format(curr, old_parent, grandparent)
            vr[old_parent].children.remove(curr)
            if grandparent is not None:
                vr[grandparent].children.remove(old_parent)
            self.update_parent_child(grandparent, curr)
            # if possible, all the children of curr should stay children of
            # curr.  this is only impossible if a descendant of that child is
            # neighbors with old_parent.
            curr_children = vr[curr].children
            old_parent_children = vr[old_parent].children
            vr[curr].children = []
            vr[old_parent].children = []
            curr_ancestor_N = vr[curr].ancestor_N
            old_parent_ancestor_N = vr[old_parent].ancestor_N
            vr[curr].ancestor_N = old_parent_ancestor_N
            vr[old_parent].ancestor_N = old_parent_ancestor_N &\
                set(self.neighbours(old_parent))
            self.update_parent_child(curr, old_parent)
            for x in old_parent_children:
                print "\thandling", x
                if x == curr:
                    print "\t  x == curr"
                    continue
                if old_parent in vr[x].ancestor_N:
                    print "\t  old_parent in", vr[x].ancestor_N
                    self.update_parent_child(old_parent, x)
                    # since old_parent has moved down in levels, we need to
                    # update the levels of its children.
                    for y in self.level_order_iter(x):
                        vr[y].depth = vr[vr[y].parent].depth + 1
                    # need to make sure old_parent knows which of its
                    # descendants have neighbors above it.
                    for y in vr[x].ancestor_N:
                        if y != old_parent:
                            vr[old_parent].ancestor_N.add(y)
                else:
                    print "\t  old_parent not in", vr[x].ancestor_N
                    self.update_parent_child(curr, x)
                    # need to make sure curr knows which of its
                    # descendants have neighbors above it.
                    # for y in vr[x].ancestor_N:
                    #     if y != curr:
                    #         vr[curr].ancestor_N.add(y)
            for x in curr_children:
                print "\thandling", x
                if x == old_parent:
                    print "\t  x == curr"
                    continue
                # x becomes a child of old_parent
                if old_parent in vr[x].ancestor_N:
                    print "\t  old_parent in", vr[x].ancestor_N
                    self.update_parent_child(old_parent, x)
                    # need to make sure old_parent knows which of its
                    # descendants have neighbors above it.
                    for y in vr[x].ancestor_N:
                        if y != old_parent:
                            vr[old_parent].ancestor_N.add(y)
                # x stays a child of curr
                else:
                    print "\t  old_parent not in", vr[x].ancestor_N
                    self.update_parent_child(curr, x)
                    # since curr has moved up in levels, we need to update the
                    # levels of its children.
                    for y in self.level_order_iter(x):
                        vr[y].depth = vr[vr[y].parent].depth + 1
                    # need to make sure old_parent knows which of its
                    # descendants have neighbors above it.
                    # for y in vr[x].ancestor_N:
                    #     if y != curr:
                    #         vr[curr].ancestor_N.add(y)
            print "child results"
            print "\t old:", vr[old_parent].children, vr[old_parent].ancestor_N
            print "\tcurr:", vr[curr].children, vr[curr].ancestor_N
            assert grandparent == vr[curr].parent
            if old_parent in self.neighbours(curr):
                vr[old_parent].ancestor_N.add(curr)

            # update the subtree color and center information
            vr[curr].centers_beneath = vr[old_parent].centers_beneath
            vr[curr].colors_beneath = vr[old_parent].colors_beneath
            center_color_counter = defaultdict(int)
            center_color_counter[self.coloring[old_parent]] += 1
            centers = set([old_parent])
            colors = set([self.coloring[old_parent]])
            for w in vr[old_parent].children:
                for c in vr[w].centers_beneath:
                    center_color_counter[self.coloring[c]] += 1
                    centers.add(c)
                colors.update(vr[w].colors_beneath)
            print "color results:"
            print "\t", center_color_counter
            print "\t", colors
            vr[old_parent].centers_beneath = set(c for c in centers if
                                                 center_color_counter[
                                                    self.coloring[c]] == 1)
            vr[old_parent].colors_beneath = colors

    def __find_next_center(self, curr, centers, unique_colors):
        if self.coloring[curr] in unique_colors:
            return curr
        else:
            for new_center in centers:
                if self.coloring[new_center] in unique_colors:
                    self.swap(new_center, curr)
                    return new_center
            else:
                raise ValueError("None of the centers have unique_colors")

    def __find_removal_order(self, other, joined, v, self_N, v_root,
                             has_unrelated_neighbors):
        print "finding ordering above {}".format(v)
        assert v_root is not None
        L_ancestors = self.rootPath(v_root)
        R_ancestors = other.rootPath(v)
        L_curr = self.root
        R_curr = other.root
        if has_unrelated_neighbors:
            L_target_depth = len(L_ancestors)
        else:
            L_target_depth = len(set(L_ancestors) & self_N)
        # R_target_depth = len(R_ancestors)
        L_depth = 0
        R_depth = 0
        for x in sorted(joined):
            print "\t", x, joined.vertexRecords[x]
        removal_order = []
        center = None
        while center != v:
            if L_depth < L_target_depth:
                L_centers = joined.vertexRecords[L_curr].centers_beneath
                L_unique = set(joined.coloring[x] for x in L_centers)
                L_all_colors = joined.vertexRecords[L_curr].colors_beneath
                assert L_curr in L_centers, "{} not in {}".format(L_curr,
                                                                  L_centers)
            else:
                L_unique = set()
                L_all_colors = set()
            R_centers = joined.vertexRecords[R_curr].centers_beneath
            R_unique = set(joined.coloring[x] for x in R_centers)
            R_all_colors = joined.vertexRecords[R_curr].colors_beneath
            assert R_curr in R_centers, "{} not in {}".format(R_curr,
                                                              R_centers)
            L_only = L_unique - R_all_colors
            R_only = R_unique - L_all_colors
            print "L_only", L_only
            print "R_only", R_only
            if R_only:
                # removal_order.append(R_curr)
                center = joined.__find_next_center(R_curr, R_centers,
                                                   R_only)
                R_depth += 1
                print "Center {} has children {} and parent {}"\
                    "".format(center,
                              joined.vertexRecords[center].children,
                              joined.vertexRecords[center].parent)
                possible = []
                for x in joined.vertexRecords[center].children:
                    if x in R_ancestors:
                        possible.append(x)
                if len(possible) == 1:
                    R_curr = possible[0]
                elif len(possible) > 1:
                    for x in possible:
                        if x != R_curr:
                            R_curr = x
                else:
                    assert center == v, "No child of {} is in the "\
                                        "ancestors {}".format(R_curr,
                                                              R_ancestors)
            elif L_only:
                # removal_order.append(L_curr)
                center = joined.__find_next_center(L_curr, L_centers,
                                                   L_only)
                if has_unrelated_neighbors or center in self_N:
                    L_depth += 1
                print "Center {} has children {} and parent {}"\
                    "".format(center,
                              joined.vertexRecords[center].children,
                              joined.vertexRecords[center].parent)
                possible = []
                for x in joined.vertexRecords[center].children:
                    if x in L_ancestors:
                        possible.append(x)
                if len(possible) == 1:
                    L_curr = possible[0]
                elif len(possible) > 1:
                    for x in possible:
                        if x != L_curr:
                            L_curr = x
                else:
                    assert center == v_root, "No child of {} is in the "\
                                             "ancestors {}".format(L_curr,
                                                                   L_ancestors)
            else:
                raise ValueError("There should be a center here somewhere")
            removal_order.append(center)
        if L_depth < L_target_depth:
            removal_order.append(L_curr)
        return removal_order

    def join(self, other, v, N):
        """
        Create a new treedepth decomposition that is the union of 'self' and
        'other'.

        :param other: a TDDecomposition
        :param v: a vertex in 'other'
        :param N: the neighbors of 'v' in the graph from which this
                  TDDecomposition was derived.
        Precondition: 'v' is the only vertex in 'other' with neighbors in this
                      decomposition.
        """
        # first add all existing vertices to the new decomposition
        # we copy the vertex records so that we can modify the records in place
        joined = self.__class__(self.coloring)
        joined.add_nodes_from(self)
        for u in self:
            assert u in joined
            assert u in self
            joined.vertexRecords[u] = copy.deepcopy(self.vertexRecords[u])
        # TODO just copy directly
        for u, w in self.edges():
            joined.add_edge(u, w)
        joined.add_nodes_from(other)
        for u in other:
            assert u in joined
            assert u in other
            joined.vertexRecords[u] = copy.deepcopy(other.vertexRecords[u])
        for u, w in other.edges():
            joined.add_edge(u, w)
        self_N = N & self.nodes
        for u in self_N:
            joined.add_edge(u, v)

        # determine the maximum at which v could be placed, which is the
        # maximum depth at which there is a vertex `v_root` such that all of
        # v's neighbors are either ancestors or descendants of v_root.
        print "Finding lowest position of {}".format(v)
        # first, determine the levels of v's neighbors in self
        N_by_level = []
        for u in self_N:
            depth = self.vertexRecords[u].depth
            if len(N_by_level) <= depth:
                N_by_level.extend([list() for _ in range(len(N_by_level),
                                   depth+1)])
            N_by_level[depth].append(u)
        # the set of all vertices at a given level that are ancestors of a
        # neighbor of v.
        active = set()
        # we need to distinguish cases where we have already found a vertex
        # that could possible be `v_root` against those where we have not found
        # such a candidate.  Otherwise we would replace `v_root` with one of
        # its ancestors.
        valid = False
        has_unrelated_neighbors = False
        print "N_by_level", N_by_level
        # iterate over levels from bottom to top
        for U in reversed(N_by_level):
            print "\tLoop with U = {}".format(U)
            for u in U:
                active.add(u)
            if len(active) == 1:
                if not valid:
                    v_root = iter(active).next()
                    v_children = self.vertexRecords[v_root].children
                    valid = True
            else:
                assert len(active) > 0
                valid = False
                has_unrelated_neighbors = True
            # to go to the next level get all the parents of the active set.
            active = set(self.vertexRecords[u].parent for u in active)
        print "Found parent {} and children {}".format(v_root, v_children)

        removal_order = self.__find_removal_order(other, joined, v,
                                                  self_N, v_root,
                                                  has_unrelated_neighbors)

        # organize the vertices of the new decomposition to respect parent-
        # child relationships.
        print "Updating parent child relationships in order", removal_order
        first = removal_order[0]
        joined.update_parent_child(None, first)
        print "\t", first, joined.vertexRecords[first]
        for u, u_parent in zip(removal_order[1:], removal_order):
            print "\t", u_parent, joined.vertexRecords[u_parent]
            print "\t", u,  joined.vertexRecords[u]
            old_parent = joined.vertexRecords[u].parent
            if old_parent != u_parent:
                if old_parent is not None:
                    joined.vertexRecords[old_parent].children.remove(u)
                joined.update_parent_child(u_parent, u)
        for u, u_parent in reversed(zip(removal_order[1:], removal_order)):
            joined.update_centers(u_parent, u)
            joined.update_ancestor_N(u_parent, u)
            print "\t", u,  joined.vertexRecords[u]
            print "\t", u_parent, joined.vertexRecords[u_parent]

        # ensure the depths are properly propagated
        print "Fixing depths"
        Q = deque(removal_order)
        max_depth = max(self.maxDepth, other.maxDepth)
        while Q:
            curr = Q.popleft()
            depth = joined.vertexRecords[curr].depth
            print depth, curr, Q
            if depth > max_depth:
                max_depth = depth
            print joined.vertexRecords[curr].children
            for child in joined.vertexRecords[curr].children:
                if child not in Q and\
                        joined.vertexRecords[child].depth != depth+1:
                    print "\t\t", child, joined.vertexRecords[child].depth
                    joined.vertexRecords[child].depth = depth+1
                    Q.append(child)
        joined.maxDepth = max_depth
        # ensure the ancestor neighbors are properly propagated
        for u in joined.vertexRecords[v].children:
            joined.update_ancestor_N(v, u)
        for u in self_N:
            # u is a descendant of v
            if joined.vertexRecords[u].depth > joined.vertexRecords[v].depth:
                curr = u
                print "\tupdating", curr
                while curr != v:
                    joined.vertexRecords[curr].ancestor_N.add(v)
                    curr = joined.vertexRecords[curr].parent
            # u is an ancestor of v
            else:
                joined.vertexRecords[v].ancestor_N.add(u)
        return joined

    def verify(self):
        vr = self.vertexRecords
        print "nodes = ", self.nodes
        for v in sorted(self):
            print "vertex records for", v, ":", vr[v]
            if vr[v].parent is None:
                assert v == self.root, "{} has no parent".format(v)
                assert vr[v].depth == 0
            else:
                parent = vr[v].parent
                assert vr[v].depth == vr[parent].depth + 1, "{} {}".format(vr[v].depth, vr[parent].depth)
                assert v in vr[parent].children
        for v in self:
            v_depth = vr[v].depth
            for u in self.neighbours(v):
                assert v in self.neighbours(u)
                u_depth = vr[u].depth
                if u_depth < v_depth:
                    curr = v
                    for _ in range(v_depth - u_depth):
                        curr = vr[curr].parent
                    assert curr == u
                    assert u in vr[v].ancestor_N
                elif u_depth == v_depth:
                    raise Exception("can't have cross edges between {} and "
                                    "{}".format(u, v))
        # real_colors_beneath = defaultdict(set)
        multiple_below = defaultdict(set)
        color_counts = defaultdict(lambda: defaultdict(int))
        ancestor_N = defaultdict(set)
        for v in reversed(list(self.level_order_iter())):
            # print "Checking subtree properties of {}".format(v)
            for u in self.neighbours(v):
                if vr[v].depth > vr[u].depth:
                    ancestor_N[v].add(u)
                    # print "adding", u, ancestor_N[v]
                # else:
                    # print "{} is below".format(u)
            color_counts[v][self.coloring[v]] = 1
            for u in vr[v].children:
                for color, count in color_counts[u].iteritems():
                    color_counts[v][color] += count
                # real_colors_beneath[v].add(color)
                multiple_below[v] |= multiple_below[u]
                ancestor_N[v] |= ancestor_N[u]
            # print "ancestor_N", ancestor_N[v]
            ancestor_N[v].discard(v)
            ancestors_diff = ancestor_N[v] ^ vr[v].ancestor_N
            assert len(ancestors_diff) == 0, "{} {}".format(v, ancestor_N[v])
            real_unique = set()
            for color, count in color_counts[v].iteritems():
                if count > 1:
                    multiple_below[v].add(color)
                else:
                    assert count == 1
                    real_unique.add(color)
                # real_colors_beneath[v].add(color)
            real_colors_beneath = set(color_counts[v].keys())
            color_diff = real_colors_beneath ^ vr[v].colors_beneath
            assert len(color_diff) == 0, "{}".format(real_colors_beneath)
            assert len(real_unique) == len(vr[v].centers_beneath),\
                "{} {}".format(v, real_unique)
            for u in vr[v].centers_beneath:
                assert self.coloring[u] in real_unique
