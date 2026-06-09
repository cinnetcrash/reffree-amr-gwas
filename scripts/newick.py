"""
newick.py — tiny dependency-free Newick parser.

Supports branch lengths, named tips, optional internal labels/support values,
and quoted labels. Enough for kinship and parsimony on a midpoint-rooted core
phylogeny. Not a general-purpose library — no NHX, no comments.
"""


class Node:
    __slots__ = ("name", "length", "children", "parent")

    def __init__(self, name=None, length=0.0):
        self.name = name
        self.length = length
        self.children = []
        self.parent = None

    def is_leaf(self):
        return not self.children

    def add(self, child):
        child.parent = self
        self.children.append(child)

    def leaves(self):
        if self.is_leaf():
            yield self
        for c in self.children:
            yield from c.leaves()

    def postorder(self):
        for c in self.children:
            yield from c.postorder()
        yield self

    def preorder(self):
        yield self
        for c in self.children:
            yield from c.preorder()


def parse(text):
    """Parse a single Newick string, return the root Node."""
    text = text.strip()
    if not text.endswith(";"):
        raise ValueError("Newick string must end with ';'")
    s = text[:-1]
    pos = 0

    def read_label():
        nonlocal pos
        if pos < len(s) and s[pos] == "'":
            pos += 1
            start = pos
            while pos < len(s) and s[pos] != "'":
                pos += 1
            label = s[start:pos]
            pos += 1
            return label
        start = pos
        while pos < len(s) and s[pos] not in "(),:;":
            pos += 1
        return s[start:pos] or None

    def read_length():
        nonlocal pos
        if pos < len(s) and s[pos] == ":":
            pos += 1
            start = pos
            while pos < len(s) and s[pos] not in "(),;":
                pos += 1
            try:
                return float(s[start:pos])
            except ValueError:
                return 0.0
        return 0.0

    def read_node():
        nonlocal pos
        node = Node()
        if pos < len(s) and s[pos] == "(":
            pos += 1
            node.add(read_node())
            while pos < len(s) and s[pos] == ",":
                pos += 1
                node.add(read_node())
            if pos >= len(s) or s[pos] != ")":
                raise ValueError("unbalanced parentheses")
            pos += 1
            node.name = read_label()      # internal label / support (ignored)
        else:
            node.name = read_label()
        node.length = read_length()
        return node

    root = read_node()
    return root


def root_to_node_depths(root):
    """Cumulative branch length from root to every node."""
    depth = {root: 0.0}
    for n in root.preorder():
        if n is root:
            continue
        depth[n] = depth[n.parent] + n.length
    return depth


def leaf_names(root):
    return [lf.name for lf in root.leaves()]
