##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# This file processes the phylogenetic tree produced by IQtree, and outputs two prcessed trees.
# Output Tree 1:
#   Rerooted using the reference sequence as an out group, and ladderised.
# Output Tree 2:
#   Rerooted as above, also branches are collapsed to polytomies if their length or support
#   values are below a threshold (defaults: min length=5e-5, min support=bootstrap 75)
#
# REQUIREMENTS:
# Python 3 and the dendropy package.

def parse_args():
    import argparse
    parser = argparse.ArgumentParser('Reroot and collapse unsupported edges in the tree')
    parser.add_argument('treefile', help='The newick tree to process')
    parser.add_argument('-l', '--min-length', default=0.00005, type=float, help='Minimum length of an edge in the output tree')
    parser.add_argument('-s', '--min-support', default=75, type=float, help='Minimum bootstrap support of an edge in the output tree')
    return parser.parse_args()


if __name__ == '__main__':
    import dendropy as dpy
    import os, sys

    args = parse_args()
    input_file = args.treefile
    min_length = args.min_length
    min_support = args.min_support

    if not os.path.exists(input_file):
        sys.stderr.write('{} not found\n'.format(input_file))
        sys.exit()

    filestub, ext = os.path.splitext(input_file)
    rooted_file = filestub + '_rooted' + ext
    collapsed_file = filestub + '_rooted_collapsed' + ext

    tree = dpy.Tree.get_from_path(input_file, schema='newick', preserve_underscores = True, rooting='force-rooted')
    tree.resolve_polytomies()
    tree.update_bipartitions()
    ref_node = tree.find_node_with_taxon_label("reference_NA_NA")
    tree.reroot_at_edge(ref_node.edge, length1=ref_node.edge_length / 2.0, length2=ref_node.edge_length / 2.0, update_bipartitions=True)
    tree.ladderize(ascending=False)
    tree.write_to_path(rooted_file, schema="newick")

    tree = dpy.Tree.get_from_path(rooted_file, schema='newick', preserve_underscores = True, rooting='force-rooted')

    for edge in tree.postorder_edge_iter():
        collapse = False
        if edge.is_internal():
            # Very small edge? -> mark for collapse
            if edge.length is not None and edge.length < min_length:
                collapse = True

            # Poorly supported edge? -> mark for collapse
            elif edge.head_node.label is not None:
                if '/' in edge.head_node.label:
                    bootstrap = float(edge.head_node.label.split('/')[1])
                else:
                    bootstrap = float(edge.head_node.label)

                if bootstrap < min_support:
                    collapse = True

            # Do the collapsing...
            if collapse:
                for child in edge.head_node.child_nodes():
                    child.edge.length += edge.length
                edge.collapse()
    tree.ladderize(ascending=False)
    tree.write_to_path(collapsed_file, schema="newick")
