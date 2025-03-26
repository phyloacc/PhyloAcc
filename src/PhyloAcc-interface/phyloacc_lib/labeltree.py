#############################################################################
# Helper function to transfer labels from a source tree to an input tree.
# Necessary for the trees generated from ASTRAL to be input into PhyloAcc
# since the labels need to match the tree in the mod file.
#
# Gregg Thomas, March 2025
#############################################################################

import phyloacc_lib.tree as TREE
import re

#############################################################################

def transferLabels(source_tree_file, input_tree_file, output_tree_file):
    with open(source_tree_file, 'r') as source_tree_stream, open(input_tree_file, 'r') as input_tree_stream:
        source_tree = TREE.Tree(source_tree_stream.read());
        input_tree = TREE.Tree(input_tree_stream.read());
    # Read the source and input trees from their respective files and
    # create Tree objects for both trees

    input_tree.label.update(source_tree.label);
    # Transfer labels from source_tree to input_tree

    label_pattern = r'(\)<\d+>)(?!,)(?:<\d+>|[^:<>\(\)]+)+(?=:)';
    output_tree_str = re.sub(label_pattern, r'\1', input_tree.tree_str);
    # Make a new version of the tree string, removing any labels from the input tree (but 
    # preserving our own <#> labels)

    #print("Input tree string:", input_tree.tree_str);
    #print("Input tree string with labels removed:", output_tree_str);
    # Print the original and modified input tree strings for debugging

    for orig_label in input_tree.label:
        if orig_label not in input_tree.tips:
            output_tree_str = output_tree_str.replace(orig_label, input_tree.label[orig_label]);
    # Replace internal labels with their corresponding values

    output_tree_str += ";";
    # Add a semicolon to the end of the tree string

    #print("Output tree string:", output_tree_str);

    with open(output_tree_file, 'w') as output_tree_stream:
        output_tree_stream.write(output_tree_str);
    # Write the modified tree string to the output file

#############################################################################