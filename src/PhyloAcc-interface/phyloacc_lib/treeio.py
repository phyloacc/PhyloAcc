#############################################################################
# Handles the reading of the input trees
# Gregg Thomas
#############################################################################

import sys
import os
import phyloacc_lib.core as CORE
import phyloacc_lib.tree as TREE
# import phyloacc_lib.tree_old as TREE

#############################################################################

def readST(globs, tree_type="species", rf_to_st=False):

    step = "Reading input " + tree_type + " tree";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    try:
        if tree_type == "species":
            for line in open(globs['mod-file']):
                if line.startswith("TREE: "):
                    tree_str = line.strip().replace("TREE: ", "");
            # Read the tree string from the MOD file
        else:
            tree_str = open(globs['coal-tree-file'], "r").read();

        tree = TREE.Tree(tree_str);

        # tree_dict, tree, root = TREE.treeParse(tree_str);
        # Parse the tree string to the tree class
    except:
        if tree_type == "species":
            CORE.errorOut("IO3", "Error reading species tree from .mod file!");
        else:
            CORE.errorOut("IO3", "Error reading tree from file!");
    # Read the tree

    tree.subtree = tree.genSubtrees();
    tree.tree_str = tree.subtree[tree.root];
    del tree.subtree;

    root_str = "unrooted";
    add_root_node = 0;
    if tree.rooted:
       root_str = "rooted";
       add_root_node = 1;
    # Node counting for log

    # tips, internals = [], [];
    # for node in tree_dict:
    #     if tree_dict[node][2] == "tip":
    #         tips.append(node);
    #     else:
    #         internals.append(node);

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + root_str + " tree read");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Tree has "+ str(tree.num_tips) + " tips and " + str(tree.num_internals + add_root_node) + " internal nodes");
    
    # step_start_time = CORE.report_step(globs, step, step_start_time, "Success: tree read");
    # CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Tree has "+ str(len(tips)) + " tips and " + str(len(internals)) + " internal nodes");
    # Status updates

    ## TODO: Compare tree topologies

    return tree_str, tree
    # return tree_str, tree_dict, tree, root, tips, internals;

#############################################################################

def readSpeciesGroups(globs):

    step = "Reading species/branch groups";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    for group in globs['groups']:
    # Check every group 

        for label in globs['input-groups'][group]:
        # Check every label in the current group for presence in tree and whether it is an internal label

            label_found = False;
            # A flag that lets us know if the label was found in the tree

            for node in globs['st'].nodes:
            # for node in globs['st']:
             # For every node read in the tree

                if node in globs['st'].tips:
                # if node in globs['tips']:
                    if label == node:
                        label_found = True;
                        globs['groups'][group].append(label);
                        break;
                # Check if the current label is a tip in the tree and if so add it to the global group and break

                else:
                    if label in [ node, globs['st'].label[node] ]:
                    # if label in [n, globs['st'][n][3]]:
                        label_found = True;
                        globs['groups'][group] += globs['st'].getClade(node);
                        # cur_clade = TREE.getClade(n, globs['st']);
                        # globs['groups'][group] += cur_clade;
                        break;
                # Check if the current label is an internal label (either from the tree input or from treeparse/--labeltree) and
                # if so get all the tips descending from it to add to the global group and break

            if not label_found:
                CORE.errorOut("IO3", "The following label was provided in a group but is not present in the tree: " + label, globs);
            # If the current label wasn't found in the tree, exit here with an error

        ## End node loop
    ## End group loop
    # A preliminary check here to make sure all provided labels are actually in the input tree

    ##########

    if not globs['groups']['conserved']:
        globs['groups']['conserved'] = [ node for node in globs['st'].tips if node not in globs['groups']['targets'] and node not in globs['groups']['outgroup'] ];
        # globs['groups']['conserved'] = [ node for node in globs['tips'] if node not in globs['groups']['targets'] and node not in globs['groups']['outgroup'] ];
    # If the conserved group isn't specified, add all remaining tip species not in the other groups to it here

    for tip in globs['st'].tips:
    # for tip in globs['tips']:
        num_in_groups = 0;
        for group in globs['groups']:
            for label in globs['groups'][group]:
                if tip == label:
                    num_in_groups += 1;
        
        if num_in_groups != 1:
            print("parsed species groups:")
            print("targets: ", globs['groups']['targets']);
            print("outgroup: ", globs['groups']['outgroup']);
            print("conserved: ", globs['groups']['conserved']);
            CORE.errorOut("IO3", "The following tip label does not appear once in only one group: " + tip + " (" + str(num_in_groups) + "). If you provided internal node labels make sure there aren't duplicate labels within that clade provided elsewhere.", globs);
    # Go through every tip label and make sure it appears once and only once in the specified groups

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");

    return globs;

#############################################################################

def writeCF(globs):

    step = "Writing: " + globs['scfstatsfile'];
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    scf_headers = ["ID", "sCF", "sCF_N", "sDF1", "sDF1_N", "sDF2", "sDF2_N", "sN", "qN" ];
    cf_headers = scf_headers + [ "label", "length" ];

    scf_labels = {};
    #print();

    globs['scfstatsfile'] = os.path.join(globs['outdir'], globs['scfstatsfile']);
    # Add the output directory to the scfstatsfile name

    with open(globs['scfstatsfile'], "w") as outfile:
        outfile.write("# Site concordance factor statistics\n");
        outfile.write("# ID:        Branch ID in the species tree\n");

        outfile.write("# sCF:       Site concordance factor averaged over qN quartets (=sCF_N/sN)\n");
        outfile.write("# sCF_N:     sCF in absolute number of sites\n");
        outfile.write("# sDF1:      Site discordance factor for alternative quartet 1 (=sDF1_N/sN)\n");
        outfile.write("# sDF1_N:    sDF1 in absolute number of sites\n");
        outfile.write("# sDF2:      Site discordance factor for alternative quartet 2 (=sDF2_N/sN)\n");
        outfile.write("# sDF2_N:    sDF2 in absolute number of sites\n");
        outfile.write("# sN:        Number of decisive sites averaged over 100 quartets\n");
        outfile.write("# qN:        Number of quartets sampled on this branch\n");

        outfile.write("# label:     The label in the original species tree string\n");
        outfile.write("# length:    The branch length in the original species tree string\n");

        ####################

        outfile.write("\t".join(cf_headers) + "\n");
        for node in globs['st'].internals:
        # for node in globs['internals']:
            outline = [];

            if node in globs['scf']:

                outline.append(round(globs['scf'][node]['scf'], 3));
                outline.append(round(globs['scf'][node]['concordant-sites'], 2));

                outline.append(round(globs['scf'][node]['sdf1'], 3));
                outline.append(round(globs['scf'][node]['disco1-sites'], 2));

                outline.append(round(globs['scf'][node]['sdf2'], 3));
                outline.append(round(globs['scf'][node]['disco2-sites'], 2));                

                outline.append(round(globs['scf'][node]['decisive-sites'], 2));

                scf_labels[node] = str(round(globs['scf'][node]['scf'], 3));

                outline.append(globs['scf'][node]['total-quartets']);

            else:
                outline += ["NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"];
                scf_labels[node] = "NA";

            outline = [node] + [ str(v) for v in outline ] + [globs['st'].label[node], globs['st'].bl[node]];
            # outline = [node] + [ str(v) for v in outline ] + [ globs['st'][node][0], globs['st'][node][3] ];
            outfile.write("\t".join(outline) + "\n");
            #print("\t".join(outline));


    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: sCF stats written");
    globs['scf-stats-written'] = True;        

    ####################

    step = "Writing: " + globs['scftreefile'];
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    globs['scftreefile'] = os.path.join(globs['outdir'], globs['scftreefile']);
    # Add the output directory to the scftreefile name

    new_labels = { node : scf_labels[node] for node in globs['scf'] };
    globs['scf-labeled-tree'] = globs['st'].addLabel(new_labels, delim="_");

    # globs['scf-labeled-tree'] = globs['labeled-tree'];
    # # Get the labeled input tree from treeParse
    # globs['scf-labeled-tree'] = TREE.addBranchLength(globs['scf-labeled-tree'], globs['st']);
    # # Add the branch lengths back onto the tree
    # for node in globs['scf']:
    #    globs['scf-labeled-tree'] = globs['scf-labeled-tree'].replace(node, node + "_" + str(round(globs['scf'][node]['avg-quartet-scf'], 2)));
    # # For every node in the tree, add the averaged scf value over all loci to the label

    with open(globs['scftreefile'], "w") as cf_tree_file:
        cf_tree_file.write(globs['scf-labeled-tree']);

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: sCF tree written");
    globs['scf-tree-written'] = True;  
    # Status update

    return globs;

#############################################################################