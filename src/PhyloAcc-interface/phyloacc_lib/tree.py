#############################################################################
# Functions to read and handle phylogenetic trees
# Gregg Thomas
#############################################################################

import os
import re
import math
import random
import itertools
import phyloacc_lib.core as PC
import multiprocessing as mp
from collections import Counter

#############################################################################
def treeParse(tree, debug=False):
# The treeParse function takes as input a rooted phylogenetic tree with or without branch lengths and labels and returns the tree with nodes
# labeled in order for a post-order traversal and a dictionary with usable info about the tree in the following format:
# treeParse node label : [ branch length, ancestral node, node type, node label ]
#
# If branch length or node label is not present, they will have "NA" as the value
# node type can be one of: 'tip', 'internal', 'root'
# treeParse node label is in the format <N>, with N being a positive integer

    tree = tree.strip();
    if tree[-1] != ";":
        tree += ";";
    # Some string handling to remove any extra lines in the tree string and add the semi-colon if not present

    nodes, bl, supports, ancs = {}, {}, {}, {};
    # Initialization of all the tracker dicts

    topology = remBranchLength(tree);

    if debug:
        print("TOPOLOGY:", topology);

    nodes = {};
    for n in topology.replace("(","").replace(")","").replace(";","").split(","):
        nodes[n] = 'tip';
    # nodes = { n : 'tip' for n in topology.replace("(","").replace(")","").replace(";","").split(",") };
    # Retrieval of the tip labels

    if debug:
        print("NODES:", nodes);

    new_tree = "";
    z = 0;
    numnodes = 1;
    while z < (len(tree)-1):
        new_tree += tree[z];
        if tree[z] == ")":
            node_label = "<" + str(numnodes) + ">";
            new_tree += node_label;
            nodes[node_label] = 'internal';
            numnodes += 1;
        z += 1;
    nodes[node_label] = 'root';
    rootnode = node_label;
    # This labels the original tree as new_tree and stores the nodes and their types in the nodes dict

    if debug:
        print("NEW TREE:", new_tree);
        print("TREE:", tree);
        print("NODES:", nodes);
        print("ROOTNODE:", rootnode);
        #sys.exit();
    topo = "";
    z = 0;
    numnodes = 1;
    while z < (len(topology)-1):
        topo += topology[z];
        if topology[z] == ")":
            node_label = "<" + str(numnodes) + ">";
            topo += node_label;
            numnodes += 1;
        z += 1;
    # This labels the topology with the same internal labels

    if debug:
        print("TOPO:", topo);
        print("----------");
        print("TOPOLOGY:", topo);

    for node in nodes:
        if node + node in new_tree:
            new_tree = new_tree.replace(node + node, node);

    for node in nodes:
    # One loop through the nodes to retrieve all other info
        if debug == 1:
            print("NODE:", node);

        if nodes[node] == 'tip':
            supports[node] = "NA";
            if node + ":" in tree:
                cur_bl = re.findall(node + ":[\d.Ee-]+", new_tree);
                cur_bl = cur_bl[0].replace(node + ":", "");
                if debug == 1:
                    print("FOUND BL:", cur_bl);
                bl[node] = cur_bl;                
            else:
                bl[node] = "NA";

        elif nodes[node] == 'internal':
            if node + node in new_tree:
                new_tree = new_tree.replace(node + node, node);

            if node + "(" in new_tree or node + "," in new_tree or node + ")" in new_tree:
                if debug == 1:
                    print("NO BL OR LABEL");
                supports[node] = "NA";
                bl[node] = "NA";

            elif node + ":" in new_tree:
                supports[node] = "NA";
                cur_bl = re.findall(node + ":[\d.Ee-]+", new_tree);
                cur_bl = cur_bl[0].replace(node + ":", "");
                if debug == 1:
                    print("FOUND BL:", cur_bl);
                bl[node] = cur_bl;                                

            else:
                cur_bsl = re.findall(node + "[\d\w<>_*+.Ee/-]+:[\d.Ee-]+", new_tree);
                if cur_bsl:
                # If the pattern above is found then the node has both support and branch length
                    cur_bs = cur_bsl[0].replace(node, "");
                    cur_bs = cur_bs[:cur_bs.index(":")];
                    cur_bl = cur_bsl[0].replace(node, "").replace(cur_bs, "").replace(":", "");
                    if debug == 1:
                        print("FOUND BL AND LABEL:", cur_bl, cur_bs);
                    supports[node] = cur_bs;
                    bl[node] = cur_bl;
                    #new_tree = new_tree.replace(cur_bs, "");
                else:
                # If it is not found then the branch only has a label
                    cur_bs = re.findall(node + "[\w*+.<> -]+", new_tree);
                    cur_bs = cur_bs[0].replace(node, "");
                    if debug == 1:
                        print("FOUND LABEL:", cur_bs);
                    supports[node] = cur_bs;
                    bl[node] = "NA";
                    #new_tree = new_tree.replace(cur_bs, "");

        elif nodes[node] == 'root':
            bl[node] = "NA";
            supports[node] = new_tree[new_tree.index(node)+len(node):];
            ancs[node] = "NA";
            continue;

        # Next we get the ancestral nodes. If the node is the root this is set to NA.
        anc_match = re.findall('[(),]' + node, new_tree);

        #if nodes[node] == 'internal':
        #    sys.exit();
        # anc_match = re.findall(node + '[\d:(),]+', new_tree);
        anc_match = re.findall(node, topo);
        if debug:
            print("ANC MATCH:", anc_match);

        anc_tree = new_tree[new_tree.index(anc_match[0]):][1:];
        # Ancestral labels are always to the right of the node label in the text of the tree, so we start our scan from the node label

        if debug:
            print("NODE:", node);
            print("ANC_MATCH:", anc_match);
            print("ANC_TREE:", anc_tree);
            
        cpar_count = 0;
        cpar_need = 1;

        for i in range(len(anc_tree)):
        # We find the ancestral label by finding the ) which matches the nesting of the number of ('s found
            if anc_tree[i] == "(":
                cpar_need = cpar_need + 1;
            if anc_tree[i] == ")" and cpar_need != cpar_count:
                cpar_count = cpar_count + 1;
            if anc_tree[i] == ")" and cpar_need == cpar_count:
                anc_tree = anc_tree[i+1:];
                ancs[node] = anc_tree[:anc_tree.index(">")+1];
                break;

        if debug:
            print("FOUND ANC:", ancs[node]);
            print("---");
    nofo = {};
    for node in nodes:
        nofo[node] = [bl[node], ancs[node], nodes[node], supports[node]];
    # Now we just restructure everything to the old format for legacy support

    if debug:
    # Debugging options to print things out
        print(("\ntree:\n" + tree + "\n"));
        print(("new_tree:\n" + new_tree + "\n"));
        print(("topology:\n" + topo + "\n"));
        print("nodes:");
        print(nodes);
        print()
        print("bl:");
        print(bl);
        print()
        print("supports:");
        print(supports);
        print()
        print("ancs:");
        print(ancs);
        print()
        print("-----------------------------------");
        print()
        print("nofo:");
        print(nofo);
        print()

    return nofo, topo, rootnode;

#############################################################################

def remBranchLength(treestring):
# Removes branch lengths from a tree.

    treestring = re.sub('[)][\d\w<>/.eE_:-]+', ')', treestring);
    treestring = re.sub(':[\d.eE-]+', '', treestring);

    return treestring;

#############################################################################

def addBranchLength(tree, treedict, no_label=False, keep_tp_label=False):
# Re-writes the branch lengths onto a topology parsed by treeParse.
	for node in treedict:
		if treedict[node][2] == 'root':
			continue;
		if treedict[node][0] != "NA":
			tree = tree.replace(node, node + ":" + treedict[node][0]);
			if treedict[node][3] != "NA" and not no_label:
				tree = tree.replace(node + ":" + treedict[node][0], node + "_" + treedict[node][3] + ":" + treedict[node][0]);
		elif treedict[node][3] != "NA" and not no_label:
			tree = tree.replace(node, node + "_" + treedict[node][3]);

	if no_label and not keep_tp_label:
		tree = re.sub("<[\d]+>", "", tree);

	return tree;

#############################################################################

def getDesc(d_node, d_treedict):
# This function takes a species in the current tree and the dictionary of the current tree
# returned by treeparse and finds the direct descendants of the species.
    d_list = [];
    for node in d_treedict:
        if d_treedict[node][1] == d_node:
            d_list.append(node);

    if d_list == []:
        return [d_node];
    else:
        return d_list;

#############################################################################

def getClade(c_node, c_treedict):
# This function takes a species in the current tree and the dictionary of the current tree
# returned by treeparse and finds all tip labels that are descendants of the current node.
# This is done by getting the direct descendants of the current node with getDesc and then
# recursively calling itself on those descendants.
    clade = [];
    c_desc = getDesc(c_node, c_treedict);
    for d in c_desc:
        if c_treedict[d][2] != 'tip':
            clade.append(getClade(d, c_treedict));
        else:
            clade.append(d);

    r_clade = [];
    for c in clade:
        if type(c) == list:
            for cc in c:
                r_clade.append(cc);
        else:
            r_clade.append(c);

    return r_clade;

#############################################################################

def sampleQuartets(globs, root_desc):

    num_quartets = 100;
    # The number of quartets to sample around each branch
    ## ADD AS INPUT OPTION

    for node in globs['tree-dict']:
        if node in globs['tree-tips'] or node in root_desc or node == globs['root-node']:
            continue;
        # Cannot calculate sCF for tips, the root, or node descendant from the root

        desc = getDesc(node, globs['tree-dict']);
        anc_desc = getDesc(globs['tree-dict'][node][1], globs['tree-dict']);
        sister_node = [ n for n in anc_desc if n != node ][0];
        # For sCF, we treat the species tree as unrooted, so for each node(*)/branch, the possible clades
        # to sample quartets from are the two clades directly descendant from the node, the clade
        # descendant from the sister node, and all other species.
        #
        #         /\
        #        /  \
        #       /\   \
        #      /  \   \
        #     *    \   \
        #    /\     \   \
        #  d1  d2   sis  other
        # ROOTED TREE
        #
        #   sister           descendant 1 (left)
        #         \         /
        #           -------*
        #         /         \
        #    other           descendant 2 (right)
        # UNROOTED TREE
        #

        cur_clades = {};
        # Keeps track of the current clades from each of the 4 relevant branches

        cur_clades['left'] = getClade(desc[0], globs['tree-dict']);
        cur_clades['right'] = getClade(desc[1], globs['tree-dict']);
        cur_clades['sister'] = getClade(sister_node, globs['tree-dict']);
        # Gets clades for the left, right, and sister descendants

        all_clade_spec = cur_clades['left'] + cur_clades['right'] + cur_clades['sister'];
        cur_clades['other'] = [ tip for tip in globs['tree-tips'] if tip not in all_clade_spec ];
        # All other tips make up the 'other' clade

        assert all(len(cur_clades[c]) > 0 for c in cur_clades), \
            " * ERROR: quartet sampling failed for node: " + node + "\n" + \
            "\tleft:   " + len(cur_clades['left']) + "\n" + \
            "\tright:  " + len(cur_clades['right']) + "\n" + \
            "\tsister: " + len(cur_clades['sister']) + "\n" + \
            "\tother:  " + len(cur_clades['other']) + "\n"
        # Make sure each clade list has species or throw an error.. this shouldn't happen

        split1_pairs = list(itertools.product(cur_clades['left'], cur_clades['right']));
        split2_pairs = list(itertools.product(cur_clades['sister'], cur_clades['other']));
        # Get every possible pair of species from each clade side of the split
        # i.e. all pairs of left-right species (split1) and all pairs of sister-other species (split2)

        quartets = list(itertools.product(split1_pairs, split2_pairs));
        # Get all pairs from the pairs of species in split1 and split2 for all possible quartets at this node

        cur_num_quartets = len(quartets);
        # Count the total number of quartets at this node

        if cur_num_quartets > num_quartets:
            random.shuffle(quartets);
            quartets = quartets[:num_quartets];
        # If there are more quartets on the current node than the number to sample, sub-sample here
        # Otherwise, use all quartets
        ## SET A SEED FOR REPRODUCIBILITY

        globs['quartets'][node] = quartets;
        # Add the current set of quartets to the global dict

    return globs;

#############################################################################

def locusSCF(locus_item):
    locus, aln, quartets, tree_dict, skip_chars = locus_item
    # Unpack the data for the current locus

    aln_len = len(aln[list(aln.keys())[0]]);
    # Get the alignment length from the first sequence

    locus_scf = {};
    # The dictionary to calculate average sCF across all nodes in the current locus
    # <locus id> : { <node id> : <scf sum>, <num quartets>, <avg scf> }

    quartet_scores = {};
    # Dictionary to hold all counts for each quartet sampled

    node_scf = [];
    # The list of sCFs in this locus

    for node in tree_dict:
    # Calculate sCF for every eligible node in the tree

        if node not in quartets:
            continue;
        # Cannot calculate sCF for tips, the root, or node descendant from the root

        locus_scf[node] = { 'scf-sum' : 0, 'num-quartets' : 0, 'avg-scf' : "NA" };
        # Initialize the scf dict for the current node

        quartet_scores[node] = {};
        # Dictionary to hold all counts for each quartet sampled

        for quartet in quartets[node]:
        # quartet is a tuple of tuples:
        # ((split1-spec1, split1-spec2), (split2-spec1, split2-spec2))

            quartet_scores[node][quartet] = { 'variable-sites' : 0, 'decisive-sites' : 0, 'concordant-sites' : 0, 'scf' : "NA" };
            # Initialize the current scf dict for the current quartet
            
            #print(quartet);

            split1 = quartet[0];
            split2 = quartet[1];
            quartet_spec_list = split1 + split2;
            # Parse out the species for easier access

            # print(split1);
            # print(split2);
            # print(quartet_spec_list);
            #print(node, desc, cur_quartet['split1'], cur_quartet['split2']);

            sub_aln = {};
            for spec in quartet_spec_list:
                sub_aln[spec] = aln[spec];
            # Get the sub-alignment for the current quartet

            #print(sub_aln);

            for j in range(aln_len):
            # Count every site in the sub-alignment

                site, split1_site, split2_site = "", "", "";
                # The full site, the site for split1 and the site for split2

                for seq in sub_aln:
                    site += sub_aln[seq][j];
                    if seq in split1:
                        split1_site += sub_aln[seq][j];
                    if seq in split2:
                        split2_site += sub_aln[seq][j];
                # Get the allele for every sequence in the sub-alignment and add it to the appropriate sites                    

                if any(state in site for state in skip_chars):
                    continue;
                # We only care about sites with full information for the quartet, so skip others here

                uniq_site = set(site);
                # The alleles at the current site

                #print(site);
                #print(uniq_site);
                #print(len(uniq_site));

                if len(uniq_site) > 1:
                    quartet_scores[node][quartet]['variable-sites'] += 1;
                # If there is more than one allele, the site is variable

                    if all(site.count(allele) == 2 for allele in uniq_site):
                        quartet_scores[node][quartet]['decisive-sites'] += 1;
                    # If all alleles have a count of 2 at this site, it is decisive

                        #print(j, site);

                        if split1_site[0] == split1_site[1] and split2_site[0] == split2_site[1] and split1_site[0] != split2_site[0]:
                        #if split1_site.count(split1_site[0]) == len(split1_site) and split2_site.count(split2_site[0]) == len(split2_site):
                            quartet_scores[node][quartet]['concordant-sites'] += 1;
                        # If the alleles from split1 match and the alleles from split2 match, the site is concordant
                        
                        ## TODO: COUNT OTHER SITE PATTERNS HERE

                        #print(split1_site, split2_site);
            ## End site loop

            if quartet_scores[node][quartet]['decisive-sites'] != 0:
                quartet_scores[node][quartet]['scf'] = quartet_scores[node][quartet]['concordant-sites'] / quartet_scores[node][quartet]['decisive-sites'];
            # Calculate the scf for the current quartet
        ## End quartet loop

        for quartet in quartet_scores[node]:
            if quartet_scores[node][quartet]['scf'] != "NA":
                locus_scf[node]['scf-sum'] += quartet_scores[node][quartet]['scf'];
                locus_scf[node]['num-quartets'] += 1;
        # For all the quartets with sCF calculated, sum the sCF for this locus to be averaged later

        if locus_scf[node]['num-quartets'] != 0:
            locus_scf[node]['avg-scf'] = locus_scf[node]['scf-sum'] / locus_scf[node]['num-quartets'];
            node_scf.append(locus_scf[node]['scf-sum'] / locus_scf[node]['num-quartets']);
        # Average all sCF from each quartet for this locus
    ## End node loop

    return node_scf, locus, quartet_scores;

#############################################################################

def scf(globs):
# A function to calculate site concordance factors for each input locus
# sCFs are calculate in two ways:
# 1. Average sCF across all nodes per LOCUS (written to aln stats file)
# 2. Average sCF across all loci per NODE (written to scf stats file)
# In both cases, a number of quartets (100) are sampled for each branch in each tree and sites are counted
# and averages are taken across quartets. See: https://doi.org/10.1093/molbev/msaa106

    root_desc = getDesc(globs['root-node'], globs['tree-dict']);
    # Get the descendants of the root node to exclude them from sCF calculations

    step = "Sampling quartets";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    globs = sampleQuartets(globs, root_desc);
    step_start_time = PC.report_step(globs, step, step_start_time, "Success");
    # Sample quartets for all nodes

    step = "Calculating per-locus sCF";
    step_start_time = PC.report_step(globs, step, False, "Processed 0 / " + str(globs['num-loci']) + " loci...", full_update=True);
    # Status update

    for node in globs['tree-dict']:
        if node in globs['tree-tips'] or node in root_desc or node == globs['root-node']:
            continue;
        # Cannot calculate sCF for tips, the root, or node descendant from the root

        globs['scf'][node] = { 'variable-sites' : 0, 'decisive-sites' : 0, 'concordant-sites' : 0, 'quartet-scf-sum' : 0,
                                    'total-quartets' : 0, 'avg-quartet-scf' : "NA" };
        # Initialize the dictionary to calculate average sCF per node across all loci

    if globs['qstats']:
        qstats_file = os.path.join(globs['outdir'], "quartet-stats.csv");
        qfile = open(qstats_file, "w");
        headers = ["locus","node","quartet","variable-sites","decisive-sites","concordant-sites"];
        qfile.write(",".join(headers) + "\n");
    # For --qstats, creates and opens a file to write site counts for each quartet within the pool loop

    with globs['scf-pool'] as pool:
        counter = 0;
        # A counter to keep track of how many loci have been completed
        for result in pool.imap_unordered(locusSCF, ((locus, globs['alns'][locus], globs['quartets'], globs['tree-dict'], globs['skip-chars']) for locus in globs['alns'])):
            # Loop over every locus in parallel to calculate sCF per node

            node_scf, locus, quartet_scores = result;
            # Unpack the current result

            globs['aln-stats'][locus]['node-scf-sum'] = sum([ n for n in node_scf ]);
            globs['aln-stats'][locus]['num-nodes'] = len(node_scf);
            globs['aln-stats'][locus]['node-scf-avg'] = "NA";
            globs['aln-stats'][locus]['perc-low-scf-nodes'] = "NA";            
            globs['aln-stats'][locus]['low-scf-nodes'] = len([ n for n in node_scf if n < globs['min-scf'] ]);
            if globs['aln-stats'][locus]['num-nodes'] != 0:
                globs['aln-stats'][locus]['node-scf-avg'] = globs['aln-stats'][locus]['node-scf-sum'] / globs['aln-stats'][locus]['num-nodes'];
                globs['aln-stats'][locus]['perc-low-scf-nodes'] = globs['aln-stats'][locus]['low-scf-nodes'] / globs['aln-stats'][locus]['num-nodes'];
            # Average sCF per locus stored with the other alignment stats

            if globs['run-mode'] == 'adaptive':
                if globs['aln-stats'][locus]['low-scf-nodes'] > math.floor(globs['aln-stats'][locus]['num-nodes'] / 3):
                    globs['aln-stats'][locus]['batch-type'] = "gt";
                    globs['gt-loci'] += 1;
                else:
                    globs['aln-stats'][locus]['batch-type'] = "st";
                    globs['st-loci'] += 1;
            elif globs['run-mode'] == 'st':
                globs['st-loci'] += 1;
            elif globs['run-mode'] == 'gt':
                globs['gt-loci'] += 1;
            # If the run mode is adaptive, partition loci with more than 1/3 of nodes with low sCF to the gene tree model and all others
            # to the species tree model
            ## Partition the sequences based on scf

            for node in quartet_scores:
                for q in quartet_scores[node]:

                    if globs['qstats']:
                        q_str = ";".join(q[0]) + ";" + ";".join(q[1]);
                        outline = [locus, node, q_str] + [str(quartet_scores[node][q][col]) for col in headers[3:]];
                        qfile.write(",".join(outline) + "\n");
                    # For --qstats, writes the quartet stats to a file for the current quartet

                    globs['scf'][node]['variable-sites'] += quartet_scores[node][q]['variable-sites'];
                    globs['scf'][node]['decisive-sites'] += quartet_scores[node][q]['decisive-sites'];
                    globs['scf'][node]['concordant-sites'] += quartet_scores[node][q]['concordant-sites'];
                    # Sum sites for this quartet for this node

                    if quartet_scores[node][q]['decisive-sites'] != 0:
                        globs['scf'][node]['quartet-scf-sum'] += quartet_scores[node][q]['concordant-sites'] / quartet_scores[node][q]['decisive-sites'];
                        globs['scf'][node]['total-quartets'] += 1;
                    # If there are decisive sites, calculate sCF for this quartet
            # For every quartet in every node, calculate sCF and sum up the number of sites in order to average
            # across nodes later

            counter += 1;
            if counter % 100 == 0:
                cur_scf_time = PC.report_step(globs, step, step_start_time, "Processed " + str(counter) + " / " + str(globs['num-loci']) + " loci...", full_update=True);
            # A counter and a status update every 100 loci
        ## End imap locus loop
    ## End pool

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + str(globs['st-loci'] ) + " st, " + str(globs['gt-loci'] ) + " gt loci.", full_update=True);
    # Status update
    
    step = "Averaging sCF per node";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    for node in globs['scf']:
        globs['scf'][node]['avg-quartet-scf'] = globs['scf'][node]['quartet-scf-sum'] / globs['scf'][node]['total-quartets'];
    step_start_time = PC.report_step(globs, step, step_start_time, "Success");
    # For every node, average the sCFs per quartet across all loci

    if globs['qstats']:
        qfile.close()
    # For --qstats, closes the quartet stats file.

    return globs;
    
#############################################################################
