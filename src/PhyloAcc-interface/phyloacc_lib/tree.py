#############################################################################
# Functions to read and handle phylogenetic trees
# Gregg Thomas
#############################################################################

import sys
import re
import copy
import random
import itertools

#############################################################################

class Tree:
# The treeParse function takes as input a rooted phylogenetic tree with branch lengths and returns the tree with node labels and a
# dictionary with usable info about the tree in the following format:
# node:[branch length (if present), ancestral node, node type, node label (if present)]

    def __init__(self, tree_string, get_subtrees=False, debug=False):

        tree = tree_string.strip();
        if tree[-1] != ";":
            tree += ";";
        # Some string handling

        tree = re.sub("#[\d.]+", "", tree);

        self.orig_tree_str = tree;
        self.topo_str = "";
        self.labeled_topo_str = "";
        self.tree_str = "";
        # Tree string types

        self.type = {};
        self.bl = {};
        self.has_bl = "";
        self.label = {};
        self.has_label = "";
        self.anc = {};
        self.desc = {};
        self.sis = {};
        self.subtree = {};
        subtrees = {};
        # Node attributes

        self.nodes = [];
        self.num_nodes = 0;

        self.tips = [];
        self.num_tips = 0;

        self.internals = [];
        self.num_internals = 0;

        self.num_polytomies = 0;

        self.rooted = "NA";
        self.root = "";
        # Node lists and counts

        ## Class attributes
        #####

        self.topo_str = re.sub('[)][\d\w<>/.eE_:-]+', ')', tree);
        self.topo_str = re.sub(':[\d.eE-]+', '', self.topo_str);
        ## Remove the branch lengths and node labels from the input tree string
        #####     

        for tip_label in self.topo_str.replace("(","").replace(")","").replace(";","").split(","):
            self.tips.append(tip_label);
            self.type[tip_label] = 'tip';
        ## Retrieval of the tip labels
        #####
       
        if debug:
            print();
            print("ORIGINAL:", self.orig_tree_str);
            print("TOPOLOGY:", self.topo_str);
        ## Node counting and root checking
        #####

        # self.bl = { node : "NA" for node in self.nodes };
        # self.label = { node : "" for node in self.internals + [self.root] };
        # self.anc = { node : "NA" for node in self.nodes };
        ## With the node lists set, initialize the other node attributes
        #####

        # z = 0;
        # numnodes = 1;
        # while z < (len(self.topo_str)-1):
        #     self.labeled_topo_str += self.topo_str[z];
        #     if self.topo_str[z] == ")":
        #         node_label = "<" + str(numnodes) + ">";
        #         self.labeled_topo_str += node_label;
        #         numnodes += 1;
        #     z += 1;
        # if debug:
        #     print("LABELED TOPOLOGY:", self.labeled_topo_str);
        #     print("----------");
        ## Add the generated internal node labels onto the topology string
        #####

        subtrees = {};

        pairs = re.findall("\([\d\w/.,:_<>-]+\)", tree);

        node_count = 1;
        while pairs:
        # Loop over all pairs of opening an closing parens at each
        # tree level

            for pair in pairs:
            # Loop over every pair in the current set

                anc_label = "<" + str(node_count) + ">";
                self.desc[anc_label] = [];
                node_count += 1;
                # The ancestral label

                node_list = pair.replace("(", "").replace(")", "").split(",");
                num_nodes = len(node_list);
                if num_nodes > 2:
                    self.num_polytomies += 1;
                # In most cases, the number of branches in a pair will be 2, but
                # for polytomies it will be more

                cur_nodes = [];
                # A list of nodes in the current pair to add
                # as descendants for the ancestral node

                for node in node_list:
                # Loop over each node in the current pair list

                    if ":" in node:
                        node, bl = node.split(":");
                    else:
                        bl = "NA"
                    # Check if the current node has a branch length

                    if node in self.tips:
                        self.label[node] = "NA";
                        self.bl[node] = bl;
                        # Tip info

                        self.desc[node] = "NA";
                        # Tips have no descendants

                        if bl != "NA":
                            subtrees[node] = node + ":" + str(bl);
                        else:
                            subtrees[node] = node;
                        # For tips, the subtree is just that node and its branch length
                        # If the current node is a tip label, add no label
                    ## Tips

                    else:
                        cur_label = re.sub("<[\d]+>", "", node, 1);
                        node = re.findall("<[\d]+>", node)[0];
                        # Parse out the provided labels in the tree from
                        # the ones generated here for all internal nodes

                        self.type[node] = "internal";
                        self.internals.append(node);
                        self.label[node] = cur_label;
                        self.bl[node] = bl;
                        # Node information

                        subtrees[node] = "(";
                        for d in self.desc[node]:
                            subtrees[node] += subtrees[d] + ",";
                        subtrees[node] = subtrees[node][:-1] + ")" + node;

                        if cur_label:
                            subtrees[node] += str(cur_label);
                        if bl != "NA":
                            subtrees[node] += ":" + str(bl);
                        # For internal nodes, the subtree is the combination of the subtrees of its
                        # descendants
                        # Add branch lengths and labels as needed
                    ## Internal nodes

                    cur_nodes.append(node);
                    self.nodes.append(node);
                    self.anc[node] = anc_label;
                    # Add the ancestral label for the current node
                ## End node loop

                self.desc[anc_label] = cur_nodes;
                # Set the descendants of the current ancestor as all nodes in the current pair

                for node in cur_nodes:
                    self.sis[node] = [ n for n in cur_nodes if n != node ];
                # Set the sister branches for each node in the current pair

                tree = tree.replace(pair, anc_label);
                # Replace the current pair with the ancestral label in the tree
            ## End pairs loop

            pairs = re.findall("\([\d\w/.,:_<>-]+\)", tree);
            # Find all new pairs in the tree
        # If there are no pairs left, we've reached the root of the tree
        ## End tree level loop

        ## Main tree parsing loop
        #####

        self.nodes.append(anc_label);
        self.internals.append(anc_label);
        self.type[anc_label] = "internal";
        self.root = anc_label;
        if tree.count(anc_label) == 1:
            orig_root_label = tree.replace(anc_label, "").replace(";", "");
            if orig_root_label == "":
                orig_root_label = "root";
            self.label[anc_label] = orig_root_label;
        else:
            self.label[anc_label] = anc_label;
        self.bl[anc_label] = "NA";
        self.anc[anc_label] = "NA";
        self.sis[anc_label] = "NA";
        # Add the final anc label as the root

        if num_nodes > 2:
            self.rooted = False;
        else:
            self.rooted = True;
        # Check if there is a polytomy at the root for rootedness

        subtrees[anc_label] = "(";
        for d in self.desc[anc_label]:
            subtrees[anc_label] += subtrees[d] + ",";
        subtrees[anc_label] = subtrees[anc_label][:-1] + ")" + anc_label;
        self.tree_str = subtrees[self.root];
        # Add the subtree for the root, which ends up being the fully labeled tree

        if debug:
            print();
            print("NODES:", self.nodes);
            print("ROOTED:", self.rooted);
            print("ROOTNODE:", self.root);
        ## Root block
        #####

        if get_subtrees:
            self.subtrees = subtrees;
        # Save the subtrees if flagged

        if all(self.bl[n] == "NA" for n in self.bl):
            self.has_bl = False;
        else:
            self.has_bl = True;
            self.bl = { n : self.bl[n] if self.bl[n] != "NA" else "0.0" for n in self.bl };
        ## Checks for branch lengths

        if any(self.label[n] in ["", "NA", ":0", "root"] for n in self.label if n in self.internals):
            self.has_label = False;
            self.label = { n : "NA" for n in self.label };
            # NOTE: does PhyloAcc need the :0 branch length at the root?
        else:
            self.has_label = True;
            self.label = { n : self.label[n] if self.label[n] != "" else "NA" for n in self.label };
        ## Checks for labels and bls
        #####

        self.num_tips = len(self.tips);
        self.num_internals = len(self.internals);
        self.num_nodes = len(self.nodes);
        ## Counts
        #####        

        if debug:
            print();
            self.showAttrib("label","anc","desc","sis");
            print("FINAL: " + self.tree_str);
            
            if get_subtrees:
                print();
                for node in self.subtrees:
                    print(node);
                    print(self.subtrees[node]);
                    print();

    #############################################################################

    ##########

    def checkRooted(self):
    # Checks if a tree object is rooted or not by counting the difference
    # in the number of tip and internal nodes

        if (self.num_internals + 1) != (self.num_tips - 1):
            return False;
        elif (self.num_internals + 1) == (self.num_tips - 1):
            return True;
        else:
            return -1;

    ##########

    def getDesc(self, node):
        # This function takes a node in the current tree object
        # and returns a list of the two direct descendant nodes of it.

        if node in self.tips:
            return [node];
        else:
            return [ n for n in self.nodes if self.anc[n] == node ];

    ##########

    def getSister(self, node):
        # This function takes a node in the current tree object
        # and returns the other direct descendant of its ancestral node

        if node == self.root:
            return "NA";
        anc_desc = self.getDesc(self.anc[node]);
        return [ n for n in anc_desc if n != node ];

    ##########

    def getClade(self, node, full=False):
    # This function takes a node in the current tree object
    # and finds all tip labels that are descendants of it.
    # This is done by getting the direct descendants of the node with getDesc and then
    # recursively calling itself on those descendants.

        clade = [];
        if node in self.tips:
            return [node];

        desc = self.desc[node];
        for d in desc:
            if d not in self.tips:
                clade += self.getClade(d, full);
                if full:
                    clade.append(d);
                # If full is true, the function will also return all internal nodes
                # descending from the given node
            else:
                clade.append(d);

        return clade;

    ##########

    def getClades(self, full=False):
    # Calls getClade() on every node in a tree object

        clades = {};
        for node in self.nodes:
            clades[node] = set(self.getClade(node, full=full));

        return clades;

    ##########

    def getSplit(self, node):
    # Returns the tips from a node that do not descend from it to define a split (with clade)

        if node == self.root:
            return "NA";

        return set(self.tips) - set(self.getClade(node));

    ##########

    def getSplits(self):
    # Calls getSplit on every node in a tree object

        splits = {};
        for node in self.nodes:
            splits[node] = self.getSplit(node);

        return splits;

    ##########

    def getQuartet(self, node):
    # Returns 4 sets of tips given an internal node:
    # 1: tips descending from one descendant node
    # 2: tips descending from the other descendant node
    # 3: tips descending from the sister node
    # 4: all other tips not in the previous 3 categories

        if self.rooted:
            if node == self.root or self.type[node] == "tip":
                return "NA";
            else:
                if self.anc[node] == self.root:
                    if self.type[self.sis[node][0]] == "tip":
                        return "NA";
                    # If the sister is a tip, quartets cannot be sampled from this branch    

                    d1 = self.desc[node][0];
                    d2 = self.desc[node][1];

                    q1 = set(self.getClade(d1));
                    q2 = set(self.getClade(d2));

                    sis_desc = self.desc[self.sis[node][0]];
                    q3 = set(self.getClade(sis_desc[0]));
                    q4 = set(self.getClade(sis_desc[1]));
                    # If the sister is not a tip, return both of its descendants as q3 and q4              
                ## For nodes that have the root as an ancestor, the sister node needs
                ## to be parsed differently

                else:
                    d1 = self.desc[node][0];
                    d2 = self.desc[node][1];
                    sis = self.sis[node][0];

                    q1 = set(self.getClade(d1));
                    q2 = set(self.getClade(d2));
                    q3 = set(self.getClade(sis));
                    q4 = set(self.tips) - q1 - q2 - q3;
                ## For most nodes, just look up the clades of the descendant and sister nodes
        ## For rooted trees

        else:
            if node == self.root or self.type[node] == "tip":
                return "NA";
            else:
                d1 = self.desc[node][0];
                d2 = self.desc[node][1];
                q1 = set(self.getClade(d1));
                q2 = set(self.getClade(d2));

                if len(self.desc[self.anc[node]]) > 2:
                    other = [ n for n in self.desc[self.anc[node]] if n != node ];
                    q3 = set(self.getClade(other[0]));
                    q4 = set(self.getClade(other[1]));    
                else:
                    q3 = set(self.getClade(self.sis[node][0]));
                    q4 = set(self.tips) - q1 - q2 - q3;
            ## If the node has more than 2 descendants, simply use the first two...
        ## For unrooted trees        

        return { 'd1' : q1, "d2" : q2, "s" : q3, "q4" : q4 };

    ##########

    def getQuartets(self, root=True):
    # Calls getQuartet on every internal node in a tree object

        quartets = {};
        if root: 
            nodes = self.internals;
        else:
            nodes = self.internals[:-1];
        for node in nodes:
            quartets[node] = self.getQuartet(node);

        return quartets;

    ##########

    def sampleQuartets(self, num_quartets=100):
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

        full_quartets = self.getQuartets();
        sampled_quartets = {};

        for node in self.internals:
            #if node in self.tips or node == self.root:
            if full_quartets[node] == "NA":
                continue;
            # Cannot calculate sCF for tips, the root, or node descendant from the root

            assert all(len(full_quartets[node][q]) > 0 for q in full_quartets[node]), \
                " * ERROR: quartet sampling failed for node: " + node + "\n" + \
                "\tleft:   " + len(full_quartets[node]['d1']) + "\n" + \
                "\tright:  " + len(full_quartets[node]['d2']) + "\n" + \
                "\tsister: " + len(full_quartets[node]['s']) + "\n" + \
                "\tother:  " + len(full_quartets[node]['q4']) + "\n"
            # Make sure each clade list has species or throw an error.. this shouldn't happen

            # if node == "<15>":
            #     print(full_quartets[node]['d1']);

            split1_pairs = list(itertools.product(full_quartets[node]['d1'], full_quartets[node]['d2']));
            split2_pairs = list(itertools.product(full_quartets[node]['s'], full_quartets[node]['q4']));
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

            sampled_quartets[node] = quartets;
            # Add the current set of quartets to the global dict

        return sampled_quartets;

    ##########

    def genSubtrees(self):
    # Generates sub-tree strings for every node in a tree object

        subtrees = {};
        # Subtree dict
        
        for node in self.tips:
            if self.has_bl:
                subtrees[node] = node + ":" + str(self.bl[node]);
            else:
                subtrees[node] = node;
        # For tips, the subtree is just that node and its branch length

        for node in self.internals + [self.root]:
            subtrees[node] = "(";
            for d in self.desc[node]:
                subtrees[node] += subtrees[d] + ",";
            subtrees[node] = subtrees[node][:-1] + ")" + node;

            if self.has_label:
                subtrees[node] += str(self.label[node]);
            if self.has_bl and node != self.root:
                subtrees[node] += ":" + str(self.bl[node]);
        # For internal nodes, the subtree is the combination of the subtrees of its
        # descendants

        return subtrees;

    ##########

    def findClades(self, tip_set):
    # A function that takes a set of tips and finds all the monophyletic
    # clades within that set, and returns the LCA of each clade found

        tip_set = set(tip_set);

        clade_nodes = [];
        for node in self.internals:
            if set(self.getClade(node)) <= tip_set:
                clade_nodes.append(node);
        # In the input tree, find all nodes that have clades that are
        # a subset of the given tip set

        nodes_to_rm = [];
        for node in clade_nodes:
            for node_check in clade_nodes:
                if node == node_check:
                    continue;

                if node in set(self.getClade(node_check, full=True)):
                #if node in self.desc[node_check]:
                    nodes_to_rm.append(node);
        clade_set = [ n for n in clade_nodes if n not in nodes_to_rm ];
        # We only want the deepest node from each possible clade, so remove
        # nodes from those found that are descendants of another node.

        nodes_to_rm = [];
        for node in tip_set:
            if any(node in set(self.getClade(n, full=True)) for n in clade_set):
                nodes_to_rm.append(node);
        clade_set += [ n for n in tip_set if n not in nodes_to_rm ];
        # From the original input set, remove any nodes that are now found
        # in a clade

        return set(clade_set);

    ##########

    def findSplits(self, tip_set, clades=False, splits=False):
    # A function that takes a set of tips and finds whether any branches
    # are defined on either side by them

        split_nodes = [];
        for node in self.nodes:
            if clades:
                clade = clades[node];
            else:
                clade = set(self.getClade(node));
            # Lookup the clade for the given node

            if splits:
                split = splits[node];
            else:
                split = self.getSplit(node);
            # Lookup the split for the given node

            if clade == tip_set or split == tip_set:
                split_nodes.append(node);
            # If the clade or the split matches the tip set, add it to the list
            # of splits

        return set(split_nodes);

    ##########

    def LCA(self, node_list):
    # Given a list of nodes, this function finds the least common ancestor of them,
    # and tells whether the nodes provided form a monophyletic clade.

        ancs = {};
        for node in node_list:
            ancs[node] = [node];
        # For each node in the input list, we make a list of the path from that node 
        # to the root of the tree, including that node

        for node in node_list:
            if node == self.root:
                continue;

            cur_anc = self.anc[node];
            ancs[node].append(cur_anc);

            while cur_anc != self.root:
                cur_anc = self.anc[cur_anc];
                ancs[node].append(cur_anc);
        # For each node, add every node between it and the root to its ancs list

        intersect_anc = set.intersection(*list(map(set, list(ancs.values()))));
        # Gets the intersect of all lists of paths to the root, unordered

        lcp = [t for t in list(ancs.values())[0] if t in intersect_anc];
        # Orders the nodes in the intersect of all paths based on their order in
        # the path of an arbitrary node (since it should be identical for all of them)

        return lcp[0];
        # Returns the first node in the common path as the LCA  

    ##########

    def Monophyletic(self, node_list):
    # Determines whether a set of nodes is within a monophyletic clade -- do
    # they all descend from the same ancestor with no other descendants?

        monophyletic = False;
        if set(self.getClade(self.LCA(node_list))) == set(node_list):
            monophyletic = True;
        return monophyletic;

    ##########

    def addBranchLength(self):
    # Re-writes the branch lengths onto a tree object topology

        tree = self.labeled_topo_str;

        for node in self.nodes:
            new_node_str = node;
            if self.has_label and node in self.internals:
                new_node_str += str(self.label[node]);

            if self.has_bl and node != self.root:
                new_node_str += ":" + str(self.bl[node]);

            tree = tree.replace(node, new_node_str);

        return tree;

    ##########

    def addLabel(self, label_dict, delim=""):
    # Given a dictionary with { node : label } format, adds those labels and the branch
    # lengths onto the given tree's topology
    # Returns: tree string

        new_tree_str = self.tree_str;
        # Extract labeled tree topology to add labels to

        for node in self.nodes:
            new_label = "";
            if self.type[node] == "tip":
                continue;
                #new_label = node;
            # If the node is a tip, always add the node as a label

            if node in label_dict and label_dict[node] != "NA":
                if delim and self.has_label:
                    new_label += self.label[node] + delim + str(label_dict[node]);
                else:
                    new_label += str(label_dict[node]);
            elif self.label[node] != "NA":
                new_label += self.label[node];
            # If the node is in the given dictionary, add the new label

            # if self.has_bl and node != self.root:
            #     new_label += ":" + str(self.bl[node]);
            # If the tree has branch lengths, add the bl

            if new_label or node == self.root:
                if self.label[node] != "NA":
                    new_tree_str = new_tree_str.replace(node + self.label[node], new_label);
                else:
                    new_tree_str = new_tree_str.replace(node, new_label);
            # If a new label was created, replace the old node label in the
            # tree with it

        return new_tree_str;

    ##########

    def labelDesc(self):
    # Adds labels in the format [descendant 1]-[descendant 2] to all internal nodes
    # of a tree
    # Descendants 1 and 2 are chosen by alphabetically sorting all tips that are 
    # children of the current node and selecting the first two tips

        label_dict = {};
        # A dictionary: [numeric label] : [desecendant label]

        tree_str = self.tree_str;
        for node in self.internals:
            selected_descs = sorted(self.getClade(node))[0:2];
            new_label = "-".join(selected_descs);
            label_dict[node] = new_label;
        # For every node, select 2 descendants as the new label
        
        tree_str = self.addLabel(label_dict);
        # Add the new labels to the tree string

        tree_str = re.sub('<[\d]+>', "", tree_str) + ";";
        # Remove the integer node labels and add a semi-colon

        if self.has_label:
            for node in self.internals:
                tree_str = tree_str.replace(self.label[node] + label_dict[node], label_dict[node]);
        # Hacky way to remove the old labels after adding the new ones... seems problematic...

        return tree_str;

    ##########

    def showAttrib(self, *attrib_list):
    # This function simply prints out the given tree attributes to the screen
    # Valid attributes are: "type", "desc", "anc", "sis", "clade", "split", "quartet"

        print("-" * 60);
        ## Seperator

        pad = 60;
        ## Width of columns

        valid_attribs = ["type", "length", "label", "desc", "anc", "sis", "clade", "split", "quartet"];
        ## List of valid attributes

        attrib_list = [ attrib for attrib in attrib_list if attrib in valid_attribs ];
        attrib_rm = [ attrib for attrib in attrib_list if attrib not in valid_attribs ];
        if attrib_rm:
            print("WARNING: The following are not valid Tree attributes to display: " + ",".join(attrib_rm));
        ## Parse passed attributes and warn if any are invalid

        headers = ["NODE"] + [ attrib.upper() for attrib in valid_attribs if attrib in attrib_list ];
        outline = [ spacedOut(header, pad) for header in headers ];
        print("".join(outline).strip());
        ## Display the attributes as headers

        header_len = (pad * len(attrib_list)) + len(attrib_list[-1]);
        print("-" * header_len);
        ## Seperator between headers and rows

        for node in self.nodes:
            outline = spacedOut(node, pad);
            for attrib in valid_attribs:
                if attrib in attrib_list:
                    if attrib == "type":
                        if node == self.root:
                            outline += spacedOut(self.type[node] + ",root", pad);
                        else:
                            outline += spacedOut(self.type[node], pad);
                    ## type

                    if attrib == "label":
                        if self.type[node] == "tip":
                            outline += spacedOut("", pad);
                        else:
                            outline += spacedOut(self.label[node], pad);
                    # label

                    if attrib == "length":
                        outline += spacedOut(self.bl[node], pad);
                    # length

                    if attrib == "desc":
                        if node in self.tips:
                            outline += spacedOut(self.desc[node], pad);
                        else:
                            outline += spacedOut(",".join(self.desc[node]), pad);
                    ## desc

                    if attrib == "anc":
                        outline += spacedOut(self.anc[node], pad);
                    ## anc

                    if attrib == "sis":
                        if node == self.root:
                            outline += spacedOut(self.sis[node], pad);
                        else:
                            outline += spacedOut(",".join(self.sis[node]), pad);
                    ## sis

                    if attrib == "clade":
                        clade = str(set(self.getClade(node)));
                        outline += spacedOut(clade, pad);
                        if len(clade) > pad:
                            outline += "\t";
                    ## clade
                        
                    if attrib == "split":
                        split = str(self.getSplit(node))
                        outline += spacedOut(split, pad);
                        if len(split) > pad:
                            outline += "\t";
                    ## split

                    if attrib == "quartet":
                        if self.type[node] == "tip" or node == self.root:
                            outline += "NA";
                        else:
                            #outline += spacedOut("", pad);
                            quartet = self.getQuartet(node);

                            outline += str(quartet["d1"]);

                            for key in ['d2', 's', 'q4']:
                                outline += "\t" + str(quartet[key]);
                    ## quartet

                ## End attrib loop
            ## End valid loop
        ## End node loop
            
            print(outline.strip());
        print("-" * 60);
        ## Seperator      

## END TREE CLASS
#############################################################################
## BEGIN TREE STRING FUNCTIONS

def categorizeBranches(globs, tree):
# For plotting, categorizes all branches by input state (rather than just tips)
    targets = [ tip for tip in globs['groups']['targets'] ];
    conserved = [ tip for tip in globs['groups']['conserved'] ];
    outgroups = [ tip for tip in globs['groups']['outgroup'] ];
    # Get the lists of tips in each category here
    # Internal nodes will be added to these lists, so we don't want to do that in globs

    for node in tree.internals:
    # Go over every internal node in the tree in a post-order traversal

        d1, d2, = tree.desc[node];
        # Get the descendants of the current node

        if d1 in outgroups or d2 in outgroups:
            outgroups.append(node);
        # If either of the descendants is an outgroup, this branch is also an outgroup
        elif d1 in targets and d2 in targets:
            targets.append(node);
        # If both of the descendants are targets, this branch is also a target
        else:
            conserved.append(node);
        # Otherwise, this is a conserved branch
    ## End internal node loop

    return targets, conserved, outgroups;

#############################################################################

def remBranchLength(tree_str):
# Removes branch lengths and labels from a tree string

    tree_str = re.sub('[)][\d\w<>/.eE_:-]+', ')', tree_str);
    tree_str = re.sub(':[\d.eE-]+', '', tree_str);

    return tree_str;

#############################################################################

def getSubtree(node, tree_str):
# Gets the subtree string at a given node from a labeled tree string from treeParse
# Much slower than genSubtrees method

    subtree = "";
    # Initialize the subtree string

    partree = tree_str[:tree_str.index(node)][::-1]
    # Slice the tree string at the index of the node label and reverse it

    cp = 0;
    op = 0;
    # Counts of closing an opening parentheses. The subtree will be complete when they
    # are equal

    for c in partree:
        if c == ")":
            cp = cp + 1;
        if c == "(":
            op = op + 1;
        subtree = subtree + c;
        if cp == op:
            break;
    # Loop through every character in the sliced and reversed tree, counting parentheses and adding
    # charcters to the subtree string one at a time. Stop when the number of closing parentheses equals
    # the number of opening

    return subtree[::-1];
    # Return the reverse of the subtree, which is the original orientation

#############################################################################

def spacedOut(string, totlen, sep=" "):
# Properly adds spaces to the end of a message to make it a given length
    spaces = sep * (totlen - len(string));
    return string + spaces;  

#############################################################################

def debugTree(globs):

    import phyloacc_lib.core as PC

    print("\n\n");
    
    tree_str = globs['tree-string']
    #tree_str = "((a,b),c,(d,e));" # Unrooted
    #tree_str = "(((a,b),c),(d,e));" # Rooted
    #tree_str = ("((a,b),c)")
    #tree_str = "(((a,b),c),(((d,e),f),h),(i,j),(k,l));";
    #tree_str = "((((a,b),c),(((d,e),f),h),i),j);";
    #tree_str = "(mspr:0.0018190933,((mwsb:0.0018060852,((mcas:0.0005908566,mpwk:0.0000009867):0.0000076480,mmus:0.0000009867):0.0000009867):0.0011628549,(((((((((((cgri:0.0152686696,((pcam:0.0012137705,psun:0.0006517732):0.0083093402,prob:0.0087755222):0.0112747594):0.0027806008,maur:0.0194531500):0.0098942322,moch:0.0330046438):0.0022180414,pman:0.0173033334):0.0096202861,((jjac:0.0617468174,itri:0.0714744806):0.0197552895,ngal:0.0425089249):0.0299773246):0.0088542882,mung:0.0286172457):0.0205036203,rnor:0.0311801045):0.0011873972,((rdil:0.0199545510,gdol:0.0117952834):0.0120430301,rsor:0.0202951611):0.0008975499):0.0043168585,(hall:0.0119687298,(pdel:0.0121964374,mnat:0.0120591970):0.0031186198):0.0070437051):0.0091272024,mpah:0.0130568501):0.0093423702,mcar:0.0073123519):0.0001651381):0.0030773597,mspi:0.0005888274);"
    #tree_str = "((jjac:0.1003165447,itri:0.1044182968)0.891:0.010544,(ngal:0.069324448,((pman:0.0406478328,(moch:0.0497980904,((cgri:0.0250402305,maur:0.0328120476)0.668:0.00352051,(prob:0.0144107446,(pcam:0.0028058377,psun:0.0031617491)0.977:0.0129734762)0.984:0.0237499666)0.93:0.0123999599)0.395:0.0020879791)0.901:0.0112615439,(mung:0.0594804999,(rnor:0.0396699258,(rsor:0.0265720262,((rdil:0.019967991,gdol:0.021105056)0.769:0.0085977828,((hall:0.0161568992,(pdel:0.0154147548,mnat:0.0156744369)0.42:0.0026038866)0.86:0.0102277821,(mpah:0.0196602678,(mcar:0.0092755425,((mmus:0.0019110393,mpwk:0.0051314407999999995)0.629:0.002015403,mspi:0.0049513531)0.808:0.0044021912)0.839:0.0074024138)0.821:0.0100953924)0.576:0.0041320338)0.387:0.0023496865)0.671:0.0051006134)0.982:0.0242160502)0.775:0.0069231995)0.98:0.0418120466)0.0:0.010544);";
    st = Tree(tree_str, get_subtrees=False, debug=False);
    # Some test trees

    #gt_str = "(mspr:0.0016324816,((mwsb:0.0000029568,mcar:0.0153129287):0.0021281260,(mcas:0.0019633923,(mmus:0.0000020851,((((hall:0.0174240612,(pdel:0.0130029362,mnat:0.0129445878):0.0029314306):0.0103158970,((((mung:0.0649440245,(((((((pcam:0.0031385524,psun:0.0033197989):0.0260127153,prob:0.0098921329):0.0196507658,cgri:0.0281001191):0.0019130533,maur:0.0355556774):0.0047041281,moch:0.0402020094):0.0057618080,pman:0.0379711174):0.0163580030,((jjac:0.0776442621,itri:0.1479591099):0.0249499942,ngal:0.0723486821):0.0508154463):0.0069335015):0.0123544338,rnor:0.0413244014):0.0021945352,(rdil:0.0203480143,gdol:0.0196458955):0.0075988551):0.0020342947,rsor:0.0270726716):0.0054596486):0.0105131536,mpah:0.0220127233):0.0095631120,mpwk:0.0009405335):0.0014971022):0.0021557143):0.0016426560):0.0028428833,mspi:0.0008185431);";

    #st = Tree(globs['orig-st-str'], debug=False);
    #gt = TREE.Tree(gt_str, debug=False);
    #st = TREE.Tree(tree_str);

    print(tree_str);    
    st.subtree = st.genSubtrees();
    st.tree_str = st.subtree[st.root];
    print(st.tree_str);


    # #print(globs['orig-st-str']);


    # # print(st.labeled_topo_str);
    # # st.showSplit();
    st.showAttrib("label", "anc", "desc");
    #print(st.rooted);
    #print(st.root);

    #new_labels = {'<1>': '0.949', '<2>': '0.896', '<4>': '0.851', '<5>': '0.933', '<6>': '0.913', '<7>': '0.276', '<8>': '0.885', '<10>': '0.374', '<3>': '0.822', '<11>': '0.792', '<9>': '0.656', '<13>': '0.78', '<12>': '0.561', '<14>': '0.535'};
    #cf_tree_str = st.addLabel(new_labels, delim="_");
    #print(cf_tree_str);


    #print(st.internals);

    #print(st.labeled_topo_str);
    #sampled_quartets = st.sampleQuartets();

    #print("\n\n\n\n\n");
    #print(sampled_quartets["<15>"])
    # print(st.tree_str);
    #pruned_tree = st.rmTips()
    #print(pruned_tree);


#############################################################################