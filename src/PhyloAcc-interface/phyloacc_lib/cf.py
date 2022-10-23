#############################################################################
# Calculates concordance factors
# Gregg Thomas
#############################################################################

import sys
import os
import time
import copy
import math
import phyloacc_lib.core as CORE
import phyloacc_lib.tree_old as TREEF

#############################################################################

def locusSCF(locus_item):
    locus, aln, quartets, st_nodes, aln_skip_chars, site_type = locus_item
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

    for node in st_nodes:
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

            quartet_scores[node][quartet] = { 'variable-sites' : 0, 'decisive-sites' : 0, 'concordant-sites' : 0, 'disco1-sites' : 0, 'disco2-sites' : 0, 'scf' : "NA" };
            # Initialize the current scf dict for the current quartet

            split1 = quartet[0];
            split2 = quartet[1];
            quartet_spec_list = split1 + split2;
            # Parse out the species for easier access

            if not all(spec in list(aln.keys()) for spec in quartet_spec_list):
                continue;
            # Skip any alignments that don't have all 4 species from the current quartet
            # Not sure if this is the best way to deal with missing data...

            ########################################
            if site_type == "zip":
                sub_aln = [ aln[spec] for spec in quartet_spec_list ];
                # Get the sub-alignment for the current quartet

                sites = zip(*sub_aln);

                for site in sites:
                    if any(char in site for char in aln_skip_chars):
                        continue;
                    # We only care about sites with full information for the quartet, so skip others here

                    uniq_site = set(site);
                    # The alleles at the current site

                    if len(uniq_site) > 1:
                        quartet_scores[node][quartet]['variable-sites'] += 1;
                    # If there is more than one allele, the site is variable                

                        if all(site.count(allele) == 2 for allele in uniq_site):
                        #if len(uniq_site) == 2:
                            quartet_scores[node][quartet]['decisive-sites'] += 1;
                        # If all alleles have a count of 2 at this site, it is decisive

                            if site[0] == site[1] and site[2] == site[3]:
                                quartet_scores[node][quartet]['concordant-sites'] += 1;
                            # If the alleles from split1 match and the alleles from split2 match, the site is concordant
                            elif site[0] == site[2] and site[1] == site[3]:
                                quartet_scores[node][quartet]['disco1-sites'] += 1;
                            elif site[0] == site[3] and site[1] == site[2]:
                                quartet_scores[node][quartet]['disco2-sites'] += 1;
                            # Otherwise they are discordant

                del sites;
                # Remove the sites from memory since it is a copy of the alignment
            ## Alternate way to loop over sites... faster and looks much nicer by but has memory leak??
            ########################################

            ########################################
            elif site_type == "loop":
            
                sub_aln = {};
                for spec in quartet_spec_list:
                    sub_aln[spec] = aln[spec];
                # Get the sub-alignment for the current quartet

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

                    if any(state in site for state in aln_skip_chars):
                        continue;
                    # We only care about sites with full information for the quartet, so skip others here
                ####################

                    uniq_site = set(site);
                    # The alleles at the current site

                    if len(uniq_site) > 1:
                        quartet_scores[node][quartet]['variable-sites'] += 1;
                    # If there is more than one allele, the site is variable                

                        if all(site.count(allele) == 2 for allele in uniq_site):
                        #if len(uniq_site) == 2:
                            quartet_scores[node][quartet]['decisive-sites'] += 1;
                        # If all alleles have a count of 2 at this site, it is decisive

                            site = split1_site + split2_site;

                            if site[0] == site[1] and site[2] == site[3]:
                                quartet_scores[node][quartet]['concordant-sites'] += 1;
                            # If the alleles from split1 match and the alleles from split2 match, the site is concordant
                            elif site[0] == site[2] and site[1] == site[3]:
                                quartet_scores[node][quartet]['disco1-sites'] += 1;
                            elif site[0] == site[3] and site[1] == site[2]:
                                quartet_scores[node][quartet]['disco2-sites'] += 1;
                            # Otherwise they are discordant


                        # if split1_site[0] == split1_site[1] and split2_site[0] == split2_site[1] and split1_site[0] != split2_site[0]:
                        # #if split1_site.count(split1_site[0]) == len(split1_site) and split2_site.count(split2_site[0]) == len(split2_site):
                        #     quartet_scores[node][quartet]['concordant-sites'] += 1;
                        # If the alleles from split1 match and the alleles from split2 match, the site is concordant
            ## End site loop
            ########################################

            if quartet_scores[node][quartet]['decisive-sites'] != 0:
                quartet_scores[node][quartet]['scf'] = quartet_scores[node][quartet]['concordant-sites'] / quartet_scores[node][quartet]['decisive-sites'];
            # Calculate the scf for the current quartet
        ## End quartet loop
        ####################

        for quartet in quartet_scores[node]:
            if quartet_scores[node][quartet]['scf'] != "NA":
                locus_scf[node]['scf-sum'] += quartet_scores[node][quartet]['scf'];
                locus_scf[node]['num-quartets'] += 1;
        # For all the quartets with sCF calculated, sum the sCF for this locus to be averaged later

        ####################

        if locus_scf[node]['num-quartets'] != 0:
            locus_scf[node]['avg-scf'] = locus_scf[node]['scf-sum'] / locus_scf[node]['num-quartets'];
            node_scf.append(locus_scf[node]['scf-sum'] / locus_scf[node]['num-quartets']);
        # Average all sCF from each quartet for this node
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

    step = "Sampling quartets";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    if globs['tree-data-type'] == 'class':
        sampled_quartets = globs['st'].sampleQuartets();

    elif globs['tree-data-type'] == 'func':
        root_desc = TREEF.getDesc(globs['root'], globs['st']);
        # Get the descendants of the root node to exclude them from sCF calculations
        sampled_quartets = TREEF.sampleQuartets(globs, root_desc);
    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # Sample quartets for all nodes

    ####################

    num_alns = len(globs['alns']);

    step = "Calculating per-locus sCF";
    step_start_time = CORE.report_step(globs, step, False, "Processed 0 / " + str(num_alns) + " alns...", full_update=True);
    # Status update

    ####################

    globs['scf'] = {};
    # The dictionary of averages over quartets

    quartet_counts = {};
    # The dictionary of site counts per quartet

    for node in sampled_quartets:
        if globs['tree-data-type'] == 'class':
            if node in globs['st'].tips or node == globs['st'].root:
                continue;
        if globs['tree-data-type'] == 'func':
            if node in globs['tips'] or node == globs['root']:
                continue;
        # Cannot calculate sCF for tips, the root, or node descendant from the root

        globs['scf'][node] = { 'variable-sites' : 0, 'decisive-sites' : 0, 
                                'concordant-sites' : 0, 'disco1-sites' : 0, 'disco2-sites' : 0,
                                 'quartet-scf-sum' : 0, 'total-quartets' : 0, 'avg-quartet-scf' : "NA",
                                 'scf' : "NA", 'sdf1' : -0.0, 'sdf2' : -0.0 };
        # The main scf dictionary with AVERAGES across quartets per node

        quartet_counts[node] = {};
        for quartet in sampled_quartets[node]:
            quartet_counts[node][quartet] = { 'variable-sites' : 0, 'decisive-sites' : 0, 
                                'concordant-sites' : 0, 'disco1-sites' : 0, 'disco2-sites' : 0 };
        # The dictionary of site COUNTS per quartet per node
    ## Initialize the dictionaries to calculate average sCF per node across all loci

    ####################

    # if globs['qstats']:
    #     qstats_file = os.path.join(globs['outdir'], "quartet-stats.csv");
    #     qfile = open(qstats_file, "w");
    #     headers = ["locus","node","quartet","variable-sites","decisive-sites","concordant-sites"];
    #     qfile.write(",".join(headers) + "\n");
    # For --qstats, creates and opens a file to write site counts for each quartet within the pool loop

    ####################

    with globs['scf-pool'] as pool:
        counter = 0;
        # A counter to keep track of how many loci have been completed

        #for locus in globs['alns']:
        #    result = locusSCF((locus, globs['alns'][locus], sampled_quartets, st, globs['aln-skip-chars']));
        # Serial version for debugging

        if globs['tree-data-type'] == 'class':
            internals = globs['st'].internals;
        elif globs['tree-data-type'] == 'func':
            internals = globs['internals'];

        # for result in pool.imap_unordered(locusSCF, ((locus, globs['alns'][locus], sampled_quartets, globs['st'].internals, globs['aln-skip-chars']) for locus in globs['alns'])):
        # for result in pool.imap_unordered(locusSCF, ((locus, globs['alns'][locus], sampled_quartets, globs['internals'], globs['aln-skip-chars']) for locus in globs['alns'])):
        for result in pool.imap_unordered(locusSCF, ((locus, globs['alns'][locus], sampled_quartets, internals, globs['aln-skip-chars'], globs['scf-site-type']) for locus in globs['alns'])):
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

            ## Summarize sCF stats for the locus
            ####################

            if globs['run-mode'] == 'adaptive':
                if globs['low-scf-branch-prop']:
                    if globs['aln-stats'][locus]['num-nodes'] != 0 and globs['aln-stats'][locus]['perc-low-scf-nodes'] >= globs['low-scf-branch-prop']:
                        globs['aln-stats'][locus]['batch-type'] = "gt";
                        globs['gt-loci'] += 1;
                    else:
                        globs['aln-stats'][locus]['batch-type'] = "st";
                        globs['st-loci'] += 1;
                # If a fraction of branches is provided with -s, check if that many branches have sCF below the min-sCF specified to put locus
                # in gene tree model

                else:
                    if globs['aln-stats'][locus]['num-nodes'] != 0 and globs['aln-stats'][locus]['node-scf-avg'] < globs['min-scf']:
                        globs['aln-stats'][locus]['batch-type'] = "gt";
                        globs['gt-loci'] += 1;
                    else:
                        globs['aln-stats'][locus]['batch-type'] = "st";
                        globs['st-loci'] += 1; 
                # If only a min scf value is provided (by default of with -scf), check if the average of all sCF values for
                # this locus is below that value to batch locus for gene tree model
            ## Adaptive run modes with sCF

            elif globs['run-mode'] == 'st':
                globs['st-loci'] += 1;
            elif globs['run-mode'] == 'gt':
                globs['gt-loci'] += 1;
            # If the run mode is adaptive, partition loci with more than 1/3 of nodes with low sCF to the gene tree model and all others
            # to the species tree model

            ## Partition the sequences based on scf            
            ####################

            for node in quartet_scores:
                for q in quartet_scores[node]:

                    # if globs['qstats']:
                    #     q_str = ";".join(q[0]) + ";" + ";".join(q[1]);
                    #     outline = [locus, node, q_str] + [str(quartet_scores[node][q][col]) for col in headers[3:]];
                    #     qfile.write(",".join(outline) + "\n");
                    # For --qstats, writes the quartet stats to a file for the current quartet

                    quartet_counts[node][q]['variable-sites'] += quartet_scores[node][q]['variable-sites'];
                    quartet_counts[node][q]['decisive-sites'] += quartet_scores[node][q]['decisive-sites'];
                    quartet_counts[node][q]['concordant-sites'] += quartet_scores[node][q]['concordant-sites'];
                    quartet_counts[node][q]['disco1-sites'] += quartet_scores[node][q]['disco1-sites'];
                    quartet_counts[node][q]['disco2-sites'] += quartet_scores[node][q]['disco2-sites'];
                    # Sum sites for this quartet for this node

            # For every quartet in every node, calculate sCF and sum up the number of sites in order to average
            # across nodes later

            ## Summarize sCF for each node to 
            ####################

            counter += 1;
            if counter % 100 == 0:
                cur_scf_time = CORE.report_step(globs, step, step_start_time, "Processed " + str(counter) + " / " + str(num_alns) + " alns...", full_update=True);
            # A counter and a status update every 100 loci
        ## End imap locus loop
    ## End pool
    ####################

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success", full_update=True);
    # Status update
    
    step = "Averaging site counts";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    for node in globs['scf']:
    # Average site counts for every node

        node_concordant = [];
        node_scfs = [];
        node_d1 = [];
        node_sdf1 = [];
        node_d2 = [];
        node_sdf2 = [];
        node_decisive = [];
        # Lists of site counts and cfs to build for each quartet

        for q in sampled_quartets[node]:
        # Append counts from every quartet on this node
            node_scfs.append(quartet_counts[node][q]['concordant-sites'] / quartet_counts[node][q]['decisive-sites']);
            node_sdf1.append(quartet_counts[node][q]['disco1-sites'] / quartet_counts[node][q]['decisive-sites']);
            node_sdf2.append(quartet_counts[node][q]['disco2-sites'] / quartet_counts[node][q]['decisive-sites']);
            node_concordant.append(quartet_counts[node][q]['concordant-sites']);
            node_d1.append(quartet_counts[node][q]['disco1-sites']);
            node_d2.append(quartet_counts[node][q]['disco2-sites']);
            node_decisive.append(quartet_counts[node][q]['decisive-sites']);

            if quartet_counts[node][q]['decisive-sites'] != 0:
                globs['scf'][node]['total-quartets'] += 1;
        ## End quartet loop

        globs['scf'][node]['scf'] = CORE.mean(node_scfs);
        globs['scf'][node]['concordant-sites'] = CORE.mean(node_concordant);
        globs['scf'][node]['sdf1'] = CORE.mean(node_sdf1);
        globs['scf'][node]['disco1-sites'] = CORE.mean(node_d1);
        globs['scf'][node]['sdf2'] = CORE.mean(node_sdf2);
        globs['scf'][node]['disco2-sites'] = CORE.mean(node_d2);

        globs['scf'][node]['decisive-sites'] = CORE.mean(node_decisive);
        # Calculate averages across quartets of sites and concordance factors
    ## End site count average loop

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");

    # if globs['qstats']:
    #     qfile.close()
    # For --qstats, closes the quartet stats file.

    return globs;
    
#############################################################################