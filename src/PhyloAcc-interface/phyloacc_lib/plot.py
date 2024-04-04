#############################################################################
# Functions to generate summary plots phyloacc datasets
# Gregg Thomas
#############################################################################

import os
import re
import math
import phyloacc_lib.core as PC
import phyloacc_lib.tree as TREE
import phyloacc_lib.tree_old as TREEF
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D as lines
import matplotlib.patches as mpatches

from Bio import Phylo
from io import StringIO

#############################################################################

def genPlots(globs):

    step = "Generating summary plots";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status updated

    mpl.rcParams["axes.spines.right"] = False;
    mpl.rcParams["axes.spines.top"] = False;

    mpl.rcParams["axes.labelsize"] = 16;
    mpl.rcParams["axes.edgecolor"] = "#595959";
    mpl.rcParams["axes.linewidth"] = 1.5;

    #mpl.rcParams["xtick.labelcolor"] = "#595959";
    mpl.rcParams["xtick.labelsize"] = 14;

    mpl.rcParams['xtick.color'] = "#595959";
    mpl.rcParams['xtick.major.size'] = 6;
    mpl.rcParams['xtick.major.width'] = 1.5;

    #mpl.rcParams["ytick.labelcolor"] = "#595959";
    mpl.rcParams["ytick.labelsize"] = 14;

    mpl.rcParams['ytick.color'] = "#595959";
    mpl.rcParams['ytick.major.size'] = 6;
    mpl.rcParams['ytick.major.width'] = 1.5;
    # Global theme settings for all matplotlib figures following
    
    ####################

    st_file = os.path.join(globs['plot-dir'], globs['input-tree-plot-file']);
    # The file to save the species tree figure

    branch_cols = PC.coreCol(pal="wilke", numcol=3);
    # The colors for the target, conserved, and outgroup branches

    if globs['tree-data-type'] == 'class':
        num_spec = globs['st'].num_tips;
        tips = globs['st'].tips;
        internals = globs['st'].internals;
    elif globs['tree-data-type'] == 'func':
        num_spec = len(globs['tips']);
        tips = globs['tips'];
        internals = globs['internals'];
    # The number of species in the input free, to adjust height of figure

    if globs['tree-data-type'] == 'class':
        targets, conserved, outgroups = TREE.categorizeBranches(globs, globs['st']);
    elif globs['tree-data-type'] == 'func':
        targets, conserved, outgroups = TREEF.categorizeBranches(globs, globs['st']);
    # Get full lists of branches for each category

    if globs['tree-data-type'] == 'class':
        tree_str = globs['st'].tree_str;
        for node in globs['st'].internals:
            tree_str = tree_str.replace(globs['st'].label[node] + ":", ":");
    elif globs['tree-data-type'] == 'func':
        tree_str = TREEF.addBranchLength(globs['labeled-tree'], globs['st'], no_label=True, keep_tp_label=True);
    # Re-add branch lengths and remove labels to the input tree for plotting, keep the treeparse label to add colors below

    handle = StringIO(tree_str);
    tree = Phylo.read(handle, "newick");
    # Parse the tree string with Bio
 
    target_lg, conserve_lg, outgroup_lg = False, False, False;
    #for clade in tree.get_terminals():
    for clade in tree.find_clades():
        if clade.name in targets:
            clade.color = branch_cols[0];
            target_lg = lines([0], [0], label="Targets", color = branch_cols[0]);
        elif clade.name in conserved:
            clade.color = branch_cols[1];
            conserve_lg = lines([0], [0], label="Conserved", color = branch_cols[1]);
        elif clade.name in outgroups:
            clade.color = branch_cols[2];
            outgroup_lg = lines([0], [0], label="Outgroup", color = branch_cols[2]);

        # if clade.name not in globs['st'].tips:
        # if clade.name not in globs['tips']:
        if clade.name not in tips:
            clade.name = None;
        # For internal nodes, remove the tree parse name so it doesn't show up on the plot
    # Color the tip branches based on their input category and specify their legend entries

    if num_spec < 20:
        fig = plt.figure(figsize=(num_spec/4, 25.4/5.08));
    elif num_spec < 60:
        fig = plt.figure(figsize=(num_spec/4, 25.4/2.54));
    else:
        #fig = plt.figure(figsize=(num_spec/4, 25.4/2.54));
        fig = plt.figure(figsize=(20, 30));
    # Specify the plot size depending on the number of species

    axes = fig.add_subplot(1, 1, 1);
    axes.axes.get_yaxis().set_visible(False);
    axes.spines['left'].set_visible(False);
    # Set the axes of the tree figure

    Phylo.draw(tree, axes=axes, show_confidence=False, do_show=False);
    # Draw the tree

    legend_handles = [ lg for lg in [target_lg, conserve_lg, outgroup_lg] if lg ];
    plt.legend(loc='upper left', handles=legend_handles);
    # Add the legend

    plt.savefig(st_file, dpi=100, bbox_inches='tight');
    # Save the figure

    # Species tree
    ####################

    aln_list = [ aln for aln in globs['aln-stats'] ];
    # A single list of alignment IDs for consistency between dictionary lookups

    aln_len_hist_file = os.path.join(globs['plot-dir'], globs['aln-len-plot-file']);
    for aln in globs['aln-stats']:
        if globs['aln-stats'][aln]['length'] == 0:
            print(aln);
    aln_lens = [ globs['aln-stats'][aln]['length'] for aln in aln_list ];
    
    plt.figure(figsize=(8,6));
    plt.hist(aln_lens, color=PC.coreCol(pal="wilke", numcol=1, offset=1)[0], bins="sturges", edgecolor="#999999");
    #plt.boxplot(locus_len)
    plt.xlabel("Alignment length");
    plt.ylabel("# loci");

    plt.savefig(aln_len_hist_file, dpi=100);
    # Locus length (hist)
    ####################

    seq_len_hist_file = os.path.join(globs['plot-dir'], globs['seq-len-plot-file']);
    seq_lens = [ globs['aln-stats'][aln]['avg-nogap-seq-len'] for aln in aln_list ];
    
    plt.figure(figsize=(8,6));
    plt.hist(seq_lens, color=PC.coreCol(pal="wilke", numcol=1, offset=2)[0], bins="sturges", edgecolor="#999999");
    #plt.boxplot(locus_len)
    plt.xlabel("Avg. sequence length without gaps (bp)");
    plt.ylabel("# loci");

    plt.savefig(seq_len_hist_file, dpi=100);
    # Avg. sequence length without gaps (hist)
    ####################

    # var_sites_hist_file = os.path.join(globs['plot-dir'], "variable-sites-hist.png");
    var_sites = [ globs['aln-stats'][aln]['variable-sites'] for aln in aln_list ];
    
    # plt.figure(figsize=(8,6));
    # plt.hist(var_sites, color=PC.coreCol(pal="wilke", numcol=1)[0], bins="sturges", edgecolor="#999999");
    # #plt.boxplot(locus_len)
    # plt.xlabel("# of variable sites");
    # plt.ylabel("# loci");

    # plt.savefig(var_sites_hist_file);
    # Variable sites (hist)
    ####################

    inf_sites_hist_file = os.path.join(globs['plot-dir'], globs['inf-sites-plot-file']);
    inf_sites = [ globs['aln-stats'][aln]['informative-sites'] for aln in aln_list ];
    
    plt.figure(figsize=(8,6));
    plt.hist(inf_sites, color=PC.coreCol(pal="wilke", numcol=1, offset=3)[0], bins="sturges", edgecolor="#999999");
    #plt.boxplot(locus_len)
    plt.xlabel("# of informative sites");
    plt.ylabel("# loci");

    plt.savefig(inf_sites_hist_file, dpi=100);
    # Informative sites (hist)
    ####################

    inf_sites_frac_hist_file = os.path.join(globs['plot-dir'], globs['inf-sites-frac-plot-file']);
    #inf_sites_frac = [ globs['aln-stats'][aln]['informative-sites'] / globs['aln-stats'][aln]['length'] for aln in aln_list ];
    inf_sites_frac = [ inf_sites[i] / aln_lens[i] for i in range(len(aln_lens)) ];

    plt.figure(figsize=(8,6));
    plt.hist(inf_sites_frac, color=PC.coreCol(pal="wilke", numcol=1, offset=4)[0], bins="sturges", edgecolor="#999999");
    #plt.boxplot(locus_len)
    plt.xlim([0, 1]);
    plt.xlabel("Fraction of sites that are informative");
    plt.ylabel("# loci");

    plt.savefig(inf_sites_frac_hist_file, dpi=100);
    # Fraction of sites that are informative (hist)
    ####################

    var_inf_sites_file = os.path.join(globs['plot-dir'], globs['var-inf-sites-plot-file']);

    slope, intercept = np.polyfit(var_sites, inf_sites, 1);
    inf_sites_pred = np.polyval([slope, intercept], var_sites);

    plt.figure(figsize=(8,6));
    plt.plot(var_sites, inf_sites_pred, color="#333333", linestyle='dashed', dashes=(5, 20));
    plt.scatter(var_sites, inf_sites, color=PC.coreCol(pal="wilke", numcol=1, offset=5)[0], alpha=0.25);
    plt.xlabel("# of variable sites");
    plt.ylabel("# of informative sites");

    plt.savefig(var_inf_sites_file, dpi=100);
    
    # Variable sites vs. informative sites (scatter w regression)
    ####################

    if globs['run-mode'] == 'adaptive':
        avg_scf_hist_file = os.path.join(globs['plot-dir'], globs['avg-scf-hist-file']);

        avg_scf = [ globs['aln-stats'][aln]['node-scf-avg'] for aln in aln_list if globs['aln-stats'][aln]['node-scf-avg'] != "NA" ];

        plt.figure(figsize=(8,6));
        plt.hist(avg_scf, color=PC.coreCol(pal="wilke", numcol=1, offset=6)[0], bins="sturges", edgecolor="#999999");
        plt.xlim([0, 1]);
        plt.xlabel("Avg. sCF across all branches per locus");
        plt.ylabel("# loci");

        plt.savefig(avg_scf_hist_file, dpi=100);        

        # Avg. scf per locus (hist)
        ####################

        low_scf_hist_file = os.path.join(globs['plot-dir'], globs['low-scf-hist-file']);

        perc_low_scf = [ globs['aln-stats'][aln]['perc-low-scf-nodes'] for aln in aln_list if globs['aln-stats'][aln]['perc-low-scf-nodes'] != "NA" ];

        plt.figure(figsize=(8,6));
        plt.hist(perc_low_scf, color=PC.coreCol(pal="wilke", numcol=1, offset=7)[0], bins="sturges", edgecolor="#999999");
        plt.xlim([0, 1]);
        plt.xlabel("% of branches with sCF below " + str(globs['min-scf']) + " per locus");
        plt.ylabel("# loci");

        plt.savefig(low_scf_hist_file, dpi=100);   

        # % of branches with low sCF per locus (hist)
        ####################

        scf_tree_file = os.path.join(globs['plot-dir'], globs['scf-tree-plot-file']);
        # The file to save the species tree figure

        # tree_str = globs['st'].tree_str;
        # tree_str = TREE.addBranchLength(globs['labeled-tree'], globs['st'], no_label=True, keep_tp_label=True);
        if globs['tree-data-type'] == 'class':
            tree_str = globs['st'].tree_str;
        elif globs['tree-data-type'] == 'func':
            tree_str = TREEF.addBranchLength(globs['labeled-tree'], globs['st'], no_label=True, keep_tp_label=True);
        # Re-add branch lengths and remove labels to the input tree for plotting, keep the treeparse label to add colors below

        for node in globs['scf']:
            tree_str = tree_str.replace(node, node + "_" + str(round(globs['scf'][node]['scf'], 2)) + "_");
        # For every node in the tree, add the averaged scf value over all loci to the label

        tree_str = re.sub("<[\d]+>[_]?", "", tree_str);
        #tree_str = re.sub("<[\d]+>", "", tree_str);

        #print(tree_str);
        # Re-add branch lengths and remove labels to the input tree for plotting

        handle = StringIO(tree_str);
        tree = Phylo.read(handle, "newick");
        # Parse the tree string with Bio
        
        # for clade in tree.get_terminals():
        #     if clade.name in globs['targets']:
        #         clade.color = branch_cols[0];
        #         target_lg = lines([0], [0], label="Targets", color = branch_cols[0]);
        #     elif clade.name in globs['conserved']:
        #         clade.color = branch_cols[1];
        #         conserve_lg = lines([0], [0], label="Conserved", color = branch_cols[1]);
        #     elif clade.name in globs['outgroup']:
        #         clade.color = branch_cols[2];
        #         outgroup_lg = lines([0], [0], label="Outgroup", color = branch_cols[2]);
        # # Color the tip branches based on their input category and specify their legend entries

        fig = plt.figure(figsize=(num_spec/2.54, 25.4/2.54));
        # Specify the plot size depending on the number of species

        axes = fig.add_subplot(1, 1, 1);
        axes.axes.get_yaxis().set_visible(False);
        axes.spines['left'].set_visible(False);
        # Set the axes of the tree figure

        Phylo.draw(tree, axes=axes, do_show=False);
        # Draw the tree

        # plt.legend(loc='upper left', handles=[target_lg, conserve_lg, outgroup_lg]);
        # Add the legend

        plt.savefig(scf_tree_file, dpi=100, bbox_inches='tight');
        # Save the figure        

        #globs['scf-labeled-tree']
        # scf tree (phylo)
        ####################

        bl_scf_file = os.path.join(globs['plot-dir'], globs['bl-scf-plot-file']);

        bls, scfs = [], [];
        # for node in globs['st'].internals:
        # for node in globs['internals']:
        for node in internals:
            if node not in globs['scf']:
                continue;

            if globs['tree-data-type'] == 'class':
                bls.append(float(globs['st'].bl[node]));
            elif globs['tree-data-type'] == 'func':
                bls.append(float(globs['st'][node][0]));
            scfs.append(globs['scf'][node]['scf']);
        # Gets the values out of their tables

        slope, intercept = np.polyfit(bls, scfs, 1);
        scfs_pred = np.polyval([slope, intercept], bls);

        plt.figure(figsize=(8,6));
        plt.plot(bls, scfs_pred, color="#333333", linestyle='dashed', dashes=(5, 20));
        plt.scatter(bls, scfs, color=PC.coreCol(numcol=1)[0], alpha=0.5);
        plt.xlabel("Branch length");
        plt.ylabel("sCF");

        plt.savefig(bl_scf_file, dpi=100);
        
        # Variable sites vs. informative sites (scatter w regression)
        ####################

    step_start_time = PC.report_step(globs, step, step_start_time, "Success");

#############################################################################

def getBFs(globs):
# A function to get the BFs into lists for plotting
# Also removes and warns about BFs that are infinite
    bf_types = [ "logBF1", "logBF2", "logBF3" ];
    bf_lists = { bf_type : [] for bf_type in bf_types };

    for locus in globs['all-loci']:
        for bf_type in bf_types:
            bf = float(globs['locus-stats']['elem_lik'][locus][bf_type]);
            if math.isinf(bf):
                print("WARNING: " + bf_type + " for locus " + locus + " is infinite. This locus may have failed to converge for this model. Excluding from plots...");
            else:
                bf_lists[bf_type].append(bf);

    return tuple(bf_lists[bf_type] for bf_type in bf_types)

#############################################################################

def genPlotsPost(globs):

    step = "Generating summary plots";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status updated

    mpl.rcParams["axes.spines.right"] = False;
    mpl.rcParams["axes.spines.top"] = False;

    mpl.rcParams["axes.labelsize"] = 16;
    mpl.rcParams["axes.edgecolor"] = "#595959";
    mpl.rcParams["axes.linewidth"] = 1.5;

    mpl.rcParams["xtick.labelcolor"] = "#595959";
    mpl.rcParams["xtick.labelsize"] = 14;

    mpl.rcParams['xtick.color'] = "#595959";
    mpl.rcParams['xtick.major.size'] = 6;
    mpl.rcParams['xtick.major.width'] = 1.5;

    mpl.rcParams["ytick.labelcolor"] = "#595959";
    mpl.rcParams["ytick.labelsize"] = 14;

    mpl.rcParams['ytick.color'] = "#595959";
    mpl.rcParams['ytick.major.size'] = 6;
    mpl.rcParams['ytick.major.width'] = 1.5;
    # Global theme settings for all matplotlib figures following
    
    ####################

    bf1s, bf2s, bf3s = getBFs(globs);
    # Just get lists of the BFs and filter out the infs for elements that failed to converge

    ####################

    bf_dist_file = os.path.join(globs['plot-dir'], globs['bf-dist-file']);
    plt.figure(figsize=(8,6));
    plt.boxplot([bf1s, bf2s, bf3s], notch=False, sym='+', vert=True, whis=1.5, labels=["BF1", "BF2", "BF3"]);
    plt.savefig(bf_dist_file, dpi=100);
    
    # All three BF dists (boxplot)
    ####################

    bf1_dist_file = os.path.join(globs['plot-dir'], globs['bf1-dist-file']);
    plt.figure(figsize=(8,6));
    plt.hist(bf1s, color=PC.coreCol(pal="wilke", numcol=1, offset=1)[0], bins="sturges", edgecolor="#999999");
    plt.savefig(bf1_dist_file, dpi=100);

    # BF1 (hist)
    ####################

    bf2_dist_file = os.path.join(globs['plot-dir'], globs['bf2-dist-file']);
    plt.figure(figsize=(8,6));
    plt.hist(bf2s, color=PC.coreCol(pal="wilke", numcol=1, offset=2)[0], bins="sturges", edgecolor="#999999");
    plt.xlabel("log BF2");
    plt.ylabel("# loci");
    plt.savefig(bf2_dist_file, dpi=100);
    
    # BF2 (hist)
    ####################

    bf3_dist_file = os.path.join(globs['plot-dir'], globs['bf3-dist-file']);
    plt.figure(figsize=(8,6));
    plt.hist(bf3s, color=PC.coreCol(pal="wilke", numcol=1, offset=2)[0], bins="sturges", edgecolor="#999999");
    plt.xlabel("log BF3");
    plt.ylabel("# loci");
    plt.savefig(bf3_dist_file, dpi=100);
    
    # BF3 (hist)
    ####################

    bf1_bf2_file = os.path.join(globs['plot-dir'], globs['bf1-bf2-plot-file']);

    #slope, intercept = np.polyfit(bf1s, bf2s, 1);
    #bf2s_pred = np.polyval([slope, intercept], bf1s);
    #bf1s_filtered = [ bf1s[i] for i in range(len(locus_list)) if bf1s[i] < 20 and bf1s[i] > -5 and bf2s[i] > -2 and bf2s[i] < 50 ]
    #bf2s_filtered = [ bf2s[i] for i in range(len(locus_list)) if bf1s[i] < 20 and bf1s[i] > -5 and bf2s[i] > -2 and bf2s[i] < 50 ]

    plt.figure(figsize=(8,6));
    #plt.plot(bf1s, bf2s_pred, color="#333333", linestyle='dashed', dashes=(5, 20));
    plt.scatter(bf1s, bf2s, color=PC.coreCol(pal="wilke", numcol=1, offset=5)[0], alpha=0.25);
    plt.axvline(x=globs['bf1-cutoff'], color='#d3d3d3', linestyle='--')
    plt.axhline(y=globs['bf2-cutoff'], color='#d3d3d3', linestyle='--')
    plt.xlabel("log BF1");
    plt.ylabel("log BF2");
    plt.savefig(bf1_bf2_file, dpi=100);

    # BF1 vs. BF2 (scatter)
    ####################

    m2_dist_file = os.path.join(globs['plot-dir'], globs['m2-locus-dist-file']);
    # Output file
    
    x_lineages = [];
    y_counts = [];
    raw_counts = [];
    # The x and y axes and the raw counts for summary calcs

    for i in range(0,globs['st'].num_nodes+1):
    # Loop over every node in the tree

        x_lineages.append(i);
        
        if i not in globs['m2-per-locus']:
            globs['m2-per-locus'][i] = 0;
            y_counts.append(0);
        else:
            y_counts.append(globs['m2-per-locus'][i]);
            raw_counts += [i] * globs['m2-per-locus'][i];
    
    plt.figure(figsize=(8,6));
    plt.bar(x_lineages, y_counts, color=PC.coreCol(pal="wilke", numcol=1, offset=2)[0], edgecolor="#999999");
    plt.xlabel("# of branches");
    plt.ylabel("# loci accelerated");
    plt.savefig(m2_dist_file, dpi=100);

    globs['summary']['m2-loci-max'] = max(raw_counts);
    globs['summary']['m2-loci-min'] = min(raw_counts);
    globs['summary']['m2-loci-med'] = PC.median(raw_counts);
    globs['summary']['m2-loci-avg'] = PC.mean(raw_counts);
    
    # Accelerated lineages per locus
    ####################

    m2_dist_file = os.path.join(globs['plot-dir'], globs['m2-lineage-dist-file']);

    sorted_lineage_counts = dict(sorted(globs['m2-per-lineage'].items(), key=lambda item: item[1], reverse=True));

    x_lineages = [];
    y_counts = [];
    num_rows = 10;
    cur_row = 1;
    for lineage in sorted_lineage_counts:
        x_lineages.append(lineage);
        y_counts.append(sorted_lineage_counts[lineage]);
        cur_row += 1;
        if cur_row > num_rows:
            break;

    plt.figure(figsize=(8,6));
    plt.barh(x_lineages, y_counts, color=PC.coreCol(pal="wilke", numcol=1, offset=1)[0], edgecolor="#999999");
    plt.xlabel("# loci accelerated (M2)");
    ax = plt.gca()
    ax.invert_yaxis()
    plt.tight_layout()
    plt.savefig(m2_dist_file, dpi=100);
    
    # Accelerated lineages per locus
    ####################

    globs['bf1'] = bf1s;
    globs['bf2'] = bf2s;
    globs['bf3'] = bf3s;
    
    step_start_time = PC.report_step(globs, step, step_start_time, "Success");
    return globs;

#############################################################################