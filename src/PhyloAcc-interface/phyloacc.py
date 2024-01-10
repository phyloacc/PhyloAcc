#!/usr/bin/env python3
#############################################################################
# This is the python front-end for PhyloAcc, a Bayesian substitution rate
# estimation program for conserved non-coding genomic elements. This script
# will handle user inputs, model selection, and batching jobs
#
# Gregg Thomas
# Summer 2021
#############################################################################

import sys
import os
import phyloacc_lib.core as PC
import phyloacc_lib.params as params
import phyloacc_lib.opt_parse as OP
import phyloacc_lib.treeio as TREEIO
import phyloacc_lib.seq as SEQ
import phyloacc_lib.tree as TREE
import phyloacc_lib.tree_old as TREEF
import phyloacc_lib.cf as CF
import phyloacc_lib.output as OUT
import phyloacc_lib.batch as BATCH
import phyloacc_lib.plot as PLOT
import phyloacc_lib.html as HTML

#############################################################################

if __name__ == '__main__':
# Main is necessary for multiprocessing to work on Windows.

    globs = params.init();
    # Get the global params as a dictionary.
    
    if any(v in sys.argv for v in ["--version", "-version", "--v", "-v"]):
        print("\n# PhyloAcc version " + globs['version'] + " released on " + globs['releasedate-patch'])
        sys.exit(0);
    # The version option to simply print the version and exit.
    # Need to get actual PhyloAcc version for this, and not just the interface version.

    if "--options" in sys.argv:
        opt_pad = 20;
        print();
        print("This is a list of other options for PhyloAcc that can be specified with -phyloacc:");
        print(PC.spacedOut("OPTION", opt_pad) + "DEFAULT");
        print("-" * 30);
        for opt in globs['phyloacc-defaults']:
            print(PC.spacedOut(opt, opt_pad) + globs['phyloacc-defaults'][opt]['default']);
        print();
        print("--options set: exiting...");
        sys.exit();
    
    globs['script-dir'] = os.path.dirname(os.path.realpath(__file__));

    if "--quiet" not in sys.argv:
        print("\n" + " ".join(sys.argv) + "\n");
        # Print the command used to call the program

        print("#");
        print("# " + "=" * 125);
        print(PC.welcome());
        if "-h" not in sys.argv:
            print("    Bayesian rate analysis of conserved");
            print("       non-coding genomic elements\n");
        # A welcome banner.
    else:
        print("# Running phyloacc_interface in quiet mode...");

    globs = OP.optParse(globs);
    # Getting the input parameters from optParse.

    #print(globs['tree-data-type']);
    #print(globs['scf-site-type']);
    # Tree debugging stuff

    if globs['info']:
        print("# --info SET. EXITING AFTER PRINTING PROGRAM INFO...\n#")
        sys.exit(0);
    if globs['norun']:
        print("# --norun SET. EXITING AFTER PRINTING OPTIONS INFO...\n#")
        sys.exit(0);
    # Early exit options

    step_start_time = PC.report_step(globs, "", "", "", start=True);
    # Initialize the step headers

    ## Initialization
    ####################

    if not globs['debug-aln']:
        if globs['tree-data-type'] == 'class':
            globs['tree-string'], globs['st'] = TREEIO.readST(globs);
        elif globs['tree-data-type'] == 'func':
            globs['tree-string'], globs['st'], globs['labeled-tree'], globs['root'], globs['tips'], globs['internals'] = TREEIO.readST(globs);

        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Tree read from mod file:", 45) + globs['tree-string']);
        # Report the tree that was read -- this line is necessary for phyloacc_post.py to read the tree used in the run

        if globs['debug-tree']:
            TREE.debugTree(globs);
            sys.exit(0);
        # Print the debug stats and exit if --debugtree is set

        if globs['label-tree']:
            if globs['tree-data-type'] == 'class':
                globs['st'].showAttrib("type", "label", "desc");

                print("\nInput tree with integer labels:");
                print(globs['st'].tree_str +"\n");

                #if not globs['st'].has_label:
                print("\nInput tree with descendant labels:");
                print(globs['st'].labelDesc());

            #print(globs['labeled-tree']);
            #print();
            sys.exit(0);
        elif not globs['st'].has_label:
            PC.errorOut("MAIN2", "Tree does not have internal nodes labeled, which is required for PhyloAcc. Use --labeltree and replace the tree in your .mod file with the labeled tree.", globs);

        #else:
        #    CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Original tree with node labels:\t" + globs['st'].tree_str);
        # Print the tree and exit if --labeltree is set

        ## Read species tree
        ####################

        globs = TREEIO.readSpeciesGroups(globs);

        ## Read species groups
        ####################

        if not globs['theta'] and globs['coal-tree-file']:
            if globs['tree-data-type'] == 'class':
                tree_str, tree = TREEIO.readST(globs, tree_type="coalescent");
            elif globs['tree-data-type'] == 'func':
                tree_str, tree_dict, tree, root, tips, internals = TREEIO.readST(globs, tree_type="coalescent");

        ## Read coalescent tree if provided
        ####################
    ## End check for --debugaln option

    globs = SEQ.readSeq(globs);

    ## Library to read input sequences
    ####################

    globs = SEQ.alnStats(globs);
    
    ## Calculate some basic alignment stats
    ####################

    if globs['debug-aln']:
        sys.exit(0);
    # Exit here if --debugaln is set    

    if globs['run-mode'] == 'adaptive':
        #globs = CF.scf(globs, globs['st'], globs['alns'], globs['scf-pool']);
        globs = CF.scf(globs);
        # globs = TREE.scf(globs);
    
    # Calculate avg. sCF per locus
    ####################

    globs = OUT.writeAlnStats(globs);
    # Write out the alignment summary stats

    if globs['run-mode'] == 'adaptive':
        #globs = OUT.writeSCFStats(globs);
        globs = TREEIO.writeCF(globs);
    # Write out the sCF summary stats

    ####################

    globs = BATCH.genJobFiles(globs);
    # Generates the locus specific job files (aln, bed, config, etc.) for phyloacc

    globs = BATCH.writeSnakemake(globs);
    # Generates the snakemake config and cluster profile

    ####################

    globs['smk-cmd'] = "snakemake -p -s " + os.path.abspath(globs['smk']);
    globs['smk-cmd'] += " --configfile " + os.path.abspath(globs['smk-config']);
    globs['smk-cmd'] += " --profile " + os.path.abspath(globs['profile-dir']);
    #globs['smk-cmd'] += " --cluster-status " + os.path.abspath(globs['status-script']);
    globs['smk-cmd'] += " --dryrun";
    # The snakemake command to run PhyloAcc

    ####################

    if globs['plot']:
       PLOT.genPlots(globs);
       globs = HTML.writeHTML(globs);

    ####################

    PC.endProg(globs);

#############################################################################

