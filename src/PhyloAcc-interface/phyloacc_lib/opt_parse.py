#############################################################################
# Parsing and printing the options and meta-info for PhyloAcc.
# Much of the error checking is done here as well.
#############################################################################

import sys
import os
import math
import argparse
import multiprocessing as mp
import phyloacc_lib.core as PC
import phyloacc_lib.tree as TREE

#############################################################################

def optParse(globs):
# This function handles the command line options and prepares the output directory and files.
# Defaults are set in params.py

    try:
        import psutil
        globs['psutil'] = True;
    except:
        globs['psutil'] = False;
    # Check if psutil is installed for memory usage stats.

    parser = argparse.ArgumentParser(description="PhyloAcc: Bayesian rate analysis of conserved non-coding genomic elements");

    parser.add_argument("-a", dest="aln_file", help="An alignment file with all loci concatenated. -b must also be specified. Expected as FASTA format for now. One of -a/-b or -d is REQUIRED.", default=False);
    parser.add_argument("-b", dest="bed_file", help="A bed file with coordinates for the loci in the concatenated alignment file. -a must also be specified. One of -a/-b or -d is REQUIRED.", default=False);
    parser.add_argument("-i", dest="id_file", help="A text file with locus names, one per line, corresponding to regions in the input bed file. -a and -b must also be specified.", default=False);
    parser.add_argument("-d", dest="aln_dir", help="A directory containing individual alignment files for each locus. Expected as FASTA format for now. One of -a/-b or -d is REQUIRED.", default=False);
    parser.add_argument("-m", dest="mod_file", help="A file with a background transition rate matrix and phylogenetic tree with branch lengths as output from PHAST. REQUIRED.", default=False);
    # Input

    parser.add_argument("-o", dest="out_dest", help="Desired output directory. This will be created for you if it doesn't exist. Default: phyloacc-[date]-[time]", default=False);
    # Output
    
    parser.add_argument("-t", dest="targets", help="Tip labels in the input tree to be used as target species. Enter multiple labels separated by semi-colons (;). REQUIRED.", default=False);
    parser.add_argument("-c", dest="conserved", help="Tip labels in the input tree to be used as conserved species. Enter multiple labels separated by semi-colons (;). Any species not specified in -t or -g will be inferred as conserved.", default=False);
    parser.add_argument("-g", dest="outgroup", help="Tip labels in the input tree to be used as outgroup species. Enter multiple labels separated by semi-colons (;).", default=False);
    # Phylo options 
    
    parser.add_argument("-burnin", dest="burnin", help="The number of steps to be discarded in the Markov chain as burnin. Default: 500", default=False);
    parser.add_argument("-mcmc", dest="mcmc", help="The total number of steps in the Markov chain. Default: 1000", default=False);
    parser.add_argument("-chain", dest="chain", help="The number of chains. Default: 1", default=False);
    # MCMC options

    parser.add_argument("-phyloacc", dest="phyloacc_opts", help="A catch-all option for other PhyloAcc parameters. Enter as a semi-colon delimited list of options: 'OPT1 value;OPT2 value'", default=False);
    # PhyloAcc options

    parser.add_argument("-path", dest="phyloacc_path", help="The path to the PhyloAcc binary. Default: PhyloAcc", default=False);
    # Dependency paths
    ## Note: For now we will likely need three dependency paths for the models, but eventually these should all be consolidated
    
    parser.add_argument("-scf", dest="low_scf", help="The value of sCF to consider as 'low' for any given branch in a locus. Default: 0.5", default=False);
    parser.add_argument("-s", dest="scf_prop", help="The proportion of branches to consider a locus for the gene tree model. Default: 0.3333, meaning if one-third of all branches for a given locus have low sCF, this locus will be run with the gene tree model.", default=False);
    parser.add_argument("-l", dest="coal_tree", help="A file containing a rooted, Newick formatted tree with the same topology as the species tree in the mod file (-m), but with branch lengths in coalescent units. When the gene tree model is used, one of -l or --theta must be set.", default=False);
    parser.add_argument("-r", dest="run_mode", help="Determines which version of PhyloAcc will be used. gt: use the gene tree model for all loci, st: use the species tree model for all loci, adaptive: use the gene tree model on loci with many branches with low sCF and species tree model on all other loci. Default: st", default=False);
    parser.add_argument("-n", dest="num_procs", help="The number of processes that this script should use. Default: 1.", type=int, default=1);
    parser.add_argument("-p", dest="procs_per_batch", help="The number of processes to use for each batch of PhyloAcc. Default: 1.", type=int, default=1);
    parser.add_argument("-j", dest="num_jobs", help="The number of jobs (batches) to run in parallel. Must be less than or equal to the total processes for PhyloAcc (-p). Default: 1.", type=int, default=1);
    # User params

    parser.add_argument("-batch", dest="batch_size", help="The number of loci to run per batch. Default: 50", default=False);
    # Batch options

    parser.add_argument("-part", dest="cluster_part", help="The partition or list of partitions (separated by commas) on which to run PhyloAcc jobs.", default=False);
    parser.add_argument("-nodes", dest="cluster_nodes", help="The number of nodes on the specified partition to submit jobs to. Default: 1.", default=False);
    parser.add_argument("-mem", dest="cluster_mem", help="The max memory for each job in GB. Default: 4.", default=False);
    parser.add_argument("-time", dest="cluster_time", help="The time in hours to give each job. Default: 1.", default=False);
    # Cluster options
    
    parser.add_argument("--theta", dest="theta", help="Set this to add gene tree estimation with IQ-tree and species estimation with ASTRAL for estimation of the theta prior. Note that a species tree with branch lengths in units of substitutions per site is still required with -m. Also note that this may add substantial runtime to the pipeline.", action="store_true", default=False);
    parser.add_argument("--labeltree", dest="labeltree", help="Simply reads the tree from the input mod file (-m), labels the internal nodes, and exits.", action="store_true", default=False);
    parser.add_argument("--overwrite", dest="ow_flag", help="Set this to overwrite existing files.", action="store_true", default=False);
    parser.add_argument("--appendlog", dest="append_log_flag", help="Set this to keep the old log file even if --overwrite is specified. New log information will instead be appended to the previous log file.", action="store_true", default=False);
    # User options
    
    #parser.add_argument("--plot", dest="plot_flag", help="Plot some summary statistics from the input data.", action="store_true", default=False);
    parser.add_argument("--summarize", dest="summarize_flag", help="Only generate the input summary plots and page. Do not write or overwrite batch job files.", action="store_true", default=False);
    parser.add_argument("--options", dest="options_flag", help="Print the full list of PhyloAcc options that can be specified with -phyloacc and exit.", action="store_true", default=False);
    parser.add_argument("--info", dest="info_flag", help="Print some meta information about the program and exit. No other options required.", action="store_true", default=False);
    parser.add_argument("--depcheck", dest="depcheck", help="Run this to check that all dependencies are installed at the provided path. No other options necessary.", action="store_true", default=False);
    parser.add_argument("--version", dest="version_flag", help="Simply print the version and exit. Can also be called as '-version', '-v', or '--v'", action="store_true", default=False);
    parser.add_argument("--quiet", dest="quiet_flag", help="Set this flag to prevent PhyloAcc from reporting detailed information about each step.", action="store_true", default=False);
    # Run options
    
    parser.add_argument("--qstats", dest="qstats", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--norun", dest="norun", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--debug", dest="debug_opt", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--debugtree", dest="debug_tree", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--nolog", dest="nolog_opt", help=argparse.SUPPRESS, action="store_true", default=False);
    # Performance tests
    
    args = parser.parse_args();
    # The input options and help messages

    ####################

    globs['call'] = " ".join(sys.argv);
    # Save the program call for later

    ####################

    if args.info_flag:
        globs['info'] = True;
        globs['log-v'] = -1;
        startProg(globs);
        return globs;
    # Parse the --info option and call startProg early if set

    ####################

    if args.norun:
        globs['norun'] = True;
        globs['log-v'] = -1;
    # Check if norun is set

    globs['overwrite'] = args.ow_flag;
    # Check if overwrite is set

    if args.debug_tree:
        globs['debug-tree'] = True;
        globs['log-v'] = -1;
    # Check if tree debugging is turned on -- will just read the input species tree for testing in the tree library

    ####################

    globs['label-tree'] = args.labeltree;
    # Parse the --labeltree option

    ####################

    if not args.mod_file:
        PC.errorOut("OP3", "A mod file must be provided with -m", globs);
    globs['mod-file'] = args.mod_file;
    # Check the mod file

    # Species tree/rate file
    ####################

    if not globs['label-tree'] and not globs['debug-tree']:
    ## Read main input options if --labeltree isn't set

        globs, deps_passed = PC.execCheck(globs, args);
        if args.depcheck:
            if deps_passed:
                print("\n# All dependencies PASSED.\n")
                sys.exit(0);
            else:
                print("\n# Some dependencies NOT FOUND. Please check your installations and provided paths.\n");
                sys.exit(1);
        # Check the dependency paths

        if (not args.aln_file and not args.aln_dir) or (args.aln_file and args.aln_dir):
            PC.errorOut("OP1", "One input method must be specified: -a or -d", globs);
        # Check that only one input type is specified

        if args.aln_file:
            if not args.bed_file:
                PC.errorOut("OP2", "A bed file with locus coordinates must be specified with -b when a concatenated alignment file is given with -a", globs);
            globs['aln-file'] = args.aln_file;
            globs['bed-file'] = args.bed_file;
            globs['id-file'] = args.id_file;
        elif args.aln_dir:
            globs['aln-dir'] = args.aln_dir;
        # Save the input type as a global param

        # Input files
        ####################

        if args.run_mode:
            if args.run_mode not in ["st", "gt", "adaptive"]:
                PC.errorOut("OP5", "Run mode (-r) must be one of: 'st', 'gt', or 'adaptive'.", globs);
            if args.run_mode in ["gt", "adaptive"] and not (args.theta or args.coal_tree):
                PC.errorOut("OP6", "When using the gene tree model with '-r gt' or '-r adaptive', a tree with branch lengths in coalescent units must also be provided with -l, or estimated with --theta.", globs);
            globs['run-mode'] = args.run_mode;

        if globs['run-mode'] == "st" and (args.theta or args.coal_tree):
            print("# WARNING: When using the species tree model '-r st', a tree with coalescent units is not required. -l and --theta will be ignored.");
        # Check the run mode option

        ## Run mode
        ####################

        job_sub_dirs = { 'job-alns' : 'alns', 'job-cfgs' : 'cfgs', 'job-bed' : 'bed', 'job-smk' : 'snakemake', 'job-out' : 'phyloacc-output' };
        # Expected job directories to create

        if globs['id-flag']:
            job_sub_dirs['job-ids'] = 'ids';

        if globs['run-mode'] != 'st':
            if args.theta and args.coal_tree:
                PC.errorOut("OP7", "Only one of -l or --theta should be specified.", globs);
            # Make sure both coalescent tree options aren't specified

            if args.theta:
                globs['theta'] = True;
                globs['coal-tree-file'] = "phyloacc-job-files/astral/astral-species-tree.treefile"
                job_sub_dirs = { 'job-alns' : 'alns', 'job-cfgs' : 'cfgs', 'job-bed' : 'bed', 'job-smk' : 'snakemake', 'job-out' : 'phyloacc-output', 'iqtree' : 'iqtree', 'astral' : 'astral' };
            # If --theta is specified, set global param to True and add directories for iqtree and astral
            elif args.coal_tree:
                if not os.path.isfile(args.coal_tree):
                    PC.errorOut("OP8", "Cannot find coalescent tree file (-l): " + args.coal_tree, globs);
                globs['coal-tree-file'] = args.coal_tree;
            # Get the file with the coalescent tree

        ## Coalescent tree/theta options
        ####################

        if not args.targets:
            PC.errorOut("OP10", "Target (-t) species must be specified.", globs);
        # Check that targets are specified. This is the only required group at this point.

        globs['input-groups']['targets'] = args.targets.replace("; ", ";").split(";");
        if args.outgroup:
            globs['input-groups']['outgroup'] = args.outgroup.replace("; ", ";").split(";");
        # Read the target and outgroup groups (if provided)

        if args.conserved:
            globs['input-groups']['conserved'] = args.conserved.replace("; ", ";").split(";");
        else:
            globs['input-groups']['conserved'] = [];
        # Read the conserved group if provided, and if not infer it from the other groups

        ## Read species group options
        ####################


        if args.low_scf:
            if not PC.isPosFloat(args.low_scf, maxval=1.0):
                PC.errorOut("OP13", "-low-scf must be a positive float, between 0 and 1.", globs);
            else:
                globs['min-scf'] = float(args.low_scf);

        if args.scf_prop:
            if not PC.isPosFloat(args.scf_prop, maxval=1.0):
                PC.errorOut("OP13", "-s must be a positive float, between 0 and 1.", globs);
            else:
                globs['low-scf-branch-prop'] = float(args.scf_prop);

        ## sCF options
        ####################

        opt_keys = {'burnin' : args.burnin, 'mcmc' : args.mcmc, 'chain' : args.chain };
        for opt in opt_keys:
            if opt_keys[opt]:
                if not PC.isPosInt(opt_keys[opt]):
                    PC.errorOut("OP13", "-b, -m, and -n must all be positive integers.", globs);
                globs[opt] = opt_keys[opt];
        # Get the MCMC options

        ## MCMC
        ####################

        if args.phyloacc_opts:
            while "; " in args.phyloacc_opts:
                args.phyloacc_opts = args.phyloacc_opts.replace("; ", ";");
            phyloacc_opts = args.phyloacc_opts.split(";");
            # A common error might be for users to put a space after each semi-colon, so we address that here and split the options

            for opt_val in phyloacc_opts:
            # For every option provided by the user

                while "  " in opt_val:
                    opt_val = opt_val.replace("  ", " ");
                opt, val = opt_val.split(" ");
                # Again a common error may be for multiple spaces to be included, so we get rid of those and split out the current
                # option and value

                opt = opt.upper();
                # PhyloAcc options are all upper case

                if opt not in globs['phyloacc-defaults']:
                    PC.errorOut("OP14", "One of the provided PhyloAcc options (-phyloacc) is invalid: " + opt, globs);
                # Check if the given option is even one of the possible ones and error out if not

                option_pass = True;
                if val == globs['phyloacc-defaults'][opt]['default']:
                    continue;
                # If the provided value for the current option is the same as the default value, just leave it out

                elif globs['phyloacc-defaults'][opt]['type'] == "POS_INT":
                    if not PC.isPosInt(val):
                        option_pass = False;                
                elif globs['phyloacc-defaults'][opt]['type'] == "POS_FLOAT":
                    if not PC.isPosFloat(val):
                        option_pass = False;           
                elif globs['phyloacc-defaults'][opt]['type'] in ["PROB", "PERC"]:
                    if not PC.isPosFloat(val, maxval=1.0):
                        option_pass = False;          
                elif globs['phyloacc-defaults'][opt]['type'] == "BOOL":
                    val = val.upper();
                    if val not in ["1", "0"]:
                        option_pass = False; 
                # Check each value against its possible values for its given type

                if not option_pass:
                    PC.errorOut("OP15", "The value of the provided PhyloAcc option " + opt + " is invalid for its type (" + globs['phyloacc-defaults'][opt]['type'] + "): " + val, globs);  
                # If the value is not valid, error out

                globs['phyloacc-opts'].append(opt + " " + val);
                # If the value for the given option passes all checks, then add this option and value to the global list
            ## End option loop

        ## Phyloacc options
        ####################

        if not args.out_dest:
            globs['outdir'] = "phyloacc-out-" + globs['startdatetime'];
        else:
            globs['outdir'] = args.out_dest;
        # Get the output directory

        if not globs['overwrite'] and os.path.exists(globs['outdir']):
            PC.errorOut("OP16", "Output directory already exists: " + globs['outdir'] + ". Specify new directory name OR set --overwrite to overwrite all files in that directory.", globs);
        if not os.path.isdir(globs['outdir']) and not globs['norun'] and not globs['info'] and not globs['debug-tree']:
            os.makedirs(globs['outdir']);
        # Main output dir

        ####################

        if args.summarize_flag:
            globs['batch'] = False;
        # Set the plot_flag to True if the plot_only_flag is True

        ## Batch option -- after checking input files check to see if the user wants to write the batch files or just get the summary
        ####################   

        if globs['plot']:
            globs['plot-dir'] = os.path.join(globs['outdir'], "plots");
            if not os.path.isdir(globs['plot-dir']) and not globs['norun']:
                os.makedirs(globs['plot-dir']);
            # Plot option

            globs['html-dir'] = os.path.join(globs['outdir'], "html");
            # if not os.path.isdir(globs['html-dir']) and not globs['norun']:
            #     os.makedirs(globs['html-dir']);
            globs['html-file'] = os.path.join(globs['outdir'], "phyloacc-pre-run-summary.html");
            # HTML directory
        
        # Parse the plot input and output locations
        ####################  

        if globs['batch']:
            globs['job-dir'] = os.path.join(globs['outdir'], "phyloacc-job-files");
            if not os.path.isdir(globs['job-dir']) and not globs['norun']:
                os.makedirs(globs['job-dir']);
            # Main job file dir

            if not globs['norun']:
                for subdir in job_sub_dirs:
                    globs[subdir] = os.path.join(globs['job-dir'], job_sub_dirs[subdir]);
                    if not os.path.isdir(globs[subdir]):
                        os.makedirs(globs[subdir]);
            # Job output subdirs

        if globs['coal-tree-file']:
            globs['coal-tree-file'] = os.path.abspath(globs['coal-tree-file']);

        globs['run-name'] = os.path.basename(os.path.normpath(globs['outdir']));
        globs['logfilename'] = os.path.join(globs['outdir'], globs['run-name'] + ".log");
        # Log file

        if not args.append_log_flag and not globs['norun'] and not globs['debug-tree']:
            logfile = open(globs['logfilename'], "w");
            logfile.write("");
            logfile.close();
        # Prep the logfile to be overwritten

        ## Output files and directories
        ####################

        if args.batch_size:
            if not PC.isPosInt(args.batch_size):
                PC.errorOut("OP17", "The number of loci per batch (-batch) must be a positive integer.", globs);
            else:
                globs['batch-size'] = int(args.batch_size);
        # Batch size

        globs['procs-per-job'] = PC.isPosInt(args.procs_per_batch, default=1);
        globs['num-jobs'] = PC.isPosInt(args.num_jobs, default=1);
        globs['total-procs'] = globs['procs-per-job'] * globs['num-jobs'];
        # Determine resource allocation for PhyloAcc

        globs['num-procs'] = PC.isPosInt(args.num_procs, default=1);
        globs['aln-pool'] = mp.Pool(processes=globs['num-procs']);
        #globs['aln-pool'] = mp.Pool(processes=48);
        globs['scf-pool'] = mp.Pool(processes=globs['num-procs']);
        # Create the pool of processes for sCF calculation here so we copy the memory profile of the parent process
        # before we've read any large data in

        # Batch size and resource allocation
        ####################

        if not args.cluster_part:
            PC.errorOut("OP18", "At least one cluster partition must be specified with -part.", globs);
        else:
            globs['partition'] = args.cluster_part;
        # Cluster partition option (required)

        if args.cluster_nodes:
            if not PC.isPosInt(args.cluster_nodes):
                PC.errorOut("OP19", "The number of nodes specified (-nodes) must be a positive integer.", globs);
            else:
                globs['num-nodes'] = args.cluster_nodes;
        # Cluster memory option

        if args.cluster_mem:
            if not PC.isPosInt(args.cluster_mem):
                PC.errorOut("OP20", "The specified cluster memory (-mem) must be a positive integer in GB.", globs);
            else:
                globs['mem'] = args.cluster_mem;
        # Cluster memory option

        if args.cluster_time:
            if not PC.isPosInt(args.cluster_time):
                PC.errorOut("OP21", "The specified cluster time (-time) must be a positive integer in hours.", globs);
            else:
                globs['time'] = args.cluster_time + ":00:00";
        # Cluster memory option

        ## Cluster options
        ####################
    ## End --labeltree check

    else:
        globs['log-v'] = -1;
    # If --labeltree or --debugtree is set, turn off log reporting

    ####################

    globs = PC.fileCheck(globs);
    # Make sure all the input file actually exist, and get their
    # full paths

    ####################

    if args.quiet_flag:
        globs['quiet'] = True;
        globs['log-v'] = 0;
    # Check the quiet option

    if args.qstats:
        globs['qstats'] = True;
    # Check for the internal quartet stats option to write to a file.

    if globs['psutil']:
        globs['pids'] = [psutil.Process(os.getpid())];
    # Get the starting process ids to calculate memory usage throughout.

    ## Internal stuff
    ####################

    startProg(globs);
    # After all the essential options have been set, call the welcome function.

    return globs;

#############################################################################

def startProg(globs):
# A nice way to start the program.
    print("#");
    PC.printWrite(globs['logfilename'], globs['log-v'], "# Welcome to PhyloAcc -- Bayesian rate analysis of conserved non-coding genomic elements.");
    PC.printWrite(globs['logfilename'], globs['log-v'], "# Version " + globs['version'] + " released on " + globs['releasedate-patch']);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# PhyloAcc was developed by " + globs['devs']);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# Citation:      " + globs['doi']);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# Website:       " + globs['http']);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# Report issues: " + globs['github']);
    PC.printWrite(globs['logfilename'], globs['log-v'], "#");
    PC.printWrite(globs['logfilename'], globs['log-v'], "# The date and time at the start is: " + PC.getDateTime());
    PC.printWrite(globs['logfilename'], globs['log-v'], "# Using Python version:              " + globs['pyver'] + "\n#");
    PC.printWrite(globs['logfilename'], globs['log-v'], "# The program was called as:         " + globs['call'] + "\n#");

    if globs['info']:
        return;
    # If --info is set, return after printing program info

    #######################

    pad = 45;
    opt_pad = 50;
    PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# INPUT/OUTPUT INFO:");

    if globs['aln-file']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Alignment file:", pad) + globs['aln-file']);
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Bed file:", pad) + globs['bed-file']);
    elif globs['aln-dir']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Alignment directory:", pad) + globs['aln-dir']);

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Tree/rate file (mod file from phyloFit):", pad) + globs['mod-file']);
    #PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Tree read from mod file:", pad) + globs['tree-string']);

    if globs['debug-tree']:
        print("\n--debugtree SET. READING TREE AND EXITING.\n");
        return;       

    if globs['label-tree']:
        print("\n--labeltree SET. READING TREE AND EXITING.\n");
        return;
    # If --labeltree is set, print a message and return. 

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Output directory:", pad) + globs['outdir']);
    if globs['batch']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# PhyloAcc run directory:", pad) + globs['job-dir']);
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Log file:", pad) + os.path.basename(globs['logfilename']));
    # Input/Output
    #######################

    PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# DEPENDENCY PATHS:");    
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Program", pad) + "Specified Path");

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# PhyloAcc", pad) + globs['phyloacc']);
    #PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# PhyloAcc-gBGC", pad) + globs['phyloacc-gbgc']);
    # Dependency paths
    #######################

    PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# SPECIES GROUPS:");    
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Group", pad) + PC.spacedOut("Species", opt_pad));

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Targets (-t)", pad) + ";".join(globs['input-groups']['targets']));
    if globs['input-groups']['conserved']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Conserved (-c)", pad) + ";".join(globs['input-groups']['conserved']));
    if globs['input-groups']['outgroup']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Outgroups (-g)", pad) + ";".join(globs['input-groups']['outgroup']));
    # Species groups
    #######################

    PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# CLUSTER OPTIONS:");    
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Option", pad) + PC.spacedOut("Setting", opt_pad));

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Partition(s)", pad) + globs['partition']);
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Number of nodes", pad) + globs['num-nodes']);
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Max mem per job (gb)", pad) + globs['mem']);
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Time per job", pad) + globs['time']); 
    # Cluster options
    #######################

    if globs['phyloacc-opts']:
        PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
        PC.printWrite(globs['logfilename'], globs['log-v'], "# PHYLOACC OPTIONS:");    
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Option", pad) + PC.spacedOut("Setting", opt_pad));

        for opt_val in globs['phyloacc-opts']:
            opt, val = opt_val.split(" ");
            PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# " + opt, pad) + val);

    # PhyloAcc options
    #######################

    PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# OPTIONS INFO:");    
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Option", pad) + PC.spacedOut("Current setting", opt_pad) + "Current action");

    if globs['aln-file']:
        if globs['id-file']:
            PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -i:", pad) + 
                PC.spacedOut(str(globs['id-file']), opt_pad) +
                "Only loci names specified in this file will be tested.");
        else:
            PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -i:", pad) + 
                PC.spacedOut(str(globs['id-file']), opt_pad) +
                "No ID file provided, all loci in input bed file will be tested.");
    # ID file option

    ####################

    if globs['run-mode'] == "st":
        info_str = "All loci will be run with the species tree model of PhyloAcc";
    elif globs['run-mode'] == "gt":
        info_str = "All loci will be run with the gene tree model of PhyloAcc";
    elif globs['run-mode'] == "adaptive":
        info_str = "Loci with many branches with low sCF will be run using the gene tree model of PhyloAcc while all others will use the species tree model";
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -r:", pad) + 
                PC.spacedOut(str(globs['run-mode']), opt_pad) +
                info_str);                  
    # Run mode option

    ####################

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -low-scf:", pad) + 
                PC.spacedOut(str(globs['min-scf']), opt_pad) +
                "The sCF value to consider as low");

    if globs['low-scf-branch-prop']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -s:", pad) + 
                    PC.spacedOut(str(globs['low-scf-branch-prop']), opt_pad) +
                    "The proportion of branches in a locus with low sCF to consider it for the gene tree model");
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -s:", pad) + 
                    PC.spacedOut("NA", opt_pad) +
                    "Branch sCF values will be averaged for each locus and if the average is below the min sCF it will be run with the gene tree model");        
    # sCF

    ####################

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -burnin:", pad) + 
                PC.spacedOut(str(globs['burnin']), opt_pad) +
                "This number of steps in the chain will discarded as burnin");         
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -mcmc:", pad) + 
                PC.spacedOut(str(globs['mcmc']), opt_pad) +
                "The number of steps in each chain");      
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -chain:", pad) + 
                PC.spacedOut(str(globs['chain']), opt_pad) +
                "The number of chains to run");      
    # MCMC options

    ####################

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Loci per batch (-batch)", pad) + 
                PC.spacedOut(str(globs['batch-size']), opt_pad) + 
                "PhyloAcc will run this many loci in a single command.");
    # Batch size

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Current processes (-n)", pad) + 
                PC.spacedOut(str(globs['num-procs']), opt_pad) + 
                "This interface will use this many processes.");
    # Interface procs

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Jobs (-j)", pad) + 
                PC.spacedOut(str(globs['num-jobs']), opt_pad) + 
                "PhyloAcc will submit this many jobs concurrently.");
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Processes per job (-p)", pad) + 
                PC.spacedOut(str(globs['procs-per-job']), opt_pad) + 
                "Each job will use this many processes.");
    # snakemake procs

    ####################

    # if globs['plot']:
    #     plot_status = "True";
    #     plot_status_str = "An HTML file summarizing the input data will be written to " + globs['html-file'];
    # else:
    #     plot_status = "False";
    #     plot_status_str = "No HTML summary file will be generated.";

    # PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --plot", pad) + 
    #             PC.spacedOut(plot_status, opt_pad) + 
    #             plot_status_str);
    # --plot option

    ####################

    if globs['batch']:
        batch_status = "False";
        batch_status_str = "PhyloAcc batch files will be generated and written to the job directory specified above.";
    else:
        batch_status = "True";
        batch_status_str = "No PhyloAcc batch files will be generated and those existing will not be overwritten.";

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --summarize", pad) + 
                PC.spacedOut(batch_status, opt_pad) + 
                batch_status_str);
    # --plotonly option

    ####################

    if globs['theta']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --theta", pad) +
                    PC.spacedOut("True", opt_pad) + 
                    "A species tree with branch lengths in coalescent units will be estimated.");
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --theta", pad) +
                    PC.spacedOut("False", opt_pad) + 
                    "A species tree with branch lengths in coalescent units will NOT be estimated.");            
    # Report the theta option.

    if globs['coal-tree-str']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -l", pad) +
                    PC.spacedOut(globs['coal-tree-file'], opt_pad) + 
                    "The species tree but with branch lengths in coalescent units.");
    # Report the -l option

    ####################

    if globs['overwrite']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --overwrite", pad) +
                    PC.spacedOut("True", opt_pad) + 
                    "PhyloAcc will OVERWRITE the existing files in the specified output directory.");
    # Reporting the overwrite option.

    ####################

    if not globs['quiet']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --quiet", pad) + 
                    PC.spacedOut("False", opt_pad) + 
                    "Time, memory, and status info will be printed to the screen while PhyloAcc is running.");
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --quiet", pad) + 
                    PC.spacedOut("True", opt_pad) + 
                    "No further information will be printed to the screen while PhyloAcc is running.");
        PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
        PC.printWrite(globs['logfilename'], globs['log-v'], "# Running...");
    # Reporting the quiet option.

    ####################

    # if globs['debug']:
    #     PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --debug", pad) + 
    #                 PC.spacedOut("True", opt_pad) + 
    #                 "Printing out a bit of debug info.");
    # Reporting the debug option.

    if globs['qstats']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --qstats", pad) + 
            PC.spacedOut("True", opt_pad) + 
            "Writing out a file with quartet site counts.");

    if globs['norun']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --norun", pad) + 
                    PC.spacedOut("True", opt_pad) + 
                    "ONLY PRINTING RUNTIME INFO.");
        PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    # Reporting the norun option.

    # Other options
    #######################

#############################################################################