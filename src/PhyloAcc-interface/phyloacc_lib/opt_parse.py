#############################################################################
# Parsing and printing the options and meta-info for PhyloAcc.
# Much of the error checking is done here as well.
#############################################################################

import sys
import os
import math
import argparse
import yaml
import multiprocessing as mp
import phyloacc_lib.core as PC
import phyloacc_lib.tree as TREE

#############################################################################

def getOpt(args_var, arg_str, arg_type, arg_default, config, flags, globs, check=True):
# This function gets the value of an option from the command line, config file, or default.
# It also checks the value for validity if check=True.

    if args_var:
        param_value = args_var;
    elif config and arg_str in config and config[arg_str]:
        param_value = config[arg_str];
    else:
        param_value = arg_default;
    # Get the value from the command line, config file, or default

    if check:
        if arg_type == str and not param_value:
            PC.errorOut("OP4", "The value for " + arg_str + " must be specified.", globs);

        if arg_type == bool and param_value not in [True, False, None, 0, 1]:
            PC.errorOut("OP4", "The value for " + arg_str + " (" + flags[arg_str] + ") must be a boolean (True/False). Current value: " + str(param_value), globs);

        if arg_type in [int, "INTSTR"] and not PC.isPosInt(param_value):
            PC.errorOut("OP4", "The value for " + arg_str + " must be a positive integer.", globs);

        if arg_type == float and not PC.isPosFloat(param_value):
            PC.errorOut("OP4", "The value for " + arg_str + " must be a positive float.", globs);

        if arg_type == "PROP" and not PC.isPosFloat(param_value, maxval=1.0):
            PC.errorOut("OP4", "The value for " + arg_str + " must be a positive float between 0 and 1.", globs);

        if arg_type == "FILE" and param_value and not os.path.isfile(param_value):
            if not optional:
                PC.errorOut("OP4", "The path provided for " + arg_str + " does not exist or is not a file. Current path: " + str(param_value), globs);

        if arg_type == "DIR" and param_value and not os.path.isdir(param_value):
            if not optional:
                PC.errorOut("OP4", "The path provided for " + arg_str + " does not exist or is not a directory. Current path: " + str(param_value), globs);

        if type(arg_type) == list and param_value not in arg_type:
            PC.errorOut("OP4", "The value for " + arg_str + " must be one of: " + ", ".join(arg_type) + ". Current value: " + str(param_value), globs);

        if arg_type == "INTSTR":
            param_value = str(param_value);

    return param_value;

#############################################################################

def addArgument(parser, flag, dest, arg_type, help_str):
# This function adds an argument to the parser with the specified flag, destination, type, and help message.

    if isinstance(flag, list):
        # If flags is a list, unpack it
        if arg_type == bool:
            parser.add_argument(*flag, dest=dest, help=help_str, action="store_true")
        else:
            parser.add_argument(*flag, dest=dest, help=help_str, type=arg_type, default=False)
        return {dest: flag[0]}
    else:
        # If flags is a single string
        if arg_type == bool:
            parser.add_argument(flag, dest=dest, help=help_str, action="store_true")
        else:
            parser.add_argument(flag, dest=dest, help=help_str, type=arg_type, default=False)
        return {dest: flag}


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
    arg_flags = {};

    ##########

    arg_flags.update(addArgument(parser, "--config", "config_file", str, 
        "A YAML formatted file with the arguments for the program to be used in lieu of command line arguments."));
    # parser.add_argument("--config", dest="config_file", 
    #     help="A YAML formatted file with the arguments for the program to be used in lieu of command line arguments.", 
    #     default=False);
    ## Config

    ##########

    arg_flags.update(addArgument(parser, "-a", "aln_file", str,
        "An alignment file in FASTA format with all loci concatenated. -b must also be specified. Expected as FASTA format for now. One of -a/-b or -d is REQUIRED."));
    #parser.add_argument("-a", dest="aln_file", help="An alignment file with all loci concatenated. -b must also be specified. Expected as FASTA format for now. One of -a/-b or -d is REQUIRED.", default=False);

    arg_flags.update(addArgument(parser, "-b", "bed_file", str,
        "A bed file with coordinates for the loci in the concatenated alignment file. -a must also be specified. One of -a/-b or -d is REQUIRED."));
    #parser.add_argument("-b", dest="bed_file", help="A bed file with coordinates for the loci in the concatenated alignment file. -a must also be specified. One of -a/-b or -d is REQUIRED.", default=False);

    arg_flags.update(addArgument(parser, "-i", "id_file", str,
        "A text file with locus names, one per line, corresponding to regions in the input bed file. -a and -b must also be specified."));    
    #parser.add_argument("-i", dest="id_file", help="A text file with locus names, one per line, corresponding to regions in the input bed file. -a and -b must also be specified.", default=False);
   
    arg_flags.update(addArgument(parser, "-d", "aln_dir", str,
        "A directory containing individual alignment files for each locus in FASTA format. One of -a/-b or -d is REQUIRED."));   
    #parser.add_argument("-d", dest="aln_dir", help="A directory containing individual alignment files for each locus. Expected as FASTA format for now. One of -a/-b or -d is REQUIRED.", default=False);

    arg_flags.update(addArgument(parser, "--softmask", "softmask", bool,
        "Treat lowercase bases in alignments as masked (converted to N) before processing."));

    arg_flags.update(addArgument(parser, "-m", "mod_file", str,
        "A file with a background transition rate matrix and phylogenetic tree with branch lengths as output from phyloFit. REQUIRED."));
    #parser.add_argument("-m", dest="mod_file", help="A file with a background transition rate matrix and phylogenetic tree with branch lengths as output from PHAST. REQUIRED.", default=False);
    ## Input

    ##########

    arg_flags.update(addArgument(parser, "-o", "out_dir", str,
        "Desired output directory. This will be created for you if it doesn't exist. Default: phyloacc-[date]-[time]"));
    #parser.add_argument("-o", dest="out_dest", help="Desired output directory. This will be created for you if it doesn't exist. Default: phyloacc-[date]-[time]", default=False);
    ## Output

    ##########

    arg_flags.update(addArgument(parser, "--filter", "filter_alns", bool,
        "Set this flag to filter out loci with no informative sites before running PhyloAcc."));
    ## Alignments

    ##########
    
    arg_flags.update(addArgument(parser, "-t", "targets", str,
        "Tip labels in the input tree to be used as target species. Enter multiple labels separated by semi-colons (;). REQUIRED."));
    #parser.add_argument("-t", dest="targets", help="Tip labels in the input tree to be used as target species. Enter multiple labels separated by semi-colons (;). REQUIRED.", default=False);
    
    arg_flags.update(addArgument(parser, "-c", "conserved", str,
        "Tip labels in the input tree to be set as non-target species. Enter multiple labels separated by semi-colons (;). Any species not specified in -t or -g will be placed in this category."));    
    #parser.add_argument("-c", dest="conserved", help="Tip labels in the input tree to be used as conserved species. Enter multiple labels separated by semi-colons (;). Any species not specified in -t or -g will be inferred as conserved.", default=False);
    
    arg_flags.update(addArgument(parser, "-g", "outgroup", str,
        "Tip labels in the input tree to be used as outgroup species. Enter multiple labels separated by semi-colons (;)."));
    #parser.add_argument("-g", dest="outgroup", help="Tip labels in the input tree to be used as outgroup species. Enter multiple labels separated by semi-colons (;).", default=False);
    ## Phylo options 

    ##########
    
    arg_flags.update(addArgument(parser, "-burnin", "burnin", int,
        "The number of steps to be discarded in the Markov chain as burnin. Default: 500"));
    #parser.add_argument("-burnin", dest="burnin", help="The number of steps to be discarded in the Markov chain as burnin. Default: 500", default=False);

    arg_flags.update(addArgument(parser, "-mcmc", "mcmc", int,
        "The total number of steps in the Markov chain. Default: 1000"));
    #parser.add_argument("-mcmc", dest="mcmc", help="The total number of steps in the Markov chain. Default: 1000", default=False);
    
    arg_flags.update(addArgument(parser, "-thin", "thin", int,
        "For the gene tree model, the number of MCMC steps between gene tree sampling. The total number of MCMC steps specified with -mcmc will be scaled by this as mcmc*thin Default: 1"));
    #parser.add_argument("-thin", dest="thin", help="For the gene tree model, the number of MCMC steps between gene tree sampling. The total number of MCMC steps specified with -mcmc will be scaled by this as mcmc*thin Default: 1", default=False);
    
    arg_flags.update(addArgument(parser, "-chain", "chain", int,
        "The number of chains. Default: 1"));
    #parser.add_argument("-chain", dest="chain", help="The number of chains. Default: 1", default=False);
    # MCMC options

    ##########

    parser.add_argument("-phyloacc", dest="phyloacc_opts", help="A catch-all option for other PhyloAcc parameters. Enter as a semi-colon delimited list of options: 'OPT1 value;OPT2 value'", default=False);
    # PhyloAcc options

    ##########

    arg_flags.update(addArgument(parser, "-st-path", "phyloacc_st_path", str,
        "The path to the PhyloAcc-ST binary. Default: PhyloAcc-ST"));
    #parser.add_argument("-st-path", dest="phyloacc_st_path", help="The path to the PhyloAcc-ST binary. Default: PhyloAcc-ST", default=False);
    
    arg_flags.update(addArgument(parser, "-gt-path", "phyloacc_gt_path", str,
        "The path to the PhyloAcc-GT binary. Default: PhyloAcc-GT"));
    #parser.add_argument("-gt-path", dest="phyloacc_gt_path", help="The path to the PhyloAcc-GT binary. Default: PhyloAcc-GT", default=False);
    
    arg_flags.update(addArgument(parser, "-iqtree-path", "iqtree_path", str,
        "The path to the IQ-Tree executable for making gene trees with --theta. Default: iqtree"));
    #parser.add_argument("-iqtree-path", dest="iqtree_path", help="The path to the IQ-Tree executable for making gene trees with --theta. Default: iqtree", default=False);
    
    arg_flags.update(addArgument(parser, "-coal-path", "coal_cmd", str,
        "The path to the program to estimate branch lengths in coalescent units with --theta (Supported programs: ASTRAL). Default: java -jar astral.jar"));
    #parser.add_argument("-coal-path", dest="coal_cmd", help="The path to the program to estimate branch lengths in coalescent units with --theta (Supported programs: ASTRAL). Default: java -jar astral.jar", default=False);
    # Dependency paths
    ## Note: For now we will likely need three dependency paths for the models, but eventually these should all be consolidated
    
    ##########

    arg_flags.update(addArgument(parser, "-scf", "scf_branch_cutoff", float,
        "The value of sCF to consider as 'low' for any given branch in a locus. Default: 0.5"));
    #parser.add_argument("-scf", dest="low_scf", help="The value of sCF to consider as 'low' for any given branch in a locus. Default: 0.5", type=float, default=0.5);

    arg_flags.update(addArgument(parser, "-s", "scf_prop", float,
        "The proportion of branches to consider a locus for the gene tree model. Default: 0.3333, meaning if one-third of all branches for a given locus have low sCF, this locus will be run with the gene tree model."));
    #parser.add_argument("-s", dest="scf_prop", help="The proportion of branches to consider a locus for the gene tree model. Default: 0.3333, meaning if one-third of all branches for a given locus have low sCF, this locus will be run with the gene tree model.", type=float, default=False);
    
    arg_flags.update(addArgument(parser, "-l", "coal_tree", str,
        "A file containing a rooted, Newick formatted tree with the same topology as the species tree in the mod file (-m), but with branch lengths in coalescent units. When the gene tree model is used, one of -l or --theta must be set."));
    #parser.add_argument("-l", dest="coal_tree", help="A file containing a rooted, Newick formatted tree with the same topology as the species tree in the mod file (-m), but with branch lengths in coalescent units. When the gene tree model is used, one of -l or --theta must be set.", default=False);
    
    arg_flags.update(addArgument(parser, "-r", "run_mode", str,
        "Determines which version of PhyloAcc will be used. gt: use the gene tree model for all loci, st: use the species tree model for all loci, adaptive: use the gene tree model on loci with many branches with low sCF and species tree model on all other loci. Default: st"));
    #parser.add_argument("-r", dest="run_mode", help="Determines which version of PhyloAcc will be used. gt: use the gene tree model for all loci, st: use the species tree model for all loci, adaptive: use the gene tree model on loci with many branches with low sCF and species tree model on all other loci. Default: st", default=False);
    
    arg_flags.update(addArgument(parser, "-n", "num_procs", int,
        "The number of processes that this script should use. Default: 1."));
    #parser.add_argument("-n", dest="num_procs", help="The number of processes that this script should use. Default: 1.", type=int, default=1);

    arg_flags.update(addArgument(parser, "-p", "procs_per_batch", int,
        "The number of processes to use for each batch of PhyloAcc. Default: 1."));
    #parser.add_argument("-p", dest="procs_per_batch", help="The number of processes to use for each batch of PhyloAcc. Default: 1.", type=int, default=1);

    arg_flags.update(addArgument(parser, "-j", "num_jobs", int,
        "The number of jobs (batches) to run in parallel. Must be less than or equal to the total processes for PhyloAcc (-p). Default: 1."));
    #parser.add_argument("-j", dest="num_jobs", help="The number of jobs (batches) to run in parallel. Must be less than or equal to the total processes for PhyloAcc (-p). Default: 1.", type=int, default=1);
    # User params

    ##########

    arg_flags.update(addArgument(parser, "-batch", "batch_size", int,
        "The number of loci to run per batch. Default: 50"));
    #parser.add_argument("-batch", dest="batch_size", help="The number of loci to run per batch. Default: 50", default=False);
    # Batch options

    ##########

    arg_flags.update(addArgument(parser, "-part", "cluster_part", str,
        "The partition or list of partitions (separated by commas) on which to run PhyloAcc jobs."));
    #parser.add_argument("-part", dest="cluster_part", help="The partition or list of partitions (separated by commas) on which to run PhyloAcc jobs.", default=False);
    
    arg_flags.update(addArgument(parser, "-nodes", "cluster_nodes", int,
        "The number of nodes on the specified partition to submit jobs to. Default: 1."));
    #parser.add_argument("-nodes", dest="cluster_nodes", help="The number of nodes on the specified partition to submit jobs to. Default: 1.", default=False);
    
    arg_flags.update(addArgument(parser, "-mem", "cluster_mem", int,
        "The max memory for each job in MB. Default: 4."));
    #parser.add_argument("-mem", dest="cluster_mem", help="The max memory for each job in GB. Default: 4.", default=False);
    
    arg_flags.update(addArgument(parser, "-time", "cluster_time", int,
        "The time in minutes to give each job. Default: 1."));
    #parser.add_argument("-time", dest="cluster_time", help="The time in hours to give each job. Default: 1.", default=False);

    arg_flags.update(addArgument(parser, "--local", "local_flag", bool,
        "Set to generate a local snakemake command to run the pipeline instead of the cluster profile and command. ONLY recommended for testing purposes."));
    # Cluster options
    
    ##########

    arg_flags.update(addArgument(parser, "--dollo", "dollo_flag", bool,
        "Set this to use the Dollo assumption in the original model (lineages descendant from an accelerated branch cannot change state)."));    
    #parser.add_argument("--dollo", dest="dollo_flag", help="Set this to use the Dollo assumption in the original model (lineages descendant from an accelerated branch cannot change state).", action="store_true", default=False);

    arg_flags.update(addArgument(parser, "--theta", "theta_flag", bool,
        "Set this to add gene tree estimation with IQ-tree and species estimation with ASTRAL for estimation of the theta prior. Note that a species tree with branch lengths in units of substitutions per site is still required with -m. Also note that this may add substantial runtime to the pipeline."));
    #parser.add_argument("--theta", dest="theta", help="Set this to add gene tree estimation with IQ-tree and species estimation with ASTRAL for estimation of the theta prior. Note that a species tree with branch lengths in units of substitutions per site is still required with -m. Also note that this may add substantial runtime to the pipeline.", action="store_true", default=False);
    
    arg_flags.update(addArgument(parser, "--labeltree", "labeltree", bool,
        "Reads the tree from the input mod file (-m), labels the internal nodes, and exits."));        
    #parser.add_argument("--labeltree", dest="labeltree", help="Simply reads the tree from the input mod file (-m), labels the internal nodes, and exits.", action="store_true", default=False);
    
    arg_flags.update(addArgument(parser, "--overwrite", "overwrite_flag", bool,
        "Set this to overwrite existing files."));       
    #parser.add_argument("--overwrite", dest="overwrite_flag", help="Set this to overwrite existing files.", action="store_true", default=False);
    
    arg_flags.update(addArgument(parser, "--appendlog", "append_log_flag", bool,
        "Set this to keep the old log file even if --overwrite is specified. New log information will instead be appended to the previous log file."));
    #parser.add_argument("--appendlog", dest="append_log_flag", help="Set this to keep the old log file even if --overwrite is specified. New log information will instead be appended to the previous log file.", action="store_true", default=False);
    # User options
    
    ##########

    #parser.add_argument("--plot", dest="plot_flag", help="Plot some summary statistics from the input data.", action="store_true", default=False);

    arg_flags.update(addArgument(parser, "--testcmd", "test_cmd_flag", bool,
        "At the end of the program, print out an example command to run a PhyloAcc batch directly."));

    arg_flags.update(addArgument(parser, "--summarize", "summarize_flag", bool,
        "Only generate the input summary plots and page. Do not write or overwrite batch job files."));
    #parser.add_argument("--summarize", dest="summarize_flag", help="Only generate the input summary plots and page. Do not write or overwrite batch job files.", action="store_true", default=False);

    arg_flags.update(addArgument(parser, "--options", "options_flag", bool,
        "Print the full list of PhyloAcc options that can be specified with -phyloacc and exit."));    
    #parser.add_argument("--options", dest="options_flag", help="Print the full list of PhyloAcc options that can be specified with -phyloacc and exit.", action="store_true", default=False);

    arg_flags.update(addArgument(parser, "--info", "info_flag", bool,
        "Print some meta information about the program and exit. No other options required."));   
    #parser.add_argument("--info", dest="info_flag", help="Print some meta information about the program and exit. No other options required.", action="store_true", default=False);
    
    arg_flags.update(addArgument(parser, "--depcheck", "depcheck", bool,
        "Run this to check that all dependencies are installed at the provided path. No other options necessary."));
    #parser.add_argument("--depcheck", dest="depcheck", help="Run this to check that all dependencies are installed at the provided path. No other options necessary.", action="store_true", default=False);
    
    arg_flags.update(addArgument(parser, ["--version", "-version", "--v", "-v"], "version_flag", bool,
        "Simply print the version and exit. Can also be called as '-version', '-v', or '--v'."));
    #parser.add_argument("--version", dest="version_flag", help="Simply print the version and exit. Can also be called as '-version', '-v', or '--v'", action="store_true", default=False);
    
    arg_flags.update(addArgument(parser, "--quiet", "quiet_flag", bool,
        "Set this flag to prevent PhyloAcc from reporting detailed information about each step."));
    #parser.add_argument("--quiet", dest="quiet_flag", help="Set this flag to prevent PhyloAcc from reporting detailed information about each step.", action="store_true", default=False);
    # Run options
    
    ##########

    parser.add_argument("-inf-frac-theta", dest="inf_frac_theta", help=argparse.SUPPRESS, type=float, default=0.2);
    # Controls the threshold for informative sites for loci to be used in the --theta estimation -- only those with a higher fraction will be used
    
    parser.add_argument("--qstats", dest="qstats", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--norun", dest="norun", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--nophyloacc", dest="no_phyloacc", help=argparse.SUPPRESS, action="store_true", default=False);
    # Generate the snakemake file but comment out the phyloacc rules, for dev
    
    parser.add_argument("-lt_source", dest="lt_source", help=argparse.SUPPRESS, type=str, default=False);
    parser.add_argument("-lt_input", dest="lt_input", help=argparse.SUPPRESS, type=str, default=False);
    parser.add_argument("-lt_output", dest="lt_output", help=argparse.SUPPRESS, type=str, default=False);
    parser.add_argument("--dev", dest="dev_opt", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--debugtree", dest="debug_tree", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--debugaln", dest="debug_aln", help=argparse.SUPPRESS, action="store_true", default=False);
    ## Dev options (hidden)
    
    args = parser.parse_args();
    # The input options and help messages

    ####################

    if any([args.lt_source, args.lt_input, args.lt_output]):
        if not all ([args.lt_source, args.lt_input, args.lt_output]):
            PC.errorOut("OP-LT", "All three options -lt_source, -lt_input, and -lt_output must be specified together.", globs);
        import phyloacc_lib.labeltree as LT
        LT.transferLabels(args.lt_source, args.lt_input, args.lt_output);
        sys.exit(0);

    ####################

    warnings = [];
    # List of warnings to print after logfile is created

    globs['call'] = " ".join(sys.argv);
    # Save the program call for later

    ####################

    if args.config_file:
        with open(args.config_file, 'r') as f:
            config = yaml.safe_load(f);
    else:
        config = {};

    ####################

    globs['version-flag'] = getOpt(args.version_flag, "version_flag", bool, False, config, arg_flags, globs);
    if globs['version-flag']:
        print("\n# PhyloAcc version " + globs['version'] + " released on " + globs['releasedate-patch'])
        sys.exit(0);
    # The version option to simply print the version and exit.        

    ####################

    globs['quiet'] = getOpt(args.quiet_flag, "quiet_flag", bool, False, config, arg_flags, globs);
    if not globs['quiet']:
        print("\n" + globs['call'] + "\n");
        print("#");
        print("# " + "=" * 125);
        print(PC.welcome());
        print("    Bayesian rate analysis of conserved");
        print("       non-coding genomic elements\n");
        # A welcome banner.
    else:
        globs['log-v'] = 0;
        print("# Running phyloacc_interface in quiet mode...");
    # Check the --quiet option and print a welcome banner if it isn't set.

    ####################

    globs['options-flag'] = getOpt(args.options_flag, "options_flag", bool, False, config, arg_flags, globs);
    if globs['options-flag']:
        PC.printOptions(globs);
    # Check if --options is set to print the PhyloAcc options and exit

    ####################

    globs['info'] = getOpt(args.info_flag, "info_flag", bool, globs['info'], config, arg_flags, globs);
    if globs['info']:
        globs['log-v'] = -1;
        startProg(globs);
        return globs;
    # Parse the --info option and call startProg early if set

    ####################

    if args.norun:
        globs['norun'] = True;
        globs['log-v'] = -1;
    # Check if norun is set

    if args.no_phyloacc:
        globs['no-phyloacc'] = True;
    # Check if --nophyloacc is set to prevent the PhyloAcc rules from being executed in 
    # the snakemake workflow (useful for debugging/testing --theta)

    if args.debug_aln:
        globs['debug-aln'] = True;
        globs['log-v'] = -1;

    if args.debug_tree:
        globs['debug-tree'] = True;
        globs['log-v'] = -1;
    # Check if tree debugging is turned on -- will just read the input species tree for testing in the tree library

    ## These options are all hidden and are not in the config file
    ####################

    #globs['overwrite'] = args.overwrite_flag;
    args.depcheck = getOpt(args.depcheck, "depcheck", bool, False, config, arg_flags, globs);
    # Check if --depcheck is set
    # It doesn't get stored in globs because it is only used here

    #globs['overwrite'] = args.overwrite_flag;
    globs['overwrite'] = getOpt(args.overwrite_flag, "overwrite_flag", bool, globs['overwrite'], config, arg_flags, globs);
    # Check if --overwrite is set

    #globs['label-tree'] = args.labeltree;
    globs['label-tree'] = getOpt(args.labeltree, "labeltree", bool, globs['label-tree'], config, arg_flags, globs);
    # Parse the --labeltree option

    globs['test-cmd-flag'] = getOpt(args.test_cmd_flag, "test_cmd_flag", bool, False, config, arg_flags, globs);
    # Parse the --testcmd option

    ####################

    if not args.depcheck and not globs['debug-aln']:
        globs['mod-file'] = getOpt(args.mod_file, "mod_file", "FILE", globs['mod-file'], config, arg_flags, globs);
        # if not globs['mod-file']:
        #     PC.errorOut("OP3", "A mod file must be provided with -m", globs);
        #globs['mod-file'] = args.mod_file;
    # Check the mod file

    # Species tree/rate file
    ####################

    if not globs['label-tree'] and not globs['debug-tree']:
    ## Read main input options if --labeltree isn't set

        globs['run-mode'] = getOpt(args.run_mode, "run_mode", ["st", "gt", "adaptive"], globs['run-mode'], config, arg_flags, globs);
        globs['coal-tree-file'] = getOpt(args.coal_tree, "coal_tree", "FILE", globs['coal-tree-file'], config, arg_flags, globs);
        globs['theta'] = getOpt(args.theta_flag, "theta_flag", bool, globs['theta'], config, arg_flags, globs);
        
        if globs['run-mode'] == "st" and (globs['theta'] or globs['coal-tree-file']):
            warnings.append("# WARNING: When using the species tree model '-r st', a tree with coalescent units is not required. -l and --theta will be ignored.");
            globs['theta'] = False;
            globs['coal-tree-file'] = None;
        if globs['run-mode'] != 'st':
            if (not globs['theta'] and not globs['coal-tree-file']) or (globs['theta'] and globs['coal-tree-file']):
                PC.errorOut("OP7", "When using the gene tree model with '-r gt' or '-r adaptive', a tree with branch lengths in coalescent units must also be provided with -l, or estimated with --theta.", globs);
        ## Run mode and coalescent tree/theta options
        ####################
        
        globs['phyloacc'] = getOpt(args.phyloacc_st_path, "phyloacc_st_path", str, globs['phyloacc'], config, arg_flags, globs, check=False);
        globs['phyloacc-gt'] = getOpt(args.phyloacc_gt_path, "phyloacc_gt_path", str, globs['phyloacc-gt'], config, arg_flags, globs, check=False);
        phyloacc_opt_input = getOpt(args.phyloacc_opts, "phyloacc_opts", str, None, config, arg_flags, globs, check=False);
        globs['iqtree-path'] = getOpt(args.iqtree_path, "iqtree_path", str, globs['iqtree-path'], config, arg_flags, globs, check=False);
        globs['coal-cmd'] = getOpt(args.coal_cmd, "coal_cmd", str, globs['coal-cmd'], config, arg_flags, globs, check=False);

        dep_check = getOpt(args.depcheck, "depcheck", bool, False, config, arg_flags, globs);

        globs, deps_passed = PC.execCheck(globs, dep_check, args.dev_opt);
        if dep_check:
            if deps_passed:
                print("\n# All dependencies PASSED.\n")
                sys.exit(0);
            else:
                print("\n# Some dependencies NOT FOUND. Please check your installations and provided paths.\n");
                sys.exit(1);
        # Check the dependency paths

        globs['aln-file'] = getOpt(args.aln_file, "aln_file", "FILE", globs['aln-file'], config, arg_flags, globs, check=False);
        globs['bed-file'] = getOpt(args.bed_file, "bed_file", "FILE", globs['bed-file'], config, arg_flags, globs, check=False);
        globs['id-file'] = getOpt(args.id_file, "id_file", "FILE", globs['id-file'], config, arg_flags, globs, check=False);
        globs['aln-dir'] = getOpt(args.aln_dir, "aln_dir", "DIR", globs['aln-dir'], config, arg_flags, globs, check=False);

        if (not globs['aln-file'] and not globs['aln-dir']) or (globs['aln-file'] and globs['aln-dir']):
            PC.errorOut("OP1", "Only one input method must be specified: -a or -d", globs);
        # Check that only one input type is specified

        if globs['aln-file'] and not globs['bed-file']:
            PC.errorOut("OP2", "A bed file with locus coordinates must be specified with -b when a concatenated alignment file is given with -a", globs);
        # Save the input type as a global param

        globs = PC.fileCheck(globs);
        # Make sure all the input files actually exist, and get their
        # full paths

        globs['filter-alns'] = getOpt(args.filter_alns, "filter_alns", bool, globs['filter-alns'], config, arg_flags, globs);
        # Alignment filtering option

        globs['softmask'] = getOpt(args.softmask, "softmask", bool, globs['softmask'], config, arg_flags, globs);
        # Soft-masked bases option

        # Input files
        ####################

        if not globs['debug-aln']:

            globs['dollo'] = getOpt(args.dollo_flag, "dollo_flag", bool, globs['dollo'], config, arg_flags, globs);
            ## Dollo option
            ####################

            if globs['theta']:
                job_sub_dirs = { 'job-alns' : 'alns', 'job-cfgs' : 'cfgs', 'job-bed' : 'bed', 'job-smk' : 'snakemake', 'job-out' : 'phyloacc-output', 'iqtree' : 'iqtree', 'astral' : 'astral' };
            else:
                job_sub_dirs = { 'job-alns' : 'alns', 'job-cfgs' : 'cfgs', 'job-bed' : 'bed', 'job-smk' : 'snakemake', 'job-out' : 'phyloacc-output' };
            # Expected job directories to create

            if globs['id-flag']:
                job_sub_dirs['job-ids'] = 'ids';

            ## Run mode and coalescent tree/theta options
            ####################

            globs['input-groups']['targets'] = getOpt(args.targets, "targets", str, globs['input-groups']['targets'], config, arg_flags, globs);
            globs['input-groups']['targets'] = globs['input-groups']['targets'].replace("; ", ";").split(";");

            globs['input-groups']['outgroup'] = getOpt(args.outgroup, "outgroup", str, globs['input-groups']['outgroup'], config, arg_flags, globs, check=False);
            if isinstance(globs['input-groups']['outgroup'], str):
                globs['input-groups']['outgroup'] = globs['input-groups']['outgroup'].replace("; ", ";").split(";");

            globs['input-groups']['conserved'] = getOpt(args.conserved, "conserved", str, globs['input-groups']['conserved'], config, arg_flags, globs, check=False);
            if isinstance(globs['input-groups']['conserved'], str):
                globs['input-groups']['conserved'] = globs['input-groups']['conserved'].replace("; ", ";").split(";");

            ## Read species group options
            ####################

            globs['min-scf'] = getOpt(args.scf_branch_cutoff, "scf_branch_cutoff", "PROP", globs['min-scf'], config, arg_flags, globs);
            globs['low-scf-branch-prop'] = getOpt(args.scf_prop, "scf_prop", "PROP", globs['low-scf-branch-prop'], config, arg_flags, globs);

            ## sCF options
            ####################

            globs['mcmc'] = getOpt(args.mcmc, "mcmc", int, globs['mcmc'], config, arg_flags, globs);
            globs['burnin'] = getOpt(args.burnin, "burnin", int, globs['burnin'], config, arg_flags, globs);
            globs['chain'] = getOpt(args.chain, "chain", int, globs['chain'], config, arg_flags, globs);
            globs['thin'] = getOpt(args.thin, "thin", int, globs['thin'], config, arg_flags, globs);

            ## MCMC options
            ####################
            
            if phyloacc_opt_input and phyloacc_opt_input != "None":
                while "; " in phyloacc_opt_input:
                    phyloacc_opt_input = phyloacc_opt_input.replace("; ", ";");
                phyloacc_opts = phyloacc_opt_input.split(";");
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

                    if globs['dollo'] and opt in ["INIT_LRATE2", "HYPER_LRATE2_A", "HYPER_LRATE2_B"]:
                        warnings.append("# WARNING: With --dollo set, LRATE2 parameters are not used. Ignoring your input for: " + opt);
                        continue;
                    # A warning in case the LRATE2 parameters are set when --dollo is assumed

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

            globs['outdir'] = getOpt(args.out_dir, "out_dir", "DIR", globs['outdir'], config, arg_flags, globs, check=False);

            if not globs['outdir']:
                globs['outdir'] = "phyloacc-out-" + globs['startdatetime'];
            # Get the output directory

            globs['outdir'] = os.path.abspath(globs['outdir']);

            if not globs['overwrite'] and os.path.exists(globs['outdir']):
                PC.errorOut("OP16", "Output directory already exists: " + globs['outdir'] + ". Specify new directory name OR set --overwrite to overwrite all files in that directory.", globs);
            if not os.path.isdir(globs['outdir']) and not globs['norun'] and not globs['info'] and not globs['debug-tree']:
                os.makedirs(globs['outdir']);
            # Main output dir

            ####################

            summarize_opt = getOpt(args.summarize_flag, "summarize_flag", bool, False, config, arg_flags, globs);
            if summarize_opt:
                globs['batch'] = False;
            else:
                globs['batch'] = True;
            #globs['batch'] = getOpt(args.summarize_flag, "summarize_flag", bool, globs['batch'], config, arg_flags, globs);

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
            elif globs['theta']:
                globs['label-coal-tree-script'] = os.path.join(globs['astral'], "label_astral_tree.sh");
                globs['coal-tree-input'] = os.path.join(globs['astral'], "input-species-tree.treefile");
                globs['coal-tree-file-unlabeled'] = os.path.join(globs['astral'], "astral-species-tree.treefile");
                globs['coal-tree-file'] = os.path.join(globs['astral'], "astral-species-tree-labeled.treefile");

            globs['run-name'] = os.path.basename(os.path.normpath(globs['outdir']));
            globs['logfilename'] = os.path.join(globs['outdir'], globs['run-name'] + ".log");
            # Log file

            append_log = getOpt(args.append_log_flag, "append_log_flag", bool, False, config, arg_flags, globs);
            if not append_log and not globs['norun'] and not globs['debug-tree']:
                logfile = open(globs['logfilename'], "w");
                logfile.write("");
                logfile.close();
            # Prep the logfile to be overwritten

            ## Output files and directories
            ####################

            globs['batch-size'] = getOpt(args.batch_size, "batch_size", int, globs['batch-size'], config, arg_flags, globs);
            # Batch size

            globs['procs-per-job'] = getOpt(args.procs_per_batch, "procs_per_batch", int, globs['procs-per-job'], config, arg_flags, globs);
            globs['num-jobs'] = getOpt(args.num_jobs, "num_jobs", int, globs['num-jobs'], config, arg_flags, globs);
            globs['total-procs'] = globs['procs-per-job'] * globs['num-jobs'];
            # Determine resource allocation for PhyloAcc

            globs['num-procs'] = getOpt(args.num_procs, "num_procs", int, globs['num-procs'], config, arg_flags, globs);
            #globs['num-procs'] = PC.isPosInt(args.num_procs, default=1);
            globs['aln-pool'] = mp.Pool(processes=globs['num-procs']);
            globs['scf-pool'] = mp.Pool(processes=globs['num-procs']);
            # Create the pool of processes for sCF calculation here so we copy the memory profile of the parent process
            # before we've read any large data in

            # Batch size and resource allocation
            ####################

            globs['local'] = getOpt(args.local_flag, "local_flag", bool, False, config, arg_flags, globs);
            if globs['local']:
                warnings.append("# WARNING: Using --local mode. Cluster options will be ignored and profile will not be generated. This is only recommended for testing purposes.");

            if not globs['local']:
                globs['partition'] = getOpt(args.cluster_part, "cluster_part", str, globs['partition'], config, arg_flags, globs, check=False);
                if not globs['partition']:
                    PC.errorOut("OP18", "At least one cluster partition must be specified with -part.", globs);
                # Cluster partition option (required)

                globs['num-nodes'] = str(getOpt(args.cluster_nodes, "cluster_nodes", int, globs['num-nodes'], config, arg_flags, globs));
                # Cluster node option

                globs['mem'] = str(getOpt(args.cluster_mem, "cluster_mem", int, globs['mem'], config, arg_flags, globs));
                #globs['mem'] = str(globs['mem'] * 1000);
                # Cluster memory option

                globs['time'] = str(getOpt(args.cluster_time, "cluster_time", int, globs['time'], config, arg_flags, globs));
                #globs['time'] = globs['time'] + ":00:00";
                # Cluster time option

            ## Cluster options
            ####################
        ## End --debugaln check
    ## End --labeltree check

    else:
        globs['log-v'] = -1;
    # If --labeltree or --debugtree is set, turn off log reporting

    ####################

    # globs = PC.fileCheck(globs);
    # # Make sure all the input file actually exist, and get their
    # # full paths

    ####################

    if args.qstats:
        globs['qstats'] = True;
    # Check for the internal quartet stats option to write to a file.

    if globs['psutil']:
        globs['pids'] = [psutil.Process(os.getpid())];
    # Get the starting process ids to calculate memory usage throughout.

    if not PC.isPosFloat(args.inf_frac_theta, maxval=1.0):
        PC.errorOut("OP22", "-inf-frac-theta must be a positive float, between 0 and 1.", globs);       
    if args.inf_frac_theta != 0.2:
        globs['inf-frac-theta'] = args.inf_frac_theta;
    # The threshold for informative sites for a locus to be used in --theta estimation
    # Useful for simulations which may have few informative sites

    #print(globs['inf-frac-theta']);
    ## Internal stuff
    ####################

    startProg(globs);
    # After all the essential options have been set, call the welcome function.

    if warnings:
        for warning in warnings:
            PC.printWrite(globs['logfilename'], globs['log-v'], warning);
            globs['warnings'] += 1; 
        PC.printWrite(globs['logfilename'], globs['log-v'], "#");
    # Print any warnings here if there were any before the logfile was created   

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

    if globs['debug-aln']:
        print("\n--debugaln SET. READING ALIGNMENTS AND EXITING.\n");
        return;
    # Check for the --debugaln option   

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
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# PhyloAcc-ST", pad) + globs['phyloacc']);
    
    if globs['run-mode'] != "st":
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# PhyloAcc-GT", pad) + globs['phyloacc-gt']);

        if globs['theta']:
            PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# IQTree", pad) + globs['iqtree-path']);
            # The path to IQ-tree

            PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Coalescent method", pad) + globs['coal-cmd']);                                        
        # The path to the coalescent species tree method (ASTRAL)
        
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

    if globs['local']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --local", pad) + str(globs['local']));
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Partition(s)", pad) + globs['partition']);
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Number of nodes", pad) + globs['num-nodes']);
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Max mem per job (mb)", pad) + globs['mem']);
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Time per (min)", pad) + globs['time']); 
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

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -gt:", pad) + 
                PC.spacedOut(str(globs['thin']), opt_pad) +
                "The number of MCMC steps between gene tree sampling for the gene tree model. -mcmc will be scaled up by this number.");   
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -burnin:", pad) + 
                PC.spacedOut(str(globs['burnin']), opt_pad) +
                "This number of steps in the chain will discarded as burnin (-burnin*-gt)");         
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -mcmc:", pad) + 
                PC.spacedOut(str(globs['mcmc']), opt_pad) +
                "The number of steps in each chain  (-burnin*-mcmc)");
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

    if globs['dollo']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --dollo", pad) +
                    PC.spacedOut("True", opt_pad) + 
                    "The Dollo assumption will be used in this run.");
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --dollo", pad) +
                    PC.spacedOut("False", opt_pad) + 
                    "he Dollo assumption will NOT be used in this run.");            
    # Report the dollo option.    

    ####################

    if globs['theta']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --theta", pad) +
                    PC.spacedOut("True", opt_pad) + 
                    "A species tree with branch lengths in coalescent units will be estimated.");
        # The path to the coalescent tree to be estimated
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

    if globs['filter-alns']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --filter", pad) +
                    PC.spacedOut("True", opt_pad) + 
                    "Low quality alignments will be filtered out before PhyloAcc is run.");
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --filter", pad) +
                    PC.spacedOut("False", opt_pad) + 
                    "All alignments with informative sites will be run with PhyloAcc.");       
    # Reporting the filter option.

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

    if globs['test-cmd-flag']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --testcmd", pad) + 
            PC.spacedOut("True", opt_pad) + 
            "Also printing single batch command.");

    if globs['inf-frac-theta'] != 0.2:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -inf-frac-theta", pad) + 
            PC.spacedOut(str(globs['inf-frac-theta']), opt_pad) + 
            "ONLY LOCI WITH A FRACTION OF INFORMATIVE SITES ABOVE THIS VALUE WILL BE USED TO ESTIMATE --theta.");

    if globs['no-phyloacc']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --nophyloacc", pad) + 
            PC.spacedOut("True", opt_pad) + 
            "PHYLOACC RULES WILL NOT BE EXECUTED IN SNAKEMAKE WORKFLOW.");

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
