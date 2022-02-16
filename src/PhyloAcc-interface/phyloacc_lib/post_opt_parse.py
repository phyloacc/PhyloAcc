#############################################################################
# Parsing and printing the options and meta-info for PhyloAcc.
# Much of the error checking is done here as well.
#############################################################################

import sys
import os
import math
import argparse
import phyloacc_lib.core as CORE
import phyloacc_lib.tree as TREE

#############################################################################

def inputPathCheck(path, path_type, req, opt, globs, glob_key):
# A function to check input files and directories

    if path_type == "file":
        checker = os.path.isfile;
    elif path_type == "dir":
        checker = os.path.isdir;
    # Get the right type of path to check

    if req and not path:
        CORE.errorOut("OCORE.", opt + " must be provided.", globs);
    # If the path is required and not provided, error out here

    if path:
        if not checker(path):
            CORE.errorOut("OCORE.", "Cannot find " + opt + ": " + path, globs);
        # If the path doesn't exist, error out here

        else:
            globs[glob_key] = os.path.abspath(path);
        # Otherwise add the provided path to the global params dict

    return globs;

#############################################################################

def optParse(globs):
# This function handles the command line options and prepares the output directory and files.
# Defaults are set in params.py

    # try:
    #     import psutil
    #     globs['psutil'] = True;
    # except:
    #     globs['psutil'] = False;
    # Check if psutil is installed for memory usage stats.

    parser = argparse.ArgumentParser(description="phyloacc_post.py combines and summarizes the output from batched PhyloAcc runs.");

    parser.add_argument("-i", dest="input_dir", help="The directory created from a phyloacc_interface run.", default=False);
    #parser.add_argument("-m", dest="mod_file", help="A file with a background transition rate matrix and phylogenetic tree with branch lengths as output from PHAST. REQUIRED.", default=False);
    # Input

    parser.add_argument("-o", dest="out_dest", help="Desired output directory. This will be created for you if it doesn't exist. Default: <input directory>/results/", default=False);
    # Output
    
    parser.add_argument("-n", dest="num_procs", help="The number of processes that this script should use. Default: 1.", type=int, default=1);
    # User params

    #parser.add_argument("--labeltree", dest="labeltree", help="Simply reads the tree from the input mod file (-m), labels the internal nodes, and exits.", action="store_true", default=False);
    parser.add_argument("--overwrite", dest="ow_flag", help="Set this to overwrite existing files.", action="store_true", default=False);
    # User options
    
    parser.add_argument("--plot", dest="plot_flag", help="Plot some summaries of the results.", action="store_true", default=False);
    parser.add_argument("--info", dest="info_flag", help="Print some meta information about the program and exit. No other options required.", action="store_true", default=False);
    #parser.add_argument("--dryrun", dest="dryrun", help="With all options provided, set this to run through the whole pseudo-it pipeline without executing external commands.", action="store_true", default=False);
    parser.add_argument("--version", dest="version_flag", help="Simply print the version and exit. Can also be called as '-version', '-v', or '--v'", action="store_true", default=False);
    parser.add_argument("--quiet", dest="quiet_flag", help="Set this flag to prevent PhyloAcc from reporting detailed information about each step.", action="store_true", default=False);
    # Run options
    
    parser.add_argument("--norun", dest="norun", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--debug", dest="debug_opt", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--nolog", dest="nolog_opt", help=argparse.SUPPRESS, action="store_true", default=False);
    # Performance tests
    
    args = parser.parse_args();
    # The input options and help messages

    ####################

    globs['call'] = " ".join(sys.argv);
    # Save the program call for later

    if args.info_flag:
        globs['info'] = True;
        globs['log-v'] = -1;
        startProg(globs);
        return globs;
    # Parse the --info option and call startProg early if set

    if args.norun:
        globs['norun'] = True;
        globs['log-v'] = -1;
    globs['overwrite'] = args.ow_flag;
    # Check run mode options.

    ## Run option checks
    ####################

    globs = inputPathCheck(args.input_dir, "dir", True, "-i", globs, "interface-run-dir");
    
    ## Interface input directory check
    ####################

    globs['run-name'] = os.path.basename(os.path.normpath(globs['interface-run-dir']));
    globs['logfilename'] = os.path.join(globs['interface-run-dir'], "final-results.log");

    # Log file
    ####################

    interface_logfiles = [ f for f in os.listdir(globs['interface-run-dir']) if f.endswith(".log") and f != os.path.basename(globs['logfilename']) ];
    if len(interface_logfiles) > 1:
        CORE.errorOut("OP9", "Found multiple logfiles in the input directory. Make sure there is only one coinciding with the run you want to analyze.", globs);
    if len(interface_logfiles) < 1:
        CORE.errorOut("OP9", "Cannot find interface logfile in the input directory.", globs);
    else:
        globs['interface-logfile'] = os.path.abspath(os.path.join(globs['interface-run-dir'], interface_logfiles[0]));

    ## Interface log file check
    ####################

    for line in open(globs['interface-logfile']):
        if line.startswith("# Tree read from mod file:"):
            globs['tree-string'] = line.strip().replace("# Tree read from mod file:", "");
            while globs['tree-string'][0] == " ":
                globs['tree-string'] = globs['tree-string'][1:];
        # Read the tree string from the MOD file

        if line.startswith("# Loci per batch (-batch)"):
            globs['batch-size'] = line.strip().replace("# Loci per batch (-batch)", "");
            globs['batch-size'] = list(filter(None, globs['batch-size'].split(" ")))[0];

        if line.startswith("# Processes per job (-p)"):
            globs['procs-per-batch'] = line.strip().replace("# Processes per job (-p)", "");
            globs['procs-per-batch'] = list(filter(None, globs['procs-per-batch'].split(" ")))[0];

    try:
        globs['tree-dict'], globs['labeled-tree'], globs['root-node'] = TREE.treeParse(globs['tree-string']);
        globs['tree-tips'] = [ n for n in globs['tree-dict'] if globs['tree-dict'][n][2] == "tip" ];
    except:
        CORE.errorOut("OP5", "Error reading tree from interface log file!", globs);
    # Read the tree from the interface log file

    ## Read info from interface log file
    ####################

    globs['phyloacc-out-dir'] = os.path.join(globs['interface-run-dir'], "phyloacc-job-files", "phyloacc-output");
    
    # The directory with all the PhyloAcc output files
    ####################

    if not args.out_dest:
        globs['outdir'] = os.path.join(globs['interface-run-dir'], "results");
    else:
        globs['outdir'] = args.out_dest;

    if not globs['overwrite'] and os.path.exists(globs['outdir']):
        CORE.errorOut("OP9", "Output directory already exists: " + globs['outdir'] + ". Specify new directory name OR set --overwrite to overwrite all files in that directory.", globs);

    if not os.path.isdir(globs['outdir']) and not globs['norun']:
        os.makedirs(globs['outdir']);
    # Main output dir

    ####################

    if args.plot_flag:
        globs['plot-dir'] = os.path.join(globs['interface-run-dir'], "plots");
        if not os.path.isdir(globs['plot-dir']) and not globs['norun']:
            os.makedirs(globs['plot-dir']);
        globs['plot'] = True;
        # Plot option

        globs['html-file'] = os.path.join(globs['interface-run-dir'], globs['html-file']);
        # HTML directory
    # Parse the --plot option

    ####################

    globs['num-procs'] = CORE.isPosInt(args.num_procs, default=1);
    # Num procs option

    ####################

    if args.quiet_flag:
        globs['quiet'] = True;
    # Checking the quiet flag

    ####################

    startProg(globs);
    return globs;

#############################################################################

def startProg(globs):
# A nice way to start the program.
    print("#");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# phyloacc_post.py combines and summarizes the output from batched PhyloAcc runs.");
    #CORE.printWrite(globs['logfilename'], globs['log-v'], "# Version (interface) " + globs['interface-version'] + " released on " + globs['releasedate']);
    #CORE.printWrite(globs['logfilename'], globs['log-v'], "# PhyloAcc was developed by Zhirui Hu, Han Yan, Taehee Lee, Gregg Thomas, Tim Sackton, Scott Edwards, and Jun Liu");
    #CORE.printWrite(globs['logfilename'], globs['log-v'], "# Citation:      " + globs['doi']);
    #CORE.printWrite(globs['logfilename'], globs['log-v'], "# Website:       " + globs['http']);
    #CORE.printWrite(globs['logfilename'], globs['log-v'], "# Report issues: " + globs['github']);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "#");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# The date and time at the start is: " + CORE.getDateTime());
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# Using Python version:              " + globs['pyver'] + "\n#");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# The program was called as:         " + globs['call'] + "\n#");

    if globs['info']:
        return;
    # If --info is set, return after printing program info

    #######################

    pad = 35;
    opt_pad = 30;
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 150);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# INPUT/OUTPUT INFO:");
    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# PhyloAcc interface dir:", pad) + globs['interface-run-dir']);
    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# PhyloAcc interface log file:", pad) + globs['interface-logfile']);
    #CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Tree/rate file (mod file from PHAST):", pad) + globs['mod-file']);
    #CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Tree read from mod file:", pad) + globs['tree-string']);
    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# PhyloAcc run directory:", pad) + globs['phyloacc-out-dir']);
    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Output directory:", pad) + globs['outdir']);
    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Log file:", pad) + os.path.basename(globs['logfilename']));
    # Input/Output
    #######################

    CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 150);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# OPTIONS INFO:");    
    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Option", pad) + CORE.spacedOut("Current setting", opt_pad) + "Current action");

    if globs['plot']:
        plot_status = "True";
        plot_status_str = "An HTML file summarizing the results will be written to " + globs['html-file'];
    else:
        plot_status = "False";
        plot_status_str = "No HTML summary file will be generated.";

    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --plot", pad) + 
                CORE.spacedOut(plot_status, opt_pad) + 
                plot_status_str);
    # --plot option

    ####################

    # CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Processes (-p)", pad) + 
    #             CORE.spacedOut(str(globs['num-procs']), opt_pad) + 
    #             "PhyloAcc and this interface will use this many total processes.");

    if globs['overwrite']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --overwrite", pad) +
                    CORE.spacedOut("True", opt_pad) + 
                    "PhyloAcc_post will OVERWRITE the existing files in the specified output directory.");
    # Reporting the overwrite option.

    ####################

    if not globs['quiet']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --quiet", pad) + 
                    CORE.spacedOut("False", opt_pad) + 
                    "Time, memory, and status info will be printed to the screen while PhyloAcc is running.");
    else:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --quiet", pad) + 
                    CORE.spacedOut("True", opt_pad) + 
                    "No further information will be printed to the screen while PhyloAcc is running.");
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 150);
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# Running...");
    # Reporting the quiet option.

    ####################

    # if globs['debug']:
    #     CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --debug", pad) + 
    #                 CORE.spacedOut("True", opt_pad) + 
    #                 "Printing out a bit of debug info.");
    # Reporting the debug option.

    if globs['norun']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --norun", pad) + 
                    CORE.spacedOut("True", opt_pad) + 
                    "ONLY PRINTING RUNTIME INFO.");
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 150);
    # Reporting the norun option.

    # Other options
    #######################

#############################################################################