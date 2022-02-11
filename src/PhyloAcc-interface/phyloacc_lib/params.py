#############################################################################
# This file holds some global variables for some of the input options.
# These global parameters should be read only -- they are not modified anywhere 
# else in the code except when reading the input options.
#############################################################################

import sys
import timeit
import phyloacc_lib.core as PC

#############################################################################

class StrictDict(dict):
# This prevents additional keys from being added to the global params dict in
# any other part of the code, just to help me limit it's scope
# https://stackoverflow.com/questions/32258706/how-to-prevent-key-creation-through-dkey-val
    def __setitem__(self, key, value):
        if key not in self:
            raise KeyError("{} is not a legal key of this StrictDict".format(repr(key)));
        dict.__setitem__(self, key, value);

#############################################################################

def init():
    globs_init = {
        'phyloacc-versions' : '',
        'interface-version' : 'Beta 1.0',
        'releasedate' : 'NA',
        'doi' : 'https://doi.org/10.1093/molbev/msz049',
        'http' : 'https://xyz111131.github.io/PhyloAcc/',
        'github' : 'https://github.com/xyz111131/PhyloAcc/issues',
        'starttime' : timeit.default_timer(),
        'startdatetime' : PC.getOutTime(),
        'startdatetimenice' : PC.getRunTimeNice(),
        # Meta info

        'pyver' :  ".".join(map(str, sys.version_info[:3])),
        # System info

        'call' : "",
        'script-dir' : "",
        # Script call info

        'aln-file' : False,
        'bed-file' : False,
        'id-file' : False,
        # Input with concatenated alignment and partitions by bed file

        'aln-dir' : False,
        # Input by a directory full of individual alignments

        'mod-file' : False,
        # Input rate and tree file from PHAST

        'seq-compression' : 'none',
        'bed-compression' : 'none',
        # The type of compression used for input sequence files

        'outdir' : '',
        'run-name' : 'phyloacc',
        'logfilename' : 'phyloacc.errlog',
        'alnstatsfile' : 'phyloacc-aln-stats.csv',
        'no-inf-loci-file' : 'no-informative-sites-loci.txt',
        'scfstatsfile' : 'phyloacc-scf-stats.csv',
        'scftreefile' : 'phyloacc-scf.tree',
        'logdir' : '',
        'tmpdir' : 'System default.',
        'overwrite' : False,
        # I/O options
        
        'run-mode' : 'st',
        # Run mode option

        'plot' : False,
        # Option to output plots/html

        'theta' : False,
        'coal-tree-file' : False,
        'coal-tree-str' : False,
        # Theta estimation

        'tree-string' : False,
        'tree-dict' : False,
        'labeled-tree' : False,
        'scf-labeled-tree' : False,
        'root-node' : False,
        'tree-tips' : False,
        # Tree variables

        'in-seqs' : {},
        'in-bed' : {},
        'locus-ids' : [],
        'alns' : {},
        'aln-stats' : {},
        'num-loci' : False,
        # Sequence variables

        'avg-aln-len' : "NA",
        'med-aln-len' : "NA",
        'avg-nogap-seq-len' : "NA",
        'med-nogap-seq-len' : "NA",
        'no-inf-sites-loci' : [],
        'gappy-loci' : [],
        # Alignment summary stats

        'targets' : [],
        'conserved' : [],
        'outgroup' : [],
        # Phylo options

        'scf' : {},
        # Tracks sCF values BY NODE over all loci

        'quartets' : {},
        # All qurtets for all nodes/branches in the species tree:
        # <node> : [ <quartets> ]
        # <quartet> is a tuple of tuples:
        # ((split1-spec1, split1-spec2), (split2-spec1, split2-spec2))

        'min-scf' : 0.5,
        # sCF params

        'burnin' : 500,
        'mcmc' : 1000,
        'chain' : 1,
        # MCMC options

        #'phyloacc-defaults-file' : "phyloacc-opts.txt",
        'phyloacc-defaults' : {},
        'phyloacc-opts' : [],
        # All other PhyloAcc options as a list

        'phyloacc' : 'PhyloAcc',
        #'phyloacc-gbgc' : 'PhyloAcc/V2_GBGC/PhyloAcc_gBGC',
        #'phyloacc-gt' : 'PhyloAcc-GT2/SRC/PhyloAcc-GT_piQ',
        'phyloacc-gt' : 'testSRC_debug4_2tree/PhyloAcc-GT_piQ',
        # Dependency paths

        'batch-size' : 50,
        'num-batches' : 0,
        'st-loci' : 0,
        'gt-loci' : 0,
        'st-batches' : [],
        'gt-batches' : [],
        # Batch variables

        'num-procs' : 1,
        # Number of procs for this script to use
        
        'num-jobs' : 1000,
        'procs-per-job' : 1,
        'total-procs' : 1,
        # Number of jobs/procs for PhyloAcc to use

        'partition' : False,
        'num-nodes' : "1",
        'mem' : "4",
        'time' : "1:00:00",
        # Cluster options

        'aln-pool' : False,
        'scf-pool' : False,
        # Process pools

        'smk' : False,
        'smk-config' : False,
        # Job files

        'id-flag' : True,

        'iqtree' : '',
        'astral' : '',
        'job-dir' : '',
        'job-alns' : '',
        'job-cfgs' : '',
        'job-bed' : '',
        'job-ids' : '',
        'job-smk' : '',
        'job-out' : '',
        'profile-dir' : False,
        # Job directories

        'plot-dir' : '',
        'input-tree-plot-file' : 'input-species-tree.png',
        'aln-len-plot-file' : 'aln-len-hist.png',
        'seq-len-plot-file' : 'avg-seq-len-nogap-hist.png',
        'inf-sites-plot-file' : 'informative-sites-hist.png',
        'inf-sites-frac-plot-file' : 'informative-sites-frac-hist.png',
        'var-inf-sites-plot-file' : 'variable-informative-sites.png',
        'avg-scf-hist-file' : 'avg-scf-per-locus.png',
        'low-scf-hist-file' : 'perc-low-scf-branchers-per-locus.png',
        'scf-tree-plot-file' : 'scf-species-tree.png', 
        'bl-scf-plot-file' : 'bl-scf.png', 
        'html-dir' : '',
        'html-file' : '',
        # Plot and HTML summary files

        'status-script' : 'slurm_status.py',
        'smk-cmd' : '',
        # The final snakemake command to report

        'label-tree' : False,
        'info' : False,
        'dryrun' : False,
        'quiet' : False,
        # Other user options

        'skip-chars' : ["-", "N"],
        'aln-stats-written' : False,
        'scf-stats-written' : False,
        'scf-tree-written' : False,
        'html-summary-written' : False,

        'pad' : 82,
        'endprog' : False,
        'exit-code' : 0,
        'log-v' : 1,
        'stats' : True,
        'progstarttime' : 0,
        'stepstarttime' : 0,
        'pids' : "",
        'psutil' : False,
        'qstats' : False,
        'norun' : False,
        'debug' : False,
        'nolog' : False,
        # Internal stuff
    }

    globs_init['logfilename'] = "phyloacc-" + globs_init['startdatetime'] + ".errlog";
    # Add the runtime to the error log file name.

    globs_init['phyloacc-defaults'] = {
        'SIMULATE'            : { 'type' : "BOOL",      'default' : "FALSE" },
        'SEED'                : { 'type' : "POS_INT",   'default' : "1" },
        'INIT_CONSERVE_RATE'  : { 'type' : "POS_FLOAT", 'default' : "0.5" },
        'INIT_ACC_RATE'       : { 'type' : "POS_FLOAT", 'default' : "1" },
        'CONSERVE_PRIOR_A'    : { 'type' : "POS_FLOAT", 'default' : "5" },
        'CONSERVE_PRIOR_B'    : { 'type' : "POS_FLOAT", 'default' : "0.04" },
        'ACCE_PRIOR_A'        : { 'type' : "POS_FLOAT", 'default' : "10" },
        'ACCE_PRIOR_B'        : { 'type' : "POS_FLOAT", 'default' : "0.2" },
        'ROPT'                : { 'type' : "POS_FLOAT", 'default' : "1" },
        'CUB'                 : { 'type' : "POS_FLOAT", 'default' : "1" },
        'NLB'                 : { 'type' : "POS_FLOAT", 'default' : "0.6" },
        'THIN'                : { 'type' : "POS_INT",   'default' : "1" },
        'INIT_LRATE'          : { 'type' : "PROB",      'default' : "0.8" },
        'INIT_GRATE'          : { 'type' : "PROB",      'default' : "0.5" },
        'HYPER_LRATE_A'       : { 'type' : "PROB",      'default' : "1" },
        'HYPER_LRATE_B'       : { 'type' : "PROB",      'default' : "1" },
        'HYPER_GRATE_A'       : { 'type' : "PROB",      'default' : "1" },
        'HYPER_GRATE_B'       : { 'type' : "PROB",      'default' : "1" },
        'WL'                  : { 'type' : "BOOL",      'default' : "FALSE" },
        'BL_WL'               : { 'type' : "POS_INT",   'default' : "15" },
        'CONSERVE_PROP'       : { 'type' : "PERC",      'default' : "0.8" },
        'CONSERVE_RATE'       : { 'type' : "PROP",      'default' : "NA" },
        'GAP_PROP'            : { 'type' : "PERC",      'default' : "0.8" },
        'CONSTOMIS'           : { 'type' : "PROB",      'default' : "0.5" },
        'BR_SAMPLE_THRESHOLD' : { 'type' : "POS_FLOAT", 'default' : "10" },
        'GAPCHAR'             : { 'type' : "CHAR",      'default' : "-" },
        'PRUNE_TREE'          : { 'type' : "BOOL",      'default' : "FALSE" },
        'TRIM_GAP_PERCENT'    : { 'type' : "PERC",      'default' : "0.9" },
        'MIN_LEN'             : { 'type' : "POS_INT",   'default' : "50" },
        'INDEL'               : { 'type' : "POS_INT",   'default' : "0" },
        'SAMPLE_HYPER'        : { 'type' : "POS_INT",   'default' : "0" },
        'VERBOSE'             : { 'type' : "POS_INT",   'default' : "FALSE" },
        'NUM_THREAD'          : { 'type' : "POS_INT",   'default' : "1" },
        'THETA_CUTOFF'        : { 'type' : "POS_FLOAT", 'default' : "NA" }
    }
    # The default values for PhyloAcc options NOT included in the interface options

    # for line in open(globs_init['phyloacc-defaults-file']):
    #     line = line.strip().split("\t");
    #     globs_init['phyloacc-defaults'][line[0]] = {'type' : line[2], 'default' : line[1]};

    globs = StrictDict(globs_init);
    # Restrict the dict from having keys added to it after this

    return globs;

#############################################################################