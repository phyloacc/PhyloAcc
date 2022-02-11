#############################################################################
# This file holds some global variables for some of the input options.
# These global parameters should be read only -- they are not modified anywhere 
# else in the code except when reading the input options.
#############################################################################

import sys
import timeit
import phyloacc_lib.core as CORE

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
        'startdatetime' : CORE.getOutTime(),
        'startdatetimenice' : CORE.getRunTimeNice(),
        # Meta info

        'pyver' :  ".".join(map(str, sys.version_info[:3])),
        # System info

        'call' : "",
        # Script call info

        'interface-run-dir' : False,
        'interface-logfile' : False,
        'phyloacc-out-dir' : False,
        'id-file' : False,
        # Input with concatenated alignment and partitions by bed file

        'bf1-cutoff' : 5,
        'bf2-cutoff' : 5,

        'locus-stats' : {},
        'accelerated-loci' : [],

        'complete-batches' : [],
        'complete-batches-st' : [],
        'complete-batches-gt' : [],
        'incomplete-batches' : [],
        'incomplete-batches-st' : [],
        'incomplete-batches-gt' : [],
        'complete-loci' : [],
        'accelerated-loci' : [],
        'batch-size' : 0,
        'procs-per-batch' : 0,
        'batch-runtimes' : [],

        'outdir' : '',
        'run-name' : 'phyloacc-post',
        'logfilename' : 'phyloacc-post.errlog',
        'tmpdir' : 'System default.',
        'overwrite' : False,
        # I/O options
        
        'plot' : False,
        # Option to output plots/html

        'tree-string' : False,
        'tree-dict' : False,
        'labeled-tree' : False,
        'root-node' : False,
        'tree-tips' : False,
        # Tree variables

        'plot-dir' : '',
        'input-tree-plot-file' : 'input-species-tree.png',
        'bf1-dist-file' : 'bf1-hist.png',
        'bf2-dist-file' : 'bf2-hist.png',
        'bf1-bf2-plot-file' : 'bf1-v-bf2.png',
        'inf-sites-frac-plot-file' : 'informative-sites-frac-hist.png',
        'var-inf-sites-plot-file' : 'variable-informative-sites.png',
        'avg-scf-hist-file' : 'avg-scf-per-locus.png',
        'low-scf-hist-file' : 'perc-low-scf-branchers-per-locus.png',
        'scf-tree-plot-file' : 'scf-species-tree.png', 
        'bl-scf-plot-file' : 'bl-scf.png', 
        'html-file' : 'phyloacc-results.html',
        # Plot and HTML summary files

        'html-summary-written' : False,
        # Output checks

        'num-procs' : 1,
        # Number of procs for this script to use

        'info' : False,
        'dryrun' : False,
        'quiet' : False,
        # Other user options

        'skip-chars' : ["-", "N"],
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

    globs_init['logfilename'] = "phyloacc-post-" + globs_init['startdatetime'] + ".errlog";
    # Add the runtime to the error log file name.

    globs = StrictDict(globs_init);
    # Restrict the dict from having keys added to it after this

    return globs;

#############################################################################