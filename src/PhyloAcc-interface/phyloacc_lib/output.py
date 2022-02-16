#############################################################################
# Output functions for the PhyloAcc interface
# Gregg Thomas
#############################################################################

import os
import phyloacc_lib.core as PC
import phyloacc_lib.tree as TREE

#############################################################################

def writeAlnStats(globs):
# Writes alignment stats to output file

    step = "Writing: " + globs['alnstatsfile'];
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status updated

    globs['alnstatsfile'] = os.path.join(globs['outdir'], globs['alnstatsfile']);
    # Updates the alignment stats file to include the output directory

    loci_sorted = sorted(globs['aln-stats']);
    # Sorts the loci

    with open(globs['alnstatsfile'], "w") as outfile:
    # Open the alignment stats file for writing

        first = True;
        # First flace

        for locus in loci_sorted:
        # Write every locus

            if first:
            # For the first locus, extract and write the headers

                keys = list(globs['aln-stats'][locus].keys())
                # The headers will be the keys in the locus dict

                headers = ['locus'] + keys;
                # Add 'locus' to the headers for the locus ID

                outfile.write(",".join(headers) + "\n");
                first = False;
                # Write the headers and set the first flag to False
            ###

            outline = [locus] + [ str(globs['aln-stats'][locus][key]) for key in keys ];
            outfile.write(",".join(outline) + "\n");
            # Extract the stats for the current locus and write to the file
        ## End locus loop
    ## Close file

    if globs['no-inf-sites-loci']:
        globs['no-inf-loci-file'] = os.path.join(globs['outdir'], globs['no-inf-loci-file']);
        with open(globs['no-inf-loci-file'], "w") as outfile:
            for locus in globs['no-inf-sites-loci']:
                outfile.write(locus + "\n");
    # Write the loci with no informative sites to a file

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: align stats written");
    globs['aln-stats-written'] = True;
    # Status update

    return globs;

#############################################################################

def writeSCFStats(globs):

    step = "Writing: " + globs['scfstatsfile'];
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status update

    headers = ["node","variable-sites","decisive-sites","concordant-sites","quartet-scf-sum","num-quartets","total-quartets","avg-quartet-scf"];
    # The headers to include in the output, must be keys from the globs['scf'] dict

    globs['scfstatsfile'] = os.path.join(globs['outdir'], globs['scfstatsfile']);
    # Add the output directory to the scfstatsfile name

    with open(globs['scfstatsfile'], "w") as outfile:
    # Open the scf stats file for writing

        outfile.write(",".join(headers) + "\n");
        # Write the headers to the output

        for node in globs['scf']:
        # Write stats for every node
            
            globs['scf'][node]['num-quartets'] = len(globs['quartets'][node])
            # Retrieve the number of quartets for this node from the global quartets dict

            outline = [node] + [ str(globs['scf'][node][header]) for header in headers if header != "node" ];
            outfile.write(",".join(outline) + "\n");
            # Extract the stats for the current node and write to the file
        ## End node loop
    ## Close file

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: sCF stats written");
    globs['scf-stats-written'] = True;
    # Status update    

    ####################

    step = "Writing: " + globs['scftreefile'];
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status update

    globs['scftreefile'] = os.path.join(globs['outdir'], globs['scftreefile']);
    # Add the output directory to the scftreefile name

    globs['scf-labeled-tree'] = globs['labeled-tree'];
    # Get the labeled input tree from treeParse

    globs['scf-labeled-tree'] = TREE.addBranchLength(globs['scf-labeled-tree'], globs['tree-dict']);
    # Add the branche lengths back onto the tree

    for node in globs['scf']:
        globs['scf-labeled-tree'] = globs['scf-labeled-tree'].replace(node, node + "_" + str(round(globs['scf'][node]['avg-quartet-scf'], 2)));
    # For every node in the tree, add the averaged scf value over all loci to the label

    with open(globs['scftreefile'], "w") as outfile:
        outfile.write(globs['scf-labeled-tree']);
    # Write the scf labeled tree to a file

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: sCF tree written");
    globs['scf-tree-written'] = True;   
    # Status update

    return globs;

#############################################################################