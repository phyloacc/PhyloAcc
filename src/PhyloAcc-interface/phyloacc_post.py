#!/usr/bin/env python3
#############################################################################
# A script to combine and summarize the outputs of a batched PhyloAcc run 
# with the interface
#
# Gregg Thomas
# Fall 2021
#############################################################################

import sys
import os
import phyloacc_lib.core as CORE
import phyloacc_lib.post_params as params
import phyloacc_lib.post_opt_parse as OP
import phyloacc_lib.treeio as TREEIO
import phyloacc_lib.plot as PLOT
import phyloacc_lib.html as HTML

from collections import defaultdict

#############################################################################

if __name__ == '__main__':
# Main is necessary for multiprocessing to work on Windows.

    globs = params.init();
    # Get the global params as a dictionary.

    print("\n" + " ".join(sys.argv) + "\n");

    if any(v in sys.argv for v in ["--version", "-version", "--v", "-v"]):
        print("# PhyloAcc interface version " + globs['version'] + " released on " + globs['releasedate'])
        sys.exit(0);
    # The version option to simply print the version and exit.
    # Need to get actual PhyloAcc version for this, and not just the interface version.

    print("#");
    print("# " + "=" * 150);

    globs = OP.optParse(globs);
    # Getting the input parameters from optParse.

    if globs['info']:
        print("# --info SET. EXITING AFTER PRINTING PROGRAM INFO...\n#")
        sys.exit(0);
    if globs['norun']:
        print("# --norun SET. EXITING AFTER PRINTING OPTIONS INFO...\n#")
        sys.exit(0);
    # Early exit options

    step_start_time = CORE.report_step(globs, "", "", "", start=True);
    # Initialize the step headers

    ####################

    globs['tree-string'], globs['st'], globs['batch-size'], globs['procs-per-batch'] = TREEIO.readST(globs, tree_type="from-log");
    # Read the tree

    for node in globs['st'].nodes:
    # Loop over every node in the tree

        # globs['m2-per-lineage'][node] = 0;
        # # Initialize that node with a count of 0

        if globs['st'].type[node] == "tip":
            globs['m2-per-lineage'][node] = 0;
        else:
            if globs['st'].has_label:
                globs['m2-per-lineage'][globs['st'].label[node]] = 0;
            else:
                globs['m2-per-lineage'][node] = 0;

        # if globs['st'].has_label:
        #     if globs['st'].type[node] == "tip":
        #         globs['st-rev-labels'][node] = node;    
        #     else:
        #         globs['st-rev-labels'][globs['st'].label[node]] = node;
        # If the tree has its own internal labels, associate them with the treeparse label here
        # globs['st-rev-labels'] = <original tree label> : <treeparse tree label>
        # For tips, original and treeparse labels are the same

    # print(globs['st-rev-labels']);
    # sys.exit();

    # Initialize the dictionary for number of times a lineage is accelerated under M2
    # globs['m2-per-lineage'] = <tp label> : <count>

    ####################

    step = "Preparing output files";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status updated

    file_suffixes = ["_elem_lik.txt", "_M0_elem_Z.txt", "_M1_elem_Z.txt", "_M2_elem_Z.txt", "_rate_postZ_M0.txt", "_rate_postZ_M1.txt", "_rate_postZ_M2.txt"];
    data_headers = ["elem_lik", "M0_elem_Z", "M1_elem_Z", "M2_elem_Z", "rate_postZ_M0", "rate_postZ_M1", "rate_postZ_M2"];
    # The files created by the various PhyloAcc runs that we care about have these suffixes
    # _elem_lik.txt MUST be first to get the locus IDs and associated numbers for each batch
        
    outfiles = [];
    # Combined output files will be created and stored in this list

    for f in file_suffixes:
        cur_outfile = os.path.join(globs['results-dir'], f[1:]);
        open(cur_outfile, "w").close()
        outfiles.append(cur_outfile);
    # Create the combined output file for each suffix and add it to the list

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: files created");
    # Status update

    ####################

    step = "Gathering batch outputs";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status updated

    batch_dirs = os.listdir(globs['phyloacc-out-dir']);
    # The directories within the PhyloAcc job directory, 1 for each batch

    id_keys = {};
    # Most output files only have locus numbers assigned by PhyloAcc. This dict stores the association with the actual
    # locus ID found in the _elem_lik.txt file
    # key:value format: <batch> : { <locus number> : <locus id> }

    first_file = True;
    # The first file flag - must be _elem_lik.txt

    for f in range(len(file_suffixes)):
    ## Loop through all file types. Not sure if I mentioned this before, but _elem_lik.txt MUST be first

        first_batch = True;
        # A first batch flag - we only write the headers in the output file during the first batch

        cur_data = data_headers[f];
        globs['locus-stats'][cur_data] = defaultdict(dict);
        # Look up the current file name and start the dict for it

        with open(outfiles[f], "a") as outfile:
        ## Open the current output file

            for batch_dir in batch_dirs:
            ## Go through every batch for the current file

                if batch_dir.endswith("-st-out"):
                    batch_type = "st";
                elif batch_dir.endswith("-gt-out"):
                    batch_type = "gt";
                # Get the batch type based on the directory suffix

                batch = batch_dir[0:batch_dir.index("-")];
                # Parse the batch string (just a number)

                if first_file:
                    id_keys[batch] = {};
                # If this is the first file (_elem_lik.txt), initialize the sub-dict in id_keys

                full_batch_dir = os.path.join(globs['phyloacc-out-dir'], batch_dir);
                cur_batch_file = os.path.join(full_batch_dir, batch + file_suffixes[f]);
                # Get the full batch directory and current batch file

                if first_file:
                    if not os.path.isfile(cur_batch_file):
                        globs['incomplete-batches'].append(batch);
                        globs['incomplete-batches-' + batch_type].append(batch);
                    else:
                        globs['complete-batches'].append(batch);
                        globs['complete-batches-' + batch_type].append(batch);
                # Check if the file exists for the current batch and if not add to unfinished_batches
                # and skip

                if batch in globs['incomplete-batches']:
                    continue;
                # Skip the batch if the file is incomplete

                batch_logfile = os.path.join(globs['phyloacc-out-dir'], batch_dir, batch + "-phyloacc.log");
                for line in open(batch_logfile):
                    if line.startswith("time used:"):
                        line = list(filter(None, line.strip().split(" ")));
                        globs['batch-runtimes'].append(int(line[2]));
                # Lookup the runtime of the current batch

                first_line = True;
                # The first line flag - the first line of each file contains the headers

                for line in open(cur_batch_file):
                ## Loop through every line in the current batch file

                    line = line.strip().split("\t");
                    # Parse the line into a list

                    if first_line:
                    ## If this is the first line (headers), we need to do some stuff, or just skip

                        cur_headers = line;
                        if first_file:
                            cur_headers = cur_headers[1:];
                        # Get the headers, and if this is the first file (_elem_lik.txt), remove the "No." header

                        cur_headers[0] = "Locus ID";
                        # Replace the "No." header with the "Locus ID" header

                        if first_batch and not first_file:
                            outfile.write("\t".join(cur_headers) + "\n");
                        # Write out the headers for all files except elem_lik

                        first_line = False;
                        continue;
                        # Skip to the next line
                    ## End first line block

                    if first_file:
                        locus_id = batch + "-" + line[1];

                        id_keys[batch][line[0]] = locus_id;
                        globs['all-loci'].append(locus_id);
                        # Get the current locus ID based on the batch number and locus number within batch
                        
                        globs['locus-stats'][cur_data][locus_id] = defaultdict(dict);
                        # Initialize the dict for the current locus

                        line = [line[0]] + line[2:];
                        # Remove the ID column

                        for col in range(len(cur_headers)):
                            globs['locus-stats'][cur_data][locus_id][cur_headers[col]] = line[col];
                        # Add the data from the current locus for the current file to the dict

                        if batch_type == "st":
                            logbf3 = float(line[-6]) - float(line[-8]);
                        elif batch_type == "gt":
                            logbf3 = float(line[-3]) - float(line[-5]);
                        else:
                            logbf3 = "NA";
                        globs['locus-stats'][cur_data][locus_id]["logBF3"] = str(round(logbf3, 6));
                        # Add BF3 to the locus stats

                        globs['locus-stats'][cur_data][locus_id]["m0-lineages-cons"] = [];
                        globs['locus-stats'][cur_data][locus_id]["m0-lineages-acc"] = [];
                        globs['locus-stats'][cur_data][locus_id]["m1-lineages-cons"] = [];
                        globs['locus-stats'][cur_data][locus_id]["m1-lineages-acc"] = [];
                        globs['locus-stats'][cur_data][locus_id]["m2-lineages-cons"] = [];
                        globs['locus-stats'][cur_data][locus_id]["m2-lineages-acc"] = [];
                        # Add a column that will store a list of the accelerated lineages from the z-matrix

                        bf1, bf2, bf3 = False, False, False;

                        if float(globs['locus-stats']['elem_lik'][locus_id]["logBF1"]) > globs['bf1-cutoff']:
                            globs['bf1-loci'].append(locus_id);
                            bf1 = True;
                        # Check the BF1 cut-off

                        if float(globs['locus-stats']['elem_lik'][locus_id]["logBF2"]) > globs['bf2-cutoff']:
                            globs['bf2-loci'].append(locus_id);
                            bf2 = True;
                        # Check the BF2 cut-off

                        if float(globs['locus-stats']['elem_lik'][locus_id]["logBF3"]) > globs['bf3-cutoff']:
                            globs['bf3-loci'].append(locus_id);
                            bf3 = True;
                        # Check the BF3 cut-off

                        if bf1 and bf2:
                            globs['m1-loci'].append(locus_id);
                        elif bf3:
                            globs['m2-loci'].append(locus_id);
                        else:
                            globs['m0-loci'].append(locus_id);
                        # Assess model fit based on BFs

                    else:
                        locus_id = id_keys[batch][line[0]];
                        outline = [locus_id] + line[1:];
                        # For all other files, replace the locus number with the locus ID based on the entry in id_keys

                        if "_elem_Z.txt" in cur_batch_file:
                        # Parse the Z matrix files to count accelerated lineages for each model

                            cur_model = file_suffixes[f][1:3];

                            z_scores = line[1:];
                            if batch_type == "gt":
                                z_scores = z_scores[:-1];
                            # Parse the z-scores for the current line
                            # For gene tree models, remove the gene tree at the end of the line

                            num_accel = 0;

                            for z in range(len(z_scores)):
                            # Loop over every z score
                                cur_label = cur_headers[z+1];
                                # lookup the lineage label from the headers

                                if z_scores[z] == "1":
                                    globs['locus-stats']["elem_lik"][locus_id][cur_model.lower() + "-lineages-cons"].append(cur_label);
                                ## End conserved block
                                if z_scores[z] == "2":
                                # If the current z score is accelerated (2)

                                    num_accel += 1;                                    
                                    # Increment the number accelerated for this locus and 

                                    # if globs['st'].has_label:
                                    #     cur_label = globs['st-rev-labels'][cur_label];
                                    # If the original tree had labels, switch to that label here

                                    if cur_model == "M2":
                                        globs['m2-per-lineage'][cur_label] += 1;
                                    # Increment the count of acclerated elements for this branch/label

                                    globs['locus-stats']["elem_lik"][locus_id][cur_model.lower() + "-lineages-acc"].append(cur_label);        
                                ## End accelerated block
                            ## End z score loop
                            
                            if cur_model == "M2":
                                globs['m2-per-locus'][num_accel] += 1;
                            # Increment the bin for the number of branches accelerated per locus 
                        ## End Z matrix file block

                        # if "_rate_postZ_" in cur_batch_file:
                        #      print(cur_batch_file);
                        #      sys.exit();
                                            

                        outfile.write("\t".join(outline) + "\n");
                        # Write the line to the output file

                        for col in range(len(cur_headers)):
                            globs['locus-stats'][cur_data][locus_id][cur_headers[col]] = line[col];
                        # Add the data from the current locus for the current file to the dict
                    ## End other file block
                ## End line loop

                first_batch = False;
                # Switch off the batch flag
            ## End batch_dir loop

            first_file = False;
            # Switch off the file flag
        ## Close batch output file
    ## End file loop

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: files combined");
    # Status update

    # print();
    # print(globs['locus-stats']['rate_postZ_M0']['50-4'].keys());
    # sys.exit();

    # print();
    # print(len(globs['all-loci']), len(globs['m0-loci']), len(globs['m1-loci']), len(globs['m2-loci']));
    # print(sum([len(globs['m0-loci']), len(globs['m1-loci']), len(globs['m2-loci'])]));
    #sys.exit();
    ####################

    if globs['incomplete-batches']:
        globs['incomplete-batches'] = list(set(globs['incomplete-batches']));
        num_unfinished = str(len(globs['incomplete-batches']));
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# WARNING: " + num_unfinished + " batches are unfinished: " + ",".join(globs['incomplete-batches']));
    # Print warnings for incomplete batches

    ####################

    step = "Writing main output file (elem_lik.txt)";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status updated

    id_key_file = os.path.join(globs['results-dir'], "elem_lik.txt");
    with open(id_key_file, "w") as idout:
    # The batch IDs assigned by the interface are different from the originals, so this block
    # associates them and writes them to a file

        idout.write("# PhyloAcc element likelihoods and ID key\n");
        idout.write("# SUMMARY: " + str(len(globs['locus-stats']['elem_lik'])) + " loci\n");
        idout.write("# SUMMARY: " + str(len(globs['complete-batches'])) + " batches completed\n");
        idout.write("# SUMMARY: " + str(len(globs['incomplete-batches'])) + " batches failed\n");
        idout.write("# SUMMARY: " + str(len(globs['m0-loci'])) + " loci best fit M0\n");
        idout.write("# SUMMARY: " + str(len(globs['m1-loci'])) + " loci best fit M1\n");
        idout.write("# SUMMARY: " + str(len(globs['m2-loci'])) + " loci best fit M2\n");
        idout.write("#\n");
        idout.write("# COLUMN DEFINITIONS\n")
        idout.write("# phyloacc.id:               The ID for the locus as assigned by PhyloAcc\n");
        idout.write("# original.id:               The ID provided in the input bed file\n");
        idout.write("# best.fit.model:            The model (M0, M1, M2) with the best fit based on Bayes Factors\n");
        idout.write("# marginal.likelihood.m0:    The marginal likelihood from M0\n");
        idout.write("# marginal.likelihood.m1:    The marginal likelihood from M1\n");
        idout.write("# marginal.likelihood.m2:    The marginal likelihood from M2\n");
        idout.write("# logbf1:                    The Bayes factor between M1 and M0 (logbf1 = marginal.likelihood.m1 - marginal.likelihood.m0, specified cutoff was " + str(globs['bf1-cutoff']) +")\n");
        idout.write("# logbf2:                    The Bayes factor between M1 and M2 (logbf2 = marginal.likelihood.m1 - marginal.likelihood.m2, specified cutoff was " + str(globs['bf2-cutoff']) +")\n");
        idout.write("# logbf3:                    The Bayes factor between M2 and M1 (logbf3 = marginal.likelihood.m2 - marginal.likelihood.m0, specified cutoff was " + str(globs['bf3-cutoff']) +")\n");
        idout.write("# conserved.rate.m0:         The posterior median of the conserved substitution rate under M0\n");
        idout.write("# accel.rate.m0:             The posterior median of the accelerated substitution rate under M0\n");          
        idout.write("# conserved.rate.m1:         The posterior median of the conserved substitution rate under M1\n");
        idout.write("# accel.rate.m1:             The posterior median of the accelerated substitution rate under M1\n");  
        idout.write("# conserved.rate.m2:         The posterior median of the conserved substitution rate under M2\n");
        idout.write("# accel.rate.m2:             The posterior median of the accelerated substitution rate under M2\n");    
        idout.write("# num.accel.m1:              The number of lineages inferred to be accelerated under M1\n");        
        idout.write("# num.accel.m2:              The number of lineages inferred to be accelerated under M2\n");        
        idout.write("# accel.lineages.m1:         A comma separated list of the accelerated lineages under M1\n");
        idout.write("# conserved.lineages.m2:     A comma separated list of the conserved lineages under M2\n");
        idout.write("# accel.lineages.m2:         A comma separated list of the accelerated lineages under M2\n");
        idout.write("#\n");        
        idout.write("# NOTE: M0 = null model   (no acceleration, classifed when BFs are below cutoffs for M1 and M2)\n");
        idout.write("# NOTE: M1 = target model (acceleration on target lineages only, classifed when logbf1 and logbf2 > cutoffs)\n");
        idout.write("# NOTE: M2 = full model   (acceleration on any branch, classifed when BFs are below cutoffs for M1 and logbf3 > cutoff)\n");
        idout.write("#\n");
        if globs['st'].has_label:
            idout.write("# TREE: " + globs['tree-string'] + "\n");
        else:
            idout.write("# TREE: " + globs['st'].tree_str + "\n");            
        idout.write("#\n");  
        id_headers = ["phyloacc.id", "original.id", "best.fit.model", "marginal.likelihood.m0", "marginal.likelihood.m1", "marginal.likelihood.m2", 
                        "logbf1", "logbf2", "logbf3", 
                        "conserved.rate.m0", "accel.rate.m0", "conserved.rate.m1", "accel.rate.m1", "conserved.rate.m2", "accel.rate.m2", 
                        "num.accel.m1", "num.accel.m2", "conserved.lineages.m1", "accel.lineages.m1", "conserved.lineages.m2", "accel.lineages.m2"];
        idout.write("\t".join(id_headers) + "\n");
        # Write some headers and explanatory info at the top of the file

        for batch_dir in batch_dirs:
        ## Go through every batch

            batch_dir = batch_dir.split("-");
            if batch_dir[0] in globs['incomplete-batches']:
                continue;
            batch = batch_dir[0] + "-" + batch_dir[2];
            # Parse the batch string and skip if the batch was incomplete

            batch_bed_file = os.path.join(globs['interface-run-dir'], "phyloacc-job-files", "bed", batch + ".bed");
            # Look up the bed file for the current batch

            for line in open(batch_bed_file):
            # Read over every line in the bed file to get the original IDs

                line = line.strip().split("\t");
                locus_id = batch[0:batch.index("-")] + "-" + line[0];
                # Parse the line and get the interface ID (original ID is third column)

                if "gt" in batch:
                    cols = ['loglik_Null_W', 'loglik_Acc_W', 'loglik_Full_W'];
                elif "st" in batch:
                    cols = ['loglik_Null', 'loglik_Acc', 'loglik_Full'];
                cols += ['logBF1', 'logBF2', 'logBF3'];
                # Get the columns to parse from the elem_lik output

                best_model = 'M0';
                if locus_id in globs['m1-loci']:
                    best_model = 'M1';
                elif locus_id in globs['m2-loci']:
                    best_model = 'M2';
                # Look up the best model for the current locus

                #cols += [best_model.lower() + "-lineages"];

                outline = [locus_id, line[3], best_model];
                for col in cols:
                    outline.append(globs['locus-stats']['elem_lik'][locus_id][col]);
                # Compile the output list based on the desired columns

                for model in ["M0", "M1", "M2"]:
                    post_z_key = "rate_postZ_" + model;
                    for rate in ["c_rate", "n_rate"]:
                        if rate in globs['locus-stats'][post_z_key][locus_id]:
                            outline.append(str(globs['locus-stats'][post_z_key][locus_id][rate]));
                        else:
                            outline.append("NA");

                # outline.append(str(globs['locus-stats']['rate_postZ_M0'][locus_id]['c_rate']));
                # outline.append(str(globs['locus-stats']['rate_postZ_M0'][locus_id]['n_rate']));
                # outline.append(str(globs['locus-stats']['rate_postZ_M1'][locus_id]['c_rate']));
                # outline.append(str(globs['locus-stats']['rate_postZ_M1'][locus_id]['n_rate']));
                # outline.append(str(globs['locus-stats']['rate_postZ_M2'][locus_id]['c_rate']));
                # outline.append(str(globs['locus-stats']['rate_postZ_M2'][locus_id]['n_rate']));

                for model in ["m1", "m2"]:
                    outline.append(str(len(globs['locus-stats']['elem_lik'][locus_id][model + "-lineages-acc"])));

                outline.append(",".join(globs['locus-stats']['elem_lik'][locus_id]["m1-lineages-cons"]));
                outline.append(",".join(globs['locus-stats']['elem_lik'][locus_id]["m1-lineages-acc"]));
                outline.append(",".join(globs['locus-stats']['elem_lik'][locus_id]["m2-lineages-cons"]));
                outline.append(",".join(globs['locus-stats']['elem_lik'][locus_id]["m2-lineages-acc"]));

                idout.write("\t".join(outline) + "\n");
                # Write the output to the id key file:

            ## End bed loop
        ## End batch_dir loop
    ## Close ID key file

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # Status update
    #sys.exit();
    ####################

    if globs['plot']:
        globs = PLOT.genPlotsPost(globs);
        globs = HTML.writeHTMLPost(globs);
    # Generate plots and write HTML file

    ####################

    CORE.endProg(globs, interface=False);

#############################################################################