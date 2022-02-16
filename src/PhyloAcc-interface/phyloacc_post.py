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
import phyloacc_lib.plot as PLOT

from collections import defaultdict

#############################################################################

if __name__ == '__main__':
# Main is necessary for multiprocessing to work on Windows.

    globs = params.init();
    # Get the global params as a dictionary.

    print("\n" + " ".join(sys.argv) + "\n");

    if any(v in sys.argv for v in ["--version", "-version", "--v", "-v"]):
        print("# PhyloAcc interface version " + globs['interface-version'] + " released on " + globs['releasedate'])
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

    ####################

    step_start_time = CORE.report_step(globs, "", "", "", start=True);
    # Initialize the step headers

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
        cur_outfile = os.path.join(globs['outdir'], f[1:]);
        open(cur_outfile, "w").close()
        outfiles.append(cur_outfile);
    # Create the combined output file for each suffix and add it to the list

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: files created");
    # Status update

    ####################

    step = "Combining batch outputs";
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

        with open(outfiles[f], "a") as outfile:
        ## Open the current output file

            for batch_dir in batch_dirs:
            ## Go through every batch for the current file

                if batch_dir.endswith("-st-out"):
                    batch_type = "st";
                elif batch_dir.endswith("-gt-out"):
                    batch_type = "gt";

                batch = batch_dir[0:batch_dir.index("-")];
                # Parse the batch string (just a number)

                if first_file:
                    id_keys[batch] = {};
                # If this is the first file (_elem_lik.txt), initialize the sub-dict in id_keys

                full_batch_dir = os.path.join(globs['phyloacc-out-dir'], batch_dir);
                cur_batch_file = os.path.join(full_batch_dir, batch + file_suffixes[f]);
                # Get the full batch directory and current batch file

                if file_suffixes[f] == "_elem_lik.txt":
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
                        # Replace the "No." header with the "Locus ID header"

                        if first_batch:
                            outfile.write("\t".join(cur_headers) + "\n");
                            # Write out the headers

                        # if first_batch:
                        # ## If this is also the first batch, we need to write the headers to the output file
        
                        #     if first_file:
                        #         cur_headers = cur_headers[1:];
                        #     # Get the headers, and if this is the first file (_elem_lik.txt), remove the "No." header

                        #     cur_headers[0] = "Locus ID"
                        #     # Replace the "No." header with the "Locus ID header"


                        # ## End first batch block

                        first_line = False;
                        continue;
                        # Skip to the next line
                    ## End first line block

                    if first_file:
                        cur_id_key = batch + "-" + line[1];
                        id_keys[batch][line[0]] = cur_id_key;
                        globs['locus-stats'][cur_data][cur_id_key] = defaultdict(dict);
                        line = line[:1] + line[2:];
                        outline = [cur_id_key] + line[1:];
                        
                    # If this is the first file (_elem_lik.txt), save the locus number:locus id in the id_keys dict and
                    # remove the "No." entry from the line
                    else:
                        cur_id_key = id_keys[batch][line[0]];
                        #print(cur_id_key);
                        outline = [cur_id_key] + line[1:];
                    # For all other files, replace the locus number with the locus ID based on the entry in id_keys

                    #print(batch);
                    #print(cur_batch_file);
                    #print(id_keys[batch][line[0]]);
                    #print(cur_headers);
                    #print(line);

                    for col in range(len(cur_headers)):
                        globs['locus-stats'][cur_data][cur_id_key][cur_headers[col]] = line[col];

                    #globs['locus-stats'][id_keys[batch][line[0]]] = { cur_headers[col] : line[col] for col in range(len(cur_headers)) };

                    #print(globs['locus-stats'][id_keys[batch][line[0]]]);
                    #print("-----------------")
                    #sys.exit();

                    outfile.write("\t".join(outline) + "\n");
                    # Write the line to the output file
                ## End line loop

                first_batch = False;
                # Switch off the batch flag
            ## End batch loop

            first_file = False;
            # Switch off the file flag
        ## Close batch output file
    ## End file loop

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: files combined");
    # Status update

    print();

    print(globs['locus-stats']['elem_lik']['12-0']);
    print(globs['locus-stats']['elem_lik']['44-0']);
    #sys.exit();
    if globs['incomplete-batches']:
        globs['incomplete-batches'] = list(set(globs['incomplete-batches']));
        num_unfinished = str(len(globs['incomplete-batches']));
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# WARNING: " + num_unfinished + " batches are unfinished: " + ",".join(globs['incomplete-batches']));

    ####################

    step = "Getting original IDs";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status updated

    id_key_file = os.path.join(globs['outdir'], "id-key.txt");
    with open(id_key_file, "w") as idout:
        for batch_dir in batch_dirs:
        ## Go through every batch for the current file

            #print(batch_dir);

            batch_dir = batch_dir.split("-");

            if batch_dir[0] in globs['incomplete-batches']:
                continue;

            batch = batch_dir[0] + "-" + batch_dir[2];
            # Parse the batch string

            batch_bed_file = os.path.join(globs['interface-run-dir'], "phyloacc-job-files", "bed", batch + ".bed");

            # print(batch);
            # print(globs['locus-stats']['M0_elem_Z'][batch].keys());
            # print(globs['locus-stats']['M1_elem_Z'][batch].keys());
            # print(globs['locus-stats']['M2_elem_Z'][batch].keys());

            for line in open(batch_bed_file):
                line = line.strip().split("\t");
                locus_id = batch[0:batch.index("-")] + "-" + line[0];

                #models = ['loglik_Null_W', 'loglik_Acc_W', 'loglik_Full_W'];
                models = ['loglik_Null', 'loglik_Acc', 'loglik_Full'];
                marginal_liks = [ float(globs['locus-stats']['elem_lik'][locus_id][model]) for model in models ];

                max_lik = max(marginal_liks);
                max_lik_ind = marginal_liks.index(max_lik);
                max_model = models[max_lik_ind];
                elem_str = "M" + str(max_lik_ind) + "_elem_Z";

                if batch == "44-gt":
                    print(max_lik, max_lik_ind, max_model, elem_str);

                outline = [ locus_id, line[3], max_model, str(max_lik)] #, globs['locus-stats'][elem_str][locus_id]['genetree'] ];
                idout.write("\t".join(outline) + "\n");

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # Status update

    ####################

    step = "Getting batch runtimes";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status updated

    for batch_dir in batch_dirs:
    ## Go through every batch for the current file

        batch = batch_dir[0:batch_dir.index("-")];
        # Parse the batch string (just a number)

        if batch in globs['incomplete-batches']:
            continue;

        # print(batch_dir);
        # print(batch);

        batch_logfile = os.path.join(globs['phyloacc-out-dir'], batch_dir, batch + "-phyloacc.log");
        if not os.path.isfile(batch_logfile):
            print(batch_logfile)

        for line in open(batch_logfile):
            if line.startswith("time used:"):
                line = list(filter(None, line.strip().split(" ")));
                globs['batch-runtimes'].append(int(line[2]));

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # Status update

    ####################

    if globs['plot']:
        globs = PLOT.genPlotsPost(globs);
        globs = PLOT.writeHTMLPost(globs);

    ####################

    CORE.endProg(globs, interface=False);

#############################################################################