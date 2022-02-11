#############################################################################
# Functions to generate files for batch jobs of PhyloAcc (with snakemake)
# Gregg Thomas
#############################################################################

import os
import phyloacc_lib.core as PC
import phyloacc_lib.tree as TREE
import phyloacc_lib.templates_2 as TEMPLATES

#############################################################################

def genJobFiles(globs):
# This function generates all files for phyloacc to run on a per locus basis with snakemake

    if globs['theta']:
        step = "Writing alignments for gene trees";
        step_start_time = PC.report_step(globs, step, False, "In progress...");
        # Status update

        iqtree_aln_dir = os.path.join(globs['iqtree'], "alns");
        os.makedirs(iqtree_aln_dir, exist_ok=True);
        iqtree_gt_dir = os.path.join(globs['iqtree'], "gene-trees");
        os.makedirs(iqtree_gt_dir, exist_ok=True);
        # Make the iqtree subdirectories

        alns_written = 0;
        # A counter for the number of alignments written

        len_sorted_alns = [ aln for aln, value in sorted(globs['aln-stats'].items(), key=lambda item: item[1]['length'], reverse=True) ];
        # Sort the alignments by length

        for aln in len_sorted_alns:
        # Go over every locus

            if globs['aln-stats'][aln]['length'] > 100 and (globs['aln-stats'][aln]['informative-sites'] / globs['aln-stats'][aln]['length']) >= 0.2:
            # Check if the current locus is long enough and has enough informative sites to make a tree

                if globs['aln-stats'][aln]['low-qual']:
                    continue;
                # Skip alignments that are low quality

                cur_aln_file = os.path.join(iqtree_aln_dir, aln + ".fa");
                # The filename for the current alignment

                with open(cur_aln_file, "w") as alnfile:
                    for spec in globs['alns'][aln]:
                        alnfile.write(">" + spec + "\n");
                        alnfile.write(globs['alns'][aln][spec] + "\n");
                alns_written += 1;
                # Write the alignment
            ## End alignment writing block

            if alns_written > 5000:
                break;
            # Don't need to make gene trees of every element, just enough to estimate branch lengths
        ## End locus alignment loop

        step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + str(alns_written) + " alignments written");
        # Status update

    ####################

        step = "Writing species tree";
        step_start_time = PC.report_step(globs, step, False, "In progress...");
        # Status update

        input_tree_no_bl = TREE.remBranchLength(globs['tree-string']);
        # Remove branch lengths from input tree

        species_tree_file = os.path.join(globs['astral'], "input-species-tree.treefile");
        with open(species_tree_file, "w") as treefile:
            treefile.write(input_tree_no_bl);
        # Write the input tree to its own file for input tp ASTRAL
        
        step_start_time = PC.report_step(globs, step, step_start_time, "Success");
        # Status update        
    ## End --theta block

    ####################

    step = "Writing PhyloAcc job files";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status update

    st_alns = { aln : globs['alns'][aln] for aln in globs['alns'] if globs['aln-stats'][aln]['batch-type'] == "st" and aln not in globs['no-inf-sites-loci'] };
    gt_alns = { aln : globs['alns'][aln] for aln in globs['alns'] if globs['aln-stats'][aln]['batch-type'] == "gt" and aln not in globs['no-inf-sites-loci'] };
    # Subset alignments based on model type. This is either specified by the user with '-r st' or '-r gt', or is decided
    # by sCF with '-r adaptive'.

    globs['st-loci'] = len(st_alns);
    globs['gt-loci'] = len(gt_alns);
    # Adjust the counts after removing alignments with no informative sites

    model_num = 0;
    batch_num = 0;
    # Trakers for the model number and batch number
    # model_num is just the index for the list [st_alns, gt_alns]

    for model_partition in [st_alns, gt_alns]:
    # Go over the partitions for both models

        if model_num == 0:
            model_type = "st";
        elif model_num == 1:
            model_type = "gt";
        model_num += 1;
        # Decide the model type string based in the model_num

        for batch in PC.dictChunk(model_partition, globs['batch-size']):
        # Split the sequence dictionary by batches and yield the result here

            batch_num += 1;
            batch_num_str = str(batch_num);
            if model_type == "st":
                globs['st-batches'].append(batch_num_str);
                coal_tree_line = "";
            elif model_type == "gt":
                globs['gt-batches'].append(batch_num_str);
                coal_tree_line = "\nTREE_IN_COALESCENT_UNIT " + globs['coal-tree-file'];
            # Batch counting

            cur_out_dir = os.path.join(globs['job-out'], batch_num_str + "-phyloacc-" + model_type + "-out");
            if not os.path.isdir(cur_out_dir):
                os.makedirs(cur_out_dir);
            # Make the phyloacc output directory for the current batch

            batch_concat = { label : "" for label in globs['tree-tips'] };
            # Dictionary to concatenate alignments for the current batch

            batch_aln_list = [];
            # As we concatenate, add the alignment name to this list so we can go through
            # in the same order as we make the bed file

            for aln in batch:
                #print(aln);
                batch_aln_list.append(aln);
                for label in globs['tree-tips']:
                    batch_concat[label] += batch[aln][label];
                # May need a check here for missing sequences
            # Concatenate alignments in the current batch together

            cur_aln_file = os.path.join(globs['job-alns'], batch_num_str + "-" + model_type + ".fa");
            with open(cur_aln_file, "w") as alnfile:
                for spec in globs['tree-tips']:
                    alnfile.write(">" + spec + "\n");
                    alnfile.write(batch_concat[spec] + "\n");
            # Write the current concatenated alignment to file

            len_sum = 0;
            # Keeps track of the last coordinate written in the bed file

            cur_bed_file = os.path.join(globs['job-bed'], batch_num_str + "-" + model_type + ".bed");
            with open(cur_bed_file, "w") as bedfile:
                batch_aln_id = 0;
                ## NOTE: Right phyloacc requires element IDs to be integers starting from 0. I think this should be changed.
                for aln in batch_aln_list:

                    aln_len = globs['aln-stats'][aln]['length'];
                    # Get the alignment length for the bed file

                    end_coord = len_sum + aln_len;
                    # Get the end coordinate of the current locus in the concatenated alignment by
                    # adding the length to the previous length sum

                    outline = [str(batch_aln_id), str(len_sum), str(end_coord), aln];
                    bedfile.write("\t".join(outline) + "\n");
                    # Write the info for the current locus

                    len_sum += aln_len;

                    batch_aln_id += 1;
                    # Add the length to the length sum as the start coordinate for the next locus
            # Write a bed file that contains all coordinates for the current alignment

            if globs['id-flag']:
                cur_id_file = os.path.join(globs['job-ids'], batch_num_str + "-" + model_type + ".id");
                cur_id_file = os.path.abspath(cur_id_file);
                with open(cur_id_file, "w") as idfile:
                    batch_aln_id = 0;
                    for aln in batch_aln_list:
                        idfile.write(str(batch_aln_id) + "\n");
                        batch_aln_id += 1;
            # Write an ID file

            phyloacc_opt_str = "";
            if globs['phyloacc-opts']:
                phyloacc_opt_str = "\n".join(globs['phyloacc-opts']);

            id_line_str = "";
            if globs['id-flag']:
                id_line_str = "\nID_FILE " + cur_id_file;

            cur_cfg_file = os.path.join(globs['job-cfgs'], batch_num_str + "-" + model_type + ".cfg");
            with open(cur_cfg_file, "w") as cfgfile:
                cfgfile.write(TEMPLATES.phyloaccConfig().format(mod_file=os.path.abspath(globs['mod-file']),
                                                                        bed_file=os.path.abspath(cur_bed_file),
                                                                        id_line=id_line_str,
                                                                        aln_file=os.path.abspath(cur_aln_file),
                                                                        coal_tree_line=coal_tree_line,
                                                                        outdir=os.path.abspath(cur_out_dir),
                                                                        batch=batch_num_str,
                                                                        burnin=str(globs['burnin']),
                                                                        mcmc=str(globs['mcmc']),
                                                                        chain=str(globs['chain']),
                                                                        targets= ";".join(globs['targets']),
                                                                        outgroup=";".join(globs['outgroup']),
                                                                        conserved=";".join(globs['conserved']),
                                                                        procs_per_job=str(globs['procs-per-job']),
                                                                        phyloacc_opts=phyloacc_opt_str
                                                                        ))
            # Write the phyloacc config file for the current concatenated alignment
        ## End batch loop
    ## End model partition loop

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + str(batch_num) + " jobs written");
    # Status update

    globs['num-batches'] = batch_num;
    # Can only compute this here after we've split the number of input alignments by batch size

    return(globs);

############################################################################# 

def writeSnakemake(globs):
# A function to write the various snakemake files

    step = "Writing Snakemake file";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status update
    
    globs['smk'] = os.path.join(globs['job-smk'], "run_phyloacc.smk");

    if globs['theta']:
        run_char = "";
    else:
        run_char = "#";
    # If --theta isn't set, we don't want to run some of the rules in the snakemake file, so we comment out their
    # expected outputs in rule all here

    with open(globs['smk'], "w") as smkfile:
        smkfile.write(TEMPLATES.snakemake().format(cmd=globs['call'],
                                                       dt=PC.getDateTime(),
                                                       run_char=run_char,
                                                       st_path=os.path.abspath(globs['phyloacc']),
                                                       gt_path=os.path.abspath(globs['phyloacc-gt']),
                                                       coal_tree_path=globs['coal-tree-file']
                                                       ))

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: Snakemake file written");
    # Status update  

    ####################

    step = "Writing Snakemake config file";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status update

    globs['smk-config'] = os.path.join(globs['job-smk'], "phyloacc-config.yaml");

    with open(globs['smk-config'], "w") as configfile:
        configfile.write(TEMPLATES.snakemakeConfig().format(indir=os.path.abspath(globs['job-cfgs']),
                                                            outdir=os.path.abspath(globs['job-out']),
                                                            st_batches=str(globs['st-batches']),
                                                            gt_batches=str(globs['gt-batches']),
                                                            iqtree=os.path.abspath(globs['iqtree']),
                                                            astral=os.path.abspath(globs['astral'])
                                                            ))

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: Snakemake config written");
    # Status update  

    ####################

    step = "Writing Snakemake cluster profile";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status update

    globs['profile-dir'] = os.path.join(globs['job-smk'], "profiles", "slurm_profile");
    if not os.path.isdir(globs['profile-dir']):
        os.makedirs(globs['profile-dir']);

    profile_file = os.path.join(globs['profile-dir'], "config.yaml");
    # This is a snakemake profile for SLURM
    ## TODO: Will likely need templates for different job schedulers

    cluster_logdir = os.path.abspath(os.path.join(globs['job-dir'], "slurm-logs"));
    # A directory to save the log files from the cluster

    with open(profile_file, "w") as profile:
        profile.write(TEMPLATES.snakemakeProfile().format(num_jobs=str(globs['num-jobs']),
                                                          cluster_logdir=cluster_logdir,
                                                          procs_per_job=str(globs['procs-per-job']),
                                                          part=globs['partition'],
                                                          num_nodes=globs['num-nodes'],
                                                          mem=globs['mem'],
                                                          time=globs['time']
                                                        ))

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: Snakemake profile written");
    # Status update     

    return globs; 

#############################################################################