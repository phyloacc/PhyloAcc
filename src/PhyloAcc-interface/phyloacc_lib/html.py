#############################################################################
# Functions to generate html files for easy visualization of 
# phyloacc datasets
# Gregg Thomas
#############################################################################

import os
import phyloacc_lib.core as PC
import phyloacc_lib.templates as TEMPLATES
import phyloacc_lib.templates_post as TEMPLATES_POST

#############################################################################

def writeHTML(globs):
    step = "Writing HTML summary file";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status updated

    # if not os.path.isdir(globs['html-dir']):
    #     shutil.copytree( os.path.join(globs['script-dir'], "html"), globs['html-dir'] );
    # Copy the associated html files (stylesheets, images) from the provided template folders

    if globs['tree-data-type'] == 'class':
        num_spec = globs['st'].num_tips;
        tips = globs['st'].tips;
        internals = globs['st'].internals;
    elif globs['tree-data-type'] == 'func':
        num_spec = len(globs['tips']);
        tips = globs['tips'];
        internals = globs['internals'];

    if globs['run-mode'] == 'adaptive':
        comment_start = "";
        comment_end = "";
    else:
        comment_start = "<!-- This block is only used when run mode (-r) is adaptive";
        comment_end = "-->";
    # Some HTML blocks will be commented out depending on the run mode

    if globs['theta']:
        theta_comment_start = "";
        theta_comment_end = "";
    else:
        theta_comment_start = "<!-- This block is only displayed when --theta is specified";
        theta_comment_end = "-->";        

    if globs['coal-tree-file']:
        coal_tree_comment_start = "";
        coal_tree_comment_end = "";
    else:
        coal_tree_comment_start = "<!-- This block is only displayed when -l is specified";
        coal_tree_comment_end = "-->";            


    if globs['batch']:
        batch_comment_start = "<!-- This block is only displayed with the --plotonly option";
        batch_comment_end = "-->";
    else:
        batch_comment_start = "";
        batch_comment_end = "";

    with open(globs['html-file'], "w") as htmlfile:
        htmlfile.write(TEMPLATES.htmlSummary().format(
            # mod_file=os.path.abspath(globs['mod-file']),
            run_name=globs['run-name'],
            run_time=globs['startdatetimenice'],
            host_name=os.uname()[1],
            script_call=globs['call'],
            num_aln=str(globs['num-loci']),
            num_no_inf_loci=str(len(globs['no-inf-sites-loci'])),
            num_st_loci=str(globs['st-loci']),
            num_gt_loci=str(globs['gt-loci']),
            num_batches=str(globs['num-batches']),
            batch_size=str(globs['batch-size']),
            procs_per_batch=str(globs['procs-per-job']),
            num_st_batches=str(len(globs['st-batches'])),
            num_gt_batches=str(len(globs['gt-batches'])),
            num_jobs=str(globs['num-jobs']),
            # num_spec=str(len(globs['st'].tips)),
            # num_spec=str(len(globs['tips'])),
            num_spec=str(num_spec),
            num_targets=str(len(globs['groups']['targets'])),
            num_conserved=str(len(globs['groups']['conserved'])),
            num_outgroups=str(len(globs['groups']['outgroup'])),
            log_file=globs['logfilename'],
            aln_stats_file=globs['alnstatsfile'],
            batch_comment_start=batch_comment_start,
            batch_comment_end=batch_comment_end,
            snakemake_cmd=globs['smk-cmd'],
            theta_comment_start=theta_comment_start,
            theta_comment_end=theta_comment_end,
            coal_tree_comment_start=coal_tree_comment_start,
            coal_tree_file=globs['coal-tree-file'],
            coal_tree_comment_end=coal_tree_comment_end,
            input_tree_plot=os.path.join("plots", globs['input-tree-plot-file']),
            avg_aln_len=str(round(globs['avg-aln-len'], 3)),
            median_aln_len=str(round(globs['med-aln-len'], 3)),
            avg_seq_len_nogap=str(round(globs['avg-nogap-seq-len'], 3)),
            med_seq_len_nogap=str(round(globs['med-nogap-seq-len'], 3)),
            aln_len_hist=os.path.join("plots", globs['aln-len-plot-file']),
            seq_len_hist=os.path.join("plots", globs['seq-len-plot-file']),
            informative_sites_hist=os.path.join("plots", globs['inf-sites-plot-file']),
            informative_sites_frac_hist=os.path.join("plots", globs['inf-sites-frac-plot-file']),
            variable_informative_sites_plot=os.path.join("plots", globs['var-inf-sites-plot-file']),
            comment_start=comment_start,
            comment_end=comment_end,
            avg_scf_hist=os.path.join("plots", globs['avg-scf-hist-file']),
            low_scf_hist=os.path.join("plots", globs['low-scf-hist-file']),
            scf_tree_plot=os.path.join("plots", globs['scf-tree-plot-file']),
            bl_scf_plot=os.path.join("plots", globs['bl-scf-plot-file']),
            date_time=PC.getFooterDateTime(),
        ))
    # Write the HTML summary file using the template

    step_start_time = PC.report_step(globs, step, step_start_time, "Success");
    globs['html-summary-written'] = True;
    # Status update

    return globs;

#############################################################################

def writeHTMLPost(globs):
    step = "Writing HTML summary file";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status updated

    # if not os.path.isdir(globs['html-dir']):
    #     shutil.copytree( os.path.join(globs['script-dir'], "html"), globs['html-dir'] );
    # Copy the associated html files (stylesheets, images) from the provided template folders

    # if globs['run-mode'] == 'adaptive':
    #     comment_start = "";
    #     comment_end = "";
    # else:
    #     comment_start = "<!-- This block is only used when run mode (-r) is adaptive";
    #     comment_end = "-->";
    # # Some HTML blocks will be commented out depending on the run mode

    # if globs['theta']:
    #     theta_comment_start = "";
    #     theta_comment_end = "";
    # else:
    #     theta_comment_start = "<!-- This block is only displayed when --theta is specified";
    #     theta_comment_end = "-->";        

    # if globs['coal-tree-file']:
    #     coal_tree_comment_start = "";
    #     coal_tree_comment_end = "";
    # else:
    #     coal_tree_comment_start = "<!-- This block is only displayed when -l is specified";
    #     coal_tree_comment_end = "-->";       
     
    comment_start = "<!--";
    comment_end = "-->"     

    if globs['incomplete-batches']:
        batch_comment_start = "";
        batch_comment_end = "";
    else:
        batch_comment_start = "<!-- This block is only displayed if incomplete batches were detected";
        batch_comment_end = "-->";

    placeholder = "";

    with open(globs['html-file'], "w") as htmlfile:
        htmlfile.write(TEMPLATES_POST.htmlSummary().format(
            # mod_file=os.path.abspath(globs['mod-file']),
            run_name=globs['run-name'],
            run_time=globs['startdatetimenice'],
            host_name=os.uname()[1],
            script_call=globs['call'],
            num_batches_complete_st=str(len(globs['complete-batches-st'])),
            num_batches_complete_gt=str(len(globs['complete-batches-gt'])),
            num_batches_incomplete_st=str(len(globs['incomplete-batches-st'])),
            num_batches_incomplete_gt=str(len(globs['incomplete-batches-gt'])),
            batch_size=str(globs['batch-size']),
            procs_per_batch=str(globs['procs-per-batch']),
            avg_runtime=str(round(PC.mean(globs['batch-runtimes']))),

            batch_comment_start=batch_comment_start,
            incomplete_batches=", ".join(globs['incomplete-batches']),
            batch_comment_end=batch_comment_end,

            bf1_cutoff=str(globs['bf1-cutoff']),
            bf1_loci=str(len(globs['bf1-loci'])),
            bf2_cutoff=str(globs['bf2-cutoff']),
            bf2_loci=str(len(globs['bf2-loci'])),
            bf3_cutoff=str(globs['bf3-cutoff']),
            bf3_loci=str(len(globs['bf3-loci'])),

            total_loci=str(len(globs['locus-stats']['elem_lik'])),
            target_loci=str(len(globs['m1-loci'])),
            full_loci=str(len(globs['m2-loci'])),

            log_file=os.path.basename(globs['logfilename']),
            results_folder=os.path.basename(globs['outdir']) + "/",

            med_bf1=str(round(PC.median(globs['bf1']), 3)),
            med_bf2=str(round(PC.median(globs['bf2']), 3)),
            med_bf3=str(round(PC.median(globs['bf3']), 3)),

            bf_boxplot=os.path.join("plots", globs['bf-dist-file']),
            #bf1_hist=os.path.join("plots", globs['bf1-dist-file']),
            #bf2_hist=os.path.join("plots", globs['bf2-dist-file']),
            bf1_bf2_plot=os.path.join("plots", globs['bf1-bf2-plot-file']),

            avg_m2_locus=str(round(PC.mean(list(globs['m2-per-locus'].values())))),
            med_m2_locus=str(round(PC.median(list(globs['m2-per-locus'].values())))),

            m2_dist=os.path.join("plots", globs['m2-locus-dist-file']),

            avg_m2_lineages=str(round(PC.mean(list(globs['m2-per-lineage'].values())))),
            med_m2_lineages=str(round(PC.median(list(globs['m2-per-lineage'].values())))),

            m2_counts=os.path.join("plots", globs['m2-lineage-dist-file']),

            comment_start=comment_start,
            comment_end=comment_end,

            num_spec=placeholder,
            num_targets=placeholder,
            num_conserved=placeholder,
            num_outgroups=placeholder,
            input_tree_plot=placeholder,

            date_time=PC.getFooterDateTime(),
        ))
    # Write the HTML summary file using the template

    step_start_time = PC.report_step(globs, step, step_start_time, "Success");
    globs['html-summary-written'] = True;
    # Status update

    return globs;

#############################################################################