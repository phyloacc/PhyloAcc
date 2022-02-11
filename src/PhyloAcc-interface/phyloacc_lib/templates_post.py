#############################################################################
# Templates for various files written by the phyloacc post-processing script
#############################################################################

def htmlSummary():

    html_template = """
<!doctype html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1">	
<link href="https://fonts.googleapis.com/css?family=Lato:300" rel="stylesheet"> 
<title>PhyloAcc results</title>
</head>

<body>
    <div class="row" id="top_grid">
        <div class="col-2-24" id="margin"></div>
        <div class="col-20-24" id="main_header">PhyloAcc results</div>
        <div class="col-2-24" id="margin"></div>
    </div>

    <div class="row" id="main-container">
        <div class="col-4-24" id="side-nav">

            <div class="row" id="side-nav-row">
                <div class="col-24-24" id="side-nav-container">
                    <div class="row" id="nav-header-container">
                        <div class="col-24-24" id="nav-header">Navigation</div>
                    </div>
                
                    <div class="row" id="nav-content-container">
                        <div class="col-2-24" id="side-nav-margin"></div>
                        <div class="col-22-24" id="nav-links">
                            <a class="nav-link" href="#run-info">1. Run info</a>
                            <a class="nav-link" href="#bf">2. Bayes factors</a>
                            {comment_start}
                            <a class="nav-link" href="#locus-lens">3. Distribution of locus lengths</a>
                            <a class="nav-link" href="#informative-sites">4. Distribution of informative sites per locus</a>
                            <a class="nav-link" href="#variable-informative-sites">5. Correlation between variable and informative sites</a>
                            <a class="nav-link" href="#scf">6. Site concordance factors (sCF) per locus</a>
                            <a class="nav-link" href="#scf-tree">7. sCF across all loci</a>
                            <a class="nav-link" href="#bl-scf">8. Branch length vs. sCF for the species tree</a> 
                            {comment_end}
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <div class="col-20-24" id="main-content">

            <div class="row content-row" id="run-info-row">
                <div class="col-24-24 content-col" id="run-info">
                    <div class="section-header">1. Run info ({run_name})</div>
                    <p>The PhyloAcc post-processing script was run on <code>{run_time}</code> on <code>{host_name}</code> as follows:</p>
                    
                    <div class="code-focus">
                        <pre>
{script_call}
                        </pre>
                    </div>

                    <div class="line"></div>

                    <div class="sub-header">Batch summary</div>
                    <div class="table-container">
                        <table class="table-content">
                            <tr>
                                <td>Complete ST batches</td>
                                <td>{num_batches_complete_st}</td>
                            </tr>

                            <tr>
                                <td>Complete GT batches</td>
                                <td>{num_batches_complete_gt}</td>
                            </tr>

                            <tr>
                                <td>Incomplete ST batches</td>
                                <td>{num_batches_incomplete_st}</td>
                            </tr>

                            <tr>
                                <td>Incomplete GT batches</td>
                                <td>{num_batches_incomplete_gt}</td>
                            </tr>

                            <tr>
                                <td>Loci per batch</td>
                                <td>{batch_size}</td>
                            </tr>

                            <tr>
                                <td>Threads per batch</td>
                                <td>{procs_per_batch}</td>
                            </tr>                        

                            <tr>
                                <td>Average runtime per batch (min.)</td>
                                <td>{avg_runtime}</td>
                            </tr>   
                        </table>
                    </div>

                    {batch_comment_start}

                    <p>
                        The following batches were not summarized and may be incomplete:
                    <p>

                    <div class="code-focus">
                        <pre>
{incomplete_batches}
                        </pre>
                    </div>

                    <p>
                        To determine the problem, navigate to the <code>phyloacc-job-files/phyloacc-output/</code> directory in the output directory for this run.
                        Find and list the contents for or open the folder corresponding to the batch. If there is no file <code>BATCH_elem_lik.txt</code>,
                        where BATCH corresponds to the batch id, then this run did not finish. You can also check the logfile, <code>BATCH-phyloacc.log</code> to see
                        if there are any error messages.
                    </p>

                    {batch_comment_end}

                    <div class="sub-header">Result summary</div>
                    <div class="table-container">
                        <table class="table-content">
                            <tr>
                                <td>Total loci</td>
                                <td>{total_loci}</td>
                            </tr>

                            <tr>
                                <td>Accelerated loci</td>
                                <td>{accelerated_loci}</td>
                            </tr>
                        </table>
                    </div>

                    <p>
                        See the <a href="{log_file}">full log file</a> for more info.
                    </p>

                    <p>
                        Below are some summary plots. Raw data is also available in in the <a href="{results_folder}">results folder</a>.
                    </p>

                    <div class="line"></div>

                </div>
            </div>

            <div class="sep-div"></div>

            {comment_start}
            <div class="row content-row" id="input-tree-row">
                <div class="col-24-24 content-col" id="input-tree">
                    <div class="section-header">2. Input species tree</div>

                    <div class="table-container">
                        <table class="table-content">
                            <tr>
                                <td>Species</td>
                                <td>{num_spec}</td>
                            </tr>

                            <tr>
                                <td>Target species</td>
                                <td>{num_targets}</td>
                            </tr>                        

                            <tr>
                                <td>Conserved species</td>
                                <td>{num_conserved}</td>
                            </tr>   

                            <tr>
                                <td>Outgroup species</td>
                                <td>{num_outgroups}</td>
                            </tr>  
                        </table>
                    </div>

                    <div class="row img-row" id="input-tree-container-row">
                        <div class="col-4-24 margin"></div>
                        <div class="col-16-24 img-container" id="input-tree-container">
                            <img id="input-tree-img" src="{input_tree_plot}">
                        </div>
                        <div class="col-4-24 margin"></div>
                    </div>

                </div>
            </div>
            {comment_end}

            <div class="sep-div"></div>

            <div class="row content-row" id="bf-row">
                <div class="col-24-24 content-col" id="bf">
                    <div class="section-header">2. Bayes factors</div>

                    {comment_start}
                    <div class="table-container">
                        <table class="table-content">
                            <tr>
                                <td>Average alignment length</td>
                                <td>{avg_aln_len}</td>
                            </tr>

                            <tr>
                                <td>Median alignment length</td>
                                <td>{median_aln_len}</td>
                            </tr>                        

                            <tr>
                                <td>Average sequence length (no gaps) per alignment</td>
                                <td>{avg_seq_len_nogap}</td>
                            </tr>   

                           <tr>
                                <td>Median sequence length (no gaps) per alignment</td>
                                <td>{med_seq_len_nogap}</td>
                            </tr>   

                        </table>
                    </div>
                    {comment_end}

                    <div class="row img-row" id="bf-container-row">
                        <div class="col-2-24 margin"></div>
                        <div class="col-9-24 img-container" id="bf1-container">
                            <img id="bf1-img" src="{bf1_hist}">
                        </div>
                        <div class="col-2-24 margin"></div>
                        <div class="col-9-24 img-container" id="bf2-container">
                            <img id="bf2-img" src="{bf2_hist}">
                        </div>
                        <div class="col-2-24 margin"></div>
                    </div>

                    <div class="small-sep-div"></div>

                    <div class="row img-row" id="bf-container-row">
                        <div class="col-7-24 margin"></div>
                        <div class="col-10-24 img-container" id="bf1-bf2-container">
                            <img id="bf1-bf2-img" src="{bf1_bf2_plot}">
                        </div>
                        <div class="col-7-24 margin"></div>
                    </div>


                </div>
            </div>

            <div class="sep-div"></div>

        </div>

    </div>


    <div class="row" id="footer">
    <div class="col-1">
        <div id="footer_text">
            <center>Page generated by the <a href="https://github.com/harvardinformatics/PhyloAcc-interface" target="_blank">PhyloAcc interface</a> | Page built: {date_time}</center>
        </div>
    </div>

</body>

<style type="text/css" media="screen">
    /*------------------------------------------------------*/
    /* Global page styles */
    body {{
        color:#1a1a1a;
        font-family: 'Helvetica', sans-serif;
        font-size: 1.1em;
        background-color: #ececec;
        overflow-y: scroll;
        margin: 0;
    }}
    a {{
        color:#006ddb;
    }}

    /*------------------------------------------------------*/
    /* Grid styles */

    .row {{ 
        width:100%;
        display: flex;
        justify-content: center;
        align-items: center;
    }}

    .col-1-24 {{ width:4.166667%; display:inline-block; }}
    .col-2-24 {{ width:8.333333%; display:inline-block; }}
    .col-3-24 {{ width:12.5%; display:inline-block; }}
    .col-4-24 {{ width:16.666667%; display:inline-block; }}
    .col-5-24 {{ width:20.833333%; display:inline-block; }}
    .col-6-24 {{ width:25%; display:inline-block; }}
    .col-7-24 {{ width:29.166667%; display:inline-block; }}
    .col-8-24 {{ width:33.333333%; display:inline-block; }}
    .col-9-24 {{ width:37.5%; display:inline-block; }}
    .col-10-24 {{ width:41.666667%; display:inline-block; }}
    .col-11-24 {{ width:45.833333%; display:inline-block; }}
    .col-12-24 {{ width:50%; display:inline-block; }}
    .col-13-24 {{ width:54.166667%; display:inline-block; }}
    .col-14-24 {{ width:58.333333%; display:inline-block; }}
    .col-15-24 {{ width:62.5%; display:inline-block; }}
    .col-16-24 {{ width:66.666667%; display:inline-block; }}
    .col-17-24 {{ width:70.833333%; display:inline-block; }}
    .col-18-24 {{ width:75%; display:inline-block; }}
    .col-19-24 {{ width:79.1666667%; display:inline-block; }}
    .col-20-24 {{ width:83.333333%; display:inline-block; }}
    .col-21-24 {{ width:87.5%; display:inline-block; }}
    .col-22-24 {{ width:91.666667%; display:inline-block; }}
    .col-23-24 {{ width:95.833333%; display:inline-block; }}
    .col-24-24 {{ width:100%; display:inline-block; }}
    /* 24 */

    .col-1-12 {{ width:8.333333%; display:inline-block; }}
    .col-2-12 {{ width:16.666667%; display:inline-block; }}
    .col-3-12 {{ width:25%; display:inline-block; }}
    .col-4-12 {{ width:33.333333%; display:inline-block; }}
    .col-5-12 {{ width:41.666667%; display:inline-block; }}
    .col-6-12 {{ width:50%; display:inline-block; }}
    .col-7-12 {{ width:58.333333%; display:inline-block; }}
    .col-8-12 {{ width:66.666667%; display:inline-block; }}
    .col-9-12 {{ width:75%; display:inline-block; }}
    .col-10-12 {{ width:83.333333%; display:inline-block; }}
    .col-11-12 {{ width:91.666667%; display:inline-block; }}
    .col-12-12 {{ width:100%; display:inline-block; }}
    /* 12 */

    .col-1-8 {{ width:12.5%; display:inline-block; }}
    .col-2-8 {{ width:25%; display:inline-block; }}
    .col-3-8 {{ width:37.5%; display:inline-block; }}
    .col-4-8 {{ width:50%; display:inline-block; }}
    .col-5-8 {{ width:62.5%; display:inline-block; }}
    .col-6-8 {{ width:75%; display:inline-block; }}
    .col-7-8 {{ width:87.5%; display:inline-block; }}
    .col-8-8 {{ width:100%; display:inline-block; }}
    /* 8 */

    .col-1-6 {{ width:16.666667%; display:inline-block; }}
    .col-2-6 {{ width:33.333333%; display:inline-block; }}
    .col-3-6 {{ width:50%; display:inline-block; }}
    .col-4-6 {{ width:66.666667%; display:inline-block; }}
    .col-5-6 {{ width:83.333333%; display:inline-block; }}
    .col-6-6 {{ width:100%; display:inline-block; }}
    /* 6 */

    .col-1-4 {{ width:25%; display:inline-block; }}
    .col-2-4 {{ width:50%; display:inline-block; }}
    .col-3-4 {{ width:75%; display:inline-block; }}
    .col-4-4 {{ width:100%; display:inline-block; }}
    /* 4 */

    .col-1-3 {{ width:33.333333%; display:inline-block; }}
    .col-2-3 {{ width:66.666667%; display:inline-block; }}
    .col-3-3 {{ width:100%; display:inline-block; }}
    /* 3 */

    .col-1-2 {{ width:50%; display:inline-block; }}
    .col-2-2 {{ width:100%; display:inline-block; }}
    /* 2 */

    .col-1-1 {{ width:100%; display:inline-block; }}
    /* 1 */

    /*------------------------------------------------------*/
    /* Banner styles */
    #top_grid {{
        height:15vh;

        /* background-image: url(html/img/bg1.png); */
        /* https://pixabay.com/illustrations/color-triangle-geometric-textured-2174049/ */

        background: linear-gradient(25deg, rgba(244,251,193,1) 0%, rgba(150,211,196,1) 50%, rgba(50,119,174,1) 100%);
        /* https://cssgradient.io/ */
        
        background-repeat: no-repeat;
        background-clip: border-box;
        background-origin: padding-box;
        -moz-background-size: cover;
        background-size: cover;
    }}
    #main_header {{
        color:#999999;
        -webkit-text-stroke-width: 0.5px;
        -webkit-text-stroke-color: black;
        font-size:3em;
        padding-top:10px;
        padding-bottom:10px;
        text-align: left;
        vertical-align: middle;
    }}

    /*------------------------------------------------------*/
    /* Header styles */

    .section-header {{
        background-color: #749abe;
        color: #fff;
        padding: 15px;
        font-size:2em;
        font-weight:bold;
        position: sticky;
        top: 0px;
        /* border-bottom: 2px solid #294157;
        border-left: 2px solid #294157; */
    }}
    @media (max-width: 1025px) {{
        .section-header  {{ 		
            text-align:center; 
        }}
    }}
    .sub-header {{
        /* background-color: #749abe; */
        /* color: #fff; */
        padding: 15px;
        font-size:1.5em;
        font-weight:bold;
    }}
    .sub-header-2 {{
        /* background-color: #749abe; */
        padding: 30px;
        font-size:1.25em;
        font-weight:bold;
    }}

    /*------------------------------------------------------*/
    /* Main container styles */

    #main-container {{
        display: flex;
        align-items: stretch;
        height: auto;
        min-height: 85vh;
    }}

    /*------------------------------------------------------*/
    /* Side nav bar styles */

    #side-nav-row {{
        height:97vh;
        background-color: #ffffff;
        border-right:1px solid black;
        position: sticky;
        top: 0px;
    }}
    #side-nav-container {{
        height:100%;
    }}
    #nav-header {{
        padding: 15px 0px;
        text-indent: 15px;
        font-size: 2em;
        overflow-wrap: break-word;
        background-color: #294157;
        color: #ffffff;
    }}
    .nav-link {{
        display: block;
        padding: 10px 10px 10px 0px;
        font-size: 1em;
        text-decoration: none;
        transition: 0.1s linear;

    }}
    .nav-link:hover {{
        color: #ececec;
        background-color: #56b4e9;
    }}
    @media (max-width: 1025px) {{
        #side-nav {{ height:94vh; }}
    }}

    /*------------------------------------------------------*/
    /* Content styles */

    .content-row p {{
        padding: 0 40px;
    }}
    .list-container {{
        padding-left: 40px;
    }}
    .code-focus {{
        display: flex;
        align-items: center;
        justify-content: center;
        font-family:'Courier New', Courier, monospace;
        font-size: 1.2em;
    }}
    .code-focus pre {{
        width: 75%;
        overflow: auto;
        background-color:#ffffff;
    }}
    .img-container {{
        display: flex;
        align-items: center;
        justify-content: center;
    }}
    .img-container img {{
        display: inline-block;
        width: 100%;
        height: auto;
    }}

    /*------------------------------------------------------*/
    /* Table styles */
    .table-container {{
        padding:20px 50px 20px 50px;
        overflow-y:auto;
        display:flex;
        align-items: center;
        justify-content: center;

    }}
    .table-content {{
        background-color:#ffffff;
        border: 1px solid #333333;
        border-radius: 10px;
        width: 25%;
    }}
    .table-content  th {{
        background-color:#000000;
        /* background-color: #A9C686; */
        /* background-color: #6B8B39; */
        padding:5px;
        width:12.5%;
        color: #ececec;
        /* border:1px solid red; */
    }}
    .table-content  tr:nth-child(even) {{
        background-color:#E0DEDE;
    }}

    /*------------------------------------------------------*/
    /* Footer grid styles */
    #footer {{
        font-family: "Courier New", Courier, monospace;
        background-color: #fff;
        color:#333;
        font-size:0.6em;
        height:3vh;
        display: flex;
        justify-content: center;
        align-content: center;
        flex-direction: column;
        -webkit-box-sizing: border-box;
        -moz-box-sizing: border-box;
        box-sizing: border-box;
        border-top:1px solid #333;
        width:100%;
    }}
    @media (max-width: 1025px) {{
        #footer {{ height:6vh; }}
        #footer_text {{ top: 25%; }}
    }}

    /*------------------------------------------------------*/
    /* Grid spacing styles */
    .small-sep-div {{
        width: 100%;
        height: 1vh;
    }}
    .sep-div {{
        width: 100%;
        height: 2vh;
    }}
    #vert_line{{
        width:1px;
        height:100%;
        background-color:#d3d3d3;
    }}
    .line {{
        height:1px;
        background-color:#d3d3d3;
        margin:20px;
    }}    
</style>

    """

    return html_template

#############################################################################