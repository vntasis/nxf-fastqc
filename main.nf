// FASTQC pipeline

/*
NXF ver 19.08+ needed because of the use of tuple instead of set
*/
if( !nextflow.version.matches('>=19.08') ) {
    println "This workflow requires Nextflow version 19.08 or greater and you are running version $nextflow.version"
    exit 1
}


// date needed to prefix results dir
DATE = new java.util.Date()
sdf = new java.text.SimpleDateFormat("yyyy-MM-dd")
fdate = sdf.format(DATE)
//println sdf.format(DATE)

// needed to pretty print read/bases counts
import java.text.DecimalFormat
df = new DecimalFormat("###,###") //TODO add symbols to fix US locale, http://tutorials.jenkov.com/java-internationalization/decimalformat.html#creating-a-decimalformat-for-a-specific-locale
//println df.format( 5.13e+23.toBigInteger() )

/*

in bash this works OK:
LC_NUMERIC=en_GB.UTF-8 printf "%'.0f\n" 1.345e+12

or in ksh (only separately, not as part of a pipe or with xargs)
printf '%#d\n' 105000000
gives
105M directly! works for G, T, P ...why not in bash?

even this works!!!
printf "%#d" 12e+12

so no need for java classes
*/

/*
 * pipeline input parameters
 */
params.readsdir = "fastq"
params.outdir = "${workflow.launchDir}/results-fastqc"
params.fqpattern = "*R{1,2}*.fastq.gz"
params.ontreads = false
//params.threads = 2 //makes no sense I think, to be removed
params.multiqc_config = "$baseDir/multiqc_config.yml" //custom config mainly for sample names
params.title = "Summarized nxf-fastqc report"
params.help = ""
params.dsrcDecompress = false

mqc_config = file(params.multiqc_config) // this is needed, otherwise the multiqc config file is not available in the docker image

if (params.help) {
    helpMessage()
    exit(0)
}

log.info """
        ===========================================
         F A S T Q C  P I P E L I N E

         Used parameters:
        -------------------------------------------
         --readsdir         : ${params.readsdir}
         --fqpattern        : ${params.fqpattern}
         --outdir           : ${params.outdir}
         --ontreads         : ${params.ontreads}
         --multiqc_config   : ${params.multiqc_config}
         --title            : ${params.title}

         Runtime data:
        -------------------------------------------
         Running with profile:   ${workflow.profile}
         Running as user:        ${workflow.userName}
         Launch dir:             ${workflow.launchDir}
         Base dir:               ${baseDir}
         """
         .stripIndent()

def helpMessage() {
log.info """
        ===========================================
         F A S T Q C   P I P E L I N E

         Usage:
        -------------------------------------------
         --readsdir         : directory with fastq files, default is "fastq"
         --fqpattern        : regex pattern to match fastq files, default is "*R{1,2}*.fastq.gz"
         --outdir           : where results will be saved, default is "results-fastp"
         --dsrcDecompress   : decompress input files with DSRC DNA Sequence Reads Compressor
         --ontreads         : use this parameter for Nanopore reads
         --multiqc_config   : path to config file for MultiQC, default is "multiqc_config.yml"
         --title            : MultiQC report title, default is "Summarized fastp report"
        ===========================================
         """
         .stripIndent()

}



//just in case trailing slash in readsdir not provided...
readsdir_repaired = "${params.readsdir}".replaceFirst(/$/, "/")
//println(readsdir_repaired)

// build search pattern for fastq files in input dir
reads = readsdir_repaired + params.fqpattern

// get counts of found fastq files
readcounts = file(reads)
println " Reads found:            ${readcounts.size()}"

// channel for read pairs --> fastp
Channel
    .fromFilePairs( reads, checkIfExists: true, size: -1 ) // default is 2, so set to -1 to allow any number of files
    .ifEmpty { error "Can not find any reads matching ${reads}" }
    .set{ read_pairs_ch }

// channel for reads --> seqtools
Channel
    .fromPath( reads, checkIfExists: true )
    .set { reads_ch }


// fastp trimmed files are published, json are only sent in the channel and used only by multiqc
process fastp {

    tag "fastp on $sample_id"
    //echo true
    publishDir params.outdir, pattern: 'fastp_trimmed/*' // publish only trimmed fastq files

    input:
        tuple sample_id, file(x) from read_pairs_ch
        val dsrc from params.dsrcDecompress

    output:
        file("${sample_id}_fastp.json") into fastp_ch
        file('fastp_trimmed/trim_*')
        file("${sample_id}_fastp.json") into fastp_ch2
        val seqmode into seqmode_ch


    script:
    def single = x instanceof Path // this is from Paolo: https://groups.google.com/forum/#!topic/nextflow/_ygESaTlCXg
    def qscore_cutoff = params.ontreads ? 7 : 15 //here ontreads matters
    if ( !single ) {
        seqmode = "PE"
        """
        function finish {
          $dsrc && rm read1.fastq read2.fastq || exit 0
        }

        if $dsrc
        then
          dsrc d -t ${task.cpus} ${x[0]} read1.fastq && dsrc d -t ${task.cpus} ${x[1]} read2.fastq
          read1='read1.fastq' && read2='read2.fastq'
        else
          read1=${x[0]} && read2=${x[1]}
        fi

        mkdir fastp_trimmed
        fastp \
        -q $qscore_cutoff \
        -i \$read1 -I \$read2 \
        -o fastp_trimmed/trim_${x[0]} -O fastp_trimmed/trim_${x[1]} \
        -j ${sample_id}_fastp.json --thread ${task.cpus}

        trap finish EXIT
        """
    }
    else {
        seqmode = "SE"
        """
        function finish {
          $dsrc && rm reads.fastq || exit 0
        }

        if $dsrc
        then
          dsrc d -t ${task.cpus} ${x} reads.fastq && reads='reads.fastq'
        else
          reads=${x}
        fi

        mkdir fastp_trimmed
        fastp \
        -q $qscore_cutoff \
        -i \$reads -o fastp_trimmed/trim_${x} \
        -j ${sample_id}_fastp.json --thread ${task.cpus}

        trap finish EXIT
        """
    }

}

//=========================

/*
* This process reads the json files from fastp (using jq)
* sums the reads/bases from all read files!
* and sends the stdout in the channels total_reads and total_bases. The stdout is in this case a single value.
* These are used later in multiqc (via the --cl_config parameter)
* to add the numbers as section comments
*/

process summary {

    input:
        file x from fastp_ch2.collect()

    output:
        stdout into total_reads //the channel contains 4 values, sep by new line

    script:
    """
    #!/usr/bin/env ksh

    # check if jq installed, relevant only in local envs
    command -v jq >/dev/null 2>&1 || \
    { echo >&2 "jq required but not installed, install it or use -profile docker or conda. Aborting."; exit 1;}

    # prefer to keep the tmp files in the scratch for reference

    jq '.summary.before_filtering.total_reads' $x | awk '{sum+=\$0} END{print sum}' > treads-before.tmp
    jq '.summary.after_filtering.total_reads' $x | awk '{sum+=\$0} END{print sum}' > treads-after.tmp
    jq '.summary.before_filtering.total_bases' $x | awk '{sum+=\$0} END{print sum}' > tbases-before.tmp
    jq '.summary.after_filtering.total_bases' $x | awk '{sum+=\$0} END{print sum}' > tbases-after.tmp

    cat treads-before.tmp treads-after.tmp tbases-before.tmp tbases-after.tmp | xargs siformat.sh
    """
}

//=========================
process multiqc {
    publishDir params.outdir
    when:
        !params.ontreads // multiqc fails on fastp json from ont files! too big

    input:
        file x from fastp_ch.collect()
        file mqc_config
        val y from total_reads // y is a string with 4 values sep by new line now, but still length 1
        val seqmode from seqmode_ch.last() // PE or SE, see process fastp, last() because it emits one entry for each fastq file

    output:
        file('multiqc_report.html')
        file('multiqc_report_data/multiqc_general_stats.txt')

    // when using --title, make sure that the --filename is explicit, otherwise
    // multiqc uses the title string as output filename
    script:
    // currently y is of length 1, the split gets the values in a list (see string split() method in groovy)
    def splitstring = y.split()

    """
    multiqc --force --interactive \
    --title "${params.title}" \
    --filename "multiqc_report.html" \
    --config $mqc_config \
    --cl_config "section_comments:
                    { fastp: '*This is ${ seqmode } data *<br>
                              Total reads before filter: ** ${ splitstring[0] } ** <br>
                              Total reads    after filter: ** ${ splitstring[1] } ** <br>
                              Total bases before filter: ** ${ splitstring[2] } ** <br>
                              Total bases    after filter: ** ${ splitstring[3] } ** <br> <br>
                              See also the report:
                              <a href="fastq-stats-report.html" target="_blank">fastq-stats-report.html</a>,
                              which contains per-read data about Phred-score and k-mer distribution, N50 etc. <br>
                              The raw data from both reports (as *.csv files) can be found in the results folder.'
                    }
                " \
    ${x}
    """
}

//=============================
Channel.fromPath("${baseDir}/bin/fastq-stats-report.Rmd").set{ fastq_stats_report_ch }
Channel.fromPath("${baseDir}/bin/fastq-stats-report-ont.Rmd").set{ fastq_stats_report_ont_ch }

if (!params.ontreads) {
    process fastq_stats_ilmn {
    publishDir params.outdir

    input:
        file x from reads_ch.collect()
        file 'fastq-stats-report.Rmd' from fastq_stats_report_ch

    output:
        file 'fastq-stats-report.html'
        file "fastq-stats.csv"
        file "fastq-stats.xlsx"

    script:
    """
    seqtools.R $x
    """
    }
} else {
    process fastq_stats_ont {
    publishDir params.outdir

    input:
        file x from reads_ch.collect()
        file 'fastq-stats-report-ont.Rmd' from fastq_stats_report_ont_ch

    output:
        file 'fastq-stats-report-ont.html'
        file "fastq-stats.csv"
        file "fastq-stats.xlsx"

    script:
    """
    seqtools-ont.R $x
    """
    }
}
//=============================
workflow.onComplete {
    if (workflow.success) {
        log.info """
            ===========================================
            Finished in ${workflow.duration}
            Processed ${ readcounts.size() } fastq files
            See the report here ==> $params.outdir/multiqc_report.html
            """
            .stripIndent()
    } else {
        log.info """
            ===========================================
            Finished with errors!
            """
            .stripIndent()
    }
}
