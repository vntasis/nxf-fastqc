# nextflow-fastp
A simple fastp-MultiQC pipeline, written in [nextflow](https://www.nextflow.io/).
For a bunch of fastq files in a directory (PE or SE), run it with:

```
nextflow run main.nf --readsdir path/to/fastqfiles/
```

The pipeline executes [fastp](https://github.com/OpenGene/fastp), saves the filtered files in `results/fastp_trimmed`, and generates a [MultiQC](https://multiqc.info/) report in `results`. That's it!

For all available options, try

```
nextflow run main.nf --help
```

If you have `conda` you can run the pipeline with:
```
nextflow run main.nf --readsdir path/to/fastqfiles/ -profile conda
```

If you don't have conda or nextflow, go get them!