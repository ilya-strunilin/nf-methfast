nextflow.enable.dsl = 2

// Quantify methylation in a shared set of genomic intervals, once per sample.

def required = ['samples', 'target_bed', 'outdir', 'postfix']
required.each { name ->
    if (!params[name]) {
        error "Missing --${name}. See README.md for usage."
    }
}

def samplesTsv = file(params.samples, checkIfExists: true)
def targetBed = file(params.target_bed, checkIfExists: true)

if (!(targetBed.name ==~ /.+\.bed(\.gz)?/)) {
    error "Target file must end in .bed or .bed.gz: ${targetBed}"
}

def postfix = params.postfix.toString()
if (!(postfix ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
    error "Unsafe postfix '${postfix}'. Use letters, digits, dots, underscores or hyphens."
}

process METHFAST {
    tag "${sample_id}"

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(methylation_bed)
    path target_bed
    val postfix

    output:
    tuple val(sample_id), path("${sample_id}.${postfix}.bed")

    script:
    """
    ml methfast/0.3.0
    
    methfast '${methylation_bed}' '${target_bed}' -o '${sample_id}.${postfix}.bed'
    """
}

process MERGE_BEDS {
    tag "${postfix}"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path beds
    val postfix

    output:
    path "merged.${postfix}.tsv.gz"

    script:
    """
    set -o pipefail

    {
        printf 'chr\tstart\tend\tn_cpg\tcoverage\tbeta\tsample_id\n'

        for bed in *.'${postfix}'.bed; do
            sample_id="\${bed%.${postfix}.bed}"

            awk -v sample_id="\$sample_id" '
                BEGIN { FS = OFS = "\t" }
                NR == 1 && (\$1 == "chr" || \$1 == "#chr") && \$2 == "start" { next }
                NF != 6 {
                    print "Expected 6 columns in " FILENAME " at line " FNR "; found " NF > "/dev/stderr"
                    exit 1
                }
                { print \$1, \$2, \$3, \$4, \$5, \$6, sample_id }
            ' "\$bed" || exit 1
        done
    } | gzip -c > merged.${postfix}.tsv.gz
    """
}

workflow {
    def sampleRows = Channel.fromPath(samplesTsv)
        .splitCsv(header: true, sep: '\t')

    def samples = sampleRows.map { row ->
        if (!row.containsKey('wgbs_directory_name')) {
            error "Missing 'wgbs_directory_name' column in ${samplesTsv}"
        }

        def id = row.wgbs_directory_name?.toString()?.trim()
        if (!id || id.equalsIgnoreCase('NA')) {
            error "Missing wgbs_directory_name in ${samplesTsv}"
        }
        if (!(id ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
            error "Unsafe sample directory name: '${id}'"
        }

        def methylationBed = file("${params.wgbs_root}/${id}/${id}.meth.bed.gz", checkIfExists: true)
        tuple(id, methylationBed)
    }.unique { sample -> sample[0] }

    METHFAST(samples, targetBed, postfix)
    MERGE_BEDS(METHFAST.out.map { sample_id, bed -> bed }.collect(), postfix)
}
