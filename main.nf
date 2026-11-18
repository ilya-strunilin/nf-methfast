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
}
