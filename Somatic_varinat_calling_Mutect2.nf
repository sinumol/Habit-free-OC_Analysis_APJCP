#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
=========================================================
 Somatic Variant Calling Pipeline (Mutect2)
 FFPE-aware with contamination estimation
 Author: Sinumol George
 Created: 12-05-2023
=========================================================
*/

params.tumor_bam     = ""
params.sample_id    = ""
params.sample_type  = "NON_FFPE"   // FFPE or NON_FFPE
params.reference    = ""
params.pon          = ""
params.germline     = ""
params.outdir       = "results"
params.gatk         = "gatk"
params.threads      = 8

workflow {

    MUTECT2(
        params.tumor_bam,
        params.sample_id,
        params.reference,
        params.pon,
        params.germline
    )

    GET_PILEUPS(
        params.tumor_bam,
        params.sample_id,
        params.reference,
        params.germline
    )

    CALCULATE_CONTAMINATION(
        GET_PILEUPS.out
    )

    if (params.sample_type == "FFPE") {

        LEARN_ORIENTATION(
            MUTECT2.out.f1r2
        )

        FILTER_FFPE(
            MUTECT2.out.vcf,
            CALCULATE_CONTAMINATION.out,
            LEARN_ORIENTATION.out,
            params.sample_id
        )

    } else {

        FILTER_STANDARD(
            MUTECT2.out.vcf,
            CALCULATE_CONTAMINATION.out,
            params.sample_id
        )
    }
}

/* ===================== PROCESSES ===================== */

process MUTECT2 {

    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/mutect2", mode: 'copy'

    input:
    path tumor_bam
    val sample_id
    path reference
    path pon
    path germline

    output:
    path "${sample_id}.unfiltered.vcf.gz", emit: vcf
    path "${sample_id}.f1r2.tar.gz", emit: f1r2

    script:
    """
    ${params.gatk} Mutect2 \
      -R ${reference} \
      -I ${tumor_bam} \
      -tumor ${sample_id} \
      --panel-of-normals ${pon} \
      --germline-resource ${germline} \
      --f1r2-tar-gz ${sample_id}.f1r2.tar.gz \
      -O ${sample_id}.unfiltered.vcf.gz
    """
}

process GET_PILEUPS {

    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/contamination", mode: 'copy'

    input:
    path tumor_bam
    val sample_id
    path reference
    path germline

    output:
    path "${sample_id}.pileups.table"

    script:
    """
    ${params.gatk} GetPileupSummaries \
      -I ${tumor_bam} \
      -V ${germline} \
      -L ${germline} \
      -R ${reference} \
      -O ${sample_id}.pileups.table
    """
}

process CALCULATE_CONTAMINATION {

    publishDir "${params.outdir}/contamination", mode: 'copy'

    input:
    path pileups

    output:
    path "contamination.table"

    script:
    """
    ${params.gatk} CalculateContamination \
      -I ${pileups} \
      -O contamination.table
    """
}

process LEARN_ORIENTATION {

    publishDir "${params.outdir}/ffpe", mode: 'copy'

    input:
    path f1r2

    output:
    path "read-orientation-model.tar.gz"

    script:
    """
    ${params.gatk} LearnReadOrientationModel \
      -I ${f1r2} \
      -O read-orientation-model.tar.gz
    """
}

process FILTER_FFPE {

    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/filtered", mode: 'copy'

    input:
    path vcf
    path contamination
    path orientation
    val sample_id

    output:
    path "${sample_id}.filtered.vcf.gz"

    script:
    """
    ${params.gatk} FilterMutectCalls \
      -V ${vcf} \
      --contamination-table ${contamination} \
      --ob-priors ${orientation} \
      -O ${sample_id}.filtered.vcf.gz
    """
}

process FILTER_STANDARD {

    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/filtered", mode: 'copy'

    input:
    path vcf
    path contamination
    val sample_id

    output:
    path "${sample_id}.filtered.vcf.gz"

    script:
    """
    ${params.gatk} FilterMutectCalls \
      -V ${vcf} \
      --contamination-table ${contamination} \
      -O ${sample_id}.filtered.vcf.gz
    """
}
