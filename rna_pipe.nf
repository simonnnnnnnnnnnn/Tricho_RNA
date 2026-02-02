#!/bin/bash

// used processes --> all are separarted in their own files
include { fastqc } from './modules/fastqc.nf'
include {trim_galore} from "./modules/trim_galore"
include {kallisto} from "./modules/kallisto"
include {multiqc} from "./modules/multiqc"
include {kallisto_index} from "./modules/kallisto_index"
include {salmon_index} from "./modules/salmon_index"
include { salmon } from './modules/salmon'
include {bowtie2_index} from './modules/bowtie2_index'
include {bowtie2} from './modules/bowtie2'
include {samtools} from './modules/samtools'
//include {stringtie} from './modules/stringtie'
include {star_index} from './modules/star_index'
include { star } from './modules/star'
include {picard_alignment_summary_bowtie2} from './modules/picard_alignment_summary_bowtie2'
include {picard_alignment_summary_star} from './modules/picard_alignment_summary_star'
/*include {picard_sort_sam} from './modules/picard_sort_sam'
include {picard_mark_duplicates} from './modules/picard_mark_duplicates'
include {picard_validate} from './modules/picard_validate'*/
include {featurecounts_bowtie2} from './modules/featurecounts_bowtie2'
include {featurecounts_star} from './modules/featurecounts_star'
//include {picard_seq_dict} from './modules/picard_seq_dict'
include {finish_pipeline} from './modules/finish_pipeline'
include {benchmarks} from './modules/benchmarks'
include {samtools_bowtie2} from './modules/samtools_bowtie2'
include {samtools_star} from './modules/samtools_star'
include { calculate_tpms_bowtie2 } from './modules/calculate_tpms_bowtie2'
include { calculate_tpms_star } from './modules/calculate_tpms_star'
include {sam_to_bam} from './modules/sam_to_bam'
include {rename_trace} from './modules/rename_trace'
include {rename_report} from './modules/rename_report'
include {create_ref_transcriptome} from './modules/create_ref_transcriptome'
include {create_gtf} from './modules/create_gtf'
include {diff_gene_expr_edger_kallisto} from './modules/diff_gene_expr_edger_kallisto'
include {diff_gene_expr_edger_salmon} from './modules/diff_gene_expr_edger_salmon'
include {diff_gene_expr_edger_bowtie2} from './modules/diff_gene_expr_edger_bowtie2'
include {diff_gene_expr_edger_star} from './modules/diff_gene_expr_edger_star'
include {make_tx2gene_from_gff} from './modules/make_tx2gene_from_gff'
include {clean_abundance} from './modules/clean_abundance'
include {clean_quant} from './modules/clean_quant'
include {clean_features_star} from './modules/clean_features_star'
include { clean_features_bowtie2 } from './modules/clean_features_bowtie2'
include {diff_gene_expr_kallisto} from './modules/diff_gene_expr_kallisto'
include {diff_gene_expr_salmon} from './modules/diff_gene_expr_salmon'
include {diff_gene_expr_bowtie2} from './modules/diff_gene_expr_bowtie2'
include {diff_gene_expr_star} from './modules/diff_gene_expr_star'

// visualization modules
include {star_stats} from './modules/star_stats'
include {bowtie2_stats} from './modules/bowtie2_stats'
include {star_calc_genome_cov} from './modules/star_calc_genome_cov'
include {bowtie_calc_genome_cov} from './modules/bowtie_calc_genome_cov'
include {star_calc_gene_cov} from './modules/star_calc_gene_cov'
include {bowtie_calc_gene_cov} from './modules/bowtie_calc_gene_cov'
include {bedtools_bowtie2} from './modules/bedtools_bowtie2'
include {star_genome_coverage} from './modules/star_genome_coverage'
include {bowtie_genome_coverage} from './modules/bowtie_genome_coverage'
include {star_gene_coverage} from './modules/star_gene_coverage'
include {bowtie_gene_coverage} from './modules/bowtie_gene_coverage'
include {star_mq} from './modules/star_mq'
include {bowtie_mq} from './modules/bowtie_mq'
include {pseudo_coverage} from './modules/pseudo_coverage'
include {compare_tpms} from './modules/compare_tpms'

// dge plotting
include {plot_volcano_edger_kallisto} from './modules/plot_volcano_edger_kallisto'
include {plot_volcano_edger_salmon} from './modules/plot_volcano_edger_salmon'
include {plot_volcano_edger_bowtie2} from './modules/plot_volcano_edger_bowtie2'
include {plot_volcano_edger_star} from './modules/plot_volcano_edger_star'
include {plot_volcano_limma_kallisto} from './modules/plot_volcano_limma_kallisto'
include {plot_volcano_limma_salmon} from './modules/plot_volcano_limma_salmon'
include {plot_volcano_limma_bowtie2} from './modules/plot_volcano_limma_bowtie2'
include {plot_volcano_limma_star} from './modules/plot_volcano_limma_star'
include {plot_volcano_deseq2_kallisto} from './modules/plot_volcano_deseq2_kallisto'
include {plot_volcano_deseq2_salmon} from './modules/plot_volcano_deseq2_salmon'
include {plot_volcano_deseq2_bowtie2} from './modules/plot_volcano_deseq2_bowtie2'
include {plot_volcano_deseq2_star} from './modules/plot_volcano_deseq2_star'

// new dge calc modules
include {dge_deseq2_bowtie2} from './modules/dge_deseq2_bowtie2'
include {dge_deseq2_kallisto} from './modules/dge_deseq2_kallisto'
include {dge_deseq2_salmon} from './modules/dge_deseq2_salmon'
include {dge_deseq2_star} from './modules/dge_deseq2_star'
include {dge_edger_bowtie2} from './modules/dge_edger_bowtie2'
include {dge_edger_kallisto} from './modules/dge_edger_kallisto'
include {dge_edger_salmon} from './modules/dge_edger_salmon'
include {dge_edger_star} from './modules/dge_edger_star'
include {dge_limma_bowtie2} from './modules/dge_limma_bowtie2'
include {dge_limma_kallisto} from './modules/dge_limma_kallisto'
include {dge_limma_salmon} from './modules/dge_limma_salmon'
include {dge_limma_star} from './modules/dge_limma_star'
//include {make_genome_txt} from './modules/make_genome_txt'
include {sort_gtf} from './modules/sort_gtf'
//include {genome_txt_from_bam} from './modules/genome_txt_from_bam'
include {get_scaffold_order} from './modules/get_scaffold_order'
include {get_scaffold_order_bowtie} from './modules/get_scaffold_order_bowtie'
include {sort_gtf_by_scaffold} from './modules/sort_gtf_by_scaffold'
include {bedtools_sort_gtf} from './modules/bedtools_sort_gtf'
include {bedtools_sort_gtf_bowtie} from './modules/bedtools_sort_gtf_bowtie'
include {mapping_rates} from './modules/mapping_rates'


// this parameter defines the experiment --> important for all paths and the structure of the results
params.experiment = "tricho_test"


// paramters --> here are the defaults, specified in params.yml
params.infile = "human_chr_21/material/read_paths.csv"
params.report_id = "pipeline_report_kallisto" // needs to be changed accordingly --> can be done in the params file
params.fasta = "human_chr_21/material/Homo_sapiens.GRCh38.dna.chromosome.21.fa" // genome fasta
//params.gtf = "human_chr_21/material/Homo_sapiens.GRCh38.91.chr21.gtf"
params.visualize = "scripts/visualize_benchmarks.r"
params.tpm_calc_script = "scripts/do_calculate_tpms.py"
params.annotation = "tricho_test/material/abc.gff3"
//params.gff = "${params.experiment}/material/Trichoplax_adhaerens.ASM15027v1.61.gff3.gz"
params.number_samples = 2
//params.aligner = "kallisto" --> better to just pass a string, that variable would have different 4 values
params.tx2gene_script = "abc"
params.edger = "aaaaaaaaaaaaaaaa"
params.all = "yeast/results/"
params.limma = "bbbb"
params.deseq2 = "cccc"
params.test = 'ffd'

params.vis_genome_coverage = "123"
params.vis_gene_coverage = "222"
params.vis_mq = "234"
params.pseudo_coverage = "dsds"
params.compare_tpms = "dsds"
params.plot_edger = "dsds"
params.plot_limma = "fdfdf"
params.plot_deseq2 = "dsdss"
params.test_gtf = "fdfd"
params.mapping_rates_script = "dsd"



// ToDo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// prüfen welchen input die verschiedenen tximport für die aligner brauchen --> also welche files
// dementsprechend muss auch das skript geändert werden!!!




//workflow --> here the processes are all chained together

workflow  {

    // workflow is dependent on input variables --> sometimes you get gff, sometimes its gtf
    def annotationfile = file(params.annotation)

    def annotation_type = annotationfile.name.endsWith(".gtf") ? "gtf"
                        : annotationfile.name.endsWith(".gff3") ? "gff"
                        : annotationfile.name.endsWith(".gff") ? "gff"
                        : "unknown"
    println("annotation type: ${annotation_type}")



    //input
    read_ch = Channel.fromPath(params.infile)
                        .splitCsv(header:true)
                        .map{row ->
                        println "processing row: ${row}"
                        [file(row.forward), file(row.reverse)]}


    // now that the input is established --> can call the channel
    fastqc(read_ch, params.experiment)// works

    trim_galore(read_ch, params.experiment) // works

    // put all quality reports together
    multiqc(
        fastqc.out.zip.mix(
            fastqc.out.html,
            trim_galore.out.trimming_reports,
            trim_galore.out.fastqc_reports_1,
            trim_galore.out.fastqc_reports_2 // i hope this works
        ).collect(),
        params.report_id, params.experiment
    )

    if (annotation_type == "gtf"){

        // when annotation is gtf
        create_ref_transcriptome(params.experiment, file(params.annotation), file(params.fasta))
        make_tx2gene_from_gff(params.experiment, file(params.annotation), file(params.tx2gene_script))


        // kallisto
        kallisto_index(create_ref_transcriptome.out.reference_transcriptome, params.experiment)
        kallisto(trim_galore.out.trimmed_reads, kallisto_index.out.idx, params.experiment)
        clean_abundance(kallisto.out.tsv, params.experiment)
        all_abundances_kallisto_ch = clean_abundance.out.collect()
        //diff_gene_expr_kallisto(all_abundances_kallisto_ch, file(params.edger), file(params.limma), file(params.deseq2), params.experiment, make_tx2gene_from_gff.out.tx2gene, "kallisto", file(params.infile))
        dge_edger_kallisto(all_abundances_kallisto_ch, file(params.edger), params.experiment, make_tx2gene_from_gff.out.tx2gene, file(params.infile))
        dge_limma_kallisto(all_abundances_kallisto_ch, file(params.limma), params.experiment, make_tx2gene_from_gff.out.tx2gene, file(params.infile))
        dge_deseq2_kallisto(all_abundances_kallisto_ch, file(params.deseq2), params.experiment, make_tx2gene_from_gff.out.tx2gene, file(params.infile))
        //plot_volcano_edger_kallisto()

        // salmon
        salmon_index(create_ref_transcriptome.out.reference_transcriptome, params.experiment)
        salmon(trim_galore.out.trimmed_reads, salmon_index.out.idx, params.experiment)
        clean_quant(salmon.out.sf, params.experiment)
        all_quants_salmon_ch = clean_quant.out.collect()
        diff_gene_expr_salmon(all_quants_salmon_ch, file(params.edger), file(params.limma), file(params.deseq2), params.experiment, make_tx2gene_from_gff.out.tx2gene, "salmon", file(params.infile))


        // Bowtie2
        bowtie2_index(file(params.fasta), params.experiment)
        bowtie2(trim_galore.out.trimmed_reads, bowtie2_index.out.idx, params.experiment)
        sam_to_bam(bowtie2.out.sam, params.experiment)
        samtools_bowtie2(sam_to_bam.out.bam, params.experiment) // this step primarily marks duplicates
        picard_alignment_summary_bowtie2(samtools_bowtie2.out.marked_bam, file(params.fasta), params.experiment)
        featurecounts_bowtie2(samtools_bowtie2.out.marked_bam, file(params.annotation), params.experiment)
        calculate_tpms_bowtie2(featurecounts_bowtie2.out.counts, file(params.tpm_calc_script), params.experiment)
        clean_features_bowtie2(calculate_tpms_bowtie2.out.complete_table, params.experiment)
        all_feature_counts_bowtie2_ch = clean_features_bowtie2.out.collect()
        diff_gene_expr_bowtie2(all_feature_counts_bowtie2_ch, file(params.edger), file(params.limma), file(params.deseq2), params.experiment, make_tx2gene_from_gff.out.tx2gene, "bowtie2", file(params.infile))


        // STAR and all related processes
        star_index(file(params.fasta), create_gtf.out.gtf, params.experiment)
        star(trim_galore.out.trimmed_reads, star_index.out.idx, params.experiment)
        samtools_star(star.out.bam, params.experiment)
        picard_alignment_summary_star(samtools_star.out.marked_bam, file(params.fasta), params.experiment) // ok until here everything is fine
        featurecounts_star(samtools_star.out.marked_bam, file(params.annotation), params.experiment)


        calculate_tpms_star(featurecounts_star.out.counts, file(params.tpm_calc_script), params.experiment)
        clean_features_star(calculate_tpms_star.out.complete_table, params.experiment)
        all_feature_counts_star_ch = clean_features_star.out.collect()
        diff_gene_expr_star(all_feature_counts_star_ch, file(params.edger), file(params.limma), file(params.deseq2), params.experiment, make_tx2gene_from_gff.out.tx2gene, "star", file(params.infile))


    } else if(annotation_type == "gff"){

        // this part is configured for gff3
        create_gtf(file(params.annotation), params.experiment)
        create_ref_transcriptome(params.experiment, file(params.annotation), file(params.fasta))
        make_tx2gene_from_gff(params.experiment, file(params.annotation), file(params.tx2gene_script))
        //make_genome_txt(file(params.annotation), params.experiment)
        //sort_gtf(create_gtf.out.gtf, params.experiment)


        // kallisto calc
        kallisto_index(create_ref_transcriptome.out.reference_transcriptome, params.experiment) //this should create the index file // works
        kallisto(trim_galore.out.trimmed_reads, kallisto_index.out.idx, params.experiment) //braucht index und fastq als input
        clean_abundance(kallisto.out.tsv, params.experiment)

        // dge kallisto
        /*
        all_abundances_kallisto_ch = clean_abundance.out.collect()
        dge_edger_kallisto(all_abundances_kallisto_ch, file(params.edger), params.experiment, make_tx2gene_from_gff.out.tx2gene, file(params.infile))
        dge_limma_kallisto(all_abundances_kallisto_ch, file(params.limma), params.experiment, make_tx2gene_from_gff.out.tx2gene, file(params.infile))
        dge_deseq2_kallisto(all_abundances_kallisto_ch, file(params.deseq2), params.experiment, make_tx2gene_from_gff.out.tx2gene, file(params.infile))
        plot_volcano_edger_kallisto(dge_edger_kallisto.out.csvs, file(params.plot_edger), params.experiment)
        plot_volcano_limma_kallisto(dge_limma_kallisto.out.csvs, file(params.plot_limma), params.experiment)
        plot_volcano_deseq2_kallisto(dge_deseq2_kallisto.out.csvs, file(params.plot_deseq2), params.experiment)*/

        // salmon
        salmon_index(create_ref_transcriptome.out.reference_transcriptome, params.experiment)
        salmon(trim_galore.out.trimmed_reads, salmon_index.out.idx, params.experiment)
        clean_quant(salmon.out.sf, params.experiment)

        // diff gene expression
        /*
        all_quants_salmon_ch = clean_quant.out.collect()
        diff_gene_expr_salmon(all_quants_salmon_ch, file(params.edger), file(params.limma), file(params.deseq2), params.experiment, make_tx2gene_from_gff.out.tx2gene, "salmon", file(params.infile))*/
        // dge salmon
        /*
        all_abundances_salmon_ch = clean_abundance.out.collect()
        dge_edger_salmon(all_abundances_salmon_ch, file(params.edger), params.experiment, make_tx2gene_from_gff.out.tx2gene, file(params.infile))
        dge_limma_salmon(all_abundances_salmon_ch, file(params.limma), params.experiment, make_tx2gene_from_gff.out.tx2gene, file(params.infile))
        dge_deseq2_salmon(all_abundances_salmon_ch, file(params.deseq2), params.experiment, make_tx2gene_from_gff.out.tx2gene, file(params.infile))
        plot_volcano_edger_salmon(dge_edger_salmon.out.csvs, file(params.plot_edger), params.experiment)
        plot_volcano_limma_salmon(dge_limma_salmon.out.csvs, file(params.plot_limma), params.experiment)
        plot_volcano_deseq2_salmon(dge_deseq2_salmon.out.csvs, file(params.plot_deseq2), params.experiment)*/


        // get coverage of pseudoaligners
        salmon_results_ch = salmon.out.json.collect()
        kallisto_results_ch = kallisto.out.json.collect()
        //pseudo_results_ch = merge(salmon_results_ch, kallisto_results_ch)
        //pseudo_results_ch = merge(kallisto.out.json.collect(), salmon.out.json.collect())
        pseudo_coverage(kallisto_results_ch, salmon_results_ch, file(params.pseudo_coverage), file(params.infile), params.experiment)

        // new Bowtie2 processes
        bowtie2_index(file(params.fasta), params.experiment)
        bowtie2(trim_galore.out.trimmed_reads, bowtie2_index.out.idx, params.experiment) 
        sam_to_bam(bowtie2.out.sam, params.experiment)
        samtools_bowtie2(sam_to_bam.out.bam, params.experiment)
        picard_alignment_summary_bowtie2(samtools_bowtie2.out.marked_bam, file(params.fasta), params.experiment)
        featurecounts_bowtie2(samtools_bowtie2.out.marked_bam, create_gtf.out.gtf, params.experiment)

        //evaluation
        bowtie2_stats(sam_to_bam.out.bam, params.experiment)
        bowtie_calc_genome_cov(params.experiment, sam_to_bam.out.bam, create_gtf.out.gtf)
        bowtie_genome_coverage(bowtie_calc_genome_cov.out.bedgraph, file(params.vis_genome_coverage), params.experiment)
        get_scaffold_order_bowtie(samtools_bowtie2.out.marked_bam, params.experiment)

        bedtools_sort_gtf_bowtie(create_gtf.out.gtf, get_scaffold_order_bowtie.out.genome_txt, params.experiment)
        bowtie_calc_gene_cov(params.experiment, samtools_bowtie2.out.marked_bam, bedtools_sort_gtf_bowtie.out.sorted_gtf, get_scaffold_order_bowtie.out.genome_txt)
        bowtie_gene_coverage(bowtie_calc_gene_cov.out.gene_coverage, file(params.vis_gene_coverage), params.experiment)
        bowtie_mq(bowtie2_stats.out.stats, file(params.vis_mq), params.experiment)

        // small processing steps
        calculate_tpms_bowtie2(featurecounts_bowtie2.out.counts, file(params.tpm_calc_script), params.experiment)
        clean_features_bowtie2(calculate_tpms_bowtie2.out.complete_table, params.experiment)

        // diff gene expression
        /*
        all_feature_counts_bowtie2_ch = clean_features_bowtie2.out.collect()
        diff_gene_expr_bowtie2(all_feature_counts_bowtie2_ch, file(params.edger), file(params.limma), file(params.deseq2), params.experiment, make_tx2gene_from_gff.out.tx2gene, "bowtie2", file(params.infile))*/
        // dge bowtie2
        /*
        all_abundances_bowtie2_ch = clean_features_bowtie2.out.collect()
        dge_edger_bowtie2(all_abundances_bowtie2_ch, file(params.edger), params.experiment, make_tx2gene_from_gff.out.tx2gene, file(params.infile))
        dge_limma_bowtie2(all_abundances_bowtie2_ch, file(params.limma), params.experiment, make_tx2gene_from_gff.out.tx2gene, file(params.infile))
        dge_deseq2_bowtie2(all_abundances_bowtie2_ch, file(params.deseq2), params.experiment, make_tx2gene_from_gff.out.tx2gene, file(params.infile))
        plot_volcano_edger_bowtie2(dge_edger_bowtie2.out.csvs, file(params.plot_edger), params.experiment)
        plot_volcano_limma_bowtie2(dge_limma_bowtie2.out.csvs, file(params.plot_limma), params.experiment)
        plot_volcano_deseq2_bowtie2(dge_deseq2_bowtie2.out.csvs, file(params.plot_deseq2), params.experiment)*/


        // STAR and all related processes
        star_index(file(params.fasta), create_gtf.out.gtf, params.experiment)
        star(trim_galore.out.trimmed_reads, star_index.out.idx, params.experiment)
        samtools_star(star.out.bam, params.experiment)
        picard_alignment_summary_star(samtools_star.out.marked_bam, file(params.fasta), params.experiment) // ok until here everything is fine
        featurecounts_star(samtools_star.out.marked_bam, create_gtf.out.gtf, params.experiment)

        // evaluation of previous steps
        star_stats(star.out.bam, params.experiment)
        star_calc_genome_cov(params.experiment, star.out.bam, create_gtf.out.gtf)
        star_genome_coverage(star_calc_genome_cov.out.bedgraph, file(params.vis_genome_coverage), params.experiment)
        //genome_txt_from_bam(star.out.bam, params.experiment)
        get_scaffold_order(samtools_star.out.marked_bam, params.experiment)
        //sort_gtf_by_scaffold(create_gtf.out.gtf, get_scaffold_order.out.scaffold_order, params.experiment)
        bedtools_sort_gtf(create_gtf.out.gtf, get_scaffold_order.out.genome_txt, params.experiment)
        star_calc_gene_cov(params.experiment, samtools_star.out.marked_bam, bedtools_sort_gtf.out.sorted_gtf, get_scaffold_order.out.genome_txt)
        star_gene_coverage(star_calc_gene_cov.out.gene_coverage, file(params.vis_gene_coverage), params.experiment)
        star_mq(star_stats.out.stats, file(params.vis_mq), params.experiment)

        // smaller processing steps
        calculate_tpms_star(featurecounts_star.out.counts, file(params.tpm_calc_script), params.experiment)
        clean_features_star(calculate_tpms_star.out.complete_table, params.experiment)

        // differential gene expression analysis
        /*
        all_feature_counts_star_ch = clean_features_star.out.collect()
        diff_gene_expr_star(all_feature_counts_star_ch, file(params.edger), file(params.limma), file(params.deseq2), params.experiment, make_tx2gene_from_gff.out.tx2gene, "star", file(params.infile))*/
        // dge star
        /*
        all_abundances_star_ch = clean_features_star.out.collect()
        dge_edger_star(all_abundances_star_ch, file(params.edger), params.experiment, make_tx2gene_from_gff.out.tx2gene, file(params.infile))
        dge_limma_star(all_abundances_star_ch, file(params.limma), params.experiment, make_tx2gene_from_gff.out.tx2gene, file(params.infile))
        dge_deseq2_star(all_abundances_star_ch, file(params.deseq2), params.experiment, make_tx2gene_from_gff.out.tx2gene, file(params.infile))
        plot_volcano_edger_star(dge_edger_star.out.csvs, file(params.plot_edger), params.experiment)
        plot_volcano_limma_star(dge_limma_star.out.csvs, file(params.plot_limma), params.experiment)
        plot_volcano_deseq2_star(dge_deseq2_star.out.csvs, file(params.plot_deseq2), params.experiment)*/


        // tpm comparison for all
        star_tpm_ch = clean_features_star.out.collect()
        bowtie2_tpm_ch = clean_features_bowtie2.out.collect()
        kallisto_tpm_ch = clean_abundance.out.cleaned_tsv.collect()
        salmon_tpm_ch = clean_quant.out.cleaned_sf.collect()
        //all_tpms_ch = merge(kallisto_results_ch, salmon_results_ch, star_results_ch, bowtie2_results_ch)
        compare_tpms(kallisto_tpm_ch, salmon_tpm_ch, bowtie2_tpm_ch, star_tpm_ch, file(params.compare_tpms), file(params.infile), params.experiment)
        kallisto_json_ch = kallisto.out.json.collect()
        salmon_jason_ch = salmon.out.json.collect()
        bowtie2_summary_ch = picard_alignment_summary_bowtie2.out.summary.collect()
        star_summary_ch = picard_alignment_summary_star.out.summary.collect()
        mapping_rates(file(params.mapping_rates_script), params.number_samples, kallisto_json_ch, salmon_jason_ch, bowtie2_summary_ch, star_summary_ch, params.experiment)


    } else {
        println("No fitting annotation file found. Please use a GTF or GFF file")
    }



    // this one works fine
    workflow.onComplete{

        // timestamp traces --> this way ONLY the most recent run is analyzed
        def timer = new Date().format("yyyy-MM-dd_HH-mm-ss")

        // get the latest tracefile...probably easier to just wrap the whole pipeline call in a script and the use bash...
        def tracefiles = file('.').listFiles().findAll {
            it.name.startsWith('trace-') && it.name.endsWith('.txt')
        }
        if(tracefiles.size() > 0){
            def latest = tracefiles.sort {-it.lastModified()}[0]
            def chosen_trace = file("trace_${timer}.txt")
            println("renaming latest trace --> include timestamp")
            "mv ${latest.name} ${chosen_trace}".execute().waitFor()
        }

        // same treatment for the report
        def reportfiles = file('.').listFiles().findAll{
            it.name.startsWith('report-') && it.name.endsWith('.html')
        }
        if (reportfiles.size > 0){
            def latest_report = reportfiles.sort {-it.lastModified()}[0]
            def chosen_report = file("report_${timer}.html")
            println("renaming latest report --> include timestamp")
            "mv ${latest_report} ${chosen_report}".execute().waitFor()
        }
        
    }

}




