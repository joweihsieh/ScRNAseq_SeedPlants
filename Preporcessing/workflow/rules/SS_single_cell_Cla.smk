# Set cellranger bin path
CELLRANGER_BIN = "/home/f06b22037/SSD2/utility/cellranger-7.1.0/cellranger"
# CELLRANGER_BIN = "/home/f06b22037/SSD2/utility/cellranger-5.0.1/cellranger"


# NOTE: Single species (single sample) single cell analysis --------------------------
# Rename fastq of PalChen2021 (for cellranger input)
## e.g. CRR299248_f1.fastq.gz -> CRR299248_S1_L001_R1_001.fastq.gz
##      CRR299248_r2.fastq.gz -> CRR299248_S1_L001_R2_001.fastq.gz
rule rename_PalChen2021_for_cellragner:
    input:
        lambda wildcards: config["OTHER_DATA_PATH"] + "/20220000_Single_cell_rawdata/20211201_scseq_Pal_2021Chen"
    output:
        bark_folder = directory("results/Fastq_renamed_PalChen2021/bark"),
        wood_folder = directory("results/Fastq_renamed_PalChen2021/wood")
    shell:
        """
        mkdir --parents {output.bark_folder}
        mkdir --parents {output.wood_folder}

        for i in {{299240..299247}}
        do
            cp {input}/CRR${{i}}_f1.fastq.gz {output.bark_folder}/CRR${{i}}_S1_L001_R1_001.fastq.gz
            cp {input}/CRR${{i}}_r2.fastq.gz {output.bark_folder}/CRR${{i}}_S1_L001_R2_001.fastq.gz
        done

        for i in {{299248..299255}}
        do
            cp {input}/CRR${{i}}_f1.fastq.gz {output.wood_folder}/CRR${{i}}_S1_L001_R1_001.fastq.gz
            cp {input}/CRR${{i}}_r2.fastq.gz {output.wood_folder}/CRR${{i}}_S1_L001_R2_001.fastq.gz
        done
        """


# Analyze SCseq data from TenX Genomics
rule cellranger_mkref_TenX:
    input:
        genome_fasta = lambda wildcards: config["REF_GENOME_FASTA"][wildcards.species],
        genome_gtf = lambda wildcards: config["REF_GENOME_GTF"][wildcards.species]
    output:
        directory("results/Single_species_analysis/cellranger_ref_TenX_{species}")
    log:
        "logs/" + DATE + "_cellranger_mkref_TenX_{species}.log"
    params:
        main_bin = CELLRANGER_BIN,
        output_folder_name = "cellranger_ref_TenX_{species}"
    threads: 8
    shell:
        """
        mkdir --parents $(dirname {output})
        {params.main_bin} mkref \
            --nthreads={threads} \
            --memgb=128 \
	    --fasta={input.genome_fasta} \
            --genes={input.genome_gtf} \
            --genome={params.output_folder_name} \
            > {log}
        mv {params.output_folder_name} {output}
        """

rule cellranger_count_TenX:
    input:
        cellranger_ref_folder = lambda wildcards: "results/Single_species_analysis/cellranger_ref_TenX_" + config["REF_SPECIES"][wildcards.batch],
        fastq_folder = lambda wildcards: [x.replace("RAW_DATA_PATH", config["RAW_DATA_PATH"]) for x in config["FASTQ_FOLDER"][wildcards.batch]]
    output:
        directory("results/Single_species_analysis/cellranger_count_TenX_{batch}")
    log:
        "logs/" + DATE + "_cellranger_count_TenX_{batch}.log"
    params:
        main_bin = CELLRANGER_BIN,
        sample_name = lambda wildcards: config["CELLRANGER_SAMPLE_NAME"][wildcards.batch],
        output_folder_name = "cellranger_count_TenX_{batch}"
    threads:
        workflow.cores
    shell:
        """
        input_fastqs_dir="{input.fastq_folder}"
        input_fastqs_dir_concated=${{input_fastqs_dir/ /,}}
        {params.main_bin} count \
            --localcores {threads} \
            --localmem 128 \
            --localvmem 128 \
            --transcriptome={input.cellranger_ref_folder} \
            --fastqs=$input_fastqs_dir_concated \
            --sample={params.sample_name} \
            --id={params.output_folder_name} \
            > {log}
        mv {params.output_folder_name} {output}
        """

# `cellranger aggr` does not output normalized UMI count matrix, so we would not use it
# Notice that it also cannot filter barcodes with total UMI counts before aggregation.
# rule cellranger_aggr_Ten:
#     input:
#         cellranger_count_folder = lambda wildcards: expand(
#             "results/Single_species_analysis/cellranger_count_{sample}",
#             sample = config["CELLRANGER_BATCH"][wildcards.aggrsample]
#         )
#     output:
#         directory("results/Single_species_analysis/cellranger_aggr_{aggrsample}")
#     params:
#         main_bin = CELLRANGER_BIN,
#         kit_version = lambda wildcards: config["CELLRANGER_BATCH_KIT"][wildcards.aggrsample],
#         output_folder_name = "cellranger_aggr_{aggrsample}"
#     threads:
#         workflow.cores
#     shell:
#         """
#             bash ../scripts/build_cellranger_aggr_csv.sh \
#                 {input.cellranger_count_folder} \
#                 {params.kit_version} \
#                 > {wildcards.aggrsample}.csv
#             {params.main_bin} aggr \
#                 --localcores {threads} \
#                 --localmem 128 \
#                 --localvmem 128 \
#                 --normalize=none \
#                 --id={params.output_folder_name} \
#                 --csv={wildcards.aggrsample}.csv
#             mv {params.output_folder_name} {output}
#         """

rule cellranger_reanalysis_TenX:
    input:
        cellranger_count_folder = "results/Single_species_analysis/cellranger_count_TenX_{batch}"
    output:
        directory("results/Single_species_analysis/cellranger_reanalysis_TenX_{batch}")
    log:
        "logs/" + DATE + "_cellranger_reanalysis_TenX_{batch}.log"
    params:
        main_bin = CELLRANGER_BIN,
        pars = lambda wildcards: config["CELLRANGER_PARAMS"][wildcards.batch],
        output_folder_name = "cellranger_reanalysis_TenX_{batch}"
    threads:
        workflow.cores
    shell:
        """
        {params.main_bin} reanalyze \
            --localcores {threads} \
            --localmem 128 \
            --localvmem 128 \
            --id {params.output_folder_name} \
            --params {params.pars} \
            --matrix {input.cellranger_count_folder}/outs/filtered_feature_bc_matrix.h5 \
             > {log}
        mv {params.output_folder_name} {output}
        """

# Analyze SCseq data from MARSseq
rule turn_MARSseq_data_into_hd5:
    conda:
        "../envs/PY311.yaml"
    input:
        umi = lambda wildcards: config["MARS_UMI_DIR"][wildcards.batch].replace("RAW_DATA_PATH", config["RAW_DATA_PATH"]),
        cds = lambda wildcards: config["MARS_CDS_PATH"][wildcards.batch].replace("RAW_DATA_PATH", config["RAW_DATA_PATH"])
    output:
        matrix_h5 = "results/h5_for_cellranger/UMI_filtered100_{batch}.h5"
    log:
        "logs/" + DATE + "_turn_MARSseq_data_into_hd5_{batch}.log"
    shell:
        """
        mkdir --parents $(dirname {output}) > {log}
        python workflow/scripts/HD5_generation.py \
            --umi_dir {input.umi} \
            --cds {input.cds} \
            --output_hd5 {output.matrix_h5} \
            --umi_criteria 100 \
            >> {log}
        """

rule cellranger_reanalysis_MARSseq:
    input:
        matrix_h5 = "results/h5_for_cellranger/UMI_filtered100_{batch}.h5"
    output:
        directory("results/Single_species_analysis/cellranger_reanalysis_MARSseq_{batch}")
    log:
        "logs/" + DATE + "_cellranger_reanalysis_MARSseq_{batch}.log"
    params:
        main_bin = CELLRANGER_BIN,
        pars = lambda wildcards: config["CELLRANGER_PARAMS"][wildcards.batch],
        output_folder_name = "cellranger_reanalysis_MARSseq_{batch}"
    threads:
        workflow.cores
    shell:
        """
        mkdir --parents $(dirname {output})
        {params.main_bin} reanalyze \
            --localcores {workflow.cores} \
            --localmem 128 \
            --localvmem 128 \
            --id {params.output_folder_name} \
            --params {params.pars} \
            --matrix {input.matrix_h5} \
             > {log}
        mv {params.output_folder_name} {output}
        """

# Turn cellranger results into csv:
## 1. geneUMI_csv: gene UMI count of each barcode (columns: Barcode, Gene1, Gene2, ...)
## 2. plotting_csv: UMAP and cluster of each barcode (columns: Barcode, UMAP.1, UMAP.2, Cluster, Color)
## 3. norm_factor_csv: normalization factor of each barcode (columns: Barcode, norm_factor)
## p.s. normalized_gene_UMI_count = norm_factor * gene_UMI_count
rule turn_cellranger_results_into_csv:
    input:
        folder = "results/Single_species_analysis/cellranger_reanalysis_{sample}"
    output:
        geneUMI_csv = "results/Single_species_analysis/all_UMI_tables/geneUMI_{sample}.csv",
        plotting_csv = "results/Single_species_analysis/all_plotting_tables/plotting_{sample}.csv",
        norm_factor_csv = "results/Single_species_analysis/all_norm_factor_tables/norm_factor_{sample}.csv"
    params:
        total_UMI_count_cutoff_on_barcode = 500,
        n_cluster = lambda wildcards: config["CELLRANGER_N_CLUSTER"][wildcards.sample],
        cluster_color = lambda wildcards: config["CELLRANGER_CLUSTER_COLOR"][wildcards.sample]
    script:
        "../scripts/turn_cellranger_results_into_csv.R"

# We have tried scMerge for batch-combined normalized UMI count matrix after batch effect correction, but result is not good.
# We haven't found a good way to combine different batch, so we will deal with each sample individually.
# rule integrate_platform_batch:
#     input:
#         geneUMI_csv_list = lambda wildcards: expand(
#             "results/Single_species_analysis/all_UMI_tables/geneUMI_{sample}.csv",
#             sample = config["MULTIBATCH"][wildcards.samplegroup]
#         )
#     output:
#         geneUMI_csv = "results/Single_species_analysis/all_UMI_tables/geneUMI_{samplegroup}.csv",
#         SEG_cutoff = 0.5
#     script:
#         "integrate_scseq_data_with_scMerge.R"

#  Calculate correlation between MQD (LCM-seq) and SC-seq
#rule get_correlation_between_MQD_and_SC:
#    input:
#        geneUMI_csv = "results/Single_species_analysis/all_UMI_tables/geneUMI_{SC_sample}.csv",
#        plotting_csv = "results/Single_species_analysis/all_plotting_tables/plotting_{SC_sample}.csv",
#        gene_abundance_csv = "results/MQD_analysis/quantification/{MQD_group}_gene_abundances.csv"
#    output:
#        plotting_csv = "results/Single_species_analysis/all_plotting_tables_cor_MQD/plotting_{SC_sample}_cor_{MQD_group}.csv"
#    params:
#        cor_method = "pearson"
#    script:
#        "../scripts/get_correlation_between_MQD_and_SC.R"

# Plot single sample data with correlation between MQD and SC data
#rule plot_single_sample_data_with_correlation:
#    input:
#        SS_MQD_plotting_csv = "results/Single_species_analysis/all_plotting_tables_cor_MQD/plotting_{SC_sample}_cor_{MQD_group}.csv"
#    output:
#        figure_folder = directory("results/Single_species_analysis/UMAP_with_correlation/{SC_sample}_cor_{MQD_group}")
#    script:
#        "../scripts/plot_single_sample_data_with_correlation.R"

# # TODO: Revise get_quatification_summary
# # Compute descriptive statistics of the quantification result
# rule get_quatification_summary:
#     input:
#         input_dirs = lambda wildcards: ["results/SCseq_" + s for s in config["MULTI_SAMPLES"][wildcards.Multi]]
#     output:
#         output_csv = "results/Summary_Multi_{Multi}.csv"
#     params:
#         species_list = config["GROUP_For_SUMMARY"]
#     script:
#         "../scripts/get_quatification_summary.R"

# # TODO: Revise plot_cellranger_umap_by_cellranger_cluster
# # Plot single sample umap by single sample cluster
# rule plot_cellranger_umap_by_cellranger_cluster:
#     input:
#         input_reanalysis_dirs = "results/SCseq_{sample}"
#     output:
#         output_dir = directory("results/Cellranger_figures/{sample}")
#     params:
#         cluster_number = lambda wildcards: config["CELLRANGER_N_CLUSTER"][wildcards.sample],
#         color = lambda wildcards: config["CELLRANGER_CLUSTER_COLOR"][wildcards.sample]
#     script:
#         "../scripts/plot_cellranger_data.R"

# # TODO: Revise get_expressed_protein_id
# # Extract protein id of which gene expressed in SCseq data
# rule get_expressed_protein_id:
#     input:
#         protein2gene_csv = "results/protein2gene.csv",
#         input_dirs = lambda wildcards: ["results/SCseq_" + s for s in config["MULTI_SAMPLES"][wildcards.Multi]]
#     output:
#         expressed_protein_id_txt = "results/Expressed_proteinID/Multi_{Multi}.txt"
#     script:
#         "../scripts/get_proteinID_of_expressed_genes.R"