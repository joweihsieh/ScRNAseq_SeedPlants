# NOTE: Multi-sample single cell analysis -------------------------------------------
# Integrate SCseq samples
def get_UMI_csv_for_integration(wildcards) -> list:
    species_set = set([
        config["REF_SPECIES"][s.split("_")[1]]
        for s in config["MULTISAMPLES"][wildcards.multisample]
    ])
    if len(species_set) == 1:
        out = expand(
            "results/Single_species_analysis/all_UMI_tables/geneUMI_{sample}.csv",
            sample = config["MULTISAMPLES"][wildcards.multisample]
        )
    else:
        out = expand(
            "results/Ortholog_analysis/all_UMI_tables/orthologUMI_{sample}.csv",
            sample = config["MULTISAMPLES"][wildcards.multisample]
        )
    return out

rule integrate_sample:
    # Integration order is the same as the sample order within multi_UMI_csv
    input:
        multi_UMI_csv = get_UMI_csv_for_integration
    output:
        integration_rds = "results/Multi_species_analysis/all_data_rds/integration_{multisample}_an_{anchor}_ft_{filters}_kan_{k_anchor}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}_base.rds"
    params:
        integrating_order = lambda wildcards: config["SEURAT_INTEGRATING_ORDER"][wildcards.multisample],
        k_params = lambda wildcards: config["SEURAT_K_PARAM"][wildcards.multisample],
        n_feature_cutoff_on_sample = 2000,
        n_barcode_cutoff_on_sample = 100
    script:
        "../scripts/integrate_scseq_data_with_CCA_anchor.R"

rule integrate_sample_runUMAP:
    conda:
        "../envs/PY311.yaml"
    input:
        integration_rds = "results/Multi_species_analysis/all_data_rds/integration_{multisample}_an_{anchor}_ft_{filters}_kan_{k_anchor}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}_base.rds"
    output:
        integrated_umap_csv = "results/Multi_species_analysis/all_plotting_tables/plotting_onlyUMAP_{multisample}_an_{anchor}_ft_{filters}_kan_{k_anchor}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}.csv"
    script:
        "../scripts/integrate_scseq_data_runUMAP.R"

# Turn Seurat CCA results into csv:
## 1. plotting_csv: plotting information of each barcode
## (columns: Barcode, Sample, Species, UMAP.1, UMAP.2, Cluster)
rule turn_CCA_results_into_csv:
    input:
        integrated_rds = "results/Multi_species_analysis/all_data_rds/integration_{multisample}_an_{anchor}_ft_{filters}_kan_{k_anchor}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}_base.rds",
        integrated_umap_csv = "results/Multi_species_analysis/all_plotting_tables/plotting_onlyUMAP_{multisample}_an_{anchor}_ft_{filters}_kan_{k_anchor}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}.csv"
    output:
        plotting_csv = "results/Multi_species_analysis/all_plotting_tables/plotting_{multisample}_an_{anchor}_ft_{filters}_kan_{k_anchor}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}.csv"
    params:
        sample = lambda wildcards: config["MULTISAMPLES"][wildcards.multisample],
        species = lambda wildcards: [
            config["SPECIES"][sample]
            for sample in config["MULTISAMPLES"][wildcards.multisample]
        ]
    script:
        "../scripts/turn_seuratCCA_results_into_csv.R"

# Merge single sample (SS) and multi-sample (MS) csv
rule merge_SS_and_MS_plotting_csv:
    input:
        MS_plotting_csv = "results/Multi_species_analysis/all_plotting_tables/plotting_{multisample}_an_{anchor}_ft_{filters}_kan_{k_anchor}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}.csv",
        SS_plotting_csv_list = lambda wildcards: expand(
            "results/Single_species_analysis/all_plotting_tables/plotting_{sample}.csv",
            sample = config["MULTISAMPLES"][wildcards.multisample]
        )
    output:
        MS_plotting_csv = "results/Multi_species_analysis/all_plotting_tables_addSS/plotting_{multisample}_an_{anchor}_ft_{filters}_kan_{k_anchor}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}.csv"
    script:
        "../scripts/merge_SS_and_MS_plotting_csv.R"

# Special case: plot integrated data with TenX_Ptr single species cluster color
rule plot_integrated_data_with_TenX_Ptr_SScolor:
    input:
        MS_plotting_csv = "results/Multi_species_analysis/all_plotting_tables_addSS/plotting_{multisample}_an_{anchor}_ft_{filters}_kan_{k_anchor}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}.csv"
    output:
        figure_folder = directory("results/Multi_species_analysis/UMAP_with_TenX_Ptr_SScolor/integration_{multisample}_an_{anchor}_ft_{filters}_kan_{k_anchor}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}")
    script:
        "../scripts/plot_integrated_data_with_TenX_Ptr_SScolor.R"

# Plot integrated data by sample
rule plot_integrated_data_by_sample:
    input:
        MS_plotting_csv = "results/Multi_species_analysis/all_plotting_tables_addSS/plotting_{multisample}_an_{anchor}_ft_{filters}_kan_{k_anchor}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}.csv"
    output:
        figure_folder = directory("results/Multi_species_analysis/UMAP_by_sample/integration_{multisample}_an_{anchor}_ft_{filters}_kan_{k_anchor}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}")
    script:
        "../scripts/plot_integrated_data_by_samples.R"

# Plot integrated data by species
rule plot_integrated_data_by_species:
    input:
        MS_plotting_csv = "results/Multi_species_analysis/all_plotting_tables_addSS/plotting_{multisample}_an_{anchor}_ft_{filters}_kan_{k_anchor}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}.csv"
    output:
        figure_folder = directory("results/Multi_species_analysis/UMAP_by_species/integration_{multisample}_an_{anchor}_ft_{filters}_kan_{k_anchor}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}")
    script:
        "../scripts/plot_integrated_data_by_species.R"

# Plot integrated data with correlation between MQD (LCM-seq) and SC-seq data
rule plot_integrated_data_with_correlation:
    input:
        MS_plotting_csv = "results/Multi_species_analysis/all_plotting_tables_addSS/plotting_{multisample}_an_{anchor}_ft_{filters}_kan_{k_anchor}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}.csv",
        SS_MQD_plotting_csv = "results/Single_species_analysis/all_plotting_tables_cor_MQD/plotting_{SC_sample}_cor_{MQD_group}.csv"
    output:
        figure_folder = directory("results/Multi_species_analysis/UMAP_with_correlation/integration_{multisample}_an_{anchor}_ft_{filters}_kan_{k_anchor}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}/{SC_sample}_cor_{MQD_group}")
    script:
        "../scripts/plot_integrated_data_with_correlation.R"


# # FIXME: Finish revision below 
# ## Plot integration results with overlapping heatmap
# rule plot_integration_umap_overlap_heatmap_pair:
#     input:
#         input_integrated_rds = "results/Integrated_rds/{Ref}_and_{Other}.rds"
#     output:
#         output_dir = directory("results/Integrated_figures/Overlap_heatmap/{Ref}_and_{Other}")
#     params:
#         ref_sample = "{Ref}"
#     script:
#         "../scripts/plot_integrated_data_as_overlap_heatmap.R"

# rule plot_integration_umap_overlap_heatmap_multi:
#     input:
#         input_integrated_rds = "results/Integrated_rds/Multi_{Multi}.rds"
#     output:
#         output_dir = directory("results/Integrated_figures/Overlap_heatmap/Multi_{Multi}")
#     params:
#         ref_sample = lambda wildcards: config["MULTI_SAMPLES"][wildcards.Multi][0]
#     script:
#         "../scripts/plot_integrated_data_as_overlap_heatmap.R"

# ## Plot integration results by cellranger cluster on sample
# rule plot_integration_umap_by_cellranger_cluster_on_sample_pair:
#     input:
#         input_integrated_rds = "results/Integrated_rds/{Ref}_and_{Other}.rds",
#         input_dirs = [
#             "results/SCseq_{Ref}",
#             "results/SCseq_{Other}"
#         ]
#     output:
#         output_dir = directory("results/Integrated_figures/UMAP_by_cellranger_cluster_on_sample/{Ref}_and_{Other}")
#     params:
#         cluster_number_list = config["CELLRANGER_N_CLUSTER"],
#         color_list = config["CELLRANGER_CLUSTER_COLOR"]
#     script:
#         "../scripts/plot_integrated_data_by_cellranger_cluster_on_sample.R"

# rule plot_integration_umap_by_cellranger_cluster_on_sample_multi:
#     input:
#         input_integrated_rds = "results/Integrated_rds/Multi_{Multi}.rds",
#         input_dirs = lambda wildcards: ["results/SCseq_" + s for s in config["MULTI_SAMPLES"][wildcards.Multi]]
#     output:
#         output_dir = directory("results/Integrated_figures/UMAP_by_cellranger_cluster_on_sample/Multi_{Multi}")
#     params:
#         cluster_number_list = config["CELLRANGER_N_CLUSTER"],
#         color_list = config["CELLRANGER_CLUSTER_COLOR"]
#     script:
#         "../scripts/plot_integrated_data_by_cellranger_cluster_on_sample.R"

# ## Plot integration results by seurat cluster
# rule plot_integration_umap_by_seurat_cluster_pair:
#     input:
#         input_integrated_rds = "results/Integrated_rds/{Ref}_and_{Other}.rds",
#         input_dirs = [
#             "results/SCseq_{Ref}",
#             "results/SCseq_{Other}"
#         ]
#     output:
#         output_dir = directory("results/Integrated_figures/UMAP_by_seurat_cluster/{Ref}_and_{Other}")
#     params:
#         cluster_number_list = config["CELLRANGER_N_CLUSTER"],
#         color_list = config["CELLRANGER_CLUSTER_COLOR"]
#     script:
#         "../scripts/plot_integrated_data_by_suerat_cluster.R"

# rule plot_integration_umap_by_seurat_cluster_multi:
#     input:
#         input_integrated_rds = "results/Integrated_rds/Multi_{Multi}.rds",
#         input_dirs = lambda wildcards: ["results/SCseq_" + s for s in config["MULTI_SAMPLES"][wildcards.Multi]]
#     output:
#         output_dir = directory("results/Integrated_figures/UMAP_by_seurat_cluster/Multi_{Multi}")
#     params:
#         cluster_number_list = config["CELLRANGER_N_CLUSTER"],
#         color_list = config["CELLRANGER_CLUSTER_COLOR"]
#     script:
#         "../scripts/plot_integrated_data_by_suerat_cluster.R"

# ## Plot integration results by cellranger cluster with lineage
# rule plot_integration_umap_by_cellranger_cluster_with_lineage_multi:
#     input:
#         input_integrated_rds = "results/Integrated_rds/Multi_{Multi}.rds",
#         input_dirs = lambda wildcards: ["results/SCseq_" + s for s in config["MULTI_SAMPLES"][wildcards.Multi]]
#     output:
#         output_dir = directory("results/Integrated_figures/UMAP_by_cellranger_cluster_with_lineage/Multi_{Multi}")
#     params:
#         cluster_number_list = config["CELLRANGER_N_CLUSTER"],
#         color_list = config["CELLRANGER_CLUSTER_COLOR"]
#     script:
#         "../scripts/plot_integrated_data_with_lineage.R"

# ## Plot ortholog distributions of given ortholog cluster id
# rule plot_ortholog_distribution_multi:
#     input:
#         ortholog_long_csv = "results/Ortholog_table/primary_group_long.csv",
#         plot_ortholog_csv = lambda wildcards: config["RAW_DATA_PATH"] + "/20220000_Single_cell_rawdata/20220331_orthogroup_list_revised.csv"
#     output:
#         output_dir = directory("results/Ortholog_figures/Ortholog_distribution")
#     params:
#         species_order = [
#             "PhP", "MaP", "SeM", "PiT", "GnM", "AmT", "OrS",
#             "LiC", "TrA", "PoT", "EuG", "ArT", "CoC", "SoL"
#         ]
#     script:
#         "../scripts/plot_ortholog_distribution.R"

# # Calculate log2-transformed normalized ortholog UMI for each sample
# rule compute_log2_norm_ortholog_expression_each:
#     input:
#         ortholog_long_csv = "results/Ortholog_table/all_group_long_convertedID.csv",
#         input_dirs = lambda wildcards: ["results/SCseq_" + s for s in config["MULTI_SAMPLES"][wildcards.Multi]]
#     output:
#         output_dir = directory("results/Ortholog_UMI_matrix/Each_log2_norm/Multi_{Multi}")
#     script:
#         "../scripts/calculate_ortholog_log2_norm_UMI_counts.R"

# # Plot ortholog expression along lineage
# rule plot_ortholog_expression_along_lineage_multi:
#     input:
#         input_dirs = lambda wildcards: ["results/SCseq_" + s for s in config["MULTI_SAMPLES"][wildcards.Multi]],
#         input_integrated_rds = "results/Integrated_rds/Multi_{Multi}.rds",
#         log2_norm_ortho_dirs = "results/Ortholog_UMI_matrix/Each_log2_norm/Multi_{Multi}",
#         plot_ortholog_csv = lambda wildcards: config["RAW_DATA_PATH"] + "/20220000_Single_cell_rawdata/20220307_orthogroup_list.csv"
#     output:
#         output_dir = directory("results/Ortholog_figures/Ortholog_UMI_along_lineage/Multi_{Multi}")
#     params:
#         cluster_number_list = config["CELLRANGER_N_CLUSTER"],
#         color_list = config["CELLRANGER_CLUSTER_COLOR"]
#     script:
#         "../scripts/plot_ortholog_expression_along_lineage.R"

# # Plot ortholog expression plot of each sample (only keep two zero expression values for each sample)
# rule plot_ortholog_expression_plot_multi:
#     input:
#         log2_norm_ortho_dirs = "results/Ortholog_UMI_matrix/Each_log2_norm/Multi_{Multi}",
#         plot_ortholog_csv = lambda wildcards: config["RAW_DATA_PATH"] + "/20220000_Single_cell_rawdata/20220331_orthogroup_list_revised.csv"
#     output:
#         output_dir = directory("results/Ortholog_figures/Ortholog_expression_plot/Multi_{Multi}")
#     params:
#         sample_list = lambda wildcards: config["MULTI_SAMPLES"][wildcards.Multi]
#     script:
#         "../scripts/plot_ortholog_expression_pie_and_violin.R"