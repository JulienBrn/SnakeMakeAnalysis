#
# Overall parameters: samples
# Overall parameters: version 
# Overall parameters: figure_version, figure_sampling, figure_size
# Overall parameters: view_sampling
# Overall parameters: fs for each sample

# Folder structure:
#   raw:             Input/electrophisiology/{sample}.mat
#   cleaned:         Clean/{version}/signal/{sample}_signal.npy
#   summary:         Clean/{version}/summary/{sample}_summary.tsv
#   png:             Clean/{version}/figures/{sample}_{figure_params}.png
#   overallsummary : Clean/{version}/overallsummary.tsv


rule rm_artefacts:
  input:
    raw
  output:
    cleaned,
    summary
  params:
    version, fs
  run:
    toolbox.rm_artefacts

rule create_fig: (optional)
  input:
    raw,
    cleaned,
    summary
  params:
    figure_sampling,
    figure_size, 
    figure_version,
    version,
    fs
  output:
    png
  run:
    scripts/artefact_figure

rule overallsummary:
  input:
    summary[samples]
  ouput: shadow (output deleted at begining of execution)
    overallsummary
  params: 
    version
  run:
    compute_features_from_individual_summaries

rule all:
  input:
    samples
    overallsummary
    summary[samples]
    cleaned[samples]
  params:
    version

rule complete:
  input:
    samples
    overallsummary
    summary[samples]
    cleaned[samples]
    png[samples]
  params:
    version
    figure.params

rule view:
  input:
    None
  output: 
    None
  params:
    fs, version
  run:
    show_iter_with_refresh with refresh doing the following:
    get_new_samples_with_summary_and_clean_computed
    compute_features_from_individual_summaries
    create_visualization_order (adding the new ones at end)
    update_figure