
samples_df= config["samples"]
samples_df["sample_id"] = samples_df["file_id"]
samples_df["sample_path"] = samples_df["full_path"]
samples_df.set_index("sample_id", drop=False, inplace=True)

mconfig = config["config"]


def get_value_from_id(sample_id, col):
  return samples_df.at[sample_id, col]

samples_path = "{result_folder}/{clean_params}/run_info_{date}/samples.txt".format(**mconfig)
raw_path = lambda wildcards: get_value_from_id(wildcards["sample_id"], "sample_path")
cleaned_path = "{result_folder}/{clean_params}/results/".format(**mconfig)+ "{sample_id}_signal.npy"
summary_path = "{result_folder}/{clean_params}/results/".format(**mconfig)+ "{sample_id}_summary.tsv"
figure_path = "{result_folder}/{clean_params}/figures/".format(**mconfig)+"{sample_id}"+"_fp={figure_params}.png".format(**mconfig)
overallsummary_path = "{result_folder}/{clean_params}/run_info_{date}/overallsummary.tsv".format(**mconfig)

all_summaries_paths = [summary_path.format(sample_id=sample_id) for sample_id in samples_df["sample_id"]]
all_cleaned_paths = [cleaned_path.format(sample_id=sample_id) for sample_id in samples_df["sample_id"]]
all_figures_paths = [figure_path.format(sample_id=sample_id) for sample_id in samples_df["sample_id"]]

fs_value = lambda wildcards: get_value_from_id(wildcards["sample_id"], "fs"),
clean_params={"version": 1, "deviation": 7}
figure_params = {"fig_version": 1, "wdith": 10, "height": 100}


rule make_sample_file:
  output:
    samples = samples_path
  run:
    with open(output["sample_file_path"], "w") as f:
      for sample in samples_df["sample_id"]:
        f.write(sample+"\n")

rule rm_artefacts:
  input:
    raw = raw_path
  params:
    clean_params,
    fs = fs_value,
  output:
    cleaned = cleaned_path,
    summary = summary_path,
  run:
    pass #toolbox.rm_artefacts


rule create_fig:
  input:
    raw = raw_path,
    cleaned = cleaned_path,
    summary = summary_path
  params:
    figure_params,
    fs = fs_value
  output:
    figure = figure_path
  run:
    scripts/artefact_figure

rule overallsummary:
  input:
    samples = samples_path,
    summaries = all_summaries_paths
  output:
    overallsummary = overallsummary_path
  run:
      pass #compute_features_from_individual_summaries

rule all:
  input:
    samples = samples_path,
    overallsummary = overallsummary_path,
    summaries = all_summaries_paths,
    all_cleaned = all_cleaned_paths
  default_target: True

rule complete:
  input:
    samples = samples_path,
    overallsummary = overallsummary_path,
    summaries = all_summaries_paths,
    all_cleaned = all_cleaned_paths,
    figures = all_figures_paths

rule view:
  input:
    samples = samples_path,
  params:
    clean_params
  run:
    pass
    # show_iter_with_refresh with refresh doing the following:
    # get_new_samples_with_summary_and_clean_computed that are within samples
    # compute_features_from_individual_summaries
    # create_visualization_order (adding the new ones at end)
    # update_figure
  
