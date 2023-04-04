import toolbox, pathlib

samples_df= config["samples"]
mconfig = config["configuration"]
clean_params_str = toolbox.clean_filename(str(config["configuration"]["exec_params"]))
figure_params_str = toolbox.clean_filename(str(config["configuration"]["figure_params"]))

mconfig = {
  "result_folder": config["configuration"]["file_structure"]["result_folder"],
  "clean_params_str": clean_params_str,
  "date": toolbox.clean_filename(str(config["configuration"]["run_info"]["date"])),
  "figure_params_str": figure_params_str,
}

def get_value_from_id(sample_id, col):
  return samples_df.at[sample_id, col]


samples_path = "{result_folder}/{clean_params_str}/run_info_{date}/samples.tsv".format(**mconfig)
raw_path = lambda wildcards: get_value_from_id(wildcards["sample_id"], "sample_path")
cleaned_path = "{result_folder}/{clean_params_str}/results/{sample_id}_signal.npy"
summary_path = "{result_folder}/{clean_params_str}/results/{sample_id}_summary.tsv"
figure_path = "{result_folder}/{clean_params_str}/figures/{sample_id}_fp={figure_params_str}.png"
overallsummary_path = "{result_folder}/{clean_params_str}/run_info_{date}/overallsummary.tsv".format(**mconfig)

samples_df["cleaned"] = samples_df.apply(lambda row: cleaned_path.format(sample_id=row["sample_id"], **mconfig), axis=1)
samples_df["summary"] = samples_df.apply(lambda row: summary_path.format(sample_id=row["sample_id"], **mconfig), axis=1)
samples_df["figure"] = samples_df.apply(lambda row: figure_path.format(sample_id=row["sample_id"], **mconfig), axis=1)

all_summaries_paths = [summary_path.format(sample_id=sample_id, **mconfig) for sample_id in samples_df["sample_id"]]
all_cleaned_paths = [cleaned_path.format(sample_id=sample_id, **mconfig) for sample_id in samples_df["sample_id"]]
all_figures_paths = [figure_path.format(sample_id=sample_id, **mconfig) for sample_id in samples_df["sample_id"]]

fs_value = lambda wildcards: get_value_from_id(wildcards["sample_id"], "fs")
clean_params = config["configuration"]["exec_params"]
clean_params["test"] = {"nested":2}
figure_params =config["configuration"]["figure_params"]

pathlib.Path(samples_path).parent.mkdir(parents=True, exist_ok=True)  
samples_df.to_csv(samples_path, sep="\t")
# rule make_sample_file:
#   params:
#     **{samples_df.at[sample_id, "sample_path"]:sample_id for sample_id in samples_df["sample_id"]}
#   output:
#     samples = samples_path
#   script: "scripts/rules/make_sample_file.py"

rule rm_artefacts:
  input:
    raw = raw_path
  params:
    **clean_params,
    fs = fs_value
  output:
    cleaned = cleaned_path,
    summary = summary_path
  script: "scripts/rules/rm_artefacts.py"
  


rule create_fig:
  input:
    raw = raw_path,
    cleaned = cleaned_path,
    summary = summary_path
  params:
    **{**figure_params, **clean_params, "fs": fs_value, 
      "_artefact_figure_path": 
        workflow.source_path(str(pathlib.Path(workflow.snakefile).parent)
          +"/scripts/common/artefact_figure.py")
    }
  output:
    figure = figure_path
  script: "scripts/rules/create_fig.py"

rule overallsummary:
  input:
    samples = samples_path,
    summaries = all_summaries_paths,
  params:
     get_features_path = workflow.source_path(
        str(pathlib.Path(workflow.snakefile).parent)+
        "/scripts/common/get_features.py"
      )
  output:
    overallsummary = overallsummary_path
  script: "scripts/rules/overallsummary.py"

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
    **{**clean_params, 
        "_get_features_path": workflow.source_path(
            str(pathlib.Path(workflow.snakefile).parent)+
            "/scripts/common/get_features.py"
          ),
        "_artefact_figure_path": 
        workflow.source_path(str(pathlib.Path(workflow.snakefile).parent)
          +"/scripts/common/artefact_figure.py")
      }
  script: "scripts/rules/view.py"
    
    # show_iter_with_refresh with refresh doing the following:
    # get_new_samples_with_summary_and_clean_computed that are within samples
    # compute_features_from_individual_summaries
    # create_visualization_order (adding the new ones at end)
    # update_figure
  
