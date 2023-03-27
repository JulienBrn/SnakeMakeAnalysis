import toolbox

configfile: "config.yaml"

all_files = toolbox.files_as_database(config["input_info"])
run_files = {
  dfname: toolbox.database_select(df, config["run_params"][dfname]["selector"]) 
  if dfname in config["run_params"] and "selector" in config["run_params"][dfname]
  else df
  for dfname, df in all_files.items()
}

mconfig = {
  "date": "date",
  "result_folder": ".",
  "clean_params": "clean_params",
  "figure_params": "figure_params"
}

config_rm_artefact = {
  "samples": run_files["electrophisiology"],
  "config": mconfig
}



module rm_artefacts:
  snakefile: "rm_artefact.smk"
  config: config_rm_artefact


use rule * from rm_artefacts as clean_*