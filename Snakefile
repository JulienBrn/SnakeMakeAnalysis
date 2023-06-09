import toolbox, datetime

configfile: "config.yaml"

all_files = toolbox.files_as_database(config["input_info"])
run_files = {
  dfname: toolbox.database_select(df, config["run_params"]["selection"][dfname]) 
  if "selection" in config["run_params"] and dfname in config["run_params"]["selection"]
  else df
  for dfname, df in all_files.items()
}

samples_df= run_files["electrophisiology"]
samples_df.set_index("sample_id", drop=False, inplace=True)

rmartconfig = config["run_params"]["exec_params"]["rm_artefacts"]
rmartconfig["file_structure"]["result_folder"] = config["output_info"]["base_folder"] + "/"+ rmartconfig["file_structure"]["result_folder"]
rmartconfig["run_info"]["date"] = toolbox.roundTime(datetime.datetime.now(), 1)

config_rm_artefact = {
  "samples": samples_df,  
  "configuration": rmartconfig
}


module rm_artefacts:
  snakefile: "rm_artefacts/rm_artefacts2.smk"
  config: config_rm_artefact


use rule * from rm_artefacts as clean_*

print(rules.clean_all.ouput)