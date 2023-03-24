
import pathlib
import pandas as pd
import logging
import beautifullogger
import matplotlib.pyplot as plt
import sys
import numpy as np
import time
import matplotlib
import joblib
import toolbox
import scipy
import tqdm
import pickle
import sklearn

mlogger=logging.getLogger(__name__)
beautifullogger.setup()

configfile: "config.yaml"

all_files = toolbox.files_as_database(config["input_info"])

run_files = {
  dfname: toolbox.database_select(df, config["run_params"][dfname]["selector"]) 
  if dfname in config["run_params"] and "selector" in config["run_params"][dfname]
  else df
  for dfname, df in all_files.items()
}

outfolder=config["output_info"]["base_folder"]
infolder=config["input_info"]["base_folder"]

electro_db = run_files["electrophisiology"]
electro_db.set_index("file_id", drop=False, inplace=True)

mlogger.info("Preprocessing done")

rule iter_rm_artefacts_figs2:
    input:
      **{file_id: outfolder+"/Cleaned_files/Summary/{file_id}.tsv".format(file_id=file_id) for file_id in electro_db["file_id"]},
      ordered_summary__=outfolder+"/Cleaned_files/ordered_input_by_interest.tsv"
    script: "scripts/iter_artefact_figure.py"

rule iter_rm_artefacts_figs:
    input:
      **{file_id: outfolder+"/Cleaned_files/Summary/{file_id}.pltpickle".format(file_id=file_id) for file_id in electro_db["file_id"]},
      ordered_summary__=outfolder+"/Cleaned_files/ordered_input_by_interest.tsv"
    script: "scripts/pltpickle2.py"

rule create_visualization_order:
  input:
      outfolder+"/Cleaned_files/summary.tsv"
  output:
      outfolder+"/Cleaned_files/ordered_input_by_interest.tsv"
  run:
      df = pd.read_csv(input[0], sep="\t")
      start_col=2
      res=pd.DataFrame()
      scaler = sklearn.preprocessing.MinMaxScaler()
      scaler.fit(df.iloc[:, start_col:], None)
      res=pd.DataFrame(scaler.transform(df.iloc[:, start_col:]), columns=df.columns[start_col:])

      order = toolbox.order_differences(res, 0)
      df["order"] = order
      df.sort_values(by=["order"], inplace=True)
      df.to_csv(output[0], sep="\t", index=False)



rule all_rm_artefacts:
    input:
        **{file_id: outfolder+"/Cleaned_files/Summary/{file_id}.tsv".format(file_id=file_id) for file_id in electro_db["file_id"]},
    output:
        outfolder+"/Cleaned_files/summary.tsv"
    run:
        mlogger.info("Gathering results")
        df = pd.DataFrame()
        df["file_id"]=[str(id) for id, file in input.items()]
        df["file"] = [str(file) for id, file in input.items()]

        def get_features(f):
            fblock=pd.read_csv(f, sep="\t")
            fblock["time"] = fblock["end"] - fblock["start"]
            if len(fblock.index) > 0:
              res = len(fblock.index), fblock["time"].mean(), fblock["time"].std()
            else:
              res = len(fblock.index), 0, 0
              res = len(fblock.index), fblock["time"].mean(), fblock["time"].std()
            return res
        df[["nb_rm_blocks", "avg_rm_time", "rm_time_std"]] = df.apply(lambda r: get_features(r["file"]), axis=1, result_type="expand")
        df.to_csv(output[0], sep="\t", index=False)

rule rm_artefacts:
    input:
        infolder+"/{file_id}.mat"
    output:
        sig =outfolder+"/Cleaned_files/Sig/{file_id}.npy",
        summary = outfolder+"/Cleaned_files/Summary/{file_id}.tsv"
    params:
      fs=lambda wildcards: electro_db.at[wildcards["file_id"], "fs"],
      deviation = config["run_params"]["electrophisiology"]["params"]["rm_artefacts"]["deviation"]
    run:
        mlogger.info("Running rm artefacts on {}".format(input[0]))
        fs = params["fs"]
        raw = scipy.io.loadmat(input[0])["RAW"][0,]
        filtered, tight, bounds =toolbox.replace_artefacts_with_nans(raw, fs)
        cleaned=toolbox.affine_nan_replace(filtered)
        mlogger.info("Writting sig to file")
        np.save(output["sig"], cleaned, allow_pickle=False)
        mlogger.info("Writting summary to file")
        summary = pd.DataFrame([[s/fs, e/fs] for s,e in bounds], columns=["start", "end"])
        summary.to_csv(output["summary"], sep="\t", index=False)
        mlogger.info("Finished Running rm artefacts on {}".format(input[0]))

rule make_artefact_figure:
    input:
      raw = infolder+"/{file_id}.mat",
      summary = outfolder+"/Cleaned_files/Summary/{file_id}.tsv",
      sig =outfolder+"/Cleaned_files/Sig/{file_id}.npy"
    output:
      pltpickle = outfolder+"/Cleaned_files/Summary/{file_id}.pltpickle",
      pngfig = outfolder+"/Cleaned_files/Summary/{file_id}.png"
    params:
      fs=lambda wildcards: electro_db.at[wildcards["file_id"], "fs"],
      deviation = config["run_params"]["electrophisiology"]["params"]["rm_artefacts"]["deviation"]
    script: "scripts/artefact_figure.py"

