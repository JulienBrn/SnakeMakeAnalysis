import toolbox, beautifullogger, logging
import pandas as pd, numpy as np, scipy
from typing import Any
import sys, pathlib
import time


snakemake : Any
input=snakemake.input
output=snakemake.output
params = snakemake.params

sys.path.append(str(pathlib.Path(snakemake.params["get_features_path"]).parent.parent))
from common.get_features import get_features

logger=logging.getLogger(__name__)
beautifullogger.setup(style={logging.INFO: "Back.CYAN"})

logger.info("Gathering results")
df = pd.read_csv(input["samples"], sep="\t")

df["summary_path"] = [str(f) for f in input["summaries"]]
df["ok"] = df.apply(lambda r: r["sample_id"] in r["summary_path"], axis=1)
all_ok = df["ok"].all()
if not all_ok:
   logger.warning("Join key seems off when computing all summaries. Keeping all columns")
if df.empty:
   logger.error("No files on which to gather results... strange")
   raise BaseException("No files on which to gather results... strange")
else:
  df[["nb_rm_blocks", "avg_rm_time", "rm_time_std"]] = df.apply(lambda r: get_features(r["summary_path"]), axis=1, result_type="expand")
  if all_ok:
     df.drop(columns=["summary_path", "sample_path", "ok"], inplace=True)
  df.to_csv(output["overallsummary"], sep="\t", index=False)