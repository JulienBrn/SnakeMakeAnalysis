import toolbox, beautifullogger, logging
import pandas as pd, numpy as np, scipy
from typing import Any

logger=logging.getLogger(__name__)
beautifullogger.setup(style={logging.INFO: "Back.CYAN"})

snakemake : Any
input=snakemake.input
output=snakemake.output
params = snakemake.params

logger.info("Running rm artefacts on {}".format(input["raw"]))
logger.info("Params are: {}".format({k:v for k,v in params.items()}))
raw = scipy.io.loadmat(input["raw"])["RAW"][0,]
filtered, tight, bounds =toolbox.replace_artefacts_with_nans2(raw, 
  **{param: params[param] for param in ["fs", "deviation_factor", "min_length", "shoulder_width", "join_width", "recursive"]}
)

if params["replace_type"]=="affine":
  cleaned=toolbox.affine_nan_replace(filtered)
elif params["replace_type"]=="keep_nans":
  cleaned=filtered
else:
  logger.error("{} is a unknown replace_type in params. Continuing with affine replace.".format(params["replace_type"]))
  cleaned=toolbox.affine_nan_replace(filtered)

logger.info("Writting sig to file")
np.save(output["cleaned"], cleaned, allow_pickle=False)
logger.info("Writting summary to file")
summary = pd.DataFrame([[s/params.fs, e/params.fs] for s,e in bounds], columns=["start", "end"])
summary.to_csv(output["summary"], sep="\t", index=False)
logger.info("Finished Running rm artefacts on {}".format(input["raw"]))