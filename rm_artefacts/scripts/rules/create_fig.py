import toolbox, beautifullogger, logging
import pandas as pd, numpy as np, scipy
from typing import Any


logger=logging.getLogger(__name__)
beautifullogger.setup(style={logging.INFO: "Back.CYAN"})

snakemake : Any
input=snakemake.input
output=snakemake.output
params = snakemake.params
file_id = snakemake.wildcards["sample_id"]

import sys,pathlib
sys.path.append(str(pathlib.Path(snakemake.params["_artefact_figure_path"]).parent.parent))
from common.artefact_figure import create_summary_figure

logger.info("Getting summary figure for {}".format(file_id))

raw = scipy.io.loadmat(input["raw"])["RAW"][0,]
cleaned = np.load(input["cleaned"], allow_pickle=False)
summary=pd.read_csv(input["summary"], sep="\t")

f, ax = create_summary_figure(raw, cleaned, summary = summary, 
  **{param: params[param] for param in ["fs", "deviation_factor", "down_sampling"]}
)

f.suptitle('{}\n nb artefacts={}'.format(file_id, len(summary.index)))


logger.info("Saving result")

f.savefig(output["figure"])

logger.info("Summary figure for {} done".format(file_id))
