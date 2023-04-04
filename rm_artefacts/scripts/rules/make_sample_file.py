import toolbox, beautifullogger, logging
import pandas as pd, numpy as np, scipy
from typing import Any

logger=logging.getLogger(__name__)
beautifullogger.setup(style={logging.INFO: "Back.CYAN"})

snakemake : Any
input=snakemake.input
output=snakemake.output
params = snakemake.params

df = pd.DataFrame(params.items(), columns=["sample_path", "sample_id"])
df.to_csv(output["samples"], sep="\t", index=False)