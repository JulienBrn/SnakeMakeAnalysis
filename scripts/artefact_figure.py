import scipy, logging, numpy as np, matplotlib.pyplot as plt, pandas as pd, pickle
import beautifullogger, matplotlib
from typing import Any

def create_summary_figure(raw, cleaned, summary, rd):
  df=pd.DataFrame()
  df["raw"] = raw
  df["cleaned"] = cleaned
  df["t"] = df.index/snakemake.params["fs"]
  mean = raw.mean()
  rawmin= raw.min()
  rawmax= raw.max()
  deviation = snakemake.params["deviation"]
  std = raw.std()

  f = plt.figure()

  ax = f.subplots(3, sharey='all', sharex='all')

  ax[0].plot(df["t"][::rd], df["raw"][::rd], label="raw", color="C0")
  ax[0].legend(loc='upper right')
  ax[0].set_xlabel("time (s)")
  ax[0].set_ylabel("signal value")

  ax[1].plot(df["t"][::rd], df["raw"][::rd], label="raw", color="C0")
  ax[1].hlines([mean-std*deviation, mean+std*deviation], 0, df["t"].iat[-1], color="C1", label="mean +/- {}*std".format(deviation))
  ax[1].vlines(summary["start"], rawmin, rawmax, color="C2", label="start artefact")
  ax[1].vlines(summary["end"], rawmin, rawmax, color="C3", label="end artefact")
  ax[1].legend(loc='upper right')
  ax[1].set_xlabel("time (s)")
  ax[1].set_ylabel("signal value")


  ax[2].plot(df["t"][::rd], df["cleaned"][::rd], color="C1", label="cleaned")
  ax[2].legend(loc='upper right')
  ax[2].set_xlabel("time (s)")
  ax[2].set_ylabel("signal value")

if __name__ == "__main__":
  snakemake: Any

  beautifullogger.setup(logfile="logtest.txt", style={logging.INFO: "Back.CYAN"})
  matplotlib.use("tkagg")
  mlogger=logging.getLogger(__name__)
  mlogger.info("Getting summary figure for {}".format(snakemake.wildcards.file_id))

  raw = scipy.io.loadmat(snakemake.input["raw"])["RAW"][0,]
  sig = np.load(snakemake.input["sig"], allow_pickle=False)
  summary=pd.read_csv(snakemake.input["summary"], sep="\t")

  f = create_summary_figure(raw, sig, summary, 10)
  # print(summary)
  # df=pd.DataFrame()
  # df["raw"] = raw
  # df["cleaned"] = sig
  # df["t"] = df.index/snakemake.params["fs"]
  # mean = raw.mean()
  # rawmin= raw.min()
  # rawmax= raw.max()
  # deviation = snakemake.params["deviation"]
  # std = raw.std()

  # f = plt.figure()
  f.suptitle('{}\n nb artefacts={}'.format(snakemake.wildcards.file_id, len(summary.index)))

  # ax = f.subplots(3, sharey='all', sharex='all')
  # rd=10



  # ax[0].plot(df["t"][::rd], df["raw"][::rd], label="raw", color="C0")
  # ax[0].legend(loc='upper right')
  # ax[0].set_xlabel("time (s)")
  # ax[0].set_ylabel("signal value")

  # ax[1].plot(df["t"][::rd], df["raw"][::rd], label="raw", color="C0")
  # ax[1].hlines([mean-std*deviation, mean+std*deviation], 0, df["t"].iat[-1], color="C1", label="mean +/- {}*std".format(deviation))
  # ax[1].vlines(summary["start"], rawmin, rawmax, color="C2", label="start artefact")
  # ax[1].vlines(summary["end"], rawmin, rawmax, color="C3", label="end artefact")
  # ax[1].legend(loc='upper right')
  # ax[1].set_xlabel("time (s)")
  # ax[1].set_ylabel("signal value")


  # ax[2].plot(df["t"][::rd], df["cleaned"][::rd], color="C1", label="cleaned")
  # ax[2].legend(loc='upper right')
  # ax[2].set_xlabel("time (s)")
  # ax[2].set_ylabel("signal value")

  mlogger.info("Saving result")

  with open(snakemake.output["pltpickle"], 'wb') as ffig: 
    pickle.dump(f, ffig) 
  f.savefig(snakemake.output["pngfig"])

  mlogger.info("Summary figure for {} done".format(snakemake.wildcards.file_id))
