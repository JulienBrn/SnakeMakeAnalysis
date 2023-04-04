import scipy, logging, numpy as np, matplotlib.pyplot as plt, pandas as pd, pickle
import beautifullogger, matplotlib

def create_summary_figure(raw, cleaned, fs, summary, down_sampling, deviation_factor):
  df=pd.DataFrame()
  df["raw"] = raw
  df["cleaned"] = cleaned
  df["t"] = df.index/fs
  mean = raw.mean()
  rawmin= raw.min()
  rawmax= raw.max()

  cleaned_mean = cleaned.mean()
  cleaned_std = cleaned.std()
  deviation = deviation_factor
  std = raw.std()
  ds= down_sampling

  f = plt.figure()

  ax = f.subplots(3, sharey='all', sharex='all')

  ax[0].plot(df["t"][::ds], df["raw"][::ds], label="raw", color="C0")
  ax[0].legend(loc='upper right')
  ax[0].set_xlabel("time (s)")
  ax[0].set_ylabel("signal value")

  ax[1].plot(df["t"][::ds], df["raw"][::ds], label="raw", color="C0")
  ax[1].hlines([mean-std*deviation, mean+std*deviation], 0, df["t"].iat[-1], color="C1", label="raw mean +/- {}*std".format(deviation))
  ax[1].hlines([cleaned_mean-cleaned_std*deviation, cleaned_mean+cleaned_std*deviation], 0, df["t"].iat[-1], color="C4", label="cleaned mean +/- {}*std".format(deviation))
  ax[1].vlines(summary["start"], rawmin, rawmax, color="C2", label="start artefact")
  ax[1].vlines(summary["end"], rawmin, rawmax, color="C3", label="end artefact")
  ax[1].legend(loc='upper right')
  ax[1].set_xlabel("time (s)")
  ax[1].set_ylabel("signal value")


  ax[2].plot(df["t"][::ds], df["cleaned"][::ds], color="C1", label="cleaned")
  ax[2].legend(loc='upper right')
  ax[2].set_xlabel("time (s)")
  ax[2].set_ylabel("signal value")
  return f, ax