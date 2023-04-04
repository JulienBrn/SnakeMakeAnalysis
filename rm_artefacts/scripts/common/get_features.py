import pandas as pd

def get_features(f):
  fblock=pd.read_csv(f, sep="\t")
  fblock["time"] = fblock["end"] - fblock["start"]
  if len(fblock.index) > 0:
    res = len(fblock.index), fblock["time"].mean(), fblock["time"].std()
  else:
    res = len(fblock.index), fblock["time"].mean(), fblock["time"].std()
  return res