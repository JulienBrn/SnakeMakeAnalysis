import matplotlib
import matplotlib.pyplot as plt
import logging, beautifullogger
import pickle
import sys
from typing import Any
import pandas as pd
import re
import toolbox

snakemake : Any

beautifullogger.setup(logfile="logtest.txt", style={logging.INFO: "Back.CYAN"})

mlogger=logging.getLogger(__name__)
df = pd.read_csv(snakemake.input.ordered_summary__, sep="\t")

def make_fig(i):
  id= df["file_id"].iat[i]
  pltpickle = snakemake.input[id]
  with open(pltpickle, "rb") as file:
    fig = pickle.load(file)
  feature_str=", ".join([col+"="+str(df[col].iat[i]) for col in df.columns[2:-1]])
  fig.suptitle('{}\n{}'.format(id, feature_str))
  return fig

def index_ok(i):
  id= df["file_id"].iat[i]
  try:
    _ =  snakemake.input[id]
  except:
    return False
  return True


toolbox.make_figure_list_interaction(make_fig, max_len = len(df.index), index_ok = index_ok)