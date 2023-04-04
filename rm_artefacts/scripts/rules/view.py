import toolbox, beautifullogger, logging
import pandas as pd, numpy as np, scipy, pathlib, sys, sklearn, time
import matplotlib, matplotlib.pyplot as plt
from typing import Any

logger=logging.getLogger(__name__)
beautifullogger.setup(style={logging.INFO: "Back.CYAN"})

snakemake : Any
input=snakemake.input
output=snakemake.output
params = snakemake.params


sys.path.append(str(pathlib.Path(snakemake.params["_get_features_path"]).parent.parent))
sys.path.append(str(pathlib.Path(snakemake.params["_artefact_figure_path"]).parent.parent))
from common.get_features import get_features
from common.artefact_figure import create_summary_figure

clean_params = {k:p for k,p in params.items() if k != "_get_features_path" and k!= "_artefact_figure_path"}
clean_params_str= ", ".join([key+"="+str(val) for key,val in clean_params.items()])
#Initial setting

df = pd.read_csv(input["samples"], sep="\t")

df["ok"] = df.apply(lambda r: r["sample_id"] in r["summary"] and r["sample_id"] in r["cleaned"], axis=1)
all_ok = df["ok"].all()
if not all_ok:
   logger.error("Join key seems off. Exiting")
   raise BaseException("Join key seems off. Exiting")

df["Done"] = False
df["Update"] = False
df["nb_rm_blocks"] = None
df["avg_rm_time"] = None
df["rm_time_std"] = None

feature_cols = ["nb_rm_blocks", "avg_rm_time", "rm_time_std"]
view_down_sampling=100


def load_fig(i):
  mdf = df.loc[df["Done"], :]
  if not mdf["Done"].iat[i]:
    logger.error("load fig")
    return None
  logger.info("Loading figure {}".format(i))
  [raw_path, cleaned_path, fs, summary_path] = mdf[["sample_path", "cleaned", "fs", "summary"]].iloc[i, :].tolist()
  raw = scipy.io.loadmat(raw_path)["RAW"][0,]
  cleaned = np.load(cleaned_path, allow_pickle=False)
  summary=pd.read_csv(summary_path, sep="\t")
  return [raw, cleaned, fs, summary, i]

ax_dsslider = None
slider = None

ax_updatebtn = None
updatebtn = None
fl = None

def change_view_sampling(val):
  global fl,view_down_sampling
  view_down_sampling=val
  fl.reload(erase_data=False)

view_samplings = [1, 5, 10, 50, 100, 500, 1000]

def refresh_callback():
  global refresh_requested, fl
  refresh_requested=True
  fl.exit_all()

def draw_fig(data):
  global fl, ax_dsslider, slider, ax_updatebtn, updatebtn
  mdf = df.loc[df["Done"], :]
  [raw, cleaned, fs, summary, i] = data
  fig, ax = create_summary_figure(raw, cleaned, fs, summary, down_sampling=view_down_sampling, deviation_factor= clean_params["deviation_factor"])
  feature_str=", ".join([col+"="+str(mdf[col].iat[i]) for col in df[feature_cols].columns])
  fig.suptitle('Sample: {}\nExtracted features: {}\nClean params: {}'.format(mdf["sample_id"].iat[i], feature_str, clean_params_str))
  

  ax_dsslider = plt.axes([0.0, 0.85, 0.08, 0.10])
  slider = matplotlib.widgets.RadioButtons(ax_dsslider, labels=[str(mdf["fs"].iat[i]/ds) + "  view fs" for ds in view_samplings], active=view_samplings.index(view_down_sampling))
  slider.on_clicked(lambda val: change_view_sampling(int(mdf["fs"].iat[i]/float(val[:-len("view fs")]))))

  ax_updatebtn = plt.axes([0.16, 0.97, 0.02, 0.03])
  updatebtn = matplotlib.widgets.Button(ax_updatebtn, "o", color='white', hovercolor='grey')
  updatebtn.on_clicked(lambda _: refresh_callback())
  return fig

# def make_fig(i):
#   global df
#   logger.info("loading")
#   [raw_path, cleaned_path, fs, summary_path] = df[["sample_path", "cleaned", "fs", "summary"]].iloc[i, :].tolist()
#   raw = scipy.io.loadmat(raw_path)["RAW"][0,]
#   cleaned = np.load(cleaned_path, allow_pickle=False)
#   summary=pd.read_csv(summary_path, sep="\t")
#   logger.info("drawing")
#   fig, ax = create_summary_figure(raw, cleaned, fs, summary, down_sampling=view_down_sampling, deviation_factor= clean_params["deviation_factor"])
#   feature_str=", ".join([col+"="+str(df[col].iat[i]) for col in df[feature_cols].columns])
#   fig.suptitle('Sample: {}\nExtracted features: {}\nClean params: {}'.format(df["sample_id"].iat[i], feature_str, clean_params_str))
#   logger.info("drawing done")
#   return fig

# def index_ok(i):
#   return df["Done"].iat[i]


#Refresh
refresh_requested = True
def refresh():
  global fl, refresh_requested
  df["Update"] = df.apply(lambda r: not r["Done"] and pathlib.Path(r["summary"]).exists(), axis=1)
  df[feature_cols] = df.apply(lambda r: get_features(r["summary"]) if not r["Done"] and pathlib.Path(r["summary"]).exists() else r[feature_cols], axis=1, result_type="expand")
  

  #Computing order
  scaler = sklearn.preprocessing.MinMaxScaler()
  scaler.fit(df[feature_cols], None)
  res=pd.DataFrame(scaler.transform(df[feature_cols]), columns=["norm_"+f for f in feature_cols])

  df["Done"] = df["Done"] | df["Update"]
  df["Update"] = False

  order = toolbox.order_differences(res.loc[df["Done"], :], 0)
  df.loc[df["Done"], "order"] = order
  df.sort_values(by=["order"], inplace=True, ignore_index=True)
  print(df)
  #Reinitializing columns
  

  
  
  #


while(refresh_requested):
  refresh()
  fl = toolbox.FigureList(load_fig, draw_fig, len(df[df["Done"]].index))
  refresh_requested=False
  fl.show()




# toolbox.make_figure_list_interaction(make_fig, max_len = len(df.index), index_ok = index_ok)



