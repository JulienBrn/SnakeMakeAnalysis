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
# matplotlib.use("tkagg")
mlogger=logging.getLogger(__name__)
df = pd.read_csv(snakemake.input.ordered_summary__, sep="\t")
mexit=False
i=-1


def call_exit(_):
  global mexit
  plt.close()
  mexit= True

def incr_i():
  global i
  tmp_i = i+1
  mexit=False
  while tmp_i < len(df.index):
    id= df["file_id"].iat[tmp_i]
    try:
      x = snakemake.input[id]
    except:
      mlogger.warning("Problem trying to display file {}, passing".format(id))
    else:
      break
    tmp_i+=1
  if tmp_i==len(df.index):
    return False
  else:
    i = tmp_i
    return True
  
def decr_i():
  global i
  tmp_i = i-1
  while tmp_i >= 0:
    id= df["file_id"].iat[tmp_i]
    try:
      x = snakemake.input[id]
    except:
      mlogger.warning("Problem trying to display file {}, passing".format(id))
    else:
      break
    tmp_i-=1
  if tmp_i==-1:
    return False
  else:
    i = tmp_i
    return True
  
def next_figure(_):
  global i 
  if incr_i():
    plt.close()

def prev_figure(_):
  global i 
  if decr_i():
    plt.close()

rnum = re.compile("(\d+)/{}".format(len(df.index)))

def go_to_figure(t):
  global i 
  matching = rnum.fullmatch(t)
  if matching:
    print(matching.group(1))
    inext = int(str(matching.group(1)))
    i=inext-1
    plt.close()



def show():
  global i
  id= df["file_id"].iat[i]
  pltpickle = snakemake.input[id]
  with open(pltpickle, "rb") as file:
    fig = pickle.load(file)
  ax_exitbutton = plt.axes([0.90, 0.97, 0.10, 0.03])
  exitbtn = matplotlib.widgets.Button(ax_exitbutton, 'Exit all', color='white', hovercolor='grey')
  exitbtn.on_clicked(call_exit)
  
  ax_numtext = plt.axes([0.045, 0.97, 0.06, 0.03])
  textbox = matplotlib.widgets.TextBox(ax_numtext, "", initial="{}/{}".format(i+1, len(df.index)))
  textbox.on_submit(go_to_figure)

  ax_nextbutton = plt.axes([0.11, 0.97, 0.04, 0.03])
  nextbtn = matplotlib.widgets.Button(ax_nextbutton, '=>', color='white', hovercolor='grey')
  nextbtn.on_clicked(next_figure)

  ax_prevbutton = plt.axes([0.00, 0.97, 0.04, 0.03])
  prevbtn = matplotlib.widgets.Button(ax_prevbutton, '<=', color='white', hovercolor='grey')
  prevbtn.on_clicked(prev_figure)

  ax_selecttoggle = plt.axes([0.8, 0.9, 0.10, 0.03])
  selecttoggle= matplotlib.widgets.CheckButtons(ax_selecttoggle, ["Remember figure"])

  manager = plt.get_current_fig_manager()
  manager.window.showMaximized()
  mlogger.warning("Problem trying to display file {}, passing".format(id))
  plt.show()

if incr_i():
  while not mexit:
    show()
else:
  mlogger.error("No figure to display")