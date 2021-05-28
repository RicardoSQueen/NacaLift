from matplotlib import rcParams
params = {
  'backend':'ps',
  'ps.usedistiller':'xpdf',
  'axes.labelsize':22,
  'axes.titlesize':20,
  'legend.fontsize':20,
  'xtick.labelsize':20,
  'ytick.labelsize':20,
  'text.usetex':True,
  'text.latex.preamble':r'\usepackage{amsmath,amsfonts,amssymb}',
  'figure.figsize':[8,6],
  'font.family':'serif',
  'font.serif':'Times New Roman'
}
rcParams.update(params)
del params