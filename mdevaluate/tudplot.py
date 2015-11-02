import matplotlib
matplotlib.use('Agg')

# Official TUD Colors
colorsa = ('#5D85C3', '#009CDA', '#50B695', '#AFCC50', '#DDDF48', '#FFE05C',
'#F8BA3C', '#EE7A34', '#E9503E', '#C9308E', '#804597')
colorsb = ('#005AA9', '#0083CC', '#009D81', '#99C000', '#C9D400', '#FDCA00',
'#F5A300', '#EC6500', '#E6001A', '#A60084', '#721085')
colorsc = ('#004E8A', '#00689D', '#008877', '#7FAB16', '#B1BD00', '#D7AC00',
'#D28700', '#CC4C03', '#B90F22', '#951169', '#611C73')
colorsd = ('#243572', '#004E73', '#00715E', '#6A8B22', '#99A604', '#AE8E00',
'#BE6F00', '#A94913', '#961C26', '#732054', '#4C226A')


import seaborn as sns
sns.reset_defaults()
sns.set_palette(sns.color_palette(colorsb, len(colorsb)), len(colorsb))
sns.set_style({'font.family': 'Charter',
               'lines.linewidth': 1.5,
               'lines.markeredgesize': 5,
               'lines.markersize': 6,
               'figure.figsize': (6,4)})
small_palette = sns.color_palette([colorsb[0],colorsb[8], colorsb[10],
    colorsb[3], colorsb[6]], 5)
full_palette = sns.color_palette(colorsc, len(colorsb))
#matplotlib.rcParams['markeredgecolor']  = 'none'

from matplotlib import pyplot
