from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = 'Helvetica'
rcParams['font.size'] = 9
rcParams['font.size'] = 9
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'sans'
rcParams['mathtext.bf'] = 'sans:bold'
rcParams['mathtext.it'] = 'sans:italic'
rcParams['mathtext.cal'] = 'sans:italic'
rcParams['mathtext.default'] = 'rm'
rcParams['xtick.major.size'] = 4
rcParams['xtick.major.width'] = 1
rcParams['ytick.major.size'] = 4
rcParams['ytick.major.width'] = 1
rcParams['axes.grid'] = False
rcParams['axes.linewidth'] = 1
rcParams['grid.linewidth'] = 1
rcParams['grid.linestyle'] = '-'
rcParams['grid.alpha'] = .15
rcParams['savefig.dpi'] = 150
rcParams['text.color'] = 'w'
rcParams['figure.facecolor'] = 'k'
rcParams['savefig.facecolor'] = 'k'
rcParams['axes.facecolor'] = 'k'
rcParams['axes.edgecolor'] = 'w'
rcParams['axes.labelcolor'] = 'w'
rcParams['xtick.color'] = 'w'
rcParams['ytick.color'] = 'w'

from matplotlib.pyplot import errorbar as eb
def errorbar(*args, **kwargs):
    return eb(*args, **{
        'ls': 'None',
        'marker': 'None',
        'capsize': 2,
        'capthick': 1,
        'ecolor': 'w',
        'elinewidth': 1,
        'zorder': -100,
        **kwargs})
