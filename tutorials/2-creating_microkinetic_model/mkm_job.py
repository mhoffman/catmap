from catmap import ReactionModel

import matplotlib
import numpy as np
#matplotlib.use('Agg')
from matplotlib import pyplot as plt
# with golden ration and the whole shebang ...
# settings size and font for revtex stylesheet
fig_width_pt = 2*246.0  # Get this from LaTeX using \showthe\columnwidth
#fig_width_pt *= 300./72 # convert to 300 dpi
inches_per_pt = 1.0/72.27               # Convert pt to inches
#inches_per_pt = 1.0/300               # Convert pt to inches
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height =fig_width*golden_mean       # height in inches
fig_size = [fig_width,fig_height]

font_size = 10
tick_font_size = 10
xlabel_pad = 8
ylabel_pad = 8
#matplotlib.rcParams['ps.usedistiller'] = 'xpdf'
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = 'Computer Modern Times'
matplotlib.rcParams['font.sans-serif'] = 'Gill Sans'
matplotlib.rcParams['text.usetex'] = 'true'

matplotlib.rcParams['lines.linewidth'] = 1.

mkm_file = 'CO_oxidation.mkm'
model = ReactionModel(setup_file=mkm_file)
model.run()

from catmap import analyze
vm = analyze.VectorMap(model)
vm.plot_variable = 'rate' #tell the model which output to plot
vm.log_scale = True #rates should be plotted on a log-scale
vm.min = 1e-48 #minimum rate to plot
vm.max = 1e2 #maximum rate to plot
vm.plot(save='rate.pdf') #draw the plot and save it as "rate.pdf"

vm.unique_only = False
vm.plot(save='all_rates.pdf')
vm.unique_only = True

model.output_variables += ['production_rate', 'rate_constant']
model.run()
vm.production_rate_map = model.production_rate_map #attach map
vm.threshold = 1e-30 #do not plot rates below this
vm.plot_variable = 'production_rate'
vm.plot(save='production_rate.pdf')

vm.descriptor_labels = ['O reactivity [eV]', 'CO reactivity [eV]']
vm.subplots_adjust_kwargs = {'left':0.2,'right':0.8,'bottom':0.15}
vm.plot(save='pretty_production_rate.pdf')

vm.plot_variable = 'coverage'
#vm.plot_mode = 'single'
vm.log_scale = False
vm.min = 0
vm.max = 1
fig = vm.plot(save='coverage.pdf')
fig.subplots_adjust(wspace=.5)
print(dir(fig))
fig.savefig('tight_coverage.pdf', bbox_inches='tight')
plt.show()


vm.include_labels = ['CO_s']
fig = vm.plot(save='CO_coverage.pdf')
fig.savefig('CO_coverage.pdf', bbox_inches='tight')

vm.include_labels = ['O_s']
fig = vm.plot(save='O_coverage.pdf')
fig.savefig('O_coverage.pdf', bbox_inches='tight')

sa = analyze.ScalingAnalysis(model)
sa.plot(save='scaling.pdf')
