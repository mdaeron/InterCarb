#! /usr/bin/env python3

'''
Data processing for the InterCarb exercise
'''

__author__	= 'Mathieu Daëron'
__contact__   = 'daeron@lsce.ipsl.fr'
__copyright__ = 'Copyright (c) 2020 Mathieu Daëron'
__license__   = 'Modified BSD License - https://opensource.org/licenses/BSD-3-Clause'
__date__	  = '2020-10'

VERBOSE =			 False
RUN_ETH1234_VS_HEG = True
RUN_INTERCARB =		 True
SAVE_RAWDATA =		 True

import os
from glob import glob
from pylab import *
from D47crunch_snapshot import *
from scipy.stats import norm, kstest, chi2
import matplotlib.patches as patches

### 18O/16O ACID FRACTIONATION AS A FUNCTION OF ACID TEMPERATURE
### (Kim et al., 2015) <https://doi.org/10.1016/j.gca.2015.02.011>
alpha18_acid_reaction_calcite = lambda T: exp( 3.48 / (T + 273.15) - 1.47e-3 )

UNKNOWNS = ['ETH-4', 'IAEA-C1', 'IAEA-C2', 'MERCK']
ANCHORS = ['ETH-1', 'ETH-2', 'ETH-3']
SAMPLES = ANCHORS + UNKNOWNS

F95_1df = chi2.ppf(.95, 1)**.5
F95_2df = chi2.ppf(.95, 2)**.5

def save_rawdata( rawdata, filename = 'rawdata.csv', mode = 'w', sep = ',', anonymize = True ):
	'''
	Save a list of analyses to a CSV file
	'''
	fields = [f.split(':') for f in [
		'UID:s',
		'Lab:s',
		'LabSession:s:Session',
		'Instrument:s',
		'Acid T:.0f',
		'Sample:s',
		'Mass:.3f',
		'd13Cwg_VPDB:.3f',
		'd18Owg_VSMOW:.3f',
		'd45:.6f',
		'd46:.6f',
		'd47:.6f',
		'd13C_VPDB:.3f',
		'd18O_VSMOW:.3f',
		'D47raw:.6f',
		'D47:.6f',
		]]
	fieldnames = [ f[2] if len(f) == 3 else f[0] for f in fields ]
	with open(filename, mode) as fid:
		if mode != 'a':
			fid.write(sep.join([f for f in fieldnames]))
		for k,r in enumerate(rawdata) :
			fid.write('\n' + sep.join([f'{r[f[0]]:{f[1]}}' for f in fields]))
	if VERBOSE:
		print(f'Wrote {len(rawdata)} records to "{filename}".')


def create_tree(path):
	'''
	Creates nested directories, ignoring the last component of the path unless it ends with a slash.

	path = "foo/bar/baz.txt" will create "foo/bar/".
	path = "foo/bar/" will also create "foo/bar/".
	path = "foo/bar" will only create "foo/".
	'''
	splitpath = path.split('/')[:-1]
	for k,x in enumerate(splitpath):
		if not os.path.exists('/'.join(splitpath[:k+1])):
			os.makedirs('/'.join(splitpath[:k+1]))


def summary_of_sessions(labdata, path = 'output/InterCarb/Table_S1_InterCarb_summary.csv'):
	'''
	Generate Table S1
	'''
	create_tree(path)
	with open(path, 'w') as f:
		f.write(','.join([
			'Lab', 'Session',
			'N (ETH-1)', 'N (ETH-2)', 'N (ETH-3)', 'N (ETH-4)', 'N (IAEA-C1)', 'N (IAEA-C2)', 'N (MERCK)',
			'Nf',
			'WG δ13C_VPDB', 'WG δ18O_VSMOW',
			'a', 'b', 'c',
			'δ13C repeatability (ppm)', 'δ18O repeatability (ppm)', 'Δ47 repeatability (ppm)',
			]))
		for lab in sorted(labdata):
			for s in labdata[lab]:
				session = labdata[lab][s].sessions[s]
				f.write(','.join([
					f'\n{lab[-2:]}', f'{s[-2:]}',
					f'{len([r for r in session["data"] if r["Sample"] == "ETH-1"])}',
					f'{len([r for r in session["data"] if r["Sample"] == "ETH-2"])}',
					f'{len([r for r in session["data"] if r["Sample"] == "ETH-3"])}',
					f'{len([r for r in session["data"] if r["Sample"] == "ETH-4"])}',
					f'{len([r for r in session["data"] if r["Sample"] == "IAEA-C1"])}',
					f'{len([r for r in session["data"] if r["Sample"] == "IAEA-C2"])}',
					f'{len([r for r in session["data"] if r["Sample"] == "MERCK"])}',
					f'{len(session["data"]) - len({r["Sample"] for r in session["data"]})}',
					f'{session["d13Cwg_VPDB"]:.2f}', f'{session["d18Owg_VSMOW"]:.2f}',
					f'{session["a"]:.2f}',
					f'{session["b"]:.1e}' if abs(session["b"]/session["SE_b"]) > F95_1df else f'({session["b"]:.1e})',
					f'{session["c"]:.3f}',
					f'{session["r_d13C_VPDB"]*1000:.0f}', f'{session["r_d18O_VSMOW"]*1000:.0f}', f'{session["r_D47"]*1000:.1f}',
					]))
				

def compute_old_to_new_conversion(final_values):
	'''
	Compute the conversion equation in section 3.6
	'''
	a,b,c = inv(array([[1, 0.010, 0.258], [1, -28.375, 0.256], [1, 0.538, 0.691]])) @ array([[final_values['ETH-1']['D47']],[final_values['ETH-2']['D47']],[final_values['ETH-3']['D47']]])
	print(f'newΔ47 = {a[0]:.6f} - {-b[0]:.6f} δ47 + {c[0]:.6f} oldΔ47')


def MS_effects(labdata, InterCarb_results, path = 'output/InterCarb/'):
	'''
	Generate Table 5 and Figure 9
	'''

	labs_with_Thermo = [lab for lab in InterCarb_results if lab not in UNKNOWNS and 'Thermo' in labdata[lab][lab+'Session01'][0]['Instrument']]
	labs_with_Nu = [lab for lab in InterCarb_results if lab not in UNKNOWNS and 'Nu' in labdata[lab][lab+'Session01'][0]['Instrument']]
	labs_with_Isoprime = [lab for lab in InterCarb_results if lab not in UNKNOWNS and 'Isoprime' in labdata[lab][lab+'Session01'][0]['Instrument']]
	labs_by_MSmodel = {
		'Thermo': labs_with_Thermo,
		'Nu': labs_with_Nu,
		'Isoprime': labs_with_Isoprime,
		}

	N_Thermo = len(labs_with_Thermo)
	N_Nu = len(labs_with_Nu)
	N_Isoprime = len(labs_with_Isoprime)

	msg = f'''
Number of mass spectrometers = {N_Thermo + N_Nu + N_Isoprime}
Number of Thermo MS = {N_Thermo}
Number of Nu MS = {N_Nu}
Number of Isoprime MS = {N_Isoprime}
'''[1:-1]
	create_tree(path)
	with open(path + 'MS_manufacturers.txt', 'w') as f:
		f.write(msg)

	fig = figure(figsize = (8.5,3))
	subplots_adjust(.08,.2,.96,.9,.33,.33)

	out = {u:{} for u in UNKNOWNS}
	for ax, xmodel, ymodel in [
		(subplot(131), 'Isoprime', 'Thermo'),
		(subplot(132), 'Thermo', 'Nu'),
		(subplot(133), 'Nu', 'Isoprime'),
		]:
		
		chisq = 0
		sca(ax)
		
		for sample in UNKNOWNS:
			xlabs = [lab for lab in labs_by_MSmodel[xmodel] if sample in InterCarb_results[lab]]
			ylabs = [lab for lab in labs_by_MSmodel[ymodel] if sample in InterCarb_results[lab]]

			X = [InterCarb_results[lab][sample]['D47'] for lab in xlabs]
			sX = [InterCarb_results[lab][sample]['SE_D47'] for lab in xlabs]

			Y = [InterCarb_results[lab][sample]['D47'] for lab in ylabs]
			sY = [InterCarb_results[lab][sample]['SE_D47'] for lab in ylabs]
			
			Xavg, sXavg = w_avg(X, sX)
			Yavg, sYavg = w_avg(Y, sY)
			chisq += (Xavg-Yavg)**2 / (sXavg**2 + sYavg**2)
			out[sample][f'{ymodel} vs {xmodel}'] = ((Yavg-Xavg), sqrt(sXavg**2 + sYavg**2))

			gca().add_artist(
				patches.Ellipse(
					xy = (Xavg,Yavg), width = 2*F95_2df*sXavg, height = 2*F95_2df*sYavg,
					lw = 1, fc = 'none', ec = 'k', ls = '-', alpha = 1 )
					)
		xmin,xmax = .325, .775
		p = 1-chi2.cdf(chisq,4)
		text(.05, .95, f'RMSWD = {sqrt(chisq/4):.2f}\n(p = {p:.2f})', va = 'top', ha = 'left', size = 9, transform = gca().transAxes)
		plot([xmin,xmax], [xmin,xmax], 'r-', lw = .75, zorder = -100)
		axis([xmin,xmax,xmin,xmax])
		xticks([.4,.5,.6,.7])
		yticks([.4,.5,.6,.7])
		xlabel(xmodel)
		ylabel(ymodel)

	savefig(path + 'Fig_9_MS_effects.pdf')
	close(fig)
	
	with open(path + 'Table_5_InterCarb_MS_effects.csv', 'w') as f:
		f.write('Sample,MAT 253 vs Isoprime 100,Nu Perspective vs MAT 253,Isoprime 100 vs Nu Perspective')
		for u in UNKNOWNS:
			f.write(f'\n{u},{out[u]["Thermo vs Isoprime"][0]:.4f} ± {out[u]["Thermo vs Isoprime"][1]:.4f},{out[u]["Nu vs Thermo"][0]:.4f} ± {out[u]["Nu vs Thermo"][1]:.4f},{out[u]["Isoprime vs Nu"][0]:.4f} ± {out[u]["Isoprime vs Nu"][1]:.4f}')
		f.write('\naverage (all samples)')
		for k in ["Thermo vs Isoprime", "Nu vs Thermo", "Isoprime vs Nu"]:
			X = sum([out[u][k][0] for u in UNKNOWNS])/len(UNKNOWNS)
			sX = sqrt(sum([out[u][k][1]**2 for u in UNKNOWNS]))/len(UNKNOWNS)
			f.write(f',{X:.4f} ± {sX:.4f}')


def acid_T_effects(labdata, InterCarb_results, path = 'output/InterCarb/'):
	'''
	Generate Table 4 and Figure 8
	'''

	labs_with_25C_acid = [lab for lab in InterCarb_results if lab not in UNKNOWNS and labdata[lab][lab+'Session01'][0]['Acid T'] == 25]
	labs_with_70C_acid = [lab for lab in InterCarb_results if lab not in UNKNOWNS and labdata[lab][lab+'Session01'][0]['Acid T'] == 70]
	labs_with_90C_acid = [lab for lab in InterCarb_results if lab not in UNKNOWNS and labdata[lab][lab+'Session01'][0]['Acid T'] == 90]
	labs_by_acid_T = {
		25: labs_with_25C_acid,
		70: labs_with_70C_acid,
		90: labs_with_90C_acid,
		}

	N_25C = len(labs_with_25C_acid)
	N_70C = len(labs_with_70C_acid)
	N_90C = len(labs_with_90C_acid)

	msg = f'''
Number of mass spectrometers = {N_25C + N_70C + N_90C}
Number of labs using 90 ºC acid = {N_90C}
Number of labs using 70 ºC acid = {N_70C}
Number of labs using 25 ºC acid = {N_25C}
'''[1:-1]
	create_tree(path)
	with open(path + 'Acid_T.txt', 'w') as f:
		f.write(msg)

	fig = figure(figsize = (3, 3))
	subplots_adjust(.18,.15,.98,.95)

	out = {u:{} for u in UNKNOWNS}
	
	chisq = 0
	
	for sample in UNKNOWNS:
		xlabs = [lab for lab in labs_by_acid_T[70] if sample in InterCarb_results[lab]]
		ylabs = [lab for lab in labs_by_acid_T[90] if sample in InterCarb_results[lab]]

		X = [InterCarb_results[lab][sample]['D47'] for lab in xlabs]
		sX = [InterCarb_results[lab][sample]['SE_D47'] for lab in xlabs]

		Y = [InterCarb_results[lab][sample]['D47'] for lab in ylabs]
		sY = [InterCarb_results[lab][sample]['SE_D47'] for lab in ylabs]
		
		Xavg, sXavg = w_avg(X, sX)
		Yavg, sYavg = w_avg(Y, sY)
		chisq += (Xavg-Yavg)**2 / (sXavg**2 + sYavg**2)
		out[sample][70] = {'D47': Xavg, 'SE_D47': sXavg}
		out[sample][90] = {'D47': Yavg, 'SE_D47': sYavg}

		gca().add_artist(
			patches.Ellipse(
				xy = (Xavg,Yavg), width = 2*F95_2df*sXavg, height = 2*F95_2df*sYavg,
				lw = 1, fc = 'none', ec = 'k', ls = '-', alpha = 1 )
				)

	xmin,xmax = .325, .775
	p = 1-chi2.cdf(chisq,4)
	text(.05, .95, f'RMSWD = {sqrt(chisq/4):.2f}\n(p = {p:.2f})', va = 'top', ha = 'left', size = 9, transform = gca().transAxes)
	plot([xmin,xmax], [xmin,xmax], 'r-', lw = .75, zorder = -100)
	axis([xmin,xmax,xmin,xmax])
	xticks([.4,.5,.6,.7])
	yticks([.4,.5,.6,.7])
	xlabel('70 ºC acid reaction')
	ylabel('90 ºC acid reaction')

	savefig(path + 'Fig_8_Acid_T_effects.pdf')
	close(fig)

	with open(path + 'Table_4_InterCarb_Acid_T_effects.csv', 'w') as f:
		f.write('Sample,Δ47 (70 ºC acid),Δ47 (90 ºC acid),Difference')
		for u in UNKNOWNS:
			f.write(f'\n{u},{out[u][70]["D47"]:.4f} ± {out[u][70]["SE_D47"]:.4f},{out[u][90]["D47"]:.4f} ± {out[u][90]["SE_D47"]:.4f},{out[u][90]["D47"]-out[u][70]["D47"]:.4f} ± {sqrt(out[u][70]["SE_D47"]**2 + out[u][90]["SE_D47"]**2):.4f}')
		f.write('\naverage (all samples),,')
		X = sum([out[u][90]["D47"]-out[u][70]["D47"] for u in UNKNOWNS])/len(UNKNOWNS)
		sX = sqrt(sum([out[u][90]["SE_D47"]**2 + out[u][70]["SE_D47"]**2 for u in UNKNOWNS]))/len(UNKNOWNS)
		f.write(f',{X:.4f} ± {sX:.4f}')


def save_InterCarb_results(InterCarb_results, path = 'output/InterCarb/'):
	'''
	Write detailed InterCarb results to Table 3
	'''
	create_tree(path)
	with open(f'{path}Table_S3_InterCarb_results.csv', 'w') as fid:
		with open(f'{path}Table_3_InterCarb_results.csv', 'w') as fid2:
			fid.write('Lab,Session,' + ','.join(UNKNOWNS))
			fid2.write('Lab,' + ','.join(UNKNOWNS))
			for lab in sorted(InterCarb_results):
				if lab not in UNKNOWNS:
					for session in sorted(InterCarb_results[lab]):
						if session not in UNKNOWNS:
							fid.write(f'\n{lab[-2:]},{session[-2:]}')
							for u in UNKNOWNS:
								try:
									fid.write(f',{InterCarb_results[lab][session][u]["D47"]:.4f} ± {InterCarb_results[lab][session][u]["SE_D47"]:.4f}')
								except KeyError:
									fid.write(f',—')
					fid.write(f'\n{lab[-2:]},w. avg')
					fid2.write(f'\n{lab[-2:]}')
					for u in UNKNOWNS:
						try:
							fid.write(f',{InterCarb_results[lab][u]["D47"]:.4f} ± {InterCarb_results[lab][u]["SE_D47"]:.4f}')
							fid2.write(f',{InterCarb_results[lab][u]["D47"]:.4f} ± {InterCarb_results[lab][u]["SE_D47"]:.4f}')
						except KeyError:
							fid.write(f',—')
							fid2.write(f',—')
			fid.write('\nall,w. avg')
			fid2.write('\nw. avg')
			for u in UNKNOWNS:
				fid.write(f',{InterCarb_results[u]["D47"]:.4f} ± {InterCarb_results[u]["SE_D47"]:.4f}')
				fid2.write(f',{InterCarb_results[u]["D47"]:.4f} ± {InterCarb_results[u]["SE_D47"]:.4f}')
		

def ETH1234_vs_HEG():
	'''
	Process the data for ETH-1/2/3/4 vs heated/equilibrated CO2
	'''
	if RUN_ETH1234_VS_HEG:

		if VERBOSE:
			print('''
===================
ETH-1/2/3/4 vs H/EG
===================''')

		ethsamples = [f'ETH-{k+1}' for k in range(4)]
		labs = sorted([g.split('/')[-1] for g in glob('input/ETH1234_vs_HEG/*')])
	
		labinfo = {}
	
		for lab in labs:

			labdir = f'input/ETH1234_vs_HEG/{lab}'
			if glob(f'{labdir}/rawdata.csv'):
				print(f'\nReading data from {lab}')
			
				rawdata = D47data(verbose = VERBOSE)

				with open(f'{labdir}/anchors.csv') as f:
					anchors = [l.strip().split(',') for l in f.readlines()]
				rawdata.Nominal_D47 = {k: float(v) for k,v in anchors}
				rawdata.read(f'{labdir}/rawdata.csv')
				print(f'Read {len(rawdata)} records from {lab}')

				with open(f'{labdir}/acid_T.txt') as f:
					acid_T = float(f.read())
				if acid_T == 90:
					a18_acid = 1.008129 # calcite reacted at 90 °C, [Kim et al., 2007]
				elif acid_T == 70:
					a18acid = exp(3.59 / (acid_T + 273.15) - 1.79e-3)
				elif acid_T == 25:
					a18_acid = 1.01025
				rawdata.wg(a18_acid = a18_acid)

				for session in rawdata.sessions:
					rawdata.sessions[session]['d13C_standardization_method'] = 'none'
					rawdata.sessions[session]['d18O_standardization_method'] = 'none'

				rawdata.crunch()

				if glob(f'{labdir}/drifts.txt'):
					with open(f'{labdir}/drifts.txt') as f:
						drifts = [l.strip() for l in f.readlines()]
					for session in rawdata.sessions:
						for drift in drifts:
							rawdata.sessions[session][drift] = True

				if not os.path.exists('output'):
					os.makedirs('output')
				if not os.path.exists('output/ETH1234_vs_HEG'):
					os.makedirs('output/ETH1234_vs_HEG')
				if not os.path.exists(f'output/ETH1234_vs_HEG/{lab}'):
					os.makedirs(f'output/ETH1234_vs_HEG/{lab}')
				if not os.path.exists(f'output/ETH1234_vs_HEG/{lab}/session_plots'):
					os.makedirs(f'output/ETH1234_vs_HEG/{lab}/session_plots')
				if not os.path.exists(f'output/ETH1234_vs_HEG/{lab}/tables'):
					os.makedirs(f'output/ETH1234_vs_HEG/{lab}/tables')

				rawdata.standardize(method = 'pooled')
				rawdata.plot_sessions(f'output/ETH1234_vs_HEG/{lab}/session_plots')
				rawdata.table_of_sessions(f'output/ETH1234_vs_HEG/{lab}/tables')
				rawdata.table_of_samples(f'output/ETH1234_vs_HEG/{lab}/tables')
				rawdata.table_of_analyses(f'output/ETH1234_vs_HEG/{lab}/tables')

				labinfo[lab] = {
					'rawdata': rawdata,
					'N_of_ETH_analyses': len([r for r in rawdata if r['Sample'] in ethsamples]),
					'N_of_CO2_analyses': len([r for r in rawdata if r['Sample'] in rawdata.anchors]),
					'N_of_sessions': len(rawdata.sessions),
					'acid_T': acid_T,
					'rD47': rawdata.repeatability['r_D47'],
					}
			
				for sample in ethsamples:
					if sample in rawdata.unknowns:
						labinfo[lab][sample] = {
							'N': len(rawdata.unknowns[sample]['data']),
							'D47': rawdata.unknowns[sample]['D47'],
							'sD47': rawdata.unknowns[sample]['SE_D47'],
							'eD47': rawdata.unknowns[sample]['SE_D47'] * rawdata.t95,
							}
						if acid_T == 70:
							labinfo[lab][sample]['D47'] += .066 # from Petersen et al. (2019)
						elif acid_T == 90:
							labinfo[lab][sample]['D47'] += .088 # from Petersen et al. (2019)


		total_N_of_analyses = sum([labinfo[lab]['N_of_CO2_analyses'] + labinfo[lab]['N_of_ETH_analyses'] for lab in labinfo])
		total_N_of_CO2_analyses = sum([labinfo[lab]['N_of_CO2_analyses'] for lab in labinfo])
		total_N_of_ETH_analyses = sum([labinfo[lab]['N_of_ETH_analyses'] for lab in labinfo])
		print(f'\nTotal number of {total_N_of_analyses} analyses from {len(labinfo)} labs.\n')
		print(f'\nTotal number of {total_N_of_CO2_analyses} CO2 analyses.\n')
		print(f'\nTotal number of {total_N_of_ETH_analyses} ETH analyses.\n')
		for sample in ethsamples:
			sum_of_weights = sum([labinfo[lab][sample]['sD47']**-2 for lab in labinfo if sample in labinfo[lab]])
			for lab in labinfo:
				if sample in labinfo[lab]:
					labinfo[lab][sample]['wD47'] = labinfo[lab][sample]['sD47']**-2 / sum_of_weights

		for lab in labinfo:
			print(f'\n{lab}:')
			for k in labinfo[lab]:
				if k != 'rawdata':
					print(f'{k} = {labinfo[lab][k]}')

		with open('output/ETH1234_vs_HEG/Table_2_ETH1234_vs_HEG.csv', 'w') as fid:
			fid.write('Lab,Sessions,N (H/E CO2),N (ETH-1), Δ47 (ETH-1), ± 95 %,Stat. weight (ETH-1),N (ETH-2), Δ47 (ETH-2), ± 95 %,Stat. weight (ETH-2),N (ETH-3), Δ47 (ETH-3), ± 95 %,Stat. weight (ETH-3),N (ETH-4), Δ47 (ETH-4), ± 95 %,Stat. weight (ETH-4)')
			for lab in labinfo:
				fid.write(f'\n{lab},{labinfo[lab]["N_of_sessions"]},{labinfo[lab]["N_of_CO2_analyses"]},'
					+ ','.join([f'{labinfo[lab][s]["N"]},{labinfo[lab][s]["D47"]:.4f},{labinfo[lab][s]["eD47"]:.4f},{labinfo[lab][s]["wD47"]:.3f}' for s in ethsamples]))

		kw_errorbar = dict(
			ecolor = 'k',
			ls = 'None',
			marker = 'None',
			elinewidth = 1.2,
			capthick = 1.2,
			capsize = 3,
			zorder = -100,
			)

		fig = figure(figsize = (7,5))
		subplots_adjust(.02, .1, .98, .95, .1, .25)
		ax = {s: subplot(221+k) for k,s in enumerate(['ETH-4', 'ETH-3', 'ETH-2', 'ETH-1'])}

		final_values = {}
		for sk, s in enumerate(ax):

			X = [labinfo[lab][s]['D47'] for lab in labinfo]
			sX = [labinfo[lab][s]['sD47'] for lab in labinfo]
			Xo, sXo = w_avg(X, sX)
			final_values[s] = {'D47': Xo, 'sD47': sXo}
		
			sca(ax[s])
			for k, lab in enumerate(sorted(labinfo, key = lambda l: labinfo[l][s]['D47'])):
				X = labinfo[lab][s]['D47']
				Y = k
				eX = labinfo[lab][s]['eD47']
				errorbar(X, Y, xerr = eX, **kw_errorbar)
			
				eX_autogenic = labinfo[lab]['rD47'] / labinfo[lab][s]['N']**.5
				gca().add_patch(Rectangle(
					(X - F95_1df * eX_autogenic, Y - 0.125), 2 * F95_1df * eX_autogenic, .25,
					linewidth = 1.2,
					edgecolor = 'k',
					facecolor = 'w',
	# 				zorder = 100
					))
			
			
				if X <= Xo:
					text(X + eX + 0.003, Y-0.025, lab,
						size = 8,
						ha = 'left',
						va = 'center',
						)
				else:
					text(X - eX - 0.003, Y-0.025, lab,
						size = 8,
						ha = 'right',
						va = 'center',
						)
			yticks([])
			axis([None, None, -1, len(labinfo)])
			if sk > 1:
				xlabel('Δ$_{47}$ ($25\,$°C acid)')
			text(0, 1.01, s, va = 'bottom', weight = 'bold', size = 12, transform = ax[s].transAxes)
		
			axvline(Xo, color = 'r', lw = 1, zorder = -200)
			text(Xo, len(labinfo)+.1, f'{Xo:.4f}  \n', color = 'r', size = 9, va = 'center', ha = 'right', weight = 'bold')
			text(Xo, len(labinfo)+.1, '±\n', color = 'r', size = 9, va = 'center', ha = 'center', weight = 'bold')
			text(Xo, len(labinfo)+.1, f'  {sXo:.4f} (1SE)\n', color = 'r', size = 9, va = 'center', ha = 'left', weight = 'bold')

			weights = {
				T: sum([labinfo[lab][s]['wD47'] for lab in labinfo if labinfo[lab]['acid_T'] == T])
				for T in sorted({labinfo[lab]['acid_T'] for lab in labinfo})
				}
			txt = 'Statistical weights:'
			for T in weights:
				txt += f'\n{T:g}$\,$°C acid: {100*weights[T]:.0f} %'
			text(.03, .95, txt, va = 'top', transform = gca().transAxes, size = 8)
	
			ax[s].xaxis.set_major_locator(MultipleLocator(.050))
			ax[s].xaxis.set_minor_locator(MultipleLocator(.010))

		savefig(f'output/ETH1234_vs_HEG/Fig_3_ETH1234_vs_HEG.pdf')
		close(fig)


		all_sigma_values = []
		fig = figure(figsize = (6.5,3.5))
		subplots_adjust(.06, .2, .975, .95, .25, .1)
		myaxes = [subplot(240+k) for k in [3,4,7,8]]

		for k, u in enumerate(['ETH-4', 'ETH-3', 'ETH-2', 'ETH-1']):
			D47, sD47 = final_values[u]['D47'], final_values[u]['sD47']

			sigma_values = sorted(
				[
					(labinfo[l][u]['D47'] - D47) / labinfo[l][u]['sD47']
					for l in labinfo if u in labinfo[l]
					]
				)
			all_sigma_values += sigma_values
			N = len(sigma_values)

			col1, col2 = [0]*3, [.4]*3

			sca(myaxes[k])
			X = array(sigma_values)
			Y = arange(N)/(N-1)
			plot(X, Y, '-', lw = .75, color = col1)
			plot(X, Y, 'wo', mec = col1, mew = .75, ms = 3.5)
			x1,x2,y1,y2 = axis()
		
			x1 = -5.5
			x2 = 5.5
			y1 = -0.1
			y2 = 1.1
		
			x = linspace(x1, x2)
			y = norm().cdf(x)
			plot(x,y,'-', lw = 2, zorder = -10, color = [.6,.8,1])
			pvalue = kstest(X, 'norm', (0, 1)).pvalue
			text(.95, .05, f'p = {100*pvalue:{".0f" if pvalue>0.1 else ".1f"}} %', va = 'bottom', ha = 'right', size = 10, transform = gca().transAxes)
			axis([x1,x2,y1,y2])
			text(.1, .9, u, va = 'top', transform = gca().transAxes, weight = 'bold')
			if k // 2:
				xlabel('Weighted Δ$_{47}$ deviation\nof each laboratory')
			else:
				xticks([])
			if k % 2 == 0:
				ylabel('Cumulative\ndistribution', labelpad = -8)
				yticks([0,1])
			else:
				yticks([])
		
		subplot(121)
		sigma_values = sorted(all_sigma_values)
		N = len(sigma_values)

		col1, col2 = [0]*3, [.4]*3

		X = array(sigma_values)
		Y = arange(N)/(N-1)
		plot(X, Y, '-', lw = .75, color = col1)
		plot(X, Y, 'wo', mec = col1, mew = .75, ms = 3.5)
		x1,x2,y1,y2 = axis()

		x1 = -5.5
		x2 = 5.5
		y1 = -0.1
		y2 = 1.1

		x = linspace(x1, x2)
		y = norm().cdf(x)
		plot(x,y,'-', lw = 2, zorder = -10, color = [.6,.8,1])
		pvalue = kstest(X, 'norm', (0, 1)).pvalue
		text(.95, .05, f'p = {100*pvalue:{".0f" if pvalue>0.1 else ".1f"}} %', va = 'bottom', ha = 'right', size = 12, transform = gca().transAxes)
		axis([x1,x2,y1,y2])
		text(.05, .95, 'All samples', va = 'top', transform = gca().transAxes, weight = 'bold')
		xlabel('Weighted Δ$_{47}$ deviation\nof each laboratory')
		ylabel('Cumulative\ndistribution', labelpad = -8)
		yticks([0,1])

		savefig('output/ETH1234_vs_HEG/Fig_7_KS_tests_ETH1234_vs_HEG.pdf')
		close(fig)


		### NEW VS OLD VALUES
		figure(figsize = (3.5, 3.5))
		subplots_adjust(.15,.17,.95,.97)


		plot([.9196], [.9196], 'wo', ms = 5, mec = 'g', mew = 1, label = '25°C equil. gases')

		X = [ .258,  .256,  .691,  .507]
		Y = [final_values[s]['D47'] for s in ethsamples]
		a,b = polyfit(X, Y, 1)
		x = array([0,min(X)])
		plot(x, a*x+b, '-', lw = .75, color = [.5]*3, dashes = (2,3), zorder = -100)
		x = array([min(X),max(X)])
		plot(x, a*x+b, '-', lw = .75, color = [.5]*3, dashes = (6,2), zorder = -100)
		x = array([max(X),1])
		plot(x, a*x+b, '-', lw = .75, color = [.5]*3, dashes = (2,3), zorder = -100)
		x = array([0,1])

		xx = 0.5
		text(xx, a*xx+b-0.1, f'slope = {a:.3f}', va = 'center', ha = 'center', size = 10, color = [.5]*3, rotation = 43)

		plot([.0266], [.0266], 'wo', ms = 5, mec = 'r', mew = 1, label = '1000°C heated gases')
		plot([.0266], [.0266+0.05], 'wo', ms = 5, mec = 'r', mew = 1)

		plot(X, Y, 'ws', mec = 'k', ms = 5, mew = 1, label = 'ETH standards')

		style = 'Simple, tail_width=0, head_width=2, head_length=2.5'
		kw = dict(arrowstyle=style, color='r', linewidth = .75)
		arr = patches.FancyArrowPatch((.0266+0.015, .0266), (.0266+0.015, .0266+0.05), connectionstyle="arc3,rad=1.5", **kw)
		gca().add_patch(arr)
		text(.0266+0.07, 0.0266+0.025, '+0.05 ‰ (partial re-equilibration of HG?)', va = 'center', ha = 'left', color = 'r', size = 9)

		axis([*x, *x])
		xlabel('Previously determined values\n(Bernasconi et al., 2018)')
		ylabel('New values (this study)')
		legend(prop = {'size': 9})
		savefig('output/ETH1234_vs_HEG/Fig_4_ETH1234_new_vs_old.pdf')

		compute_old_to_new_conversion(final_values)

		return final_values

	else:
		return {
			'ETH-1': {'D47': 0.2920368766271781, 'sD47': 0.0018830130448527374},
			'ETH-2': {'D47': 0.2948560677529254, 'sD47': 0.0017907001951658636},
			'ETH-3': {'D47': 0.7007967024253421, 'sD47': 0.0016113573967636997},
			'ETH-4': {'D47': 0.5369344061975312, 'sD47': 0.0019727851254116363},
			}


def run_InterCarb():	
	'''
	Process the data for ETH-4, IAEA-C1/2, and MERCK vs ETH-1/2/3
	'''
	if RUN_INTERCARB:

		create_tree('output/InterCarb/')

		print('\nNOMINAL Δ47 VALUES USED IN INTERCARB:')
		try:
			D47data.Nominal_D47 = {}
			for sample in ['ETH-1', 'ETH-2', 'ETH-3']:
				D47data.Nominal_D47[sample] = round(final_values[sample]['D47'], 4)
		except NameError:
			D47data.Nominal_D47 = {'ETH-1': 0.2920, 'ETH-2': 0.2949, 'ETH-3': 0.7008}
		print(D47data.Nominal_D47)

		print('Reading raw data...')
		rawdata = D47data()
		rawdata.read('input/InterCarb/rawdata.csv')
		labs = sorted({r['Lab'] for r in rawdata})

		for r in rawdata:
			r['LabSession'] = r['Session']
			r['Session'] = r['Lab'] + r['LabSession']
# 		rawdata.refresh()


		if SAVE_RAWDATA:
			save_rawdata([], 'output/InterCarb/Table_S2_InterCarb_data.csv')
		
		labdata = {}
		InterCarb_results = {}
		for lab in labs:
			InterCarb_results[lab] = {}
			labdata[lab] = {}
			for session in sorted({r['Session'] for r in rawdata if r['Lab'] == lab}):
				InterCarb_results[lab][session] = {}
				rd = D47data([r for r in rawdata if r['Lab'] == lab and r['Session'] == session])
				rd.ALPHA_18O_ACID_REACTION = alpha18_acid_reaction_calcite(rd[0]['Acid T'])
				for session in rd.sessions:
					rd.sessions[session]['d13C_standardization_method'] = '1pt'
					rd.sessions[session]['d18O_standardization_method'] = '1pt'
				rd.wg()
				rd.crunch()
				rd.standardize(method = 'indep_sessions')
				labdata[lab][session] = rd
				if SAVE_RAWDATA:
					save_rawdata(rd, 'output/InterCarb/Table_S2_InterCarb_data.csv', 'a')
				for u in UNKNOWNS:
					if u in rd.unknowns:
						InterCarb_results[lab][session][u] = {k: rd.unknowns[u][k] for k in ['D47', 'SE_D47']}
						InterCarb_results[lab][session][u]['autogenic_SE_D47'] = rd.repeatability['r_D47'] / sqrt(rd.unknowns[u]['N'])

		for lab in labs:
			for u in UNKNOWNS:
				X = [InterCarb_results[lab][s][u]['D47'] for s in InterCarb_results[lab] if u in InterCarb_results[lab][s]]
				sX = [InterCarb_results[lab][s][u]['SE_D47'] for s in InterCarb_results[lab] if u in InterCarb_results[lab][s]]
				if X:
					InterCarb_results[lab][u] = {}
					InterCarb_results[lab][u]['D47'], InterCarb_results[lab][u]['SE_D47'] = w_avg(X, sX)
				sX = [InterCarb_results[lab][s][u]['autogenic_SE_D47'] for s in InterCarb_results[lab] if u in InterCarb_results[lab][s]]
				if X:
					InterCarb_results[lab][u]['autogenic_D47'], InterCarb_results[lab][u]['autogenic_SE_D47'] = w_avg(X, sX)

		for u in UNKNOWNS:
			InterCarb_results[u] = {}

			X = [InterCarb_results[lab][u]['D47'] for lab in InterCarb_results if u in InterCarb_results[lab]]
			sX = [InterCarb_results[lab][u]['SE_D47'] for lab in InterCarb_results if u in InterCarb_results[lab]]
			InterCarb_results[u]['D47'], InterCarb_results[u]['SE_D47'] = w_avg(X, sX)
		
			X = [InterCarb_results[lab][u]['autogenic_D47'] for lab in InterCarb_results if u in InterCarb_results[lab]]
			sX = [InterCarb_results[lab][u]['autogenic_SE_D47'] for lab in InterCarb_results if u in InterCarb_results[lab]]
			InterCarb_results[u]['autogenic_D47'], InterCarb_results[u]['autogenic_SE_D47'] = w_avg(X, sX)

		
		summary_of_sessions(labdata)
		save_InterCarb_results(InterCarb_results)
		acid_T_effects(labdata, InterCarb_results)
		MS_effects(labdata, InterCarb_results)
		interlab_plot(InterCarb_results)
		KS_tests(InterCarb_results)
		intra_lab_session_plots(labdata, InterCarb_results)
		single_session_plots(labdata)


def compute_lab_rsmwd( InterCarb_results, lab ):
	chisq, Nf = 0, 0
	for sample in UNKNOWNS :
		D47, sD47 = [], []
		sessions = sorted([s for s in InterCarb_results[lab] if s not in UNKNOWNS and sample in InterCarb_results[lab][s]])
		if len(sessions) > 1:
			Nf -= 1
			for s in sessions :
				chisq += ((InterCarb_results[lab][s][sample]['D47'] - InterCarb_results[lab][sample]['D47']) / InterCarb_results[lab][s][sample]['SE_D47'])**2
				Nf += 1
	if Nf :
		RMSWD = sqrt(chisq / Nf)#		 print(f'RMSWD for {lab} sessions is {RMSWD:.6f}')
		return RMSWD
	else :
		return 0

def intra_lab_session_plots(labdata, InterCarb_results, path = 'output/InterCarb/Session plots/'):
	create_tree(path)
	for lab in labdata:
		rmswd = compute_lab_rsmwd( InterCarb_results, lab )
		sessions = sorted([s for s in labdata[lab]])

		Ns = len(sessions)
		if Ns < 2 : continue
		e,f = .55, 1.4
		g = 2*e+Ns*f
		if Ns > 2:
			fig = figure( figsize = (g,g) )
			subplots_adjust(e/g,e/g,1-e/g,1-e/g,0.08,.08)
		else:
			fig = figure( figsize = (g,g+1) )
			subplots_adjust(e/g,e/g*2.3,1-e/g,1-e/g*0.8,0.08,.08)

		for i in range(Ns) :
			for j in range(Ns-i-1) :
				j += i+1
				session1, session2 = sessions[j], sessions[i]
				samples = sorted([s for s in labdata[lab][session1].unknowns if s in labdata[lab][session2].unknowns])
				subplot(Ns-1,Ns-1,(Ns-1)*i+j)
#				 gca().tick_params(axis='both', direction='in')
				for sample in samples:
					allX1 = [r['D47'] for r in labdata[lab][session1].samples[sample]['data']]
					X1 = labdata[lab][session1].samples[sample]['D47']
					SE_X1 = labdata[lab][session1].samples[sample]['SE_D47']
					autogenic_SE_X1 = labdata[lab][session1].repeatability['r_D47'] / len(allX1)**.5

					allX2 = [r['D47'] for r in labdata[lab][session2].samples[sample]['data']]
					X2 = labdata[lab][session2].samples[sample]['D47']
					SE_X2 = labdata[lab][session2].samples[sample]['SE_D47']
					autogenic_SE_X2 = labdata[lab][session2].repeatability['r_D47'] / len(allX2)**.5

					plot( allX1, [X2]*len(allX1), 'k+', ms = 2, mew = .5 )
					plot( [X1]*len(allX2), allX2, 'k+', ms = 2, mew = .5 )
					gca().add_artist(
						patches.Ellipse(
							 xy = (X1,X2), width = 2 * F95_2df * autogenic_SE_X1, height = 2 * F95_2df * autogenic_SE_X2,
							 lw = .5, fc = 'none', ec = 'k', ls = '-', alpha = 1/3 )
							 )
					gca().add_artist(
						 patches.Ellipse(
							 xy = (X1,X2), width = 2 * F95_2df * SE_X1, height = 2 * F95_2df * SE_X2,
							 lw = 1, fc = 'none', ec = 'k', ls = '-' )
							 )
				xmin,xmax = .325, .775
				plot([xmin,xmax], [xmin,xmax], 'r-', lw = .75, zorder = -100)
				axis([xmin,xmax,xmin,xmax])
				xticks([.4,.5,.6,.7])
				yticks([.4,.5,.6,.7])
				grid(alpha = .333)
				if j == i+1 :
					xlabel(session1[5:])
					ylabel(session2[5:])
				else :
					gca().set_xticklabels( ['']*4 )
					gca().set_yticklabels( ['']*4 )
					gca().tick_params(axis='both', length=0)

		caption = f'''
Comparison of Δ$_{{47}}$ values
for unknown samples obtained
in {Ns} sessions by Lab {lab[-2:]}.
Black markers show individual
analyses from each session.
Thick black ellipses correspond
to 95 % confidence limits on
session-averaged values with fully
propagated uncertainties. Thin
gray ellipses correspond to 95 %
confidence limits not accounting
for standardization errors.
'''[1:-1]
		if Ns > 2:
			text(
				.5,
				1-e/g*0.8,
				f'$\\mathbf{{{lab}}}$   (RMSWD = {rmswd:.3f})',
				ha = 'center',
				va = 'bottom',
				size = 10,
	#			 weight = 'bold',
				transform = fig.transFigure
				)
			text(
				.05, .05, caption,
				va = 'bottom',ha = 'left',
				transform = fig.transFigure,
				size = 8,
				)
		else:
			text(
				.5,
				.97,
				f'$\\mathbf{{{lab}}}$\n(RMSWD = {rmswd:.3f})',
				ha = 'center',
				va = 'top',
				size = 10,
	#			 weight = 'bold',
				transform = fig.transFigure
				)
			text(
				.5, 0.05, '\n'.join([f'{x} {y}' for x,y in zip(caption.split('\n')[::2], caption.split('\n')[1::2]+['']) ]),
				va = 'bottom',ha = 'center',
				transform = fig.transFigure,
				size = 8,
				)

		savefig(f'{path}{lab}.pdf')
		close(fig)


def single_session_plots(labdata, path = 'output/InterCarb/Session plots/'):
	create_tree(path)
	for lab in labdata:
		for session in labdata[lab]:

			d47max = max([r['d47'] for r in labdata[lab][session]])

			sp = labdata[lab][session].plot_single_session(
				session,
				error_contour_interval = 3e-3,
				xylimits = (d47max - 65, d47max + 5, 0.1, 0.9),
				)

			sp.fig.set_size_inches(7,7)
			subplots_adjust(left=.1, bottom=0.25, right=.95, top=.9)
			
			text(
				.17, .03, f'''
All analyses (cross markers) in Session {session[-2:]} of Lab {lab[-2:]}.
Anchors are plotted in red and unknown in blue.
Average Δ$_{{47}}$ values for each sample are plotted as thick blue and red lines.
Standardization errors are mapped as gray contours.
				'''[1:-1],
				va = 'bottom',ha = 'left',
				transform = sp.fig.transFigure,
				)

			for x in sp.anchor_avg + sp.unknown_avg:
				x.set_alpha(0.25)
				x.set_linewidth(3)

			title(f'{lab} - Session {session[-2:]}\nΔ$_{{47}}$ repeatability = {labdata[lab][session].repeatability["r_D47"]*1000:.1f} ppm')

			savefig(f'{path}/{session}.pdf')
			close(sp.fig)


def interlab_plot(InterCarb_results, path = 'output/InterCarb/Fig_5_InterCarb_results.pdf'):
	'''
	Generate Figure 5
	'''
	
	create_tree(path)
	
	fig = figure(figsize = (8,8))
	subplots_adjust( left = .11, right = .98, bottom = .04, top = .95, hspace = .11, wspace = .29)
	global_yrange = -1.
	axes = [None]*4
	for ks, sample in enumerate(UNKNOWNS):
		axes[ks] = subplot(2, 2, ks+1)
		chisq_t, Nf = 0, -1

		sorted_labs = sorted(
			[l for l in InterCarb_results if l not in UNKNOWNS and sample in InterCarb_results[l]],
			key = lambda l: InterCarb_results[l][sample]['SE_D47']
			)
		x, y_inf, y_sup, x_labels = 0, 1e3, -1e3, {}
		for lab in sorted_labs :
			x += 1/(1 + len(sorted_labs))
			chisq_t += ( (InterCarb_results[lab][sample]['D47'] - InterCarb_results[sample]['D47']) / InterCarb_results[lab][sample]['SE_D47'] )**2
			Nf += 1
			errorbar( x, InterCarb_results[lab][sample]['D47'], F95_1df * InterCarb_results[lab][sample]['SE_D47'],
				marker = 'None',
				ls = 'None',
				ecolor = 'k',
				elinewidth = 1,
				capsize = 4,
				capthick = 1,
				)
			rect_width = .015
			gca().add_patch(Rectangle(
				(x-rect_width/2,InterCarb_results[lab][sample]['autogenic_D47']-F95_1df*InterCarb_results[lab][sample]['autogenic_SE_D47']), rect_width, 3.92*InterCarb_results[lab][sample]['autogenic_SE_D47'],
				linewidth = 1,
				edgecolor = 'k',
				facecolor = 'w',
				zorder = 100,
				))
			y_inf = min(y_inf, InterCarb_results[lab][sample]['D47'] - F95_1df * InterCarb_results[lab][sample]['SE_D47'])
			y_sup = max(y_sup, InterCarb_results[lab][sample]['D47'] + F95_1df * InterCarb_results[lab][sample]['SE_D47'])
			x_labels[lab] = x
			text(
				x_labels[lab],
				0.005 + InterCarb_results[lab][sample]['D47'] + F95_1df * InterCarb_results[lab][sample]['SE_D47'],
				str(int(lab[-2:])),
				ha = 'center',
				va = 'center',
				color = (0,0,.5),
				size = 7,
				)

		axhline( InterCarb_results[sample]['D47'], c = 'r', zorder = -100, lw = 1.5 )

		x1,x2,y1,y2 = axis([0, 1, y_inf - 0.05*(y_sup-y_inf), y_sup + 0.1*(y_sup-y_inf)])
		global_yrange = max(global_yrange, y2-y1)
		xticks([])
		kw = dict(size = 10, transform = gca().transAxes, va = 'bottom', ha = 'left')
		text(0.02, 0.01,
			f"$\\mathbf{{{sample}}}$: {InterCarb_results[sample]['D47']:.4f} ± {InterCarb_results[sample]['SE_D47']:.4f} ‰ (1SE)",
			**kw)
# 		if Nf:
# 			p = 1-chi2.cdf(chisq_t,Nf)
# 			kw['ha'] = 'left'
# 			kw['va'] = 'bottom'
# 			text(0.02, 0.02,
# 				f'RMSWD = {(chisq_t/Nf)**.5:.2f}',
# 				**kw)
		ylabel('Δ$_{47}$ (‰)')

	for ks, sample in enumerate(UNKNOWNS):
		sca(axes[ks])
		x1,x2,y1,y2 = axis()
		axis([x1, x2, (y1+y2-global_yrange)/2, (y1+y2+global_yrange)/2])
		
	savefig(path)
	close(fig)


def KS_tests(InterCarb_results, path = 'output/InterCarb/Fig_6_InterCarb_KS_tests.pdf'):
	'''
	Generate Figure 6
	'''

	fig = figure(figsize = (8,3.7))
	subplots_adjust(.05, .16, .99, .93, .05, .05)
	myaxes = [subplot(2,5,1+k) for k in range(10)]
	myaxes = myaxes[-5:] + myaxes[:5]
	sdev = []
	sdevu = []
	for k, u in enumerate(UNKNOWNS):
		D47, sD47 = InterCarb_results[u]['D47'], InterCarb_results[u]['SE_D47']
		dev = sorted([
			(InterCarb_results[l][u]['D47'] - D47)/InterCarb_results[l][u]['SE_D47']
			for l in InterCarb_results if l not in UNKNOWNS and u in InterCarb_results[l]
			])
		chi2 = sum([x**2 for x in dev])
		N = len(dev)

		D47u, sD47u = InterCarb_results[u]['autogenic_D47'], InterCarb_results[u]['autogenic_SE_D47']
		devu = sorted([
			(InterCarb_results[l][u]['autogenic_D47'] - D47)/InterCarb_results[l][u]['autogenic_SE_D47']
			for l in InterCarb_results if l not in UNKNOWNS and u in InterCarb_results[l]
			])
		chi2u = sum([x**2 for x in devu])

		sdev += dev
		sdevu += devu

		col1, col2 = [0]*3, [0]*3
		sca(myaxes[k+1])
		X = array(dev)
		Y = arange(N)/(N-1)
		plot(X, Y, '-', lw = .75, color = col1)
		plot(X, Y, 'wo', mec = col1, mew = .75, ms = 3.5)
		x1,x2,y1,y2 = axis()

		sca(myaxes[k+6])
		Xu = array(devu)
		plot(Xu, Y, '-', lw = .75, color = col2)
		plot(Xu, Y, 'wo', mec = col2, mew = .75, ms = 3.5)
		x3,x4,y3,y4 = axis()
		
		x1 = -7.5
		x2 = 7.5
		y1 = -0.15
		y2 = 1.15
		
		sca(myaxes[k+1])
		x = linspace(x1, x2)
		y = norm().cdf(x)
		plot(x,y,'-', lw = 2, zorder = -10, color = [.6,.8,1])
		pvalue = kstest(X, 'norm', (0, 1), mode = 'asymp').pvalue
		text(.03, .97, f'fully propagated\nerrors', va = 'top', ha = 'left', size = 8, transform = gca().transAxes)
		text(.95, .05, f'p = {100*pvalue:{".0f" if pvalue>0.1 else ".1f"}} %', va = 'bottom', ha = 'right', size = 12, transform = gca().transAxes)
		axis([x1,x2,y1,y2])
		yticks([])
		xlabel('Sigma-deviation\nof each laboratory')

		sca(myaxes[k+6])
		plot(x,y,'-', lw = 2, zorder = -10, color = [.6,.8,1])
		pvalue = kstest(Xu, 'norm', (0, 1), mode = 'asymp').pvalue
		text(.03, .97, f'not accounting for\nstandardization\nerrors', va = 'top', ha = 'left', size = 8, transform = gca().transAxes)
		text(.95, .05, f'p = {100*pvalue:{".0f" if pvalue>0.1 else ".1f"}} %', va = 'bottom', ha = 'right', size = 12, transform = gca().transAxes)
		axis([x1,x2,y1,y2])
		yticks([])
		xticks([])
		title(u,weight = 'bold', size = 11)
		
	sdev = sorted(sdev)
	sdevu = sorted(sdevu)
	N = len(sdev)

	sca(myaxes[0])
	X = array(sdev)
	Y = arange(N)/(N-1)
	plot(X, Y, '-', lw = .75, color = col1)
	plot(X, Y, 'wo', mec = col1, mew = .5, ms = 2)
	x = linspace(x1, x2)
	y = norm().cdf(x)
	plot(x,y,'-', lw = 2, zorder = -10, color = [.6,.8,1])
	pvalue = kstest(X, 'norm', (0, 1), mode = 'asymp').pvalue
	text(.03, .97, f'fully propagated\nerrors', va = 'top', ha = 'left', size = 8, transform = gca().transAxes)
	text(.95, .05, f'p = {100*pvalue:{".0f" if pvalue>0.1 else ".1f"}} %', va = 'bottom', ha = 'right', size = 12, transform = gca().transAxes)
	axis([x1,x2,y1,y2])
	xlabel('Sigma-deviation\nof each laboratory')
	ylabel('Cumulative\ndistribution', labelpad = -8)
	yticks([0,1])

	sca(myaxes[5])
	Xu = array(sdevu)
	plot(Xu, Y, '-', lw = .75, color = col2)
	plot(Xu, Y, 'wo', mec = col2, mew = .5, ms = 2)
	plot(x,y,'-', lw = 2, zorder = -10, color = [.6,.8,1])
	pvalue = kstest(Xu, 'norm', (0, 1), mode = 'asymp').pvalue
	text(.03, .97, f'not accounting for\nstandardization\nerrors', va = 'top', ha = 'left', size = 8, transform = gca().transAxes)
	text(.95, .05, f'p = {100*pvalue:{".0f" if pvalue>0.1 else ".1f"}} %', va = 'bottom', ha = 'right', size = 12, transform = gca().transAxes)
	axis([x1,x2,y1,y2])
	ylabel('Cumulative\ndistribution', labelpad = -8)
	xticks([])
	yticks([0,1])
	title('All samples',weight = 'bold', size = 11)
		
	savefig(path)
	close(fig)
	
	
if __name__ == '__main__':

	final_values = ETH1234_vs_HEG()
	run_InterCarb()		
