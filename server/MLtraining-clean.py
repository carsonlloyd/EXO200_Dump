#!/usr/bin/env python
import cPickle
import numpy
import sys

## At SLAC
import ROOT
ROOT.gSystem.Load("libEXOUtilities")

## At home
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#from sklearn import tree #DecisionTreeRegressor
#from sklearn import cross_validation

norm = numpy.array([1.,1.,1.,200./50000.])
algorithm = "DTR"

def main(argv):
	# This will only work at SLAC, or a system at which ROOT and EXO Analysis are installed.
	if 'create' in argv[1]:
		# ./MLtraining.py create cut opt
		create_training_data(argv[2], argv[3]) 
#		raw_input("Press enter to quit")

	if 'plot' in argv[1]: 
		# ./MLtraining.py plot cut opt ttip
		# These are E vs. z plots and 3d position plots.
		plot_pickle(argv[2], argv[3], argv[4])
		raw_input("Press enter to quit")

	if 'posdiff' in argv[1]:
		# ./MLtraining.py posdiff initial_input.pkl predictions_input.pkl cut
		plot_pos_diffs(argv[2], argv[3], argv[4])
		raw_input("Press enter to quit")

	if 'predict' in argv[1]: 
		# ./MLtraining.py predict cut opt
		train_DTR(argv[2], argv[3])
		raw_input("Press enter to quit")

	if 'recon' in argv[1]: 
		# ./MLtraining.py recon recon_file cut recon_opt ttip
		recon_data(argv[2], argv[3], argv[4], argv[5])
		raw_input("Press enter to quit")		

	if 'scint' in argv[1]: # ./MLtraining.py scint recon_file opt cut ttip
		scint_fits(argv[2], argv[3], argv[4], argv[5])
		raw_input("Press enter to quit")

	if 'split' in argv[1]: #./MLtraining.py split cut opt splitsize
		split_set(argv[2], argv[3], float(argv[4]))
		raw_input("Press enter to quit")

def create_training_data(cut, opt):
	afile = ROOT.TFile.Open("alpha_runs.root","r")
	t = afile.Get("tree")
	num = t.GetEntries()
	
	count = 0 
	apd_data, position = [], []
	for i in range(num):
		if i %int(num/10) == 0:
			print "%d of %d"%(i,num)
		t.GetEntry(i)
		ed = t.EventBranch
		nsc = ed.GetNumScintillationClusters()
		if not nsc == 1: continue ## Rn-Po
		sc = ed.GetScintillationCluster(0)
		if not sc.GetNumChargeClusters() == 1: continue ## Rn-Po
#		if sc.fRawEnergy != 0: continue ## this will be 0 if the event is not fully reconstructed.
		cc = sc.GetChargeClusterAt(0) # remember to loop when not Rn-Po

		## fiducial volume cut
#		if abs(cc.fZ) < 5.: continue ## get rid of events on cathode
#		if abs(cc.fZ) > 185.: continue ## get rid of events on anode
#		if numpy.sqrt(cc.fX**2+cc.fY**2) > 172.: continue ## get rid of events on/too close to surface

		nscounts = 0
		if sc.GetNumChargeClusters() == 1:
			nscounts = sc.GetCountsSumOnAPDPlane(0) + sc.GetCountsSumOnAPDPlane(1)
			if float(nscounts)/float(cc.fPurityCorrectedEnergy) < 33.864: continue

#			if float(sc.fRawEnergy)/float(cc.fPurityCorrectedEnergy) < 33.864: continue ## cutoff for betas, gammas

		## could add more cuts here
#		position.append((cc.fX, cc.fY, cc.fZ, sc.fRawEnergy))
#		position.append((cc.fX, cc.fY, cc.fZ, nscounts))
		apd_data.append(numpy.zeros(226-150,dtype='float32'))
		count += 1
		apd_data[-1][0] = sc.GetCountsOnAPDPlane(0)
		apd_data[-1][1] = sc.GetCountsOnAPDPlane(1)
		napd = sc.GetNumAPDSignals()
		for j in range(napd):
			sig = sc.GetAPDSignalAt(j)
			if sig.fChannel < 3: continue
			apd_data[-1][sig.fChannel-150] = sig.fRawCounts
		
	del afile
	del t
	print "Data check finished - %i events survived out of %i" %(count, num)
	data = {}
	data['apd_data'] = numpy.array(apd_data)
#	data['position'] = numpy.array(position)
	cPickle.dump(data,open('%s_%s.pkl'%(cut,opt),'w'))
	print "Pickle finished"

def plot_pickle(cut, opt, ttip):
	print "Plotting", cut, opt, ttip
	data = cPickle.load(open("%s/%s/%s_%s_%s.pkl"%(algorithm,cut,cut,opt,ttip),"r"))

	print "E vs z plot"
	E_vs_z = plt.figure(figsize=(6,4))
	#H, xedges, yedges = numpy.histogram2d(data['position'][:,2]/norm[2],data['position'][:,3]/norm[3],bins=(200,200),range=((-200,200),(25000,70000)))
	H, xedges, yedges = numpy.histogram2d(data['position'][:,2]/norm[2],data['position'][:,3]/norm[3],bins=(200,200),range=((-200,200),(0,100000)))
	H = numpy.rot90(H)
	H = numpy.flipud(H)	
	H2 = numpy.ma.masked_where(H==0,H)
	H3 = numpy.log(H2)	
	plt.pcolormesh(xedges,yedges,H3)
	plt.xlabel('Z pos (mm)')
	plt.ylabel('fRawEnergy (LM corrected)')
	plt.savefig('%s/%s/%s_%s_%s_data_e_vs_z.png'%(algorithm,cut,cut,opt,ttip),dpi=300)
	plt.close(E_vs_z)
	print "Finished E vs z plot"
	
	print "3d position plot"
	pos_3d = plt.figure(figsize=(6,4))
	ax = pos_3d.add_subplot(111, projection='3d')
	ax.scatter(data['position'][:,0]/norm[0],data['position'][:,1]/norm[1],data['position'][:,2]/norm[2], s=2, lw=0)
	ax.set_xlabel('cc.fX')
	ax.set_ylabel('cc.fY')
	ax.set_zlabel('cc.fZ')
	ax.set_xlim([-200,200])
	ax.set_ylim([-200,200])
	ax.set_zlim([-200,200])
	plt.savefig('%s/%s/%s_%s_%s_data_3d_pos.png'%(algorithm,cut,cut,opt,ttip),dpi=300)
	plt.close(pos_3d)	
 	print "Finished 3d position plot"
	print "End plotting", cut, opt, ttip

	print "histograms"
	nhists = 40
#	z_edges = numpy.linspace(-200.0, 200.0, nhists+1)
	z_edges = numpy.array([-200, -190, -10, 0, 10, 190, 200])

	nbins = 20
	for (z1, z2) in zip(z_edges[:-1], z_edges[1:]):
		mask = (z1 <= data['position'][:,2]/norm[2]) & (data['position'][:,2]/norm[2] < z2)
		# make a 2D histogram with event_xs[mask] and event_ys[mask]

		event_x = data['position'][:,0]/norm[0]
		event_y = data['position'][:,1]/norm[1]

		hists, xedges, yedges = numpy.histogram2d(event_x[mask], event_y[mask], bins=nbins)
		hists = numpy.rot90(hists)
		hists = numpy.flipud(hists)
		hists_masked = numpy.ma.masked_where(hists==0,hists)

		fig = plt.figure(figsize=(6,4))
		plt.pcolormesh(xedges,yedges,hists_masked)
		plt.xlabel('cc.fX')
		plt.ylabel('cc.fY')
		cbar = plt.colorbar()
		cbar.ax.set_ylabel('Counts')
		plt.savefig('%s/%s/%s_%s_%s_%3.1f_to_%3.1f.png' %(algorithm,cut,cut,opt,ttip,z1,z2),dpi=300)
		plt.close(fig)

def train_DTR(cut,opt):
	## unload data
	print "Loading pickle"
	data = cPickle.load(open("%s/%s/%s_%s.pkl"%(algorithm,cut,cut,opt),"r"))
	apd_data = data['apd_data']
	position = data['position']*norm
		
	training_initial = {}
	testing_initial = {}	
	training_predictions = {}
	testing_predictions = {}

	## split training and test sets
	apd_train, apd_test, pos_train, pos_test = cross_validation.train_test_split(apd_data, position, test_size = 0.3)
	dtr = tree.DecisionTreeRegressor()

	## Pickle these later just in case we need them
	training_initial['apd_data'] = apd_train
	training_initial['position'] = pos_train
	testing_initial['apd_data'] = apd_test
	testing_initial['position'] = pos_test

	## fit on the training set
	print "Fitting..."
	testing_predictions['dtrfit'] = dtr.fit(apd_train, pos_train)

	## predict training set - should have a score of very close to 1, if not 1
	print "Predicting training set."
	training_predictions['position'] = dtr.predict(apd_train)	

	## predict on test set
	print "Predicting test set."
	testing_predictions['position'] = dtr.predict(apd_test)	

	## get scores
	training_predictions['score'] = dtr.score(apd_train, pos_train)
	testing_predictions['score'] = dtr.score(apd_test, pos_test)

	## Making pickle of differences
	print "Pickling %s"%cut
	cPickle.dump(training_initial,open('%s/%s/%s_%s_training_initial.pkl'%(algorithm,cut,cut,opt),'w'))
	cPickle.dump(training_predictions,open('%s/%s/%s_%s_training_predictions.pkl'%(algorithm,cut,cut,opt),'w'))
	cPickle.dump(testing_initial,open('%s/%s/%s_%s_test_initial.pkl'%(algorithm,cut,cut,opt),'w'))
	cPickle.dump(testing_predictions,open('%s/%s/%s_%s_test_predictions.pkl'%(algorithm,cut,cut,opt),'w'))

	print "Sending output pickle to plotting differences function."
	plot_pos_diffs("training_initial.pkl", "training_predictions.pkl", cut, opt)
	plot_pos_diffs("test_initial.pkl", "test_predictions.pkl", cut, opt)
	print "Done with difference plots."

	print "Sending output pickle to plot E vs z and 3d plot function."
	plot_pickle(cut, opt, "training_initial")
	plot_pickle(cut, opt, "training_predictions")
	plot_pickle(cut, opt, "test_initial")
	plot_pickle(cut, opt, "test_predictions")
	print "Done with E vs z and 3d plots."

def plot_pos_diffs(pickle_name_initial, pickle_name_predictions, cut, opt):
	print "Loading pickle %s"%pickle_name_initial
	initial = cPickle.load(open("%s/%s/%s_%s_%s"%(algorithm,cut,cut,opt,pickle_name_initial),"r"))
	print "Loading pickle %s"%pickle_name_predictions
	predictions = cPickle.load(open("%s/%s/%s_%s_%s"%(algorithm,cut,cut,opt,pickle_name_predictions),"r"))

	diff = predictions['position'] - initial['position']

	if 'test' in pickle_name_initial: filename = "test"
	if 'training' in pickle_name_initial: filename = "training"

	## plot position differences
	print "Plotting differences."
	nbins = 200
	hist_x_train, bins_x_train = numpy.histogram(diff[:,0]/norm[0],nbins,(-200,200))
	hist_y_train, bins_y_train = numpy.histogram(diff[:,1]/norm[1],nbins,(-200,200))
	hist_z_train, bins_z_train = numpy.histogram(diff[:,2]/norm[2],nbins,(-200,200))
	hist_e_train, bins_e_train = numpy.histogram(diff[:,3]/norm[3],nbins,(-10000,10000))
	bincenters_pos = 0.5*(bins_x_train[1:]+bins_x_train[:-1])
	bincenters_ene = 0.5*(bins_e_train[1:]+bins_e_train[:-1])

	fig = plt.figure(figsize=(15,6))
	plt.subplot(1,4,1)
	plt.step(bincenters_pos,hist_x_train)
	locs, labels = plt.xticks()
	plt.setp(labels, rotation=90)
	plt.title("X difference (mm)")
	plt.subplot(1,4,2)
	plt.step(bincenters_pos,hist_y_train)
	locs, labels = plt.xticks()
	plt.setp(labels, rotation=90)
	plt.title("Y difference (mm)")
	plt.subplot(1,4,3)
	plt.step(bincenters_pos,hist_z_train)
	locs, labels = plt.xticks()
	plt.setp(labels, rotation=90)
	plt.title("Z difference (mm)")
	plt.subplot(1,4,4)
	plt.step(bincenters_ene,hist_e_train)
	locs, labels = plt.xticks()
	plt.setp(labels, rotation=90)
	plt.title("Energy difference (APD Counts)")
	plt.savefig("%s/%s/%s_%s_%s_difference.png"%(algorithm,cut,filename,cut,opt))

	print "Score is", predictions['score']	

def split_set(cut, opt, splitsize):
	print "Loading data pickle"
	data = cPickle.load(open("%s/%s/%s_%s.pkl"%(algorithm,cut,cut,opt),"r"))
	apd_data = data['apd_data']
	position = data['position']

	apd_traintest, apd_verify, pos_traintest, pos_verify = cross_validation.train_test_split(apd_data, position, test_size = splitsize)
	traintest = {}
	verify = {}

	traintest['apd_data'] = apd_traintest
	traintest['position'] = pos_traintest
	verify['apd_data'] = apd_verify
	verify['position'] = pos_verify
	traintest_size = 1.-splitsize

	print "Pickling data"
	cPickle.dump(traintest,open("%s/%s/%s_%2.1f.pkl"%(algorithm,cut,cut,traintest_size), "w"))
	cPickle.dump(verify,open("%s/%s/%s_%2.1f.pkl"%(algorithm,cut,cut,splitsize), "w"))	
	

def recon_data(recon_file, cut, recon_opt, ttip):
	print "Loading data pickle"
	data = cPickle.load(open("%s/%s/%s_%s.pkl"%(algorithm,recon_file,recon_file,recon_opt),"r"))
#	data = cPickle.load(open("%s/%s/%s_0.3.pkl"%(algorithm,recon_file,recon_file),"r"))
	apd_data = data['apd_data']
#	position = data['position']*norm

	print "Loading pickle with fit"
	fit = cPickle.load(open("%s/%s/%s_%s_test_predictions.pkl"%(algorithm,cut,cut,recon_opt),"r"))
	dtr = fit['dtrfit']
	MLposition = dtr.predict(apd_data)
	MLpredict = {}
	MLpredict['position'] = MLposition
	MLpredict['apd_data'] = apd_data

	print "Pickling output."
	cPickle.dump(MLpredict,open("%s/%s/%s_%s_%s_%s.pkl"%(algorithm, recon_file, recon_file, recon_opt, cut, ttip),"w"))
	print "Sending to plots"
	plot_pickle(recon_file,str(recon_opt)+"_"+str(cut),ttip)


def scint_fits(recon_file, opt, cut, ttip):
	print "Loading data pickle" 
	data = cPickle.load(open("%s_%s_%s_%s.pkl"%(recon_file, opt, cut, ttip),"r"))
	apd_data = data['apd_data']
	position = data['position'] # position data is already normalized at this point if this pickle has been made with this code


	## Brian's code
	# Here the fitting to the data happens and the number of alphas extracted
	h1 = ROOT.TH1I("h1","h1",100,30000,60000)
	h2 = ROOT.TH1I("h2","h2",100,30000,60000)
	h1.Sumw2() #Set bin error calc
	h2.Sumw2() #Set bin error calc
	data_py = {'tpc1':[],'tpc2':[]}

	## change this if you use a different fiducial cut!!!
	for i in position:
		if i[2] > 5 and i[2] < 185:
			h1.Fill(i[3]/norm[3])
			data_py['tpc1'].append(i[3]/norm[3])
		if i[2] > -185 and i[2] < -5:
			h2.Fill(i[3]/norm[3])
			data_py['tpc2'].append(i[3]/norm[3])

	py_gaus = lambda p,x: p[0]*numpy.exp(-.5*((x-p[1])/p[2])**2)
	gaus = ROOT.TF1("gaus","[0]*exp(-.5*((x-[1])/[2])^2)")
	# Triple gauss for fittting
	#py_skew_trpl_gaus = lambda p,x: p[0]*numpy.exp(-.5*((x-p[1])/p[2])**2)*erfc(p[3]*(x-p[1])/(2.**.5*p[2])) + p[4]*numpy.exp(-.5*((x-p[5])/(p[5]/p[1]*p[2]))**2)*erfc(p[3]*(x-p[5])/(2.**.5*(p[5]/p[1]*p[2]))) + p[6]*numpy.exp(-.5*((x-p[7])/(p[7]/p[1]*p[2]))**2)*erfc(p[3]*(x-p[7])/(2.**.5*(p[7]/p[1]*p[2])))
	#skew_trpl_gaus = ROOT.TF1("skewtrplgaus","[0]*exp(-.5*((x-[1])/[2])^2)*TMath::Erfc([3]*(x-[1])/(sqrt(2)*[2]))+[4]*exp(-.5*((x-[5])/([5]/[1]*[2]))^2)*TMath::Erfc([3]*(x-[5])/(sqrt(2)*[5]/[1]*[2]))+[6]*exp(-.5*((x-[7])/([7]/[1]*[2]))^2)*TMath::Erfc([3]*(x-[7])/(sqrt(2)*[7]/[1]*[2]))")
	#skew_trpl_gaus.SetParameters(400,39000,1000,1,200,43000,50,55000)
	py_trpl_gaus = lambda p,x: p[0]*numpy.exp(-.5*((x-p[1])/p[2])**2) + p[3]*numpy.exp(-.5*((x-p[4])/p[5])**2) + p[6]*numpy.exp(-.5*((x-p[7])/p[8])**2)
	trpl_gaus = ROOT.TF1("trpl_gaus","[0]*exp(-.5*((x-[1])/[2])^2)+[3]*exp(-.5*((x-[4])/[5])^2)+[6]*exp(-.5*((x-[7])/[8])^2)")
	trpl_gaus.SetParameters(1400,39000,500,1400,43000,500,100,55000,500)

	# Fit of Rn+Po218+Po214
	h1.Fit("trpl_gaus","0L","goff",30000,60000)
	h1_fit_res = h1.GetFunction("trpl_gaus")
	par_tpc1 = h1_fit_res.GetParameter(0), h1_fit_res.GetParameter(1), h1_fit_res.GetParameter(2),\
				h1_fit_res.GetParameter(3), h1_fit_res.GetParameter(4), h1_fit_res.GetParameter(5),\
				h1_fit_res.GetParameter(6), h1_fit_res.GetParameter(7), h1_fit_res.GetParameter(8),\
				h1_fit_res.GetChisquare(), h1_fit_res.GetNDF()
	h2.Fit("trpl_gaus","0L","goff",30000,60000)
	h2_fit_res = h2.GetFunction("trpl_gaus")
	par_tpc2 = h2_fit_res.GetParameter(0), h2_fit_res.GetParameter(1), h2_fit_res.GetParameter(2),\
				h2_fit_res.GetParameter(3), h2_fit_res.GetParameter(4), h2_fit_res.GetParameter(5),\
				h2_fit_res.GetParameter(6), h2_fit_res.GetParameter(7), h2_fit_res.GetParameter(8),\
				h2_fit_res.GetChisquare(), h2_fit_res.GetNDF()

	fig_trpl_gaus = plt.figure(figsize=(6,4))
	plt.hist(data_py['tpc1'],bins=100,range=(30000,60000),histtype='step',label='TPC1')
	plt.hist(data_py['tpc2'],bins=100,range=(30000,60000),histtype='step',label='TPC2')
	x = numpy.linspace(30000,60000,500)
	y1 = py_trpl_gaus(par_tpc1,x)
	y2 = py_trpl_gaus(par_tpc2,x)
	plt.plot(x,y1,label='Fit 1')
	plt.plot(x,y2,label='Fit 2')
	plt.xlim(35000,60000)
	plt.legend()
	plt.savefig('Beta_Neut_zpos_corr_scint_trpl_%s_%s.png'%(recon_file, cut),dpi=200)
	plt.close(fig_trpl_gaus)

	fig_trpl_gaus_ind,axes = plt.subplots(nrows=2,ncols=1,sharex=True,figsize=(6,6))
	axes[0].hist(data_py['tpc1'],bins=100,range=(35000,60000),histtype='step',label='TPC1')
	#plt.hist(data_py['tpc2_zcor_sce'],bins=200,range=(20000,60000),histtype='step',label='TPC2')
	x = numpy.linspace(30000,60000,500)
	y1 = py_trpl_gaus(par_tpc1,x)
	y1Rn = py_gaus((par_tpc1[0],par_tpc1[1],par_tpc1[2]),x)
	y1P8 = py_gaus((par_tpc1[3],par_tpc1[4],par_tpc1[5]),x)
	y1P4 = py_gaus((par_tpc1[6],par_tpc1[7],par_tpc1[8]),x)
	axes[0].plot(x,y1,label='Overall')
	axes[0].plot(x,y1Rn,label='Rn222')
	axes[0].plot(x,y1P4,label='Po214')
	axes[0].plot(x,y1P8,label='Po218')
	axes[0].set_xlim(30000,60000)
	axes[0].legend(prop={'size':12})
	axes[0].set_ylabel('Events')
	#plt.subplot(2,1,2)
	axes[1].hist(data_py['tpc2'],bins=100,range=(35000,60000),histtype='step',label='TPC2')
	y2 = py_trpl_gaus(par_tpc2,x)
	y2Rn = py_gaus((par_tpc2[0],par_tpc2[1],par_tpc2[2]),x)
	y2P8 = py_gaus((par_tpc2[3],par_tpc2[4],par_tpc2[5]),x)
	y2P4 = py_gaus((par_tpc2[6],par_tpc2[7],par_tpc2[8]),x)
	axes[1].plot(x,y2,label='Overall')
	axes[1].plot(x,y2Rn,label='Rn222')
	axes[1].plot(x,y2P4,label='Po214')
	axes[1].plot(x,y2P8,label='Po218')
	axes[1].set_xlim(35000,60000)
	axes[1].legend(prop={'size':12})
	axes[1].set_xlabel('Scintillation Counts')
	axes[1].set_ylabel('Events')
	plt.savefig('Beta_Neut_zpos_corr_scint_trpl_ind_tpc1_%s_%s.png'%(recon_file, cut),dpi=200)
	plt.close(fig_trpl_gaus_ind)

	# Lets print results of fits and integrals
	# TPC1
	print "TPC1 Fit Chi^2/NDF = %.2f/%d"%(par_tpc1[9],par_tpc1[10])
	int1_Rn = gaus.Integral(30000,70000,numpy.array((par_tpc1[0],par_tpc1[1],par_tpc1[2]),dtype=numpy.double))
	print 'TPC1 Rn222 = ', int1_Rn/h1.GetXaxis().GetBinWidth(0)
	int1_P8 = gaus.Integral(30000,70000,numpy.array((par_tpc1[3],par_tpc1[4],par_tpc1[5]),dtype=numpy.double))
	print 'TPC1 Po218 = ', int1_P8/h1.GetXaxis().GetBinWidth(0)
	int1_P4 = gaus.Integral(30000,70000,numpy.array((par_tpc1[6],par_tpc1[7],par_tpc1[8]),dtype=numpy.double))
	print 'TPC1 Po214 = ', int1_P4/h1.GetXaxis().GetBinWidth(0)
	# TPC2
	print "TPC2 Fit Chi^2/NDF = %.2f/%d"%(par_tpc2[9],par_tpc2[10])
	int2_Rn = gaus.Integral(30000,70000,numpy.array((par_tpc2[0],par_tpc2[1],par_tpc2[2]),dtype=numpy.double))
	print 'TPC2 Rn222 = ', int2_Rn/h2.GetXaxis().GetBinWidth(0)
	int2_P8 = gaus.Integral(30000,70000,numpy.array((par_tpc2[3],par_tpc2[4],par_tpc2[5]),dtype=numpy.double))
	print 'TPC2 Po218 = ', int2_P8/h2.GetXaxis().GetBinWidth(0)
	int2_P4 = gaus.Integral(30000,70000,numpy.array((par_tpc2[6],par_tpc2[7],par_tpc2[8]),dtype=numpy.double))
	print 'TPC2 Po214 = ', int2_P4/h2.GetXaxis().GetBinWidth(0)

if __name__ == '__main__':
    main(sys.argv)
