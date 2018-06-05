######### #!/usr/bin/env python
import cPickle
import numpy
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import tree #DecisionTreeRegressor
from sklearn import cross_validation

norm = numpy.array([1.,1.,1.,200./50000.])
algorithm = "DTR"

def main(argv):
	if 'split' in argv[1]: #./MLtraining.py split cut opt splitsize
		split_set(argv[2], argv[3], float(argv[4]))
		raw_input("Press enter to quit")
	if 'predict' in argv[1]: 
		# ./MLtraining.py predict cut opt
		train_DTR(argv[2], argv[3])
		raw_input("Press enter to quit")
	if 'recon' in argv[1]: 
		# ./MLtraining.py recon recon_file cut recon_opt ttip
		recon_data(argv[2], argv[3], argv[4], argv[5])
		raw_input("Press enter to quit")

	if 'plot' in argv[1]: 
		# ./MLtraining.py plot cut opt ttip
		# These are E vs. z plots and 3d position plots.
		plot_pickle(argv[2], argv[3], argv[4])
		raw_input("Press enter to quit")
	if 'posdiff' in argv[1]:
		# ./MLtraining.py posdiff initial_input.pkl predictions_input.pkl cut
		plot_pos_diffs(argv[2], argv[3], argv[4])
		raw_input("Press enter to quit")

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
	cPickle.dump(verify,   open("%s/%s/%s_%2.1f.pkl"%(algorithm,cut,cut,splitsize),   "w"))	
	

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

if __name__ == '__main__':
	main(sys.argv)
