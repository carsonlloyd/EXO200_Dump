#!/usr/bin/env python
import cPickle
import numpy
import sys
import ROOT
ROOT.gSystem.Load("libEXOUtilities")

norm = numpy.array([1.,1.,1.,200./50000.])
algorithm = "DTR"

def main(argv):
	if 'create' in argv[1]:
		# ./MLtraining.py create cut opt
		create_training_data(argv[2])
#		raw_input("Press enter to quit")
	else:
		print "Run stripped version with -create- -filename- arguments"

def create_training_data(fname):
	afile = ROOT.TFile.Open('%s.root'%(fname),"r")
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
		if abs(cc.fZ) < 5.: continue ## get rid of events on cathode
		if abs(cc.fZ) > 185.: continue ## get rid of events on anode
		if numpy.sqrt(cc.fX**2+cc.fY**2) > 172.: continue ## get rid of events on/too close to surface

		nscounts = 0
		if sc.GetNumChargeClusters() == 1:
			nscounts = sc.GetCountsSumOnAPDPlane(0) + sc.GetCountsSumOnAPDPlane(1)
			if float(nscounts)/float(cc.fPurityCorrectedEnergy) < 33.864: continue

#			if float(sc.fRawEnergy)/float(cc.fPurityCorrectedEnergy) < 33.864: continue ## cutoff for betas, gammas

		## could add more cuts here
		position.append((cc.fX, cc.fY, cc.fZ, sc.fRawEnergy))
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
	data['position'] = numpy.array(position)
	cPickle.dump(data,open('%s-pos1-VOLcut.pkl'%(fname),'w'))
	print "Pickle finished"

if __name__ == '__main__':
    main(sys.argv)
