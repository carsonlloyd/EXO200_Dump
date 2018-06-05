#!/usr/bin/env python
import cPickle
import numpy
import sys

## At SLAC
import ROOT
ROOT.gSystem.Load("libEXOUtilities")

norm = numpy.array([1.,1.,1.,200./50000.])
algorithm = "DTR"

def main(argv):
		create_training_data(argv[1]) 

def create_training_data(fname):
	afile = ROOT.TFile.Open("%s.root"%(fname),"r")
	t = afile.Get("tree")
	num = t.GetEntries()
	
	count = 0 
	apd_data, position = [], []
	charge_energy, scint_count = [], []
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

		if cc.fPurityCorrectedEnergy == 0.0: continue

		nscounts = 0
		if sc.GetNumChargeClusters() == 1:
			nscounts = sc.GetCountsSumOnAPDPlane(0) + sc.GetCountsSumOnAPDPlane(1)
			#if float(nscounts)/float(cc.fPurityCorrectedEnergy) < 33.864: continue
		
		scint_count.append(nscounts)
		charge_energy.append(cc.fPurityCorrectedEnergy)
		count += 1
		
	del afile
	del t
	print "Data check finished - %i events survived out of %i" %(count, num)
	data = {}
	#data['apd_data'] = numpy.array(apd_data)
	#data['position'] = numpy.array(position)
	data['scint_count'] = numpy.array(scint_count)
	data['charge_energy'] = numpy.array(charge_energy)
	cPickle.dump(data,open('%sCL.pkl'%(fname),'w'))
	print "Pickle finished"

if __name__ == '__main__':
    main(sys.argv)
