######### #!/usr/bin/env python
import cPickle
import numpy
import sys

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def main(argv):
		cplot(argv[1])

def cplot(fname):
	data = cPickle.load(open(fname,"r"))

	newdata = {}
	charge_energy, scint_count = [], []
	for i in range(len(data['charge_energy'])):
		if data['charge_energy'][i] <= 600.: # 600 to match paper
			charge_energy.append(data['charge_energy'][i])
			scint_count.append(data['scint_count'][i])

	newdata['scint_count'] = numpy.array(scint_count)
	newdata['charge_energy'] = numpy.array(charge_energy)

	print "Scintillation vs Charge plot"
	scint_v_charge = plt.figure(figsize=(10,15))
	ax = scint_v_charge.add_subplot(111)
	ax.scatter(newdata['charge_energy'],newdata['scint_count'], s = 1)
	
        linex, liney = [], []
        i = 0.0
        while i < 600.0:
            linex.append(i)
            liney.append(i*33.864)
            i += 1
        ax.plot(linex,liney,'g-')
	
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(16)

	ax.set_xlabel('Charge Energy (kev)')
	ax.set_ylabel('Scintillation Counts')
	plt.savefig('charge_energy.png')
	plt.close(scint_v_charge)
	print("Finished Scintillation vs Charge plot")

if __name__ == '__main__':
    main(sys.argv)
