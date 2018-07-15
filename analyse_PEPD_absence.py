### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###


__usage__ = """
	python analyse_PEPD_absence.py
	--in <FULL_PATH_TO_DATA_SET_REPORT>
	--group <NAME_OF_TAXONOMIC_GROUP>
				"""

from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import sys

# ---end of imports --- #

def main( arguments ):
	"""! @brief run analysis """
	
	input_file = arguments[ arguments.index( '--in' )+1 ]
	group = arguments[ arguments.index( '--group' )+1 ]
	
	pepd = []
	no_pepd = []
		
	with open( input_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			percentage = float( parts[-1].split('%')[0] )
			if parts[1] == group:
				if parts[2] == '+':
					pepd.append( percentage )
				elif parts[2] == '-':
					no_pepd.append( percentage )
			line = f.readline()
	
	print "average completion (PEPD): " + str( np.median( pepd ) )
	print "average completion (no PEPD): " + str( np.median( no_pepd ) )

	print "average sample size (PEPD): " + str( len( pepd ) )
	print "average sample size (no PEPD): " + str( len( no_pepd ) )

	print stats.mannwhitneyu( pepd, no_pepd )

	fig, ax = plt.subplots()
	ax.hist( pepd, color="green", bins=20 )
	ax.hist( no_pepd, color="red", bins=20 )
	plt.show()


if __name__ == '__main__':
	
	if '--in' in sys.argv and '--group' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
