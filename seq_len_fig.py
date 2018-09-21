### Boas Pucker ###
### Hanna Schilbert ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python seq_len_fig.py
	--in <FULL_PATH_TO_INPUT_DIR>
	--fig <FULL_PATH_TO_FIGURE_FILE>
				"""

import matplotlib.pyplot as plt
import glob, re, sys
# --- end of imports --- #


def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip()
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq.upper() } )
					header = line.strip()[1:]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq.upper() } )
	return sequences


def main( arguments ):
	"""! @brief run everything """
	
	inputdir = arguments[ arguments.index( '--in' )+1 ]
	fig_file = arguments[ arguments.index( '--fig' )+1 ]
	
	if inputdir[-1] != "/":
		inputdir += "/"
	
	filenames = glob.glob( inputdir + "*.fasta" )
	lengths = {}
	for filename in filenames:
		seqs = load_sequences( filename )
		lens = []
		for seq in seqs.values():
			lens.append( len( seq ) )
		lengths.update( { filename.split('/')[-1].split('.')[0]: lens } )


	fig, ax = plt.subplots()
	colors = [ "red", "purple", "black",  "blue", "lime" ]
	for idx, ID in enumerate( sorted( lengths.keys() ) ):
		n, b, p = ax.hist( lengths[ ID ], bins=int( ( max( lengths[ ID ] )-min( lengths[ ID ] ) )  ), color=colors[ idx ], label=ID, alpha=0.7, normed=1 )
		n = list( n )
		print n.index( max( n ) )+min( lengths[ ID ] )
		
	ax.set_xlim( 200, 700 )
	ax.set_xlabel( "PEPD length [amino acids]" )
	ax.set_ylabel( "frequency of sequence length" )
	ax.legend( bbox_to_anchor=(0.9, 0.9) )
	fig.savefig( fig_file, dpi=600 )


if __name__ == '__main__':
	if '--in' in sys.argv and '--fig' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
