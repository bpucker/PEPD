### Boas Pucker ###
### Hanna Schilbert ### 
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python construct_new_data_base.py
	--in <FULL_PATH_TO_TAXON_OUTPUTFILE>
	--out <FULL_PATH_TO_TAXON_OUTPUTFILE>
	--taxon_out <FULL_PATH_TO_TAXON_OUTPUTFILE>
					"""

import os, glob, sys

# --- end of imports --- #

def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip().split(' ')[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:].split(' ')[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	return sequences


def main( arguments ):
	"""! @brief run all parts of this script """

	input_dir = arguments[ arguments.index( '--in' )+1 ]
	output_dir = arguments[ arguments.index( '--out' )+1 ]
	taxon_file = arguments[ arguments.index( '--taxon_out' )+1 ]	#output file with IDs

	extensions = [ ".fa", ".fasta", ".faa" ]

	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )

	# --- get all subject files --- #
	subject_files = []
	for extension in extensions:
		subject_files += glob.glob( input_dir + "*" + extension )
	print "number of detected subject files: " + str( len( subject_files ) )

	taxons = []
	for filename in subject_files:
		taxons.append( filename.split('/')[-1].split('.')[0] )

	# --- consruct new databases based on selected species names --- #
	specs_done = []
	with open( taxon_file, "w" ) as doc:
		for filename in subject_files:
			taxon_ID = filename.split('/')[-1].split('.')[0]
			if len( taxon_ID ) > 6:
				mod_tax_ID = taxon_ID[:3].upper() + taxon_ID[-4:]
			else:
				mod_tax_ID = taxon_ID
			if taxon_ID in taxons:
				seqs = load_sequences( filename )
				output_file = output_dir + mod_tax_ID + ".pep.fa"
				with open( output_file, "w" ) as out:
					for idx, key in enumerate( seqs.keys() ):
						out.write( '>' + mod_tax_ID + "@" + str( idx+1 ) + '\n' + seqs[ key ].upper().replace( '*','' ) + '\n' )
				doc.write( mod_tax_ID + '\t' + taxon_ID + '\n' )
				specs_done.append( taxon_ID )

	# --- find missing species --- #
	print "missing data sets: "
	for spec in taxons:
		if spec not in specs_done:
			print spec


if __name__ == '__main__':
	
	if '--in' in sys.argv and '--out' in sys.argv and '--taxon_out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
