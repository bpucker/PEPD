### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """	
	python run_BUSCO_across_specs.py
	--in <FULL_PATH_TO_INPUT_DIR>
	--out <FULL_PATH_TO_OUTPUT_DIR>
	WARNING: installation of several tools required
				"""

import sys, os, re, glob
from operator import itemgetter

# --- end of imports --- #

def run_BUSCO( input_file, prefix, busco_path, augustus_path, augustus_config_path ):
	"""! @brief run BUSCO in genome mode on given assembly file """
	
	os.chdir( prefix )
	os.environ["PATH"] = augustus_path + ":" + os.environ["PATH"]
	print os.environ["PATH"]
	os.environ["AUGUSTUS_CONFIG_PATH"] = augustus_config_path
	print os.environ["AUGUSTUS_CONFIG_PATH"]
	cmd = "python " + busco_path + " --in " + input_file + " --cpu 10 -m proteins --out busco_run > " + prefix +"log.txt"
	os.popen( cmd )


def main( arguments ):
	"""! @brief run everything """
	
	input_dir = arguments[ arguments.index( '--in' )+1 ]
	output_dir = arguments[ arguments.index( '--out' )+1 ]
	
	if input_dir[-1] != '/':
		input_dir += "/"
	if output_dir[-1] != '/':
		output_dir += "/"
	
	output_file = output_dir + "SUMMARY.txt"
	
	augustus = "augustus"
	augustus_seqs_ex_script = "getAnnoFasta.pl"
	
	busco = "run_BUSCO.py"
	augustus_path = "augustus/bin/"
	augustus_config_path = "augustus/config/"
	
	active = True
	
	filenames = glob.glob( input_dir + "*.fasta" ) + glob.glob( input_dir + "*.faa" ) + glob.glob( input_dir + "*.fa" ) 
	results = {}
	for filename in filenames:
		ID = filename.split('/')[-1].split('.')[0]
		working_dir = output_dir + ID + '/'
		if not os.path.exists( working_dir ):
			os.makedirs( working_dir )
		
		# --- run BUSCO v3 --- #
		prefix = working_dir + "busco/"
		if not os.path.exists( prefix ):
			os.makedirs( prefix )
		if active:
			run_BUSCO( filename, prefix, busco, augustus_path, augustus_config_path )
		
		# ---- collect results --- #
		log_file = prefix + "log.txt"
		with open( log_file, "r" ) as f:
			content = f.read()
		try:
			result_string = re.findall( "C:\d+\.\d+%\[S:\d+\.\d+%,D:\d+\.\d+%\],F:\d+\.\d+%,M:[-]*\d+\.\d+%,n:\d+", content )[0]
		except IndexError:
			print content
		results.update( { ID: result_string } )
	
	with open( output_file, "w" ) as out:
		for key in results.keys():
			out.write( key + '\t' + results[ key ] + '\n' )


if __name__ == "__main__":
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	
