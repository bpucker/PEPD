### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	pyhton alignment_sensitivity_check.py
	--in <NAME_OF_INPUT_FILE>
	--out <NAME_OF_OUTPUT_DIR>
	--name <NAME_OF_REFERENCE_SEQUENCE>
					"""


import os, random, glob, sys
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter

# --- end of imports --- #

def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip().split(" ")[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:].split(" ")[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	return sequences


def calculate_results( data ):
	"""! @brief evaluate amino acid distribution per position and return values """
	
	results_per_pos = []
	for idx, each in enumerate( data ):
		aa_counter_dict_list = []
		for key in each.keys():
			aa_counter_dict_list.append( { 'id': key, 'value': each[ key ] } )
		sorted_aa_distr = sorted( aa_counter_dict_list, key=itemgetter('value') )
		try:
			gaps = 0	#each['-']
		except KeyError:
			gaps = 0
		try:
			results_per_pos.append( { 'pos': idx, '1_aa': sorted_aa_distr[-1]['id'], '1_ratio': sorted_aa_distr[-1]['value'] / float( sum( each.values() )-gaps ),
														'2_aa': sorted_aa_distr[-2]['id'], '2_ratio': sorted_aa_distr[-2]['value'] / float( sum( each.values() )-gaps )
													} )
		except IndexError:
			try:
				results_per_pos.append( { 'pos': idx, '1_aa': sorted_aa_distr[-1]['id'], '1_ratio': sorted_aa_distr[-1]['value'] / float( sum( each.values() )-gaps ) } )
			except IndexError:
				pass	#print sorted_aa_distr
	return results_per_pos


def analyze_aa_composition_per_position( alignment ):
	"""! @brief analyze amino acid composition per position """
	
	final_aa = []
	for i in range( len( alignment.values()[0] ) ):
		final_aa.append( {} )
		for key in alignment.keys():
			aa = alignment[ key ][ i ]
			try:
				final_aa[-1][ aa ] += 1
			except KeyError:
				final_aa[-1].update( { aa: 1 } )
	# --- process data --- #
	final_results = calculate_results( final_aa )
	return final_results


def get_ref_seq_pos_for_all_residues( ref_seq_name, raw_alignment, occ ):
	"""! @brief match alignment positions to refseq positions """
	
	ref_seq_pos_mapping_table = {}
	ref_seq_pos = 1
	clean_alignment_pos = 0
	pos_to_keep = []
	for idx, aa in enumerate( raw_alignment[ ref_seq_name ] ):
		occupancy = []
		for key in raw_alignment.keys():
			if raw_alignment[ key ][ idx ] == '-':
				occupancy.append( 0 )
			else:
				occupancy.append( 1 )
		if sum( occupancy ) / float( len( occupancy ) ) > occ:	#alignment column will be kept
			if aa == "-":
				ref_seq_pos_mapping_table.update( { clean_alignment_pos: "-" } )
			else:
				ref_seq_pos_mapping_table.update( { clean_alignment_pos: ref_seq_pos } )
				ref_seq_pos += 1
			clean_alignment_pos += 1
			pos_to_keep.append( idx )
		else:	#alignment column will be removed
			if aa == "-":
				pass
			else:
				ref_seq_pos += 1
	return ref_seq_pos_mapping_table, pos_to_keep


def generate_output_table( final_results, output_file, ref_seq_pos ):
	"""! @brief write alignment results into output table """
	
	with open( output_file, "w" ) as out:
		out.write( "AlignmentPosition\tRefSeqPos\tAAF\n" )
		
		for idx, entry in enumerate( final_results ):
			new_line = [ str( idx+1 ), str( ref_seq_pos[ idx ] )  ]
			try:
				entry = final_results[ idx ]
				new_line.append( entry['1_aa']+"_"+( str( entry['1_ratio']  )[:5] ) )
				try:
					new_line.append( entry['2_aa']+"_"+ (str( entry['2_ratio']  )[:5] ) )
				except KeyError:
					new_line.append( "-" )
			except IndexError:
				new_line.append( "-" )
			out.write( "\t".join( new_line ) + '\n' )


def construct_clean_alignment_file( raw_alignment, pos_to_keep, clean_alignment_file ):
	"""! @brief construct clean alignment file and return clean alignment """
	
	clean_alignment = {}
	with open( clean_alignment_file, "w" ) as out:
		for key in sorted( raw_alignment.keys() ):
			seq = raw_alignment[ key ]
			new_seq = []
			for idx in pos_to_keep:
				new_seq.append( seq[ idx ] )
			new_seq = "".join( new_seq )
			clean_alignment.update( { key: new_seq } )
			out.write( '>' + key + '\n' + new_seq + '\n' )
	return clean_alignment


def assess_conservation( pep_file, ref_seq_name, alignment_file, occ ):
	"""! @brief runs everything """
		
	peps = load_sequences( pep_file )
	raw_alignment = load_sequences( alignment_file )
	ref_pos_map_table, pos_to_keep = get_ref_seq_pos_for_all_residues( ref_seq_name, raw_alignment, occ )
	
	clean_alignment_file = alignment_file + ".cln"
	alignment = construct_clean_alignment_file( raw_alignment, pos_to_keep, clean_alignment_file )
	
	total_final_results = analyze_aa_composition_per_position( alignment )
	ref_seq = alignment[ ref_seq_name ]
	
	output_file = alignment_file + ".values.txt"
	generate_output_table( total_final_results, output_file, ref_pos_map_table )	#per residue table showing conservation in group


def load_values( value_file ):
	"""! @brief load all animal values """
	
	values = {}
	with open( value_file, "r" ) as f:
		f.readline()	#header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				values.update( { int( parts[1] ): float( parts[2].split('_')[1] ) } )
			except ValueError:
				pass
			line = f.readline()
	return values


def compare_results( all_value_files, figfile, ref_seq_name ):
	"""! @brief analyze conservation per residue differences """
	
	# --- calculating and saving values --- #
	values = [ ]
	for filename in all_value_files:
		values.append( load_values( filename ) )

	min_values = []
	average_values = []
	max_values = []
	value_file = figfile + ".txt"
	with open( value_file, "w" ) as out:
		for key in sorted( values[0].keys() ):
			conservation = []
			for each in values:
				try:
					conservation.append( each[ key ] )
				except KeyError:
					print "MissingKey: " + str( key )
			min_value = min( conservation )
			avg_value = np.median( conservation )
			max_value = max( conservation )
			
			min_values.append( min_value )
			average_values.append( avg_value )
			max_values.append( max_value )
			out.write( "\t".join( map( str, [ key, min_value, avg_value, max_value ] ) ) + '\n' )
	
	# --- constructing the figure --- #
	fig, ax = plt.subplots( figsize=(30, 5) )
	ax.plot( average_values, color="green", linewidth=3 )
	ax.plot( min_values, color="red", linewidth=1 )
	ax.plot( max_values, color="red", linewidth=1 )
	ax.set_xlabel( "position in " + ref_seq_name )
	ax.set_ylabel( "residue conservation" )
	
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	ax.set_frame_on(False)
	plt.subplots_adjust( left=0.03, right=0.9999, top=0.99, bottom=0.1 )
	
	fig.savefig( figfile, dpi=300 )


def main( arguments ):
	input_file = arguments[ arguments.index( '--in' ) ]
	output_dir = arguments[ arguments.index( '--out' ) ]
	
	ref_seq_name = arguments[ arguments.index( '--name' ) ]

	iterations = 50
	occ = 0.7

	if output_dir[-1] != '/':
		output_dir += "/"

	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	# --- moving input file into right place --- #
	os.popen( "cp " + input_file + " " + output_dir )
	input_file = output_dir + input_file.split('/')[-1]
	
	# --- generate initial MAFFT alignment --- #
	alignment_file = output_dir + "alignmennt.mafft.aln"
	if not os.path.isfile( alignment_file ):
		os.popen( "mafft " + input_file + " > " + alignment_file )
		assess_conservation( input_file, ref_seq_name, alignment_file, occ )
	
	# --- generate initial CLUSTAL omega alignment --- #
	alignment_file = output_dir + "alignmennt.clustalo.aln"
	if not os.path.isfile( alignment_file ):
		os.popen( "clustalo -i " + input_file + " -o " + alignment_file )
		assess_conservation( input_file, ref_seq_name, alignment_file, occ )
	
	# --- generate initial MUSCLE alignment --- #
	alignment_file = output_dir + "alignmennt.muscle.aln"
	if not os.path.isfile( alignment_file ):
		os.popen( "muscle -in " + input_file + " -out " + alignment_file + " -quiet" )
		assess_conservation( input_file, ref_seq_name, alignment_file, occ )
	

	seqs = load_sequences( input_file )
	all_seq_ids = seqs.keys()
	for i in range( iterations ):
		# --- generating shuffled input sequence file --- #
		output_file = output_dir + str( i ) + ".fasta"
		if not os.path.isfile( output_file ):
			random.shuffle( all_seq_ids )
			with open( output_file, "w" ) as out:
				for key in all_seq_ids:
					out.write( '>' + key + '\n' + seqs[ key ] + '\n' )
		
		# --- MAFFT alignment --- #
		alignment_file = output_file + ".mafft.aln"
		if not os.path.isfile( alignment_file ):
			os.popen( "mafft " + output_file + " > " + alignment_file )
			assess_conservation( output_file, ref_seq_name, alignment_file, occ )
		
		# --- CLUSTAL omega alignment --- #
		alignment_file = output_file + "alignmennt.clustalo.aln"
		if not os.path.isfile( alignment_file ):
			os.popen( "clustalo -i " + output_file + " -o " + alignment_file )
			assess_conservation( output_file, ref_seq_name, alignment_file, occ )
		
		# --- MUSCLE alignment --- #
		alignment_file = output_file + "alignmennt.muscle.aln"
		if not os.path.isfile( alignment_file ):
			os.popen( "muscle -in " + output_file + " -out " + alignment_file + " -quiet" )
			assess_conservation( output_file, ref_seq_name, alignment_file, occ )
	
	# --- assess alignment results on sequence basis --- #
	all_clean_alignment_files = glob.glob( output_dir + "*.cln" )
	results_per_seq = {}
	for key in all_seq_ids:
		result_seqs = []
		for alignmentfile in all_clean_alignment_files:
			result_seqs.append( load_sequences( alignmentfile )[ key ] )
		results_per_seq.update( { key: len( list( set( result_seqs ) ) ) } )

	fig_file = output_dir + "different_alignment_results_per_seq.png"
	fig, ax = plt.subplots()
	ax.hist( results_per_seq.values(), bins=10, color="green" )
	ax.set_xlabel( "number of different alignment results" )
	ax.set_ylabel( "number of sequences" )
	fig.savefig( fig_file, dpi=300 )
	
	# --- alignment sensitivity --- #
	figfile = output_dir + "impact_on_conservation_per_pos.png"
	all_value_files = glob.glob( output_dir + "*.values.txt" )
	compare_results( all_value_files, figfile, ref_seq_name )


if __name__ == '__main__':
	if '--in' in sys.argv and '--out' in sys.argv and '--name' in sys.argv:
		main()
	else:
		sys.exit( __usage__ )
