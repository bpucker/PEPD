### Boas Pucker ###
### Hanna Schilbert ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python conservation_per_pos.py
	--in <MULTIPLE_FASTA_FILE_WITH_PEPTIDE_SEQ_INPUT>
	--name <NAME_OF_REFERENCE_SEQUENCE>
	--out <FULL_PATH_TO_OUTPUT_DIR>
	--res <FULL_PATH_TO_RESIDUE_OF_INTEREST_FILE>
	example for file content: 'D276\\nD287\\nH370'
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import re, glob, os, sys
from operator import itemgetter
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import matplotlib.patches as mpatches
import numpy as np

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


def analyze_aa_composition_per_position( alignment, groups ):
	"""! @brief analyze amino acid composition per position over groups """
	
	total_final_results = {}
	
	for group in sorted( list( set( groups.values() ) ) ):
		# --- extract info from alignment and assign to group --- #
		final_aa = []
		for i in range( len( alignment.values()[0] ) ):
			final_aa.append( {} )
			for key in alignment.keys():
				if '@' in key:
					ID = key.split('@')[0]
					if groups[ ID ] == group:
						aa = alignment[ key ][ i ]
						try:
							final_aa[-1][ aa ] += 1
						except KeyError:
							final_aa[-1].update( { aa: 1 } )
		# --- process data --- #
		final_results = calculate_results( final_aa )
		total_final_results.update( { group: final_results } )
	return total_final_results


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


def generate_output_table( final_results, output_file, ref_seq_pos, groups ):
	"""! @brief write alignment results into output table """
	
	group_order = sorted( list( set( groups.values() ) ) )
	
	with open( output_file, "w" ) as out:
		string = "AlignmentPosition\tRefSeqPos"
		for group in group_order:
			string += "\tAAF(" + group + ")\tAAF(" + group + ")"
		out.write( string + "\n" )
		
		for idx, entry in enumerate( final_results.values()[0] ):
			new_line = [ str( idx+1 ), str( ref_seq_pos[ idx ] )  ]
			for group in group_order:
				try:
					entry = final_results[ group ][ idx ]
					new_line.append( entry['1_aa']+"_"+( str( entry['1_ratio']  )[:5] ) )
					try:
						new_line.append( entry['2_aa']+"_"+ (str( entry['2_ratio']  )[:5] ) )
					except KeyError:
						new_line.append( "-" )
				except IndexError:
					new_line.append( "-" )
			out.write( "\t".join( new_line ) + '\n' )


def generate_table( alignment, residues_of_interest, output_file, ref_pos_map_table ):
	"""! @brief check all given residues across all sequences and generate output table """
	
	# --- invert positional mapping dictionary --- #
	inverted_dict = {}
	for key in ref_pos_map_table.keys():
		inverted_dict.update( { ref_pos_map_table[ key ]: key } )
	
	# --- run analysis --- #
	with open( output_file, "w" ) as out:
		out.write( "SeqID\t" + "\t".join( residues_of_interest ) + '\n' )
		for key in sorted( alignment.keys() ):
			new_line = [ key ]
			for residue in residues_of_interest:
				try:
					real_pos = inverted_dict[ int( residue[1:] ) ]
					new_line.append(alignment[ key ][ real_pos ]  )
				except KeyError:
					new_line.append( '?' )
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


def load_values( input_file, groups ):
	"""! @brief load values """
	
	group_order = sorted( list( set( groups.values() ) ) )
	all_values = {}
	for group in group_order:
		all_values.update( { group: {} } )
	
	with open( input_file, "r" ) as f:
		headers = f.readline().strip().split('\t')[1:]
		for each in headers:
			for group in group_order:
				all_values[ group ].update( { each: [ 0, 0 ] } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			group = groups[ parts[0].split('@')[0] ]
			for idx, aa in enumerate( parts[1:] ):
				if aa == headers[ idx ][0]:
					all_values[ group ][ headers[ idx ] ][0] += 1
				else:
					all_values[ group ][ headers[ idx ] ][1] += 1
			line = f.readline()
	all_return_values = {}
	for group in group_order:
		all_return_values.update( { group: [] } )
	for each in headers:
		for group in group_order:
			all_return_values[ group ].append( all_values[ group ][ each ][0]/float( sum( all_values[ group ][ each ] ) ) )
	return all_return_values, headers


def construct_bar_plot( values, labels, output_file, groups ):
	"""! @brief construct stacked bar plot """
	
	group_order = [ "animal", "plant", "fungi", "bacteria", "archaea" ]
	
	fig, ax = plt.subplots( figsize=(15, 5) )
	colors = [ "red", "lime", "black", "orange", "purple" ]
	for idx in range( len( values.values()[0] ) ):
		for k, group in enumerate( group_order ):
			ax.bar( idx+k*0.1, values[ group ][ idx ], width=0.1, color=colors[k] )

	my_legend = [ ]
	for k, group in enumerate( group_order ):
		my_legend.append( mpatches.Patch( color=colors[k], label=group ) )
		
	ax.legend( handles=my_legend, fontsize=10, bbox_to_anchor=(0.95, 0.9) )	#loc=2, 
	ax.set_xlabel( "positions of interest", fontsize=12 )
	ax.set_ylabel( "frequency of conserved amino acid", fontsize=12 )
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	start, end = ax.get_xlim()
	ax.xaxis.set_ticks( np.arange(0, len( values.values()[0] ), 1) )
	ax.set_xticklabels( labels )
	plt.subplots_adjust( bottom=0.1, left=0.05, right=0.92, top=0.95 )
	fig.savefig( output_file, dpi=600 )


def load_groups( group_file ):
	"""! @brief load groups from given file """
	
	groups = {}
	with open( group_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			groups.update( { parts[0]: parts[1] } )	#.split('@')[0]
			line = f.readline()
	return groups


def main( arguments ):
	"""! @brief runs everything """
	
	pep_file = arguments[ arguments.index( '--in' )+1 ]
	ref_seq_name = arguments[ arguments.index( '--name' )+1 ]
	output_dir = arguments[ arguments.index( '--out' )+1 ]
	if '--res' in arguments:
		res_file = arguments[ arguments.index( '--res' )+1 ]
	else:
		res_file = ""
	
	group_file = arguments[ arguments.index( '--group' )+1 ]
	
	occ = 0.3
	mafft_path = "mafft"
	
	if output_dir[-1] != '/':
		output_dir += "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	peps = load_sequences( pep_file )
	groups = load_groups( group_file )
	in_file_name = pep_file.split('/')[-1]
	alignment_file = output_dir + in_file_name + ".aln"
	os.popen( mafft_path + " " + pep_file + " > " + alignment_file )
	
	raw_alignment = load_sequences( alignment_file )
	ref_pos_map_table, pos_to_keep = get_ref_seq_pos_for_all_residues( ref_seq_name, raw_alignment, occ )
	
	clean_alignment_file = output_dir + in_file_name + ".cln"
	alignment = construct_clean_alignment_file( raw_alignment, pos_to_keep, clean_alignment_file )
	
	total_final_results = analyze_aa_composition_per_position( alignment, groups )
	ref_seq = alignment[ ref_seq_name ]
	
	output_file = output_dir + "values.txt"
	generate_output_table( total_final_results, output_file, ref_pos_map_table, groups )	#per residue table showing conservation in group
	
	if len( res_file ) > 0:
		table_file = output_dir + "table.txt"
		with open( res_file, "r" ) as f:
			residues_of_interest = f.read().strip().split('\n')
		generate_table( alignment, residues_of_interest, table_file, ref_pos_map_table )	#construct table per residue of interest and species
		
		output_file = output_dir + "conserved_pos_fig.png"
		values, labels = load_values( table_file, groups )
		construct_bar_plot( values, labels,  output_file, groups )


if __name__ == '__main__':
	
	if '--in' in sys.argv and '--name' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
