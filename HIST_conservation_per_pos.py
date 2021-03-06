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
			results_per_pos.append( { 'pos': idx, '1_aa': sorted_aa_distr[-1]['id'], '1_ratio': sorted_aa_distr[-1]['value'] / float( sum( each.values() )-gaps ) } )
	return results_per_pos


def analyze_aa_composition_per_position( alignment ):
	"""! @brief analyze amino acid composition per position over groups """
	
	# --- extract info from alignment and assign to group --- #
	final_aa = []
	
	for i in range( len( alignment.values()[0] ) ):
		final_aa.append( {} )
		for key in alignment.keys():
			if '@' in key:
				ID = key.split('@')[0]
			else:
				ID = key
			aa = alignment[ key ][ i ]
			try:
				final_aa[-1][ aa ] += 1
			except KeyError:
				final_aa[-1].update( { aa: 1 } )
	
	# --- process data --- #
	final_results = calculate_results( final_aa )
	return final_results


def construct_figure( final_results, fig_file, ref_seq ):
	"""! @brief construct figure to illustrate differences in amino acid distribution per position """
	
	colors = { 	'A': "purple", 'C': "green", "D": "blue", "E": "orange", "F": "gold", "G": "yellow", "H": "brown", "I": "black", "K": "lime", "L": "darkgreen", "M": "pink", "N": "grey", 'P': "salmon",
						'Q': "mediumspringgreen", "R": "cyan", "S": "navy", "T": "crimson", "V": "lawngreen", "W": "maroon", "Y": "magenta", "-": "grey" }
	
	### ADJUST COLORS!!!
	
	x_miss_matches = []
	y_miss_matches = []
	
	x1 = []	#is not really necessary. would be enough to use index
	y1 = []
	x1_colors = []
	x1_labels = []
	x2 = []
	y2 =[]
	x2_colors = []
	x2_labels = []
	for idx, entry in enumerate( final_results ):
		x1.append( entry['pos'] )
		x2.append( entry['pos'] )
		y1.append( entry['1_ratio'] )
		if entry['1_ratio'] > 0.2:
			x1_labels.append( entry['1_aa'] )
		else:
			x1_labels.append( "" )
		x1_colors.append( colors[ entry['1_aa'] ] )
		try:
			y2.append( -entry['2_ratio'] )
			x2_labels.append( entry['2_aa'] )
			x2_colors.append( colors[ entry['2_aa'] ] )
		except KeyError:
			y2.append( 0 )
			x2_labels.append( "" )
			x2_colors.append( "white" )
	
	x_importance = []
	y_importance = []
	for idx, each in enumerate( y1 ):
		if each > 0.95:
			x_importance.append( idx )
			y_importance.append( 0 )
	
	fig, ax = plt.subplots( figsize=(30, 5) )
	for idx in range( len( x1 ) ):
		if len( x1_labels[ idx ] ) > 0:
			if idx % 2 == 0:
				ax.plot( [ idx, idx ], [ -1, 1 ], color="grey", linewidth=0.1 )	#line to associated characters to bars
			else:
				ax.plot( [ idx, idx ], [ -0.95, 0.95 ], color="grey", linewidth=0.1 )	#line to associated characters to bars
		ax.plot( [ x1[ idx ], x1[ idx ] ], [ 0, y1[ idx ] ], color=x1_colors[ idx ], zorder=1 )
		ax.plot( [ x2[ idx ], x2[ idx ] ], [ 0, y2[ idx ] ], color=x2_colors[ idx ], zorder=1 )
		#plotting characters at slightly different positions to reduce overlap
		if idx % 2 == 0:
			ax.text( x1[ idx ], 1, x1_labels[ idx ], fontsize=3, ha="center" )
			#ax.text( x2[ idx ], -1, x2_labels[ idx ], fontsize=3, ha="center" )
		else:
			ax.text( x1[ idx ], 0.95, x1_labels[ idx ], fontsize=3, ha="center" )
			#ax.text( x2[ idx ], -0.95, x2_labels[ idx ], fontsize=3, ha="center" )
		
		# --- add sequences of Ath and rice as a reference --- #
		txt = ax.text( idx, -0.025, ref_seq[ idx ], fontsize=3, ha="center" )	#Reference sequence
		txt.set_path_effects([PathEffects.withStroke( linewidth=0.5, foreground='w')] )
		
	ax.scatter( x_importance, y_importance, color="red", s=1, zorder=2 )
	
	# domains = [ { 'name': "xxx", 'start': 1, 'end': 11, 'color': 'red' } ]
	# for domain in domains:
		# ax.plot( [ domain['start'], domain['end'] ], [ 0.1, 0.1 ], color=domain['color'] )
		# ax.text( domain['start'], 0.15, domain['name'], ha="left", color=domain[ 'color' ] )
	
	ax.set_xlim( -0.5, len( x1 )-0.5 )
	ax.set_ylim( -1, 1 )
	
	ax.set_ylabel( "frequency of top (up) and 2nd (down) amino acid" )
	ax.set_xlabel( "position within alignment" )
	ax.set_yticklabels( [ 1.0, 0.75, 0.5, 0.25, 0.0, 0.25, 0.5, 0.75, 1.0 ] )
	
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	
	plt.subplots_adjust( left=0.03, right=0.99, top=0.9, bottom=0.1 )	
	
	fig.savefig( fig_file, dpi=600 )


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
		out.write( "AlignmentPosition\tRefSeqPos\tAminoAcid_Frequency\tAminoAcid_Frequency\n" )
		for idx, entry in enumerate( final_results ):
			new_line = [ str( idx+1 ), str( ref_seq_pos[ idx ] ) ]
			new_line.append( entry['1_aa']+"_"+str( entry['1_ratio']  ) )
			try:
				new_line.append( entry['2_aa']+"_"+str( entry['2_ratio']  ) )
			except KeyError:
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


def load_values( input_file ):
	"""! @brief load values """
	
	values = []
	with open( input_file, "r" ) as f:
		headers = f.readline().strip().split('\t')[1:]
		values = {}
		for each in headers:
			values.update( { each: [ 0, 0 ] } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			for idx, aa in enumerate( parts[1:] ):
				if aa == headers[ idx ][0]:
					values[ headers[ idx ] ][0] += 1
				else:
					values[ headers[ idx ] ][1] += 1
			line = f.readline()
	return_values = []
	for each in headers:
		return_values.append( [ values[ each ][0]/float( sum( values[ each ] ) ), values[ each ][1]/float( sum( values[ each ] ) ) ] )
	return return_values, headers


def construct_bar_plot( values, labels, output_file ):
	"""! @brief construct stacked bar plot """
	
	fig, ax = plt.subplots( figsize=(15, 5) )
	my_plots = []
	colors = [ "green", "grey" ]
	for idx, each in enumerate( values ):
		for i in range( 2 ):
			if i > 0:
				my_plots.append( ax.bar( idx, each[i], width=0.2, bottom=sum( each[ :i ] ), color=colors[i] ) )
			else:
				my_plots.append( ax.bar( idx, each[i], width=0.2, color=colors[i] ) )

	my_legend = [ 	mpatches.Patch( color='green', label='conserved' ),
								mpatches.Patch( color='grey', label='other' )
							]
	ax.legend( handles=my_legend, fontsize=10, bbox_to_anchor=(0.95, 0.9) )	#loc=2, 
	ax.set_xlabel( "positions of interest", fontsize=12 )
	ax.set_ylabel( "frequency of amino acid", fontsize=12 )
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	start, end = ax.get_xlim()
	ax.xaxis.set_ticks( np.arange(0, len( values ), 1) )
	ax.set_xticklabels( labels )
	plt.subplots_adjust( bottom=0.1, left=0.05, right=0.92, top=0.95 )
	fig.savefig( output_file, dpi=600 )


def main( arguments ):
	"""! @brief runs everything """
	
	pep_file = arguments[ arguments.index( '--in' )+1 ]
	ref_seq_name = arguments[ arguments.index( '--name' )+1 ]
	output_dir = arguments[ arguments.index( '--out' )+1 ]
	if '--res' in arguments:
		res_file = arguments[ arguments.index( '--res' )+1 ]
	else:
		res_file = ""
	
	occ = 0.3
	mafft_path = "mafft"
	
	if output_dir[-1] != '/':
		output_dir += "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	peps = load_sequences( pep_file )
	in_file_name = pep_file.split('/')[-1]
	alignment_file = output_dir + in_file_name + ".aln"
	os.popen( mafft_path + " " + pep_file + " > " + alignment_file )
	
	raw_alignment = load_sequences( alignment_file )
	ref_pos_map_table, pos_to_keep = get_ref_seq_pos_for_all_residues( ref_seq_name, raw_alignment, occ )
	
	clean_alignment_file = output_dir + in_file_name + ".cln"
	alignment = construct_clean_alignment_file( raw_alignment, pos_to_keep, clean_alignment_file )
	
	final_results = analyze_aa_composition_per_position( alignment )
	ref_seq = alignment[ ref_seq_name ]
	
	fig_file = output_dir + "hist.png"
	construct_figure( final_results, fig_file, ref_seq )
	output_file = output_dir + "values.txt"
	
	generate_output_table( final_results, output_file, ref_pos_map_table )
	
	if len( res_file ) > 0:
		table_file = output_dir + "table.txt"
		with open( res_file, "r" ) as f:
			residues_of_interest = f.read().strip().split('\n')
		generate_table( alignment, residues_of_interest, table_file, ref_pos_map_table )
		
		output_file = output_dir + "conserved_pos_fig.png"
		values, labels = load_values( table_file )
		construct_bar_plot( values, labels,  output_file)


if __name__ == '__main__':
	
	if '--in' in sys.argv and '--name' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
