import sys
import argparse

from acdc_nn import util
from acdc_nn import nn

import numpy as np
import pandas as pd

from pkg_resources import resource_filename
import functools


def load_nn(weight_path):
	num_H=[32,16]
	d=0.2
	model1, model2 = nn.build_acdc_3d(num_H,d,25)
	weight_path = resource_filename('acdc_nn', 'weights/full_dataset_TL')
	model1.load_weights(weight_path)
	return model1, model2

@functools.lru_cache(100) # memory cache to avoid reloading the same protein data multiple times
def	load_prot(profile_path, pdb_path, chain):
	r = dict(zip(['structure', 'pchain', 'seq', 'd_seq2pdb', 'd_pdb2seq'], util.pdb2info(pdb_path, chain)))
	r['prof'] = util.getProfile(profile_path)
	r['profile_path'] = profile_path
	r['pdb_path'] = pdb_path
	return r

def compute_nn_input(mut, prot):
	kvar = (mut[0], prot['d_pdb2seq'][mut[1:-1]], mut[-1]) # position in the sequence
	kvar_pdb = (mut[0], mut[1:-1], mut[-1]) # position in the pdb
	from_residue = prot['seq'][int(kvar[1]) - 1]
	if kvar[0] != from_residue:
		print('WARNING: mutation {} does not correspond to residue {} from {}'.format(kvar[0], from_residue, prot['pdb_path']))

	dist_neigh_3d = util.get_neigh_ps(kvar_pdb, 5, prot['d_seq2pdb'], prot['pchain']) 
	list_dist_neigh_3d = dist_neigh_3d[kvar]

	# extracting features
	codif = util.getMutCod(mut)
	all_profile = util.Unified_prof(kvar[1], prot['prof'], prot['seq'], list_dist_neigh_3d)
	#all_data=[*codif,*all_profile,*np.zeros(600-len(all_profile))]
	return pd.DataFrame([*codif, *all_profile, *np.zeros(600 - len(all_profile))]).T

def predict(model, wt_mutation, wt_profile, wt_structure, wt_chain, mt_mutation, mt_profile, mt_structure, mt_chain):
		wt_prot = load_prot(wt_profile, wt_structure, wt_chain)
		to_predict_dir = compute_nn_input(wt_mutation, wt_prot)
		if mt_mutation is None: # compute input features for the inverse mutation inverting the input features for the direct mutation
			to_predict_inv = to_predict_dir.copy()
			to_predict_inv.iloc[:, :20] = to_predict_inv.iloc[:, :20].replace([1.0, -1.0],[-1.0, 1.0])
		else: # compute input features for the inverse mutation using the mutated protein
			mt_prot = load_prot(mt_profile, mt_structure, mt_chain)
			to_predict_inv = compute_nn_input(mt_mutation, mt_prot)
	
		# Making input in the proper shape 
		Xm_d, X1D_d, X3D_d = nn.mkInp_3d(np.asarray(to_predict_dir).astype(np.float32), 500)
		Xm_i, X1D_i, X3D_i = nn.mkInp_3d(np.asarray(to_predict_inv).astype(np.float32), 500)

		prediction = model.predict([X3D_d, X1D_d, Xm_d , X3D_i, X1D_i, Xm_i])
		return prediction[0][0][0]
	

def main():
	parser = argparse.ArgumentParser(description='Predict the DDG of point mutations.', epilog='''
	Notes:
	Mutations are written as XNY, meaning that the residue X at position N changes to Y. X and Y are given as a one letter amino acid code and 
	N is 1-based and referred to the the PDB numbering of the relevant chain, and not the position on the sequence.
	PDB and profile files will be automatically decompressed (by gzip) if the paths end with ".gz".
	The batch file is a tab-separated table with a row for each mutation. Two row formats are possible, either the four field:
	MUT PROF PDB CHAIN
	or the eight field:
	MUT WT-PROF WT-PDB WT-CHAIN INV-MUT MT-PROF MT-PDB MT-CHAIN
	''')
	parser.add_argument('--output', metavar='OUTFILE', 
		type=argparse.FileType('w'), default=sys.stdout, 
		help='output file path, use standard output otherwise')
	parser.add_argument('--weights', metavar='NN-WEIGHTS', 
		default=resource_filename('acdc_nn', 'weights/full_dataset_TL'), 
		help='path to the network weights, in tensorflow format') # FIXME what is the name of the format?
	subparsers = parser.add_subparsers(metavar='COMMAND', dest='command', title='COMMAND')

	single_parser = subparsers.add_parser('single', help='predict a single mutation')
	single_parser.add_argument('wt_mutation', metavar='MUT', help='mutation (e.g. Q339N)')
	single_parser.add_argument('wt_profile', metavar='PROF', help='profile file path')
	single_parser.add_argument('wt_structure', metavar='PDB', help='PDB file path')
	single_parser.add_argument('wt_chain', metavar='CHAIN', 
		help='ID of the PDB chain where the mutation occurs')

	inverse_parser = subparsers.add_parser('inverse', 
		help='predict a single mutation using also the mutated protein information')
	inverse_parser.add_argument('wt_mutation', metavar='MUT', help='direct mutation (e.g. Q339N)')
	inverse_parser.add_argument('wt_profile', metavar='WT-PROF', 
		help='profile file path for the wild-type protein')
	inverse_parser.add_argument('wt_structure', metavar='WT-PDB', 
		help='PDB path for the wild-type protein')
	inverse_parser.add_argument('wt_chain', metavar='WT-CHAIN', 
		help='ID of the wild-type PDB chain where the mutation occurs')
	inverse_parser.add_argument('mt_mutation', metavar='INV-MUT', help='inverse mutation (e.g. N339Q)')
	inverse_parser.add_argument('mt_profile', metavar='MT-PROF', help='profile file path for the mutated protein')
	inverse_parser.add_argument('mt_structure', metavar='MT-PDB', 
		help='PDB file path for the mutated protein')
	inverse_parser.add_argument('mt_chain', metavar='MT-CHAIN', 
		help='ID of the mutated PDB chain where the mutation occurs')

	batch_parser = subparsers.add_parser('batch', help='predict a list of mutations from a file')
	batch_parser.add_argument('mutation_table', metavar='MUT-FILE', type=argparse.FileType('r'), 
		help='tab-separated mutation table, see later for format guidelines')

	#profile_parser = subparsers.add_parser('profile', help='compute the profile from multiple alignment file') # we could also simply accept a multialn file and make a profile
	#profile_parser.add_argument('multialn', metavar='ALN-FILE', type=argparse.FileType('r'))

	args = parser.parse_args()

	if args.command is None:
		parser.error('the following arguments are required: COMMAND')
	
	base_arg_list = ['wt_mutation', 'wt_profile', 'wt_structure', 'wt_chain', 'mt_mutation', 'mt_profile', 'mt_structure', 'mt_chain']

	model1, model2 = load_nn(args.weights) # needed by all command
	if args.command in {'single', 'inverse'}:
		args_dict = vars(args)
		mutation_args = [args_dict.get(a) for a in base_arg_list]
		ddg = predict(model1, *mutation_args)
		print(ddg, file=args.output)
	elif args.command == 'batch':
		muts = pd.read_csv(args.mutation_table, sep='\t', names=base_arg_list)
		for row, mut in muts.iterrows():
			# check na values
			if mut[base_arg_list[:4]].count() < 4:
				print('ERROR at line {}: missing values in the first four columns of batch file'.format(row + 1), file=sys.stderr)
				sys.exit(1)
			if mut[base_arg_list[4:]].count() not in {0, 4}:
				print('ERROR at line {}: columns 4-8 of the batch file must be either all missing or all present'.format(row + 1), file=sys.stderr)
				sys.exit(1)
					
			mutation_args = [None if pd.isnull(v) else v for v in mut.values]
			ddg = predict(model1, *mutation_args)
			print(ddg, file=args.output)
		
	#elif args.command == 'profile':
	else:
		raise NotImplementedError('Command {} not implemented'.format(args.command))

if __name__ == '__main__':
	main()
