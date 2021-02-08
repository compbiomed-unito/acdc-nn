import sys
import argparse

from acdc_nn import util
from acdc_nn import nn

import numpy as np
import pandas as pd

from pkg_resources import resource_filename


def load_nn(weight_path):
	num_H=[32,16]
	d=0.2
	model1, model2 = nn.ACDC(num_H,d,25)
	weight_path = resource_filename('acdc_nn', 'weights/full_dataset_TL')
	model1.load_weights(weight_path)
	return model1, model2

def	load_prot(profile_path, pdb_path, chain):
	r = dict(zip(['structure', 'pchain', 'seq', 'd_seq2pdb', 'd_pdb2seq'], util.pdb2info(pdb_path, chain)))
	r['prof'] = util.getProfile(profile_path)
	return r

def compute_nn_input(mut, prot):
	# FIXME mutation dependent
	kvar = (mut[0], prot['d_pdb2seq'][mut[1:-1]], mut[-1])
	kvar_pdb = (mut[0], mut[1:-1], mut[-1])

	dist_neigh_3d= util.get_neigh_ps(kvar_pdb, 5, prot['d_seq2pdb'], prot['pchain']) 
	list_dist_neigh_3d = dist_neigh_3d[kvar]

	# extracting features
	codif = util.getMutCod(mut)
	all_profile = util.Unified_prof(kvar[1], prot['prof'], prot['seq'], list_dist_neigh_3d)
	#all_data=[*codif,*all_profile,*np.zeros(600-len(all_profile))]
	to_predict_dir = pd.DataFrame([*codif, *all_profile, *np.zeros(600 - len(all_profile))]).T
	return to_predict_dir



def main():
	parser = argparse.ArgumentParser(description='Predict the DDG of one or more mutations.')
	parser.add_argument('--output', metavar='OUTFILE', type=argparse.FileType('w'), default=sys.stdout, help='output file path, use stdout otherwise')
	parser.add_argument('--weights', metavar='NN-WEIGHTS', 
		default=resource_filename('acdc_nn', 'weights/full_dataset_TL'), 
		help='path to the network weights, in tensorflow format') # FIXME what is the name of the format?
	subparsers = parser.add_subparsers(dest='command', metavar='CMD', required=True)

	single_parser = subparsers.add_parser('single', help='predict a single mutation')
	single_parser.add_argument('wt_mutation', metavar='MUT', help='mutation in the form Q339N')
	single_parser.add_argument('wt_profile', metavar='PROF', help='')
	single_parser.add_argument('wt_structure', metavar='PDB', help='')
	single_parser.add_argument('wt_chain', metavar='CHAIN', help='')

	inverse_parser = subparsers.add_parser('inverse', help='predict a single mutation using also the inverse structure')
	inverse_parser.add_argument('wt_mutation', metavar='MUT', help='mutation in the form Q339N')
	inverse_parser.add_argument('wt_profile', metavar='WT-PROF', help='')
	inverse_parser.add_argument('wt_structure', metavar='WT-PDB', help='')
	inverse_parser.add_argument('wt_chain', metavar='WT-CHAIN', help='')
	inverse_parser.add_argument('mt_mutation', metavar='INV-MUT', help='inverse mutation in the form Q339N')
	inverse_parser.add_argument('mt_profile', metavar='MT-PROF', help='')
	inverse_parser.add_argument('mt_structure', metavar='MT-PDB', help='')
	inverse_parser.add_argument('mt_chain', metavar='MT-CHAIN', help='')

	batch_parser = subparsers.add_parser('batch', help='predict a list of mutations from a file')
	batch_parser.add_argument('mutations', metavar='MUT-FILE', type=argparse.FileType('r'))

	#profile_parser = subparsers.add_parser('profile', help='compute the profile from multiple alignment file') # we could also simply accept a multialn file and make a profile
	#profile_parser.add_argument('multialn', metavar='ALN-FILE', type=argparse.FileType('r'))

	args = parser.parse_args()
	
	if args.command in {'single', 'inverse'}:
		# ACDN-NN building
		model1, model2 = load_nn(args.weights)

		# compute input features for the direct mutation using the wild-type protein
		wt_prot = load_prot(args.wt_profile, args.wt_structure, args.wt_chain)
		to_predict_dir = compute_nn_input(args.wt_mutation, wt_prot)
		if args.command == 'inverse': # compute input features for the inverse mutation using the mutated protein
			mt_prot = load_prot(args.mt_profile, args.mt_structure, args.mt_chain)
			to_predict_inv = compute_nn_input(args.mt_mutation, mt_prot)
		else: # compute input features for the inverse mutation inverting the input features for the direct mutation
			to_predict_inv = to_predict_dir.copy()
			to_predict_inv.iloc[:, :20] = to_predict_inv.iloc[:, :20].replace([1.0, -1.0],[-1.0, 1.0])
	
		# Making input in the proper shape 
		Xm_d, X1D_d, X3D_d = nn.mkInp(np.asarray(to_predict_dir).astype(np.float32), 500)
		Xm_i, X1D_i, X3D_i = nn.mkInp(np.asarray(to_predict_inv).astype(np.float32), 500)

		prediction=model1.predict([X3D_d, X1D_d, Xm_d , X3D_i, X1D_i, Xm_i])
		print(prediction[0][0][0], file=args.output)
	#elif args.command == 'batch':
	#elif args.command == 'profile':
	else:
		raise NotImplementedError('Command {} not implemented'.format(args.command))

if __name__ == '__main__':
	main()
