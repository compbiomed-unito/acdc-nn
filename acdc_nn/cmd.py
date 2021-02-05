import sys
import argparse



from acdc_nn import util
from acdc_nn import nn

import numpy as np
import pandas as pd

#path_weights="weights/Pesi_ACDCNN_postTL"
#/usr/local/lib/python3.8/site-packages/acdc_nn
from pkg_resources import resource_filename

def main():
	print(resource_filename('acdc_nn', 'weights/Pesi_ACDCNN_postTL'))
	parser = argparse.ArgumentParser(description='Predict the DDG of one or more mutations.')
	# use subcommands
	# https://docs.python.org/3/library/argparse.html#sub-commands
	
	#parser.add_argument('profile', metavar='PROF', type=argparse.FileType('r'), help='')
	#parser.add_argument('structure', metavar='PDB', type=argparse.FileType('r'), help='')
	parser.add_argument('mutation', metavar='MUT', help='mutation in the form Q339N') # FIXME better help
	parser.add_argument('profile', metavar='PROF', help='')
	parser.add_argument('structure', metavar='PDB', help='')
	parser.add_argument('chain', metavar='CHAIN', help='')

	parser.add_argument('--inverse', nargs=4, help='')
	#parser.add_argument('--from-file', type=argparse.FileType('r'), help='')
	

	args = parser.parse_args()
	#print(args)

	#if len(sys.argv)!= 5 and len(sys.argv)!=9:

	#	print ("\n usage with 1 structure: mut prof pdb chain  \n usage with 2 structures:  mut  prof_w pdb_w chain_w mut_inv prof_m pdb_m chain_m")
	#	sys.exit()

	
	# information processing
	# get structure
	structure, pchain, seq, d_seq2pdb, d_pdb2seq  = util.pdb2info(args.structure, args.chain)
	prof = util.getProfile(args.profile)

	# FIXME mutation dependent
	mut = args.mutation
	kvar = (mut[0], d_pdb2seq[mut[1:-1]], mut[-1])
	kvar_pdb = (mut[0], mut[1:-1], mut[-1])

	dist_neigh_3d= util.get_neigh_ps(kvar_pdb, 5, d_seq2pdb, pchain) 
	list_dist_neigh_3d = dist_neigh_3d[kvar]

	# extracting features
	codif = util.getMutCod(mut)
	all_profile = util.Unified_prof(kvar[1], prof,seq, list_dist_neigh_3d)
	#all_data=[*codif,*all_profile,*np.zeros(600-len(all_profile))]
	To_predict_dir = pd.DataFrame([*codif, *all_profile, *np.zeros(600 - len(all_profile))]).T

	# if clauses for the second structure

	if args.inverse:
	#if len(sys.argv)==9:
		mut_inv, prof_path_m, pdb_path_m, chain_m = args.inverse
		#mut_inv=sys.argv[5]
		#prof_path_m = sys.argv[6]
		#pdb_path_m = sys.argv[7]
		#chain_m = sys.argv[8]
		
		# information processing
		# get structure
		structure_m, pchain_m, seq_m, d_seq2pdb_m, d_pdb2seq  = util.pdb2info(pdb_path_m, chain_m)
		prof_m = util.getProfile(prof_path_m)

		kvar_m=(mut_inv[0],d_pdb2seq[mut_inv[1:-1]],mut_inv[-1])
		kvar_pdb_m=(mut_inv[0],mut_inv[1:-1],mut_inv[-1])
		
		dist_neigh_3d_m= util.get_neigh_ps(kvar_pdb_m,5,d_seq2pdb_m,pchain_m) 
		list_dist_neigh_3d_m = dist_neigh_3d_m[kvar_m]
	  
		# extracting features

		codif_m=util.getMutCod(mut_inv)
		all_profile_m = util.Unified_prof(kvar_m[1],prof_m,seq_m, list_dist_neigh_3d_m)
		To_predict_inv=pd.DataFrame([*codif_m,*all_profile_m,*np.zeros(600-len(all_profile_m))]).T

	 
	# ACDN-NN 1-input
	else:
	#if len(sys.argv)==5:
		To_predict_inv=To_predict_dir.copy()
		To_predict_inv.iloc[:,:20]=To_predict_inv.iloc[:,:20].replace([1.0,-1.0],[-1.0,1.0])	  

	
	# ACDN-NN building
	# Making input in the proper shape 

	Xm_d, X1D_d, X3D_d = nn.mkInp(np.asarray(To_predict_dir).astype(np.float32), 500)
	Xm_i, X1D_i, X3D_i = nn.mkInp(np.asarray(To_predict_inv).astype(np.float32), 500)  
  
	num_H=[32,16]
	d=0.2
	model1, model2 = nn.ACDC(num_H,d,25)
	model1.load_weights(path_weights)
	#model1.summary()
	prediction=model1.predict([X3D_d, X1D_d, Xm_d , X3D_i, X1D_i, Xm_i])
	print(prediction[0][0][0])		


if __name__ == '__main__':
	main()
