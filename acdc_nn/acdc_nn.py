import sys
import argparse

from acdc_nn import util
from acdc_nn import nn

import numpy as np
import pandas as pd

from pkg_resources import resource_filename

#@functools.lru_cache(100) # memory cache to avoid reloading the same protein data multiple times
def	load_prot(profile_path, pdb_path=None, chain=None):
	if pdb_path is None:
		r = {}
	else:
		assert chain is not None
		r = dict(zip(['structure', 'pchain', 'seq', 'd_seq2pdb', 'd_pdb2seq'], util.pdb2info(pdb_path, chain)))
	r['profile'] = util.getProfile(profile_path)
	r['pprofile'] = util.load_profile(profile_path)
	if False:
		try:
			r['pprofile'] = pd.read_csv(profile_path, sep='\t', index_col='POS')
		except ValueError:
			r['pprofile'] = pd.read_csv(profile_path, sep=' ', index_col='POS')
	r['profile_path'] = profile_path
	r['pdb_path'] = pdb_path
	return r


class ACDCSeq:
	def __init__(self, weights=True):
		self.nn = nn.build_acdc_seq(num_H=[128, 64], d=0.35)
		if weights is False:
			pass
		elif weights is True:
			self.load_weights()
		else:
			self.load_weights(weights)
		
	def load_weights(self, path=None):
		if path is None:
			path = resource_filename('acdc_nn', 'weights/acdc_nn_seq.h5')
		self.nn.load_weights(path)

	def predict(self, wt_mut, wt_prof): # TODO accept mt_prot
		codif = util.getMutCod(str(wt_mut))
		
		# local profile
		#pos = int(wt_mut[1:-1])
		wt_prof_ctx = wt_prof.get_context(wt_mut.aa_pos, context=3)
		#profile = util.profile_context(wt_prof, pos, 3)
		if False:
			profile_old = np.array(util.prof(wt_mut, wt_prot['profile'], wt_prot['profile'].keys()))
			if not np.array_equal(profile, profile_old):
				print(profile.shape)
				print(profile)
				print(profile_old.shape)
				print(profile_old)
				assert False

		# direct mutation input preparation
		To_predict_dir = pd.DataFrame([*codif,*wt_prof_ctx]).T
		Xm_d, X1D_d = nn.mkInp_seq(np.asarray(To_predict_dir).astype(np.float32))

		# inverse mutation input preparation
		To_predict_inv = To_predict_dir.copy()
		To_predict_inv.iloc[:,:20]=To_predict_inv.iloc[:,:20].replace([1.0,-1.0],[-1.0,1.0])
		Xm_i, X1D_i = nn.mkInp_seq(np.asarray(To_predict_inv).astype(np.float32))

		prediction=self.nn.predict([X1D_d, Xm_d, X1D_i, Xm_i])
		return prediction[0][0][0]

class Structure:
	def __init__(self, path, chain):
		names = ['structure', 'pchain', 'seq', 'd_seq2pdb', 'd_pdb2seq']
		values = util.pdb2info(path, chain)
		for n, v in zip(names, values):
			setattr(self, n, v)
		self.pdb_path = path
		self.pdb_chain = chain
		# FIXME maybe pdb_chain is the same as pchain
		
def prepare_3d_input(mut, prof, struct):
	kvar = (mut[0], struct.d_pdb2seq[mut[1:-1]], mut[-1]) # position in the sequence
	kvar_pdb = (mut[0], mut[1:-1], mut[-1]) # position in the pdb
	from_residue = struct.seq[int(kvar[1]) - 1]
	if kvar[0] != from_residue:
		print('WARNING: mutation {} does not correspond to residue {} from {}'.format(kvar[0], from_residue, struct.pdb_path)) # FIXME use logging

	dist_neigh_3d = util.get_neigh_ps(kvar_pdb, 5, struct.d_seq2pdb, struct.pchain)
	list_dist_neigh_3d = dist_neigh_3d[kvar]

	# extracting features
	codif = util.getMutCod(mut)
	all_profile = util.Unified_prof(kvar[1], prof, struct.seq, list_dist_neigh_3d)
	#all_data=[*codif,*all_profile,*np.zeros(600-len(all_profile))]
	return pd.DataFrame([*codif, *all_profile, *np.zeros(600 - len(all_profile))]).T

class ACDC3D:
	def __init__(self, weights=True):
		self.nn, _ = nn.build_acdc_3d(num_H=[32, 16], d=0.2, num_3d=25)
		if weights is False:
			pass
		elif weights is True:
			self.load_weights()
		else:
			self.load_weights(weights)

	def load_weights(self, path=None):
		if path is None:
			path = resource_filename('acdc_nn', 'weights/full_dataset_TL')
		self.nn.load_weights(path)

	def predict(self, wt_mut, wt_prof, wt_struct, mt_mut=None, mt_prof=None, mt_struct=None):
		direct_nn_input = prepare_3d_input(wt_mut, wt_prof, wt_struct)
		
		if mt_mut is None: # compute input features for the inverse mutation inverting the input features for the direct mutation
			inverse_nn_input = direct_nn_input.copy()
			inverse_nn_input.iloc[:, :20] = inverse_nn_input.iloc[:, :20].replace([1.0, -1.0],[-1.0, 1.0])
		else: # compute input features for the inverse mutation using the mutated protein
			inverse_nn_input = prepare_3d_input(mt_mut, mt_prof, mt_struct)
		
		# Reshaping input
		Xm_d, X1D_d, X3D_d = nn.mkInp_3d(np.asarray(direct_nn_input).astype(np.float32), 500)
		Xm_i, X1D_i, X3D_i = nn.mkInp_3d(np.asarray(inverse_nn_input).astype(np.float32), 500)

		prediction = self.nn.predict([X3D_d, X1D_d, Xm_d , X3D_i, X1D_i, Xm_i])
		return prediction[0][0][0]

	
def run_tests():
	# from acdc_nn.acdc_nn import *
	#p=load_prot('tests/profiles/2ocjA.prof.gz')
	p=load_prot('acdc-nn/tests/profiles/2ocjA.prof.gz')
	#lp = util.profile_context(p['pprofile'], 2, 3)
	#lp1 = util.prof('_3_', p['profile'], p['profile'].keys())
	print(p.keys())

	print('ACDC-NN Sequence tests')
	nn = ACDCSeq()
	nn.load_weights()
	print(nn.predict('A1C', p))

	print('ACDC-NN Structure tests')
	nn = ACDC3D()
	nn.load_weights()
	print('OK')
if __name__ == '__main__':
	run_tests()
