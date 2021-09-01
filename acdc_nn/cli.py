from acdc_nn import acdc_nn
from acdc_nn import util
import ddgun

import click
from warnings import warn
import functools

class Substitution(click.ParamType):
	'''Click parameter class for substitutions'''
	name = 'amino acid substitution'
	def convert(self, value, param, ctx):
		if isinstance(value, ddgun.Substitution):
			return value
		try:
			return ddgun.Substitution.parse(value)
		except Exception as e:
			self.fail(f"{value!r} is not a valid {self.name}", param, ctx)

help_notes = '''Notes:
Mutations are written as XNY, meaning that the residue X at position N changes to Y. X and Y are given as a one letter amino acid code and
N is 1-based and refers to the the PDB numbering of the relevant chain, and not the position on the sequence.
PDB and profile files will be automatically decompressed (by gzip) if the paths end with ".gz".
'''

@click.group(epilog=help_notes)
def cli():
	pass

@cli.command(epilog=help_notes)  # TODO add option for the weights
@click.argument("sub", type=Substitution())
@click.argument("profile", type=click.Path(exists=True, readable=True))  # FIXME use File object
def seq(sub, profile):
	'''Predict DDG of SUB from the protein PROFILE.

\b
SUB is an amino acid substitution (e.g. Q339P).
PROFILE is the path to a protein profile file.

Uses a trained ACDC-NN Seq that does not require protein structural information.'''
	wt_prof = ddgun.Profile(profile)
	net = acdc_nn.ACDCSeq()
	ddg = net.predict(sub, wt_prof)
	click.echo(ddg)

@cli.command(epilog=help_notes)
@click.argument("sub", type=Substitution())
@click.argument("profile", type=click.Path(exists=True, readable=True))
@click.argument("pdb", type=click.Path(exists=True, readable=True))
@click.argument("chain")
#@click.option('--inverse', type=(Substitution(), click.Path(exists=True, readable=True), click.Path(exists=True, readable=True), str), help=')
def struct(sub, profile, pdb, chain):  #FIXME add inverse mut 
	'''Predict DDG of SUB from the protein PROFILE and PDB structure.

\b
SUB is an amino acid substitution (e.g. Q339P).
PROFILE is the path to a protein profile file.
PDB is the path to the protein structure in PDB file.
CHAIN is the PDB chain to be used.

Uses a trained ACDC-NN that requires protein structural information.'''
	wt_prof = util.getProfile(profile)  #FIXME use Profile
	wt_struct = acdc_nn.Structure(pdb, chain)
	net = acdc_nn.ACDC3D()
	ddg = net.predict(str(sub), wt_prof, wt_struct)
	click.echo(ddg)

@cli.command(epilog=help_notes)
@click.argument("sub", type=Substitution())
@click.argument("profile", type=click.Path(exists=True, readable=True))
@click.argument("pdb", type=click.Path(exists=True, readable=True))
@click.argument("chain")
@click.argument("isub", type=Substitution())
@click.argument("iprofile", type=click.Path(exists=True, readable=True))
@click.argument("ipdb", type=click.Path(exists=True, readable=True))
@click.argument("ichain")
#@click.option('--inverse', type=(Substitution(), click.Path(exists=True, readable=True), click.Path(exists=True, readable=True), str), help=')
def istruct(sub, profile, pdb, chain, isub, iprofile, ipdb, ichain):
	'''Predict DDG using both the wild-type and mutated protein structures.

\b
SUB is an amino acid substitution (e.g. Q339P).
PROFILE is the path to a protein profile file.
PDB is the path to the protein structure in PDB file.
CHAIN is the PDB chain to be used.
ISUB, IPROFILE, IPDB and ICHAIN are the same for the mutated protein.

Uses a trained ACDC-NN that requires protein structural information.'''
	wt_prof = util.getProfile(profile)  #FIXME use Profile
	wt_struct = acdc_nn.Structure(pdb, chain)
	mt_prof = util.getProfile(iprofile)  #FIXME use Profile
	mt_struct = acdc_nn.Structure(ipdb, ichain)
	net = acdc_nn.ACDC3D()
	ddg = net.predict(str(sub), wt_prof, wt_struct, str(isub), mt_prof, mt_struct)
	click.echo(ddg)


# caching for functions
@functools.lru_cache(10)
def load_nn(seq):
	return acdc_nn.ACDCSeq() if seq else acdc_nn.ACDC3D()

@functools.lru_cache(100)
def load_prot_seq(profile):
	return ddgun.Profile(profile)

@functools.lru_cache(100)
def load_prot_3d(profile, pdb, chain):
	return util.getProfile(profile), acdc_nn.Structure(pdb, chain) #FIXME use Profile


@cli.command(epilog=help_notes)
@click.argument("subs", type=click.File())
def batch(subs):
	'''Predict DDG of SUBS using available information.

SUBS is a table containing one amino acid substitution per row and paths to protein profiles and optionally protein structure.

\b
Each row can have 2, 4 or 8 fields of tab separated values that are interpreted with the following schema. 
# \tPredictor   \tFields
2 \tACDC-NN Seq \tSUB PROFILE
4 \tACDC-NN     \tSUB PROFILE PDB CHAIN
8 \tACDC-NN     \tWT-SUB WT-PROFILE WT-PDB WT-CHAIN MT-SUB MT-PROFILE MT-PDB MT-CHAIN
For rows with 2 fields, that is, without structural information, the sequence-based ACDC-NN Seq predictor is used. For rows with 4 or 8 fields, the structure-based ACDC-NN is used. Outputs one DDG value for each table row.
'''
	for row, line in enumerate(subs, 1):
		fields = line.rstrip('\n\r').split('\t')
		
		# detect how many fields are available
		null_fields = [f in {'', '.', 'NA', 'na'} for f in fields]
		for i, nf in enumerate(null_fields):
			if nf:
				n = i
				break
		else:
			n = len(fields)

		if n == 2: # only profile
			seq = True
			pargs = ddgun.Substitution.parse(fields[0]), load_prot_seq(fields[1])
		elif n == 4 or n == 8: # also structure
			seq = False
			pargs = (
					str(ddgun.Substitution.parse(fields[0])),
					*load_prot_3d(*fields[1:4]))
			if n == 8: # also inverse substitution
				pargs = (
						*pargs,
						str(ddgun.Substitution.parse(fields[4])), 
						*load_prot_3d(*fields[5:8])) 
		else:
			raise ValueError(f"found {n} fields at line {row}: fields must be 2, 4 or 8")

		# check that there
		for i, nf in enumerate(null_fields[n:], n):
			if not nf:
				warn(f"found value at column {i} after missing value at column {n}, at line {row}")

		ddg = load_nn(seq=seq).predict(*pargs)

		click.echo(ddg)

