import gzip
import numpy as np
from Bio.PDB import *
from Bio.PDB.Polypeptide import three_to_one, is_aa

aa_modif={"TRN":"TRP", "CSO":"CYS","M3L":"LYS", "PCA":"GLU"}

def pdb2seq(pp):
    ''' pdb2seq(pp) takes a pdb_structure_chain 
    and return its sequence '''
    seq = [] # pp.get_sequence()
    reslist = []
    for ppc  in pp:
        reslist += [res for res in ppc]
        seq += [str(ppc.get_sequence())]
    return "".join(seq)

def map_pdb_pos(pp):
    ''' map_pdb_pos
    Returns two dicts seq2pdb[seq_pos], pdb2seq[pdb_pos]'''
    reslist = []
    for ppc  in pp:
        reslist += [res for res in ppc]
    seq2pdb = dict(zip( map(str,range(1,len(reslist)+1)), [str(r.get_id()[1])+r.get_id()[2].strip() for r in reslist]))
    pdb2seq = dict(zip( [str(r.get_id()[1])+r.get_id()[2].strip() for r in reslist], map(str,range(1,len(reslist)+1)) ))
    return seq2pdb, pdb2seq

def magic_open(path):
    return (gzip.open if path.endswith('.gz') else open)(path, 'rt')

def pdb2info(pdb_file, chain):
    ''' pdb2info(pdb_file) 
    Returns structure, polypeptide '''
    parser=PDBParser(QUIET=True)
    with magic_open(pdb_file) as f:
        structure = parser.get_structure('X', f)
    pchain=structure[0][chain]
    ppb=PPBuilder()
    pp = ppb.build_peptides(pchain, aa_only=False) #[0]
    return (structure, pchain, pdb2seq(pp), *map_pdb_pos(pp)) 

def mk_d_renumb_KseqVpdb(structure):
    d={}
    ppb=PPBuilder()
    s=""
    for pp in ppb.build_peptides(structure):
        s+=pp.get_sequence()
    i=1
    for residue in structure:
        t,n,lab=residue.get_id()
        if t == ' ':
            #d[(n,lab.strip())]=i
            d[str(i)]=str(n)+lab.strip()
            i+=1
    return d

def invert_keyval(din): 
    do={}
    for kin in din: do[din[kin]]=kin
    return do


def getMutCod(mut):
    code=list(np.zeros(20))
    d = {'A': 0, 'C': 1, 'D': 2,'E':3,'F':4,'G':5,'H':6,'I':7,'K':8,'L':9,
    	'M':10,'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,'Y':19}
    w=mut[0]
    m=mut[-1]
    pos_w=d[w]
    pos_m=d[m]
    cod_mut=code.copy()
    cod_mut[pos_w]=-1
    cod_mut[pos_m]=1
    return cod_mut


def getProfile(profile_path):
    prof = {}
    lines = magic_open(profile_path).readlines()
    # POS A C D E F G H I K L M N P Q R S T V W Y -
    #aa = lines[0].split()[1:-1]
    aa = lines[0].split()[1:]
    for line in lines[1:]:
         v = line.split()
         if aa[-1] == 'SEQ':
            v = v[:-1]
         prof[int(v[0])] = dict(zip(aa, map(float,v[1:])))
    return prof


def get_neighb(pchain,pos,radius_angst):
    # Next line: A for atoms 
    all_atom_list = Selection.unfold_entities(pchain, 'A')
    ns = NeighborSearch(all_atom_list)
    dv, dist={}, {} 
    dv[pos]=[]
    dist[pos] = []
    try: p,rlabel=int(pos),' '
    except: p,rlabel=int(pos[:-1]),pos[-1]
    for atomo in pchain[(' ',p,rlabel)]: # for atomo in chain[p]:
       rvicini = ns.search(atomo.get_coord(), radius_angst, "R")
       for aa in rvicini:
            l_dist = []
            if (aa.get_id()[0]==' ') and (aa.get_id()[1]!=p):
                for atoms in aa:
                    l_dist.append(atomo-atoms)
                if aa not in dv[pos]: 
                    dv[pos].append(aa)  # Returns object
                    dist[pos].append((aa, min(l_dist)))
                else:
                    for bb in dist[pos]:
                        if (bb[0] == aa) and (min(l_dist) < bb[1]):
                            dist[pos].remove(bb)
                            dist[pos].append((aa, min(l_dist)))
    return dv, dist
    
    

def get_neigh_pp(pchain,mut_pdb,Cangstroms):
    pdb_pos=mut_pdb[1]
    # Compute neighbor for mutated positions
    dneigh, dist =get_neighb(pchain,pdb_pos, Cangstroms)
    # Makes a dict: {mutatedpos_pdbnum : [ neighbor1, neighbor2, ... ], ...}
    dout={}
    dout[mut_pdb]=[]
    for aa in dist[pdb_pos]:
        pdbpos=(str(aa[0].get_id()[1])+aa[0].get_id()[2]).strip()
        dout[mut_pdb].append((three_to_one(aa[0].get_resname())+pdbpos, round(aa[1], 2)))
    return dout




def get_neigh_ps(mut_pdb,Cangstroms,d_seq2pdb_numb,pdbchain):
    '''From pdb to seq'''
    # Reads mutated positions for each pdb
    d_pdb=get_neigh_pp(pdbchain,mut_pdb,Cangstroms)
    dist_all={}
    d_pdb2seq_numb=invert_keyval(d_seq2pdb_numb)
    mut = list(d_pdb.keys())[0]
    mut_seqn=(mut[0],d_pdb2seq_numb[mut[1]],mut[2])
    dist_all[mut_seqn]=[]
    for wtpos in d_pdb[mut]:
        wt=wtpos[0][0]
        seqpos=d_pdb2seq_numb[wtpos[0][1:]]
        dist_all[mut_seqn].append((wt+seqpos, wtpos[1]))
    return dist_all
    


def Unified_prof(mut,profile,seq,l_dist_3d):
    pos_mut=int(mut)
    #sequence first
    prof_seq=profseq(pos_mut, profile, seq)
    #the rest of the 3d neigh        
    l_dist_3d.sort(key = lambda elem: (elem[1], elem[0]))
    prof_3d=prof3d(pos_mut,l_dist_3d, profile)
    return(prof_seq+prof_3d)

# context = 3 for seq, 2 for 3d
def profseq(pos_mut,profile, seq):
    prof_s=[]
    for j in range(-2, +3):
        if (pos_mut+j)>=1 and (pos_mut+j) <=len(seq):
            prof_s.append(list(profile[pos_mut+j].values())[:-1])
        else: prof_s.append(list(np.zeros(20)))
    unified_prof_s=[x for l in prof_s for x in l]
    return unified_prof_s

def prof(kprint,profile,seq):#seq       
    pos_mut=int(kprint[1:-1])
    profili=[]
    for i in range(-3,+4):
        if (pos_mut+i)>=1 and (pos_mut+i) <len(seq):
            profili.append(list(profile[pos_mut+i].values())[:-1])
        else:
            profili.append(list(np.zeros(20)))
    unified_prof=[]
    for j in profili:
        for i in j:
            unified_prof.append(i)
    return(unified_prof)

def prof3d(pos_mut,l_dist_3d, profile):
    profiles=[]
    for i in l_dist_3d:
        pos_res = int(i[0][1:])
        avoid_res=np.arange(-2,3)+pos_mut
        if pos_res not in avoid_res:  profiles.append(list(profile[pos_res].values())[:-1])       
    unified_prof3d=[x for l in profiles for x in l]
    return unified_prof3d

def bias(direct,inverse):
  bias=np.mean(direct+inverse)
  return (bias/2)

import pandas
from warnings import warn
# profile functions
aa1 = pandas.Index(list('ACDEFGHIKLMNPQRSTVWY')) # standard 20 amino acids
def load_profile(path_or_file):
	df = pandas.read_csv(path_or_file, sep=None, engine='python') # python engine needed for automatic detection of separator

	# check columns
	missing_aa = aa1.difference(df.columns)
	if len(missing_aa) > 0:
		raise ValueError('Missing AA columns in profile table: ' + ', '.join(missing_aa))
	unknown_cols = df.columns.difference(aa1).difference(['POS', 'SEQ', '-'])
	if len(unknown_cols) > 0:
		warn("Found unknown columns in profile table: " + ", ".join(unknown_cols))

	# check position column if any
	if 'POS' in df.columns:
		if 'POS' != df.columns[0]:
			warn("Found 'POS' column in profile but not in first position")
		if (df['POS'].astype(int) - df['POS'] != 0).any():
		#if df['POS'].values.kind not in ('i', 'u'):
			raise ValueError("Found 'POS' column in profile but it is not an integer")
		df = df.set_index('POS')
		
	if df.index[0] == 0: # 0-based, change into 1-based position
		warn("Found zero-based positions, adding 1")
		df.index = df.index + 1

	return df
			
def profile_context(profile, pos: int, context: int): # position is based on table index
	'''Extract local subset of profile.

	profile: protein profile as pandas DataFrame
	pos: protein position
	context: number of position to add around pos
	'''
	
	idx = list(range(pos - context, pos + context + 1))
	return profile.reindex(index=idx, columns=aa1, fill_value=0).values.flatten()
