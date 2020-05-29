import pandas as pd
import os,sys,glob,psutil
from multiprocessing import Pool
RES = sys.argv[1]
#tf = sys.argv[3]
pdbid = RES.split('_')[0]
inseq = RES.split('_')[1]
hla = RES.split('_')[2]
print pdbid,inseq,hla
bf = []
pre_rmsd = []

def ff(x):
	df = pd.read_csv(x + '_rmsd.txt',sep='\t',header=None)
	ppp = df[0].str.split('.').str[0] + '.pdb'
	df[3] = x + '_' + df[0].str.split('.').str[0]
	pre_rmsd.append(df)
	for i in ppp:
		os.system('sed \'s/HETATM/ATOM  /\' ' + i.split('.')[0] + '_het.pdb > ../pdbs/' +  '_'.join(i.split('.')[:4]) + '.pdb')
	os.chdir('../')

def pred(x):
	os.system('python ~/peptide_docking/predictFluct.py ' + x + ' XRAY')
def b_part(x):
	nni = ['B_factor-P%d'% i for i in range(1,seqlen+1)]
	print nni
	for i in sorted(glob.glob('*.pdb.csv')):
		pdbn = pdbid + '_' + i.split('_')[0] + '_' + i.split('_')[6]
		df = pd.read_csv(i,sep=',')
		df_pep = df[df['Chain index'] == ' B']
		df_pep['new_index'] = nni
		df_pep_d = df_pep.drop(['Residue index','Chain index'],axis=1).set_index('new_index').T.rename(index={'Predicted MD fluctuation value' : pdbn})
		bf.append(df_pep_d)
	pd.concat(bf).reset_index().rename(columns={'index':'PDB'}).to_csv('total_bfact.txt',sep='\t',index=False)
	os.system('cp total_bfact.txt ../' + RES + '_energy_matrix/')
	os.chdir('../')
os.chdir('/home/user1/' + RES)
if not os.path.exists('pdbs'):
	os.mkdir('pdbs')
else: pass
if not os.path.exists(RES + '_energy_matrix/'):
	os.mkdir(RES + '_energy_matrix/')
else: pass

for nn in ['1ST','2ND']:
	os.chdir('PDB_' + nn)
	ff(nn)
os.chdir('pdbs')
pool = Pool(psutil.cpu_count()-4)#psutil.cpu_count()-4)
pool.map(pred,sorted(glob.glob('*.pdb')))
pool.close()
pool.join()
seqlen = len(inseq)
b_part('total')
os.chdir(RES + '_energy_matrix/')
dd = pd.concat(pre_rmsd)
dd.to_csv('colad_rmsd.txt',sep='\t',index=False)

df_ac = pd.read_csv('total_ac.txt',sep='\t')
df_nac = pd.read_csv('total_nac.txt',sep='\t')
df_sk = pd.read_csv('total_sk.txt',sep='\t')
df_hh = pd.read_csv('total_hh.txt',sep='\t')
df_bf = pd.read_csv('total_bfact.txt',sep='\t')
rmsd = pd.read_csv('colad_rmsd.txt',sep='\t')
rmsd[4] = pdbid + '_' + rmsd['3'].str.split('_').str[0] + '_' + rmsd['3'].str.split('_').str[7]
rmsd.columns = ['pdb','refer','rmsd','pre_PDB','PDB']
rmsd.drop(['pdb','pre_PDB'],axis=1).set_index(['PDB']).reset_index().to_csv('total_rmsd1.txt',sep='\t',index=False)
df_rmsd = pd.read_csv('total_rmsd1.txt',sep='\t').drop(['refer'],axis=1)

total_df = [df_rmsd,df_hh,df_nac,df_ac,df_bf]
df_final1 = reduce(lambda left,right:pd.merge(left,right,on=['PDB'],how='outer'),total_df)
df_final2 = pd.merge(df_final1.set_index('rmsd').reset_index(),df_ac).dropna()
df_final2['PDB1'] = df_final2['PDB'].str.split('_').str[0] + '_' + inseq + '_' + df_final2['PDB'].str.split('_').str[1] + '_' + df_final2['PDB'].str.split('_').str[2]
df_final3 = df_final2.drop(['PDB'],axis=1).rename(columns={'PDB1':'PDB'}).set_index('PDB').reset_index()
df_final3.to_csv(pdbid + '_' + inseq + '_' + hla +'_total.txt',sep='\t',index=False,na_rep='-')
os.chdir('/home/user1/')
