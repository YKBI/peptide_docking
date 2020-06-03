import pandas as pd
import os,sys,glob,psutil,shutil
from multiprocessing import Pool,Process,Queue
#lists
rec_n = []
lig_n = []
ac = []
nac = []
sk = []
hh = []
rmsd = []
bf = []
rfeats = []
rfeats1 = []
pep_feat = []
#dictionaries
rec_dic = {}
lig_dic = {}
cc_dic = {}
pep_dic = {}
def reading_DB(x):
    hladb = pd.read_csv(x,sep='\t')
    indb = hladb[hladb['PDB ID'] == pdbid.upper()].dropna()
    print indb
    o_seq = indb['Peptide Sequence'].values[0]
    hlaclass = indb['HLA CLASS'].values[0]
    hlatype = indb['HLA TYPE'].values[0]

    seqlen = len(o_seq)
    hla = hlaclass.split('-')[1] + ''.join(hlatype.split(':'))
    return hla,seqlen,o_seq
def reT(x):
    nam = x.split('.')
    nam1 = '.'.join([nam[0], nam[3],'red',nam[2]])
    Olines = []
    recl = []
    pepl = []
    os.system('csplit -f \'%s_\' %s \'/TER/\' --quiet'%(x,x))
    with open(x + '_00', 'r') as F:
        for line in F.readlines():
            tline = line.strip()
            if tline.startswith('ATOM'):
                if tline[-1] != 'H':
                    recl.append(tline[:21] + 'A ' + tline[23:])
    with open(x + '_01','r') as F:
        for line in F.readlines():
            tline = line.strip()
            if tline.startswith('ATOM'):
                if tline[-1] != 'H':
                    pepl.append(tline[:21] + 'B' + str(pep_dic[tline[23:26].strip()]).rjust(4) + tline[26:])
    with open(nam1,'w') as W:
        for line in recl:
            W.write(line + '\n')
        W.write('TER\n')
        for line in pepl:
            W.write(line + '\n')

def revise_origin(x):
    rec_ser = []
    lig_ser = []
    rnn = 1
    nnn = 1
    with open(RECLIB + x + '_rec.pdb','r') as R: #receptor reading
        for line in R.readlines():
            rec_n.append(line[7:13].strip())
            if line[13:17].strip() == 'OXT' or line[13:17].strip() == '' :
                pass
            else:
                if rnn == int(line[22:26].strip()):
                    rec_ser.append(line[13:17].strip())
                    rec_dic[rnn] = rec_ser
                else:
                    rec_dic[rnn] = rec_ser
                    rnn = int(line[22:26].strip())
                    rec_ser = []
                    rec_ser.append(line[13:17].strip())
    with open(x + '_' + inseq + '_pep.pdb','r') as P: #peptide reading
        for line in P.readlines():
            lig_n.append(line[7:13].strip())
            if line[13:17].strip() == 'OXT' or line[13:17].strip() == '':
                pass
            else:
                if nnn == int(line[22:26].strip()):
                    lig_ser.append(line[13:17].strip())
                    lig_dic[nnn] = lig_ser
                else:
                    lig_dic[nnn] = lig_ser
                    nnn = int(line[22:26].strip())
                    lig_ser = []
                    lig_ser.append(line[13:17].strip())

def atom_revise(x):
    nam = x.split('.')
    revnam = '.'.join([nam[0],nam[1],'rev',nam[3]])
    rnn = 1
    nnn = 1
    olig_ser = {}
    orec_ser = {}
    lig_dic2 = {}
    rec_dic2 = {}
    lig_lines = []
    out_lines = []
    with open(x,'r') as F:
        for line in F.readlines():
            if line[21:23].strip() == 'A':
                if rnn == int(line[23:26].strip()):
                    orec_ser[line[13:17].strip()] = line
                    rec_dic2[rnn] = orec_ser
                else:
                    rec_dic2[rnn] = orec_ser
                    rnn = int(line[23:26].strip())
                    orec_ser = {}
                    orec_ser[line[13:17].strip()] = line
            elif line[21:23].strip() == 'B':
                if nnn == int(line[23:26].strip()):
                    olig_ser[line[13:17].strip()] = line
                    lig_dic2[nnn] = olig_ser
                else:
                    lig_dic2[nnn] = olig_ser
                    nnn = int(line[23:26].strip())
                    olig_ser = {}
                    olig_ser[line[13:17].strip()] = line
    feature = list(set(pep_dic.values()))
    for i in feature:
        for j in lig_dic[i]:
            lig_lines.append(lig_dic2[i][j].strip())
    for i in rec_dic:
        for j in rec_dic[i]:
            out_lines.append(rec_dic2[i][j].strip())
    with open(revnam,'w') as W:
        for i,j in zip(out_lines,rec_n):
            W.write(i[:7] + j.rjust(4) + '  ' + i[13:] + '\n')
        W.write('TER\n')
        for i,j in zip(lig_lines,lig_n):
            W.write(i[:7] + j.rjust(4) + '  ' + i[13:] + '\n')
        W.write('END\n')

def sheba(x):
    nam = x.split('.')
    num = str(nam[1])
    nana = pdbid + '_' + num
    namtrf = '.'.join([nam[0], num, nam[2]])
    nnn = '.'.join([nam[0], nam[1], nam[2]])
    outform1 = pdbid + '_' + num + '.pdb'
    outform2 = pdbid + '_' + num + '_new.pdb'
    outform3 = pdbid + '_' + num + '_het.pdb'
    os.system('clean_pdb.py %s A'%x)
    os.system('clean_pdb.py %s B'%x)
    os.system('cat ' + nnn + '_A.pdb ' + nnn + '_B.pdb > ' + nnn + '_AB.pdb')
    os.system('sheba_01 -x ' + refer + ' ' + nnn + '_AB.pdb')
    os.system('sheba_01 -t ' + namtrf + '_AB.trf ' + nnn + '_AB.pdb')
    shutil.move(nnn + '_AB.pdb.pdb',outform1)
    os.system('clean_pdb.py ' + outform1 + ' B')
    os.system('cat ' + RECLIB + '%s_rec.pdb %s_%s_B.pdb > %s'%(pdbid,pdbid,num,outform2))
    os.system('sed \'s/ATOM  /HETATM/\' %s_%s_B.pdb > %s_%s_02'%(pdbid,num,pdbid,num))
    os.system('cat ' + RECLIB + '%s_rec.pdb %s_%s_02 > %s'%(pdbid,pdbid,num,outform3))
    '''
    if not os.path.exists(nana + '_a.out'):
        os.system('enva.v2 -a ' + outform2 + ' > ' + nana + '_a.out')
    if not os.path.exists(nana + '_b.out'):
        os.system('enva.v2 -b ' + outform3 + ' > ' + nana + '_b.out')
    if not os.path.exists(nana + '_m.out'):
        os.system('enva.v2 -m ' + outform2 + ' ' + pdbid +  'B > ' + nana + '_m.out' )
    os.system('enva.v2 -e ' + outform2 + ' B')'''

def enva(x):
    nam = x.split('.')
    num = str(nam[1])
    nana = pdbid + '_' + num
    namtrf = '.'.join([nam[0], num, nam[2]])
    nnn = '.'.join([nam[0], nam[1], nam[2]])
    outform1 = pdbid + '_' + num + '.pdb'
    outform2 = pdbid + '_' + num + '_new.pdb'
    outform3 = pdbid + '_' + num + '_het.pdb'
    if not os.path.exists(nana + '_a.out'):
        os.system('enva.v2 -a ' + outform2 + ' > ' + nana + '_a.out')
    if not os.path.exists(nana + '_b.out'):
        os.system('enva.v2 -b ' + outform3 + ' > ' + nana + '_b.out')
    if not os.path.exists(nana + '_m.out'):
        os.system('enva.v2 -m ' + outform2 + ' ' + pdbid +  'B > ' + nana + '_m.out' )
    os.system('enva.v2 -e ' + outform2 + ' B')
    if not os.path.exists(nana + '_new.pdb.csv'):
        os.system('python ../../../../peptide_docking/predictFluct.py ' + outform2 + ' XRAY')

def rmsd_part(x):
    for i in sorted(glob.glob(pdbid + '_*[0-9]_B.pdb')):
        #os.system('csplit -f \'' + i + '_\'' + i + ' \'/TER/\' --quiet')
        os.system('rmsd_total ' + '../../' + pdbid + '_' + inseq + '_pep.pdb ' + i + ' bb >> ' + x + '_rmsd.txt')
def nac_sk_part(x):
    pdbl = []
    mat_lst = []
    newl = []
    tttt = []
    ave_acc = []
    tot_lig_acc = 0
    ave_lig_acc = 0
    tot_lig_num = 0
    for f, g in zip(sorted(glob.glob('*a.out')),sorted(glob.glob('*m.out'))):
        sidx = []
        pdb = '_'.join([pdbid, str(x), f.split('.')[0].split('_')[1]])
        pdbl.append(pdb)
        adf = pd.read_csv(f, sep='\s+', skiprows=[0, 0], header=None)
        mdf = pd.read_csv(g, sep='\s+', header=None)
        mdf_filter = mdf[mdf[17] == 1].iloc[:, 9:16]
        for i in mdf_filter[9]:
            tot_lig_acc = tot_lig_acc + float(i)
            tot_lig_num += 1
            if tot_lig_num > 0:
                ave_lig_acc = tot_lig_acc / float(tot_lig_num)
            else:
                ave_lig_acc = 0
        ave_acc.append(ave_lig_acc)

        new_acc = {}
        num = 0
        for i in rfeats:
            new_acc['AA_' + cc_dic[str(i)]] = sum([adf[adf[1] == int(i.split('_')[0])][adf[2] == i.split('_')[2]][adf[3] == i.split('_')[1]][11].values,1])
        for i in mdf_filter.T:
            tidx = str(mdf_filter.T[i][13]) + '_' + str(mdf_filter.T[i][14]) + '_' + str(mdf_filter.T[i][15])
            if tidx in rfeats:
                tttt.append(tidx)
                zdx = rfeats.index(tidx) + 1
                sidx.append(zdx)
                num += 1
        newl.append(pd.DataFrame(new_acc))
        if not os.path.exists(pdb.split('.')[0] + '.ser'):
            with open(pdb.split('.')[0] + '.ser', 'w') as W:
                for sid in list(set(sidx)):
                    W.write(str(sid) + '\n')
        os.system('iskew ' + pdb.split('.')[0] + '.ser >> total_' + x + '_sk.txt')
        ratio = float(num) / float(len(rfeats))
        mat_lst.append(ratio)
    new_ac = pd.concat(newl)
    new_ac['total_rec_acc'] = new_ac.sum(axis=1)
    new_ac['ave_lig_acc'] = ave_acc
    new_ac['%Match'] = mat_lst
    new_ac['PDB'] = pdbl
    new_ac.set_index('PDB').reset_index()
    new_ac.to_csv(x + '_nac.txt', sep='\t', index=False)
def hh_part(x):
    ff_list = []
    hh_list = []
    nn_list = []

    for f in sorted(glob.glob('*b.out')):
        with open(f, 'r') as F:
            n = 0
            for line in F.readlines():
                tt = line[7:13].strip() + '_' + line[17:21].strip() + '_' + line[13:17].strip()
                if tt in rfeats:
                    n += 1
            nn_list.append(n)
        with open(f, 'r') as F:
            hh_list.append(len(F.readlines()))
        ff_list.append(pdbid + '_' + str(x) + '_' + f.split('.')[0].split('_')[1])
    with open(x + '_hh.txt', 'w') as W:
        W.write('PDB\tN.of.BB_full\tN.of.BB_feat\n')
        for i, j, k in zip(ff_list, hh_list, nn_list):
            W.write(str(i) + '\t' + str(j) + '\t' + str(k) + '\n')
def ac_part(x):
    pdbl = []
    tt_acc = []
    tt_phi = []
    tt_psi = []
    acc_columns = ['P%d' % i for i in range(1, seqlen + 1)]
    phi_columns = ['PHI%d' % i for i in range(1, seqlen + 1)]
    psi_columns = ['PSI%d' % i for i in range(1, seqlen + 1)]
    for f in sorted(glob.glob('*.env')):
        pdb = '_'.join([pdbid, str(x), f.split('_')[0].split('B')[1]])
        pdbl.append(pdb)
        with open(f, 'r') as F:
            acc = []
            psi = []
            phi = []
            for line in F.readlines():
                if line.find('chain') and line.startswith('ATOM'):
                    envs = line.split()[9:]
                    acc.append(envs[2])
                    phi.append(envs[4])
                    psi.append(envs[5])
            tt_acc.append(pd.DataFrame(acc, index=acc_columns).T)
            tt_phi.append(pd.DataFrame(phi, index=phi_columns).T)
            tt_psi.append(pd.DataFrame(psi, index=psi_columns).T)
    pdbdf = pd.DataFrame(pdbl, columns=['PDB'])
    acc_df = pd.concat(tt_acc, ignore_index=True).reindex()
    phi_df = pd.concat(tt_phi, ignore_index=True).reindex()
    psi_df = pd.concat(tt_psi, ignore_index=True).reindex()
    act_df = pd.concat([pdbdf, acc_df, phi_df, psi_df], axis=1)
    act_df.to_csv(x + '_ac.txt', sep='\t', index=False)

def b_part(x):
    nni = ['B_factor-P%d'%i for i in range(1,seqlen+1)]
    for i in sorted(glob.glob('*.pdb.csv')):
        pdbn = pdbid + '_' + x + '_' + i.split('_')[1]
        df = pd.read_csv(i,sep=',')
        df_pep = df[df['Chain index'] == ' B']
        df_pep['new_index'] = nni
        df_pep_d = df_pep.drop(['Residue index','Chain index'],axis=1).set_index('new_index').T.rename(index={'Predicted MD fluctuation value':pdbn})
        bf.append(df_pep_d)
    pd.concat(bf).reset_index().rename(columns={'index':'pdb'}).to_csv('total_bfact.txt',sep='\t',index=False)
pdbid = sys.argv[1]
infile = sys.argv[3]
inseq = sys.argv[2]
(hla,seqlen,o_seq) = reading_DB(infile)
with open('feat.list','r') as F:
    for line in F.readlines():
        pep_feat.append(line.strip())
for i,j in zip(pep_feat,range(1,len(pep_feat)+1)):
    pep_dic[i] = j

os.environ['rosetta'] = '/awork08/93_hong/rosetta_src_2019.14.60699_bundle/'
os.environ['neogear'] = '/awork06-1/neoscan_gear/'
gear = os.environ['neogear']
rosetta = os.environ['rosetta']
os.environ['PATH'] += ':' + os.environ['PATH'] + ':' + rosetta + 'main/source/bin/:' + rosetta + 'tools/protein_tools/scripts/'
os.environ['PATH'] += ':' + os.environ['PATH'] + ':' + gear
os.environ['neogear'] = '/awork06-1/neoscan_gear/'
gear = os.environ['neogear']
os.environ['PATH'] += ':' + os.environ['PATH'] + ':' + gear
PDBLIB = '/awork06-1/YKLee/pdpdb/Neoscan_V2/Native/'
RECLIB = '/awork06-1/YKLee/pdpdb/Neoscan_V2/v2/revise_pdb/receptor/' + hla + '/'
PEPLIB = '/awork06-1/YKLee/pdpdb/Neoscan_V2/v2/revise_pdb/peptide/' + hla + '/'
ROSETTA_DB = rosetta + 'main/database'
'''
PDBLIB = '/awork06-1/YKLee/final_class2/' + ffreq + '/' + hla + '/minimize/'
RECLIB = '/awork06-1/YKLee/final_class2/' + ffreq + '/' + hla + '/minimize/'
PEPLIB = '/awork06-1/YKLee/final_class2/' + ffreq + '/' + hla + '/minimize/'
ROSETTA = '/awork06-1/rosetta_src_2019.40.60963_bundle/tools/protein_tools/scripts/'
'''
refer = PDBLIB + '/' + hla + '_native/'+ pdbid + '_native.pdb'
ncpu = str(psutil.cpu_count()-1)
feats = pd.read_csv(PDBLIB  + '/' + hla + '_native/' + pdbid + '.out',sep='\s+',header=None)
feats_filt = feats[feats[17] == 1].iloc[:,13:16]
for i in feats_filt.values.tolist():
    rfeats.append(str(i[0]) + '_' + '_'.join(i[1:]))
feats1_filt = feats[feats[17] == 1].iloc[:,12:16]
for i in feats1_filt.values.tolist():
    rfeats1.append(str(i[0]) + '_' + '_'.join(i[2:]))
for i,j in zip(rfeats,rfeats1):
    cc_dic[i] = j
revise_origin(pdbid)
os.chdir('traj_1/')
if not os.path.exists(pdbid + '_' + hla + '_energy_matrix'):
    os.mkdir(pdbid + '_' + hla + '_energy_matrix')
else:
    pass

if not os.path.exists('pdbs'):
    os.mkdir('pdbs')
else:
    pass
os.chdir('pdb_from_prod/')
list1 = sorted(glob.glob('*.pdb.*'))[:100]
pool = Pool(int(ncpu))
pool.map(reT,list1)
pool.close()
pool.join()

list22 = sorted(glob.glob('*.red.*'))
pool1 = Pool(int(ncpu))
pool1.map(atom_revise,list22)
pool1.close()
pool1.join()

list2 = sorted(glob.glob('*.rev.*'))
pool2 = Pool(int(ncpu))
pool2.map(sheba,list2)
pool2.close()
pool2.join()
proc_one = Process(target=rmsd_part,args=('traj_1',))
proc_one.start()
proc_one.join()
pre_rmsd = pd.read_csv('traj_1_rmsd.txt',sep='\t',header=None)
ppp = pre_rmsd[pre_rmsd[2]<=2.0][0].str.split('_').str[1]
os.system('cp traj_1_rmsd.txt ../pdbs')
for i in ppp:
    os.system('cp ' + pdbid + '_' + i + '_new.pdb ../pdbs')
    os.system('cp ' + pdbid + '_' + i + '_het.pdb ../pdbs')
    os.system('cp ' + pdbid + '_' + inseq + '.' + i + '.rev.pdb ../pdbs')
os.chdir('../pdbs')
list3 = sorted(glob.glob('*.rev.*'))
pool3 = Pool(int(ncpu))
pool3.map(enva,list2)
proc_two = Process(target=nac_sk_part,args=('traj_1',))
proc_three = Process(target=hh_part,args=('traj_1',))
proc_four = Process(target=ac_part,args=('traj_1',))
proc_five = Process(target=b_part,args=('traj_1',))

proc_two.start()
proc_three.start()
proc_four.start()
proc_five.start()

proc_two.join()
proc_three.join()
proc_four.join()
proc_five.join()

print os.getcwd()
os.system('cp *.txt ../' + pdbid + '_' + hla + '_energy_matrix/')
os.chdir('../' + pdbid + '_' + hla + '_energy_matrix/')

df_bf = pd.read_csv('total_bfact.txt',sep='\t')
df_ac = pd.read_csv('traj_1_ac.txt',sep ='\t')
df_rmsd = pd.read_csv('traj_1_rmsd.txt',sep='\t',header=None)
df_nac = pd.read_csv('traj_1_nac.txt',sep='\t')
df_sk = pd.read_csv('total_traj_1_sk.txt',sep='\t',header=None)
df_hh = pd.read_csv('traj_1_hh.txt',sep='\t')
df_sk.columns = ['PDB','skewness','Class','Decision']
df_rmsd.columns = ['pdb','reference','rmsd']

total_df = [df_hh,df_nac,df_sk]
df_final1 = reduce(lambda left,right: pd.merge(left,right,on=['PDB'],how='outer'),total_df)
df_final1['rmsd'] = df_rmsd.reset_index()['rmsd']
df_final2 = pd.merge(df_final1.set_index('rmsd').reset_index().set_index('PDB').reset_index(),df_ac)
df_final2.to_csv(pdbid + '_' + hla + '_total.txt',sep='\t',index=False,na_rep='-')
os.chdir('../')
