import pandas as pd
import os,sys,glob,psutil

def reading_DB(x):
    hladb = pd.read_csv(x,sep='\t')
    indb = hladb[hladb['PDB ID'] == pdbid].dropna()
    hlaclass = indb['HLA CLASS']
    hlatype = indb['HLA TYPE']
    o_seq = indb['Peptide Sequence']
    seqlen = len(o_seq)

def mutate(x,y):
    mtd_lines = []
    os.environ['MMTSB'] = '/lwork01/mmtsb'
    mmtsb = os.environ['MMTSB']
    os.environ['PATH'] +=  os.environ['PATH'] + ':' + mmtsb + '/perl:' + mmtsb + '/bin'
    os.system('mutate.pl -seq 1:' + inseq + ' ' + pdbid + '_B.pdb > ' + inseq + '_mtd.pdb')
    with open(inseq + '_mtd.pdb','r') as MP:
        for line in MP.readlines():
            if not line.startswith('TER\n' or 'END\n') :
                mtd_lines.append(line[:21] + 'B ' + line[23:56] + '1' + line[57:])
            else:
                mtd_lines.append(line)
    del mtd_lines[-1]
    with open(pdbid + '_' + inseq + '_pep.pdb','w') as W:
        W.write('TER\n')
        for line in mtd_lines:
            W.write(line)
    os.system('cat ' + pdbid + '_A.pdb ' + pdbid + '_' + inseq + '_pep.pdb > ' + pdbid + '_' + inseq + '_inp.pdb')

def write_exec(pdb,seq):
    prepack_lines=['-s ' + pdb + '_' + seq + '_inp.pdb',\
                   '-out:path:all ' + RES,\
                   '-out:prefix PPK_',\
                   '-out:file:scorefile score_ppk_' + pdb + '_' + seq + '.sc',\
                   '-ex1',\
                   '-ex2aro',\
                   '-user_input_sc',\
                   '-flexpep_prepack',\
                   '-nstruct 1']
    firststep_lines =['-s ' + RES + '/PPK_' + pdb + '_' + seq + '_inp_0001.pdb',\
                  '-native ' + ROOT + 'Native/' + HLAclaty + '_native/' + pdb + '_native.pdb',\
                  '-out:path:all ' + RES + '/PDB_1ST',\
                  '-out:file:scorefile score_1st_' + pdb + '_' + seq + '.sc',\
                  'nstruct ' + nout,\
                  '-flexPepDocking:pep_refine',\
                  '-flexPepDocking:flexpep_score_only',\
                  '-ex1',\
                  '-ex2aro']
    secondstep_lines=['-s ' + RES + '/PPK_' + pdb + '_' + seq + '_inp_0001.pdb',\
                  '-native ' + ROOT + 'Native/' + HLAclaty + '_native/' + pdb + '_native.pdb',\
                  '-out:path:all ' + RES + '/PDB_2ND',\
                  '-out:file:scorefile score_2nd_' + pdb + '_' + seq + '.sc',\
                  'nstruct ' + nout,\
                  '-flexPepDocking:pep_refine',\
                  '-flexPepDocking:flexpep_score_only',\
                  '-ex1',\
                  '-ex2aro']
    with open('prepack_flags','w') as W:
        for line in prepack_lines:
            W.write(line + '\n')
    with open('run_flag1','w') as W:
        for line in firststep_lines:
            W.write(line + '\n')
    with open('run_flag2','w') as W:
        for line in secondstep_lines:
            W.write(line + '\n')
def prepare_input(x):
    if inseq != o_seq:
        mutate(x,inseq)
        write_exec(pdbid,inseq)
    else:
        os.system('cat ' + x + '_A.pdb ' + x + '_B.pdb > ' + x + '_' + o_seq + '_inp.pdb')
        write_exec(x,o_seq)
pdbid = sys.argv[1]
inseq = sys.argv[2]
nout = sys.argv[3]
ncpu = str(psutil.cpu_count()-1)
if int(ncpu) >= 17 :
    os.environ['rosetta'] = '/lwork01/rosetta_src_2019.40.60963_bundle/main/source/bin/'
    os.environ['rosetta_clean'] = '/lwork01/rosetta_src_2019.40.60963_bundle/tools/protein_tools/scripts/'
    os.environ['neogear'] = '/lwork01/neoscan_gear/'
    gear = os.environ['neogear']
    rosetta = os.environ['rosetta']
    rosetta_clean = os.environ['rosetta_clean']
    os.environ['PATH'] += ':' + os.environ['PATH'] + rosetta + ':' + os.environ['PATH']
    os.environ['PATH'] += ':' + os.environ['PATH'] + rosetta_clean + ':' + os.envirion['PATH']
    os.environ['PATH'] += ':' + os.environ['PATH'] + gear + ':' + os.environ['PATH']
else:
    os.environ['rosetta'] = '/awork06-1/rosetta_src_2019.40.60963_bundle/main/source/bin/'
    os.environ['rosetta_clean'] = '/awork06-1/rosetta_src_2019.40.60963_bundle/tools/protein_tools/scripts/'
    os.environ['neogear'] = '/awork06-1/neoscan_gear/'
    gear = os.environ['neogear']
    rosetta = os.environ['rosetta']
    rosetta_clean = os.environ['rosetta_clean']
    os.environ['PATH'] += ':' + os.environ['PATH'] + rosetta + ':' + os.environ['PATH']
    os.environ['PATH'] += ':' + os.environ['PATH'] + rosetta_clean + ':' + os.envirion['PATH']
    os.environ['PATH'] += ':' + os.environ['PATH'] + gear + ':' + os.environ['PATH']


#rosetta clean py, divide by chain
os.system('clean_pdb.py ' + pdbid + '.pdb A')
os.system('clean_pdb.py ' + pdbid + '.pdb B')
prepare_input(pdbid)
os.system('FlexPepDocking.linuxgccrelease -database ' + ROSETTA_DB + ' @' + RES + '/preapack_flags > ' + RES + '/prepack.log')
for i,j in zip(['run_flag1','run_flag2'],['run1.log','run2.log']):
    os.system('mpirun -np ' + ncpu + ' FlexPepDocking.mpi.linuxgccrelease -database ' + ROSETTA_DB + ' @' + RES + '/' + i + '> ' RES + j)

