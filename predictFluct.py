#!/usr/bin/env python
# Michal Jamroz, 2011, 2015, jamroz@chem.uw.edu.pl, public domain.
from math import sqrt,pi
from os import remove, path
from re import search
from sys import argv
from subprocess import Popen


# SET THESE PARAMETERS :
svm_scale="/home/user1/peptide_docking/flexPred/libsvm-3.24/svm-scale"
svm_predict="/home/user1/peptide_docking/flexPred/libsvm-3.24/svm-predict"

#XRAY
svr_modelXRAY_path="/home/user1/peptide_docking/flexPred/Set15SVRModel"
svr_modelXRAY_range_path="/home/user1/peptide_docking/flexPred/Set15data.range"
#NMR
svr_modelNMR_path="/home/user1/peptide_docking/flexPred/Set16SVRModel"
svr_modelNMR_range_path="/home/user1/peptide_docking/flexPred/Set16data.range"

class Vector:
   def __init__(self,x,y,z,bf,idx, name):
      self.x = x
      self.y = y
      self.z = z
      self.bfactor = bf
      self.index = idx
      self.name = name
   def distance(self, v2):
      d = sqrt( (v2.x-self.x)**2 + (v2.y-self.y)**2 + (v2.z-self.z)**2)
      return d

class PDBRead:
  def __init__(self,filename):
      self.chains = {}
      self.chain_list = []
      self.parse(filename)
  def parse(self,filename):
      with open(filename) as f:
        for l in f.readlines():
            if (search("^ATOM.........CA.[A ]",l)):
                x1 = float(l[29:38])
                y1 = float(l[38:46])
                z1 = float(l[46:54])
                bf = float(l[60:66])
                idx = int(l[22:26])
                name= l[17:20]
                ch = l[21]
                v = Vector(x1,y1,z1,bf,idx,name)
                if ch in self.chains:
                    self.chains[ch].append(v)
                else:
                    self.chain_list.append(ch)
                    self.chains[ch] = [v]

            if search("^ENDMDL", l): # only first model
                break

  def chain(self, chain_id):
      if chain_id in self.chains:
          return self.chains[chain]
      else:
          return None
  def chainIndices(self):
      return self.chain_list
  def chainsNum(self):
      return len(self.chain_list)
  def allAtoms(self):
      return [i for sub in self.chains.values() for i in sub]

class Contacts:
    def __init__(self, reference_structure, cutoff=6):
        self.cutoff = cutoff
        self.ref = reference_structure
    def getContacts(self, atoms):
      length = len(atoms)
      cont = [-1 for i in range(length)]
      for i in range(length):
         for j in range(len(self.ref)):
             if atoms[i].distance(self.ref[j])  <= self.cutoff: 
                 cont[i] += 1
      return cont


def predictionForChain(chain_id, pdb):
    global svm_scale, svm_predict, svr_modelXRAY_path, svr_modelXRAY_range_path, svr_modelNMR_path, svr_modelNMR_range_path
    atoms = pdb.chain(chain_id)
    c16 = Contacts(pdb.allAtoms(),16.0).getContacts(atoms)
    c15 = Contacts(pdb.allAtoms(),15.0).getContacts(atoms)
    c18 = Contacts(pdb.allAtoms(),18.0).getContacts(atoms)
    c12 = Contacts(pdb.allAtoms(),12.0).getContacts(atoms)
    c8 = Contacts(pdb.allAtoms(),8.0).getContacts(atoms)
    c6 = Contacts(pdb.allAtoms(),6.0).getContacts(atoms)
    c20 = Contacts(pdb.allAtoms(),20.0).getContacts(atoms)
    c22 = Contacts(pdb.allAtoms(),22.0).getContacts(atoms)

    bf = [e.bfactor for e in atoms]
    bfsum = sum(bf)

    if bfsum<0.001:
        print "<br><strong>Empty B-factor column. Uses SVR model without B-factor as a feature!</strong>"
# write SVR formatted output
    sofname = argv[1]+".svrinput"
    predicted_fname = argv[1]+".svrprediction"
    with open(sofname,"w") as fsvr:
        if method == 'NMR' or bfsum < 0.001:
            rescale_cmd = svm_scale+" -r "+svr_modelNMR_range_path+" "+sofname+" > "+sofname+".rescaled"
            predict_cmd = svm_predict+" "+sofname+".rescaled "+svr_modelNMR_path+" "+predicted_fname+" >/dev/null"
            for i in range(len(c16)):
                line = "1 1:%1d 2:%1.2f 3:%1d 4:%1d 5:%1d 6:%1d 7:%1d 8:%1d\n" %(c16[i],c18[i],c12[i],c8[i],c6[i],c15[i],c20[i],c22[i])
                fsvr.write(line)
        else:
            rescale_cmd = svm_scale+" -r "+svr_modelXRAY_range_path+" "+sofname+" > "+sofname+".rescaled"
            predict_cmd = svm_predict+" "+sofname+".rescaled "+svr_modelXRAY_path+" "+predicted_fname+" >/dev/null"
            for i in range(len(c16)):
                line = "1 1:%1d 2:%1.2f 3:%1d 4:%1d 5:%1d 6:%1d 7:%1d 8:%1d 9:%1d\n" %(c16[i],bf[i],c18[i],c12[i],c8[i],c6[i],c15[i],c20[i],c22[i])
                fsvr.write(line)

    p = Popen(rescale_cmd, shell=True)
    p.communicate()
    p = Popen(predict_cmd, shell=True)
    p.communicate()

# read predicted values and save it into formatted file
    with open(predicted_fname) as pred:
        predictions = [float(v) for v in pred]

    with open(argv[1]+".csv", "a+") as fwcsv:
        if len(fwcsv.readlines()) == 0:
            fwcsv.write("Residue index,Chain index,Predicted MD fluctuation value\n")
        for i in range(len(predictions)):
            fwcsv.write("%4d,%2s,%f\n" % (atoms[i].index, chain_id, predictions[i]))
# delete temporary files  
    for e in [argv[1]+".svrinput", argv[1]+".svrinput.rescaled", argv[1]+".svrprediction"]:
        remove(e)


# MAIN PROGRAM
method = argv[2]

# parse first chain
pdb = PDBRead(argv[1])
remove(argv[1]+".csv") if path.exists(argv[1]+".csv") else None
[predictionForChain(chain, pdb) for chain in pdb.chainIndices()]
# output in argv[1].csv

