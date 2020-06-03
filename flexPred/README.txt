FlexPred takes a protein structure in PDB format as input, and predict the
fluctuation of each residue (Calpha atom) by considering  
C-alpha atom contacts number (cutoffs: 16Å, 15Å, 18Å, 12Å, 8Å, 6Å, 20Å, 22Å) with set of cutoffs of the protein with the Support Vector Regression.
As output program create text file with data and png image with fluctuation plot (and DSSP sec. str. assignment).


predictFluct.py usage:

1) Download libsvm package (http://www.csie.ntu.edu.tw/~cjlin/libsvm/) and install it
2) Edit predictFluct.py file and set path to files as above:
svm_scale="/some/path/svm-scale"        # path to previously compiled libsvm application
svm_predict="/some/path/svm-predict"    # path to previously compiled libsvm application
set paths to SVR files:
svr_modelXRAY_path
svr_modelXRAY_range_path
svr_modelNMR_path
svr_modelNMR_range_path


XRAY = Set15 (with Bfactor)
NMR = Set16 (without Bfactor)

gnuplot="/usr/bin/gnuplot" # if you have gnuplot and want to plot png with predicted data - set path. If not - comment this line by "#"
dssp_path="/path/to/dssp" # you can download dssp from swift.cmbi.ru.nl/gv/dssp/ . Script use it to draw secondary structure on fluctuation plot.

3) set chmod u+x predictFluct.py

4) in example dir you have 1BFG protein, try to predict fluctuation on this model:
  cd example
  ../predictFluct.py 1BFG.pdb XRAY
  and after few seconds you should got 1BFG.pdb.fPred and 1BFG.pdb.fPred.png files. 

5) general usage:
/path/to/predictFluct.py  model_file.pdb [NMR|XRAY]


Please cite our paper: 
"M Jamroz, A Kolinski & D. Kihara: Structural features of proteins that predict real-value fluctuation of globular proteins" 
if you publish results of this program.


This script was developed by Michal Jamroz (jamroz@chem.uw.edu.pl) in 2011.
Script license: public domain
