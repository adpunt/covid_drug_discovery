import pandas as pd 
from rdkit.Chem import Descriptors
import rdkit
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit import Chem

def test_compute_descriptors():
	# List of descriptors from rdkit
	descriptors_list = [x[0] for x in Descriptors._descList]
	# Create a dictionary representing all the descriptors output by rdkit 
	rdkit_dict = dict.fromkeys(descriptors_list)
	for key in rdkit_dict:
		rdkit_dict[key] = []
	index = 0
	# Load Covid moonshot csv file
	df = pd.read_csv(r'activity_data.csv')
	# Iterate through SMILES strings and calculate descriptors
	for smile_str in df['SMILES']:
		# Convert SMILES string to rdkit Molecule object
		mol = Chem.MolFromSmiles(smile_str)
		# Calculate rdkit descriptors and add them to the dictionary
		calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
		ds = calc.CalcDescriptors(mol)
		for i in range(len(ds)):
			rdkit_dict[descriptors_list[i]].append(ds[i])
		# Print index to show progress 
		print(index)
		index += 1
	# Convert the rdkit dictionary 
	df_dictionary = pd.DataFrame([rdkit_dict])
	# Add the 
	df = pd.concat([df, df_dictionary], ignore_index=True)
	print(df)

# 6lu7
# find non-covalent inhibitors of covid to validate docking
test_compute_descriptors()


# http://molprobity.biochem.duke.edu/index.php
# Add hydrogens and flip states for mpro

# https://github.com/forlilab/Meeko
# prepare ligand

# https://ccsb.scripps.edu/agfr/
# agfrgui for preparing the protein
# set grid box for docking

# keep pole hydrogens 

# check fregalysis for non-covalent ligands