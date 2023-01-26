import pandas as pd 
from rdkit.Chem import Descriptors
import rdkit
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit import Chem
from vina import Vina
import subprocess
import re
from rdkit.Chem import AllChem
import meeko
import os
import vina
import numpy as np
import pickle

# Ideally you should run this script and it will create a dataframe with everything we need
# If something goes wrong with creating that dataframe, I save the dictionaries containing the features, and we should be able to add that to the dataframe in another script

# To run: 
# Before running, you need to replace the files that I had set on my local device
# Search for `Users/apunt` in the code and replace the files, which I uploaded to our repo, with their locations on your machine
# Here are the files you'll need:
# 	`/Users/apunt/Downloads/7jkv.pdb`
#	`/Users/apunt/Downloads/mproFH.pdbqt`
# You'll also need to replace the following code and download autodock_vina on your machine
# 	`/Users/apunt/Downloads/autodock_vina_1_1_2_mac_catalina_64bit/bin/vina` - vina program, I sent you a zip of this, or you can download it online
# 	`/Users/apunt/Downloads/autodock_vina_1_1_2_mac_catalina_64bit/bin/vina_split` - comes with vina program

# You'll also probably need to pip install some of the imported modules

# Once you're done replacing these, run this file however you would run a python script
# I use `python3 -m create_features`, although that's because of how I set up virtual python environments on my device, it'll probably be different for you


def add_rdkit_descriptors(df):
	# List of descriptors from rdkit
	descriptors_list = [x[0] for x in Descriptors._descList]
	# Create a dictionary representing all the descriptors output by rdkit 
	rdkit_dict = dict.fromkeys(descriptors_list)
	for key in rdkit_dict:
		rdkit_dict[key] = []
	# Iterate through SMILES strings and calculate descriptors
	for smile_str in df['SMILES']:
		# Convert SMILES string to rdkit Molecule object
		mol = Chem.MolFromSmiles(smile_str)
		# Calculate rdkit descriptors and add them to the dictionary
		calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
		ds = calc.CalcDescriptors(mol)
		for i in range(len(ds)):
			rdkit_dict[descriptors_list[i]].append(ds[i])
	# Convert the rdkit dictionary 
	df_dictionary = pd.DataFrame([rdkit_dict])
	# Add the rdkit dictionary to the main dataframe
	df = pd.concat([df, df_dictionary], ignore_index=True)
	return df

def autodock_features(df):
	# Load the receptor pdb
	receptor = '/Users/apunt/Downloads/7jkv.pdb'
	# Append features to keys of dictionary
	keys = []
	keys.append('SMILE')
	for i_smile in range(1,10):
		keys.append('gauss1_' + str(i_smile))
		keys.append('gauss2_' + str(i_smile))
		keys.append('repulsion_' + str(i_smile))
		keys.append('hydrophobic_' + str(i_smile))
		keys.append('hydrogen_' + str(i_smile))
		keys.append('affinity_' + str(i_smile))
		keys.append('distl_' + str(i_smile))
		keys.append('distr_' + str(i_smile))
	docking_dict = dict.fromkeys(keys)
	for key in docking_dict:
		docking_dict[key] = []
	count = 0
	for smile_str in df['SMILES']:
		docking_dict['SMILE'] = smile_str
		print("current SMILE string: " + smile_str)
		try: 
			lig = rdkit.Chem.MolFromSmiles(smile_str)
			# Add hydrogens
			protonated_lig = rdkit.Chem.AddHs(lig)
			# Generate 3D coordinates
			rdkit.Chem.AllChem.EmbedMolecule(protonated_lig)
			# Prepare the molecule using meeko
			meeko_prep = meeko.MoleculePreparation()
			meeko_prep.prepare(protonated_lig)
			lig_pbdqt = meeko_prep.write_pdbqt_string()
			# Convert to pdbqt file
			text_file = open("output.txt", "w")
			text_file.write(lig_pbdqt)
			text_file.close()
			src = "output.txt"
			dest = "output.pdbqt"
			os.rename(src, dest)
			# Docking
			cmd = f'/Users/apunt/Downloads/autodock_vina_1_1_2_mac_catalina_64bit/bin/vina --receptor /Users/apunt/Downloads/mproFH.pdbqt --ligand output.pdbqt --center_x 8.460 --center_y -3.609 --center_z 19.031 --size_x 20 --size_y 20 --size_z 20 --out lig_test.pdbqt'
			#  Parse through output of docking
			result = subprocess.run(cmd.split(), stdout=subprocess.PIPE)
			result = str(result)
			lines = result.split('\\n')
			# Get individual modes
			index = 1
			for i_line in range(29,38):
				line = lines[i_line]
				docking_dict = save_vals(line, docking_dict, index)
				index += 1

			# Split multi-model pbdqt 
			cmd = f'/Users/apunt/Downloads/autodock_vina_1_1_2_mac_catalina_64bit/bin/vina_split --input lig_test.pdbqt'
			result = subprocess.run(cmd.split())

			# Generate scores based on docking results for each mode
			for i in range(1,10):
				try:
					saved = False # check to see if we're actually able to save the scores
					# Generate filename created by vina_split on multi-model output from docking
					filename = "lig_test_ligand_" + str(i)
					# Calculate scores for each mode
					cmd = f'/Users/apunt/Downloads/autodock_vina_1_1_2_mac_catalina_64bit/bin/vina --receptor /Users/apunt/Downloads/mproFH.pdbqt --ligand {filename}.pdbqt --center_x 8.460 --center_y -3.609 --center_z 19.031 --size_x 20 --size_y 20 --size_z 20 --score_only'
					# Parse through output of vina for the mode's scores
					result = subprocess.run(cmd.split(), stdout=subprocess.PIPE)
					result = str(result)
					lines = result.split('\\n')
					for line in lines:
						if "affinity" in line:
							val = float(re.findall("[-+]?\d+\.\d+", line)[0])
							docking_dict['affinity_' + str(i)].append(val)
							saved = True
						if "gauss 1" in line:
							val = float(re.findall("[-+]?\d+\.\d+", line)[0])
							docking_dict['gauss1_' + str(i)].append(val)
						if "gauss 2" in line:
							val = float(re.findall("[-+]?\d+\.\d+", line)[0])
							docking_dict['gauss2_' + str(i)].append(val)
						if "repulsion" in line:
							val = float(re.findall("[-+]?\d+\.\d+", line)[0])
							docking_dict['repulsion_' + str(i)].append(val)
						if "hydrophobic" in line: 
							val = float(re.findall("[-+]?\d+\.\d+", line)[0])
							docking_dict['hydrophobic_' + str(i)].append(val)
						if "Hydrogen" in line:
							val = float(re.findall("[-+]?\d+\.\d+", line)[0])
							docking_dict['hydrogen_' + str(i)].append(val)
					if not saved:
						# If scores did not generate for a given mode 
						docking_dict = add_nan(docking_dict, i, False)
				except:
					# If scores did not generate for a given mode due to an error
					print("scores failed for " + smile_str + " and mode " + str(i))
					docking_dict = add_nan(docking_dict, i, False)
		except:
			# If docking failed for given SMILE string, add NaN to entire row 
			print("docking failed for " + smile_str)
			for i in range(1,10):
				docking_dict = add_nan(docking_dict, i, True)
		# save dictionary every 100 iterations, as a backup
		if count == 5:
			count = 0
			with open('saved_docking_dictionary_intermediary.pkl', 'wb') as f:
				pickle.dump(docking_dict, f)
		count += 1
	# Save final docking dictionary
	with open('saved_docking_dictionary.pkl', 'wb') as f:
		pickle.dump(docking_dict, f)
	# Add the scores to the dataframe
	for i in range(1,10):
		df['gauss1_' + str(i)] = docking_dict['gauss1_' + str(i)]
		df['gauss2_' + str(i)] = docking_dict['gauss2_' + str(i)]
		df['repulsion_' + str(i)] = docking_dict['repulsion_' + str(i)]
		df['hydrophobic_' + str(i)] = docking_dict['hydrophobic_' + str(i)]
		df['hydrogen_' + str(i)] = docking_dict['hydrogen_' + str(i)]
		df['affinity_' + str(i)] = docking_dict['affinity_' + str(i)]
		df['distl_' + str(i)] = docking_dict['distl_' + str(i)]
		df['distr_' + str(i)] = docking_dict['distr_' + str(i)]
	return df

# Helper function to add NaN values to the row or scores associated with a given mode, depending on the value of the bool docking
def add_nan(docking_dict, i, docking):
	docking_dict['affinity_' + str(i)].append(np.nan)
	docking_dict['gauss1_' + str(i)].append(np.nan)
	docking_dict['gauss2_' + str(i)].append(np.nan)
	docking_dict['repulsion_' + str(i)].append(np.nan)
	docking_dict['hydrophobic_' + str(i)].append(np.nan)
	docking_dict['hydrogen_' + str(i)].append(np.nan)
	if docking:
		docking_dict['distl_' + str(i)].append(np.nan)
		docking_dict['distr_' + str(i)].append(np.nan)
	return docking_dict

# Helper function to save values from docking to docking_dict
def save_vals(line, docking_dict, j):
	vals = re.findall("[-+]?\d+\.\d+", line)
	# docking_dict['affinity' + str(j)].append(float(vals[0]))
	docking_dict['distl_' + str(j)].append(float(vals[1]))
	docking_dict['distr_' + str(j)].append(float(vals[2]))
	return docking_dict

# Generates features using rdkit and autovina docking
def add_all_features():
	# Load Covid moonshot csv file
	df = pd.read_csv(r'activity_data.csv')
	df = add_rdkit_descriptors(df)
	df = autodock_features(df)
	df.to_pickle("full_dataset.pkl")

add_all_features()

# Random notes, you can ignore these: 
# /Users/apunt/Downloads/autodock_vina_1_1_2_mac_catalina_64bit/bin/vina --receptor ~/Downloads/mproFH.pdbqt --ligand ~/Downloads/lig.pdbqt --center_x 8.460 --center_y -3.609 --center_z 19.031 --size_x 20 --size_y 20 --size_z 20 --out lig_test.sdf --log scoring_result.log --score_only

# 6lu7
# find non-covalent inhibitors of covid to validate docking
# add_rdkit_descriptors()


# http://molprobity.biochem.duke.edu/index.php
# Add hydrogens and flip states for mpro

# https://github.com/forlilab/Meeko
# prepare ligand

# https://ccsb.scripps.edu/agfr/
# agfrgui for preparing the protein
# set grid box for docking

# keep pole hydrogens 

# check fregalysis for non-covalent ligands

# /Users/apunt/Downloads/autodock_vina_1_1_2_mac_catalina_64bit/bin/vina --receptor ~/Downloads/mproFH.pdbqt --ligand ~/Downloads/lig.pdbqt --center_x 8.460 --center_y -3.609 --center_z 19.031 --size_x 20 --size_y 20 --size_z 20

# Validate using known, non-covalent binding
# Make sure affinity seems reasonable (-8.0)
# Find center coordinates for docking


