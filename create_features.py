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

def add_rdkit_descriptors():
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
	# Add the rdkit dictionary to the main dataframe
	df = pd.concat([df, df_dictionary], ignore_index=True)
	print(df)

def autodock_scores():
	df = pd.read_csv(r'activity_data.csv')
	receptor = '/Users/apunt/Downloads/7jkv.pdb'
	docking_dict = dict.fromkeys(['gauss1', 'gauss2', 'repulsion', 'hydrophobic', 'hydrogen'])
	for key in docking_dict:
		docking_dict[key] = []
	for smile_str in df['SMILES']:
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
			text_file = open("output.txt", "w")
			text_file.write(lig_pbdqt)
			text_file.close()
			src = "output.txt"
			dest = "output.pdbqt"
			os.rename(src, dest)
			cmd = f'/Users/apunt/Downloads/autodock_vina_1_1_2_mac_catalina_64bit/bin/vina --receptor /Users/apunt/Downloads/mproFH.pdbqt --ligand output.pdbqt --center_x 8.460 --center_y -3.609 --center_z 19.031 --size_x 20 --size_y 20 --size_z 20 --out lig_test.sdf --log scoring_result.log --score_only'
			result = subprocess.run(cmd.split(), stdout=subprocess.PIPE)
			result = str(result)
			lines = result.split('\\n')
			for line in lines:
				if "gauss 1" in line:
					val = float(re.findall("\d+\.\d+", line)[0])
					docking_dict['gauss1'].append(val)
				if "gauss 2" in line:
					val = float(re.findall("\d+\.\d+", line)[0])
					docking_dict['gauss2'].append(val)
				if "repulsion" in line:
					val = float(re.findall("\d+\.\d+", line)[0])
					docking_dict['repulsion'].append(val)
				if "hydrophobic" in line: 
					val = float(re.findall("\d+\.\d+", line)[0])
					docking_dict['hydrophobic'].append(val)
				if "Hydrogen" in line:
					val = float(re.findall("\d+\.\d+", line)[0])
					docking_dict['hydrogen'].append(val)
		except:
			print("didn't work for: ")
			print(smile_str )
	df_dictionary = pd.DataFrame([docking_dict])
	print(df_dictionary)
	df = pd.concat([df, df_dictionary], ignore_index=True)

def autodock_features():
	df = pd.read_csv(r'activity_data.csv')
	receptor = '/Users/apunt/Downloads/7jkv.pdb'
	docking_dict = dict.fromkeys(['affinity1', 'distl1', 'distr1', 'affinity2', 'distl2', 'distr2', 'affinity3', 'distl3', 'distr3', 'affinity4', 'distl4', 'distr4', 'affinity5', 'distl5', 'distr5', 'affinity6', 'distl6', 'distr6', 'affinity7', 'distl7', 'distr7', 'affinity8', 'distl8', 'distr8', 'affinity9', 'distl9', 'distr9'])
	for key in docking_dict:
		docking_dict[key] = []
	for smile_str in df['SMILES']:
		# try: 
		lig = rdkit.Chem.MolFromSmiles(smile_str)
		# Add hydrogens
		protonated_lig = rdkit.Chem.AddHs(lig)
		# Generate 3D coordinates
		rdkit.Chem.AllChem.EmbedMolecule(protonated_lig)
		# Prepare the molecule using meeko
		meeko_prep = meeko.MoleculePreparation()
		meeko_prep.prepare(protonated_lig)
		lig_pbdqt = meeko_prep.write_pdbqt_string()
		text_file = open("output.txt", "w")
		text_file.write(lig_pbdqt)
		text_file.close()
		src = "output.txt"
		dest = "output.pdbqt"
		os.rename(src, dest)
		cmd = f'/Users/apunt/Downloads/autodock_vina_1_1_2_mac_catalina_64bit/bin/vina --receptor /Users/apunt/Downloads/mproFH.pdbqt --ligand output.pdbqt --center_x 8.460 --center_y -3.609 --center_z 19.031 --size_x 20 --size_y 20 --size_z 20 --out lig_test.sdf'
		result = subprocess.run(cmd.split(), stdout=subprocess.PIPE)
		result = str(result)
		lines = result.split('\\n')
		for i in range(len(lines)):
			line = lines[i]
			if "mode" in line:
				for j in range(1,9):
					line = lines[i+j]
					docking_dict = save_vals(line, docking_dict, j)
		print(docking_dict)
		break
		# except:
		# 	print("didn't work for: ")
		# 	print(smile_str )
	df_dictionary = pd.DataFrame([docking_dict])
	print(df_dictionary)
	df = pd.concat([df, df_dictionary], ignore_index=True)

def save_vals(line, docking_dict, j):
	vals = re.findall("\d+\.\d+", line)
	print(vals)
	docking_dict['affinity' + str(j)].append(float(vals[0]))
	docking_dict['distl' + str(j)].append(float(vals[1]))
	docking_dict['distr' + str(j)].append(float(vals[2]))
	return docking_dict


autodock_features()
# try box of 25 instead of 20

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





