import pandas as pd 
from rdkit.Chem import Descriptors

df = pd.read_csv(r'activity_data.csv')
print(df)

# Helper function to compute descriptors for a single molecule
def compute_descriptors(molecule):
	descriptors = {d[0]: d[1](molecule) for d in Descriptors.descList}
	descriptors = pd.Series(descriptors)
	return descriptors