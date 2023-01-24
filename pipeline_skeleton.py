###vina script skeleton


#note: sdf and pdbqt files are interchangable interms of the workflow! but we need to be consistent

### PLAN

#use vina and rdkit to generate features for XGboost to run regression model
#can then use the outputs of these to run a graph net
#can run GNN to find a structural binding affinity which we can concat to the features and let XGBoost figure out what components it wants to use
# this will then create a model which uses both structural and ligand based features - this combined with use of XGboost should give us really good results!
#in garrets presentation he showed us the some of the highest accuracy results was achieved with a SB+LB w/ Vina feat framework



#use Vina to go from SMILES to descriptors and generate sdf files


### to change!
#save descriptors and sdf files using the CID as a name - save to 
#if we are going to train multiple models then maybe we just select the highest affinity pose?
#only working on one sdf model makes our downstream tasks much simpler! (but we can change this)




#use sdf files to generate descriptors from RDkit
#save using CID 






#use MDAnalysis to grab bond information and genarate topology which we can use to parse graph links
#node embeddings(just use everything in the sdf): element, charge, occupancy etcccc
#edge embeddings: bond type(single/double), bond length

#what kind of graph ml model do we want to use to calculate the binding affinity from this?



#concatenate all the data into one table
#can then jsut shove this into XGBoost?