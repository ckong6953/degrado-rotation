# This is a code to create fasta files per chain for Rosetta outputs
# Please place this file outside of the output folder

from Bio.PDB import *
import os
import gzip

structures_list = []
filepath = "./output"

AA_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

directory = os.fsencode(filepath)
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".pdb"):
        parser = PDBParser(QUIET=True)
        path = filepath + "/" + filename
        structure = parser.get_structure(filename.split(".pdb")[0],path)
        structures_list.append(structure)

chain_A_list = []
chain_B_list = []
chain_C_list = []
chain_D_list = []
for structure in structures_list:
    chain_count = 0
    for model in structure:
        for chain in model:
            chain_list = []
            for residue in chain:
                chain_list.append(AA_dict[residue.get_resname()])
            new_list = "".join(chain_list)
            match chain_count:
                case 0:
                    chain_A_list.append(new_list)
                case 1:
                    chain_B_list.append(new_list)
                case 2:
                    chain_C_list.append(new_list)
                case 3:
                    chain_D_list.append(new_list)
            chain_count += 1


def createFasta(header, seq):
    sequences = ""
    for i in range(0,len(header)):
        sequence = ">" + str(header[i]) + " \n" +seq[i] +"\n"
        sequences += sequence
    return sequences


chain_A_data = createFasta(structures_list,chain_A_list)
saveChainAFasta = open(r'./chainA.fasta','w+')
saveChainAFasta.write(chain_A_data)
saveChainAFasta.close()

chain_B_data = createFasta(structures_list,chain_B_list)
saveChainBFasta = open(r'./chainB.fasta','w+')
saveChainBFasta.write(chain_B_data)
saveChainBFasta.close()

chain_C_data = createFasta(structures_list,chain_C_list)
saveChainCFasta = open(r'./chainC.fasta','w+')
saveChainCFasta.write(chain_C_data)
saveChainCFasta.close()

chain_D_data = createFasta(structures_list,chain_D_list)
saveChainDFasta = open(r'./chainD.fasta','w+')
saveChainDFasta.write(chain_D_data)
saveChainDFasta.close()

print("Succesfully created fasta files for each chain!")


