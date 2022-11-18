# This is a code to create fasta files per chain for Rosetta outputs
# Please place this file outside of the output folder

from Bio.PDB import *
from Bio import SeqIO
import os
import pandas as pd
import shutil

directory_in_str = "./output/"

AA_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def gen_structures_list(directory_in_str):
    structures_list = []
    descript_list = []
    directory = os.fsencode(directory_in_str)
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".pdb"):
            parser = PDBParser(QUIET=True)
            path = directory_in_str + filename
            structure = parser.get_structure(filename.split(".pdb")[0],path)
            structures_list.append(structure)
            descript_list.append(filename[:-4])
    return structures_list, descript_list

def create_chains(structures_list):
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
    return chain_A_list, chain_B_list, chain_C_list, chain_D_list

def create_fasta(header, seq):
    sequences = ""
    for i in range(0,len(header)):
        sequence = ">" + str(header[i]) + " \n" +seq[i] +"\n"
        sequences += sequence
    return sequences

def create_seq_dict(list_A, list_B, list_C, list_D, descript_list):
    seq_dict = {}
    for i in range (0, len(structures_list)):
        seq_descript = list_A[i] + list_B[i] + list_C[i] + list_D[i] 
        seq_dict[descript_list[i]] = seq_descript
    return seq_dict

def gen_df(directory_in_str, chain_A_list, chain_B_list, chain_C_list, chain_D_list, descript_list):
    
    seq_dict = create_seq_dict(chain_A_list, chain_B_list, chain_C_list, chain_D_list, descript_list)

    df_description = []
    df_sequence = []
    df_total_score = []
    df_dslf_fa13 = []
    df_fa_atr = []
    df_fa_dun = []
    df_fa_elec = []
    df_fa_intra_rep = []
    df_fa_intra_fa_intra_sol_xover4 = []
    df_fa_rep = []
    df_fa_sol = []
    df_hbond_bb_sc = []
    df_hbond_lr_bb = []
    df_hbond_sc = []
    df_hbond_sr_bb = []
    df_lk_ball_wtd = []
    df_net_charge = []
    df_netcharge = []
    df_omega = []
    df_p_aa_pp = []
    df_pro_close = []
    df_pstat = []
    df_rama_prepro = []
    df_ref = []
    df_yhh_planarity = []

    directory = os.fsencode(directory_in_str)
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        filepath = directory_in_str + filename
        if filename.endswith(".sc"):
            file_data = open(filepath,"r")
            score_data = file_data.readlines(0)[2].split()
            df_description.append(score_data[24])
            df_sequence.append(seq_dict[score_data[24]])
            df_total_score.append(score_data[1])
            df_dslf_fa13.append(score_data[2])
            df_fa_atr.append(score_data[3])
            df_fa_dun.append(score_data[4])
            df_fa_elec.append(score_data[5])
            df_fa_intra_rep.append(score_data[6])
            df_fa_intra_fa_intra_sol_xover4.append(score_data[7])
            df_fa_rep.append(score_data[8])
            df_fa_sol.append(score_data[9])
            df_hbond_bb_sc.append(score_data[10])
            df_hbond_lr_bb.append(score_data[11])
            df_hbond_sc.append(score_data[12])
            df_hbond_sr_bb.append(score_data[13])
            df_lk_ball_wtd.append(score_data[14])
            df_net_charge.append(score_data[15])
            df_netcharge.append(score_data[16])
            df_omega.append(score_data[17])
            df_p_aa_pp.append(score_data[18]) 
            df_pro_close.append(score_data[19])
            df_pstat.append(score_data[20])
            df_rama_prepro.append(score_data[21])
            df_ref.append(score_data[22])
            df_yhh_planarity.append(score_data[23])
    d = {"Description": df_description, "Sequence": df_sequence, "Total Score": df_total_score,
        "dslf_fa13": df_dslf_fa13, "fa_atr": df_fa_atr, "fa_dun": df_fa_dun, "fa_elec": df_fa_elec, 
        "fa_intra_rep": df_fa_intra_rep, "fa_intra_fa_intra_sol_xover4": df_fa_intra_fa_intra_sol_xover4, 
        "fa_rep": df_fa_rep, "fa_sol": df_fa_sol, "hbond_bb_sc": df_hbond_bb_sc, "hbond_lr_bb": df_hbond_lr_bb,
        "hbond_sc": df_hbond_sc, "hbond_sr_bb": df_hbond_sr_bb, "lk_ball_wtd": df_lk_ball_wtd,
        "net_charge": df_net_charge, "netcharge": df_netcharge, "omega": df_omega, "p_aa_pp": df_p_aa_pp,
        "pro_close":df_pro_close, "pstat":df_pstat, "rama_prepro": df_rama_prepro, "ref": df_ref, "yhh_planarity": df_yhh_planarity}
    return d


structures_list, descript_list = gen_structures_list(directory_in_str)
chain_A_list, chain_B_list, chain_C_list, chain_D_list = create_chains(structures_list)


chain_A_data = create_fasta(structures_list,chain_A_list)
saveChainAFasta = open(r'./chainA.fasta','w+')
saveChainAFasta.write(chain_A_data)
saveChainAFasta.close()

chain_B_data = create_fasta(structures_list,chain_B_list)
saveChainBFasta = open(r'./chainB.fasta','w+')
saveChainBFasta.write(chain_B_data)
saveChainBFasta.close()

chain_C_data = create_fasta(structures_list,chain_C_list)
saveChainCFasta = open(r'./chainC.fasta','w+')
saveChainCFasta.write(chain_C_data)
saveChainCFasta.close()

chain_D_data = create_fasta(structures_list,chain_D_list)
saveChainDFasta = open(r'./chainD.fasta','w+')
saveChainDFasta.write(chain_D_data)
saveChainDFasta.close()

print("Succesfully created fasta files for each chain! Check folder \"fasta-files\" for output.")

df_dict = gen_df(directory_in_str, chain_A_list, chain_B_list, chain_C_list, chain_D_list, descript_list)
df = pd.DataFrame(data=df_dict)
df = df.sort_values(by=['Total Score'],ascending=[False])
df.to_csv('score_data.csv',index=False)

print("Successfully generated CSV file of scores!")

# appended code to filter through the current inputs to create a new

score_df = df
score_df = score_df.iloc[0:10,:]
pstat_df = df.sort_values(by=['pstat'],ascending=[False])
pstat_df = pstat_df.iloc[0:10,:]

# print("Total score data frame:\n",score_df.iloc[0:15,:])
# print("pstat data frame:\n",pstat_df.iloc[0:15,:])

def get_name(row):
    return row + ".pdb"


def filter_score(filter_list,id):
    new_dir = os.path.join(os.getcwd(),id)
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    output_path = os.path.join(os.getcwd(),"output")
    output_list = os.listdir(output_path)
    for file in output_list:
        filename = os.fsdecode(file)
        if filename in filter_list:
            file_loc = output_path + "/" + filename
            shutil.copy(file_loc,new_dir)

def filter_chains(filter_list,filename):
    sequences = ""
    for seq_id in SeqIO.parse(filename,"fasta"):
        str_id = str(seq_id.description.split('<Structure id=')[1][:-1]) + ".pdb"
        if str_id in filter_list:
            seq_string = str(repr(seq_id.seq)).split('Seq(\'')[1][:-2]
            sequence = ">" + str_id + " \n" + seq_string +"\n"
            sequences += sequence
    return sequences

def move_fasta():
    cwd = os.getcwd()
    files = os.listdir(cwd)

    fasta_dir = os.path.join(os.getcwd(),'fasta-files')
    if not os.path.exists(fasta_dir):
        os.makedirs(fasta_dir)

    for file in files:
        filename = os.fsdecode(file)
        if filename.endswith(".fasta"):
            shutil.move(file,fasta_dir)

        

print("\nGenerating filtered outputs according to total score and packing statistics. Grabbing best pdb files...")
tot_filter = [get_name(x) for x in score_df['Description']]
pstat_filter = [get_name(x) for x in pstat_df['Description']]
filter_score(tot_filter,"total-score")
print("Created a separate folder for lowest 10 total scores located in \"total-score\".")
filter_score(pstat_filter,"pstat-score")
print("Created a separate folder for top 10 packing scores located in \"pstat-score\".")
print("May this help you with further design.")


tot_A_chain = filter_chains(tot_filter,"chainA.fasta")
saveChainAFasta = open(r'./tot-chainA.fasta','w+')
saveChainAFasta.write(tot_A_chain)
saveChainAFasta.close()
tot_B_chain = filter_chains(tot_filter,"chainB.fasta")
saveChainBFasta = open(r'./tot-chainB.fasta','w+')
saveChainBFasta.write(tot_B_chain)
saveChainBFasta.close()
tot_C_chain = filter_chains(tot_filter,"chainC.fasta")
saveChainCFasta = open(r'./tot-chainC.fasta','w+')
saveChainCFasta.write(tot_C_chain)
saveChainCFasta.close()
tot_D_chain = filter_chains(tot_filter,"chainD.fasta")
saveChainDFasta = open(r'./tot-chainD.fasta','w+')
saveChainDFasta.write(tot_D_chain)
saveChainDFasta.close()

pstat_A_chain = filter_chains(pstat_filter,"chainA.fasta")
saveChainAFasta = open(r'./pstat-chainA.fasta','w+')
saveChainAFasta.write(pstat_A_chain)
saveChainAFasta.close()
pstat_B_chain = filter_chains(pstat_filter,"chainB.fasta")
saveChainBFasta = open(r'./pstat-chainB.fasta','w+')
saveChainBFasta.write(pstat_B_chain)
saveChainBFasta.close()
pstat_C_chain = filter_chains(pstat_filter,"chainC.fasta")
saveChainCFasta = open(r'./pstat-chainC.fasta','w+')
saveChainCFasta.write(pstat_C_chain)
saveChainCFasta.close()
pstat_D_chain = filter_chains(pstat_filter,"chainD.fasta")
saveChainDFasta = open(r'./pstat-chainD.fasta','w+')
saveChainDFasta.write(pstat_D_chain)
saveChainDFasta.close()

move_fasta()