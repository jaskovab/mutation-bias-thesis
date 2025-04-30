# Importing dependencies
import pandas as pd

# Genetic code dictionary with one-letter amino acid codes
genetic_code = {
    'UUU': 'F', 'UUC': 'F',  # Phenylalanine (F)
    'UUA': 'L', 'UUG': 'L',  # Leucine (L)
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',  # Leucine (L)
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I',  # Isoleucine (I)
    'AUG': 'M',  # Methionine (M) - Start codon
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',  # Valine (V)
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',  # Serine (S)
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',  # Proline (P)
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',  # Threonine (T)
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',  # Alanine (A)
    'UAU': 'Y', 'UAC': 'Y',  # Tyrosine (Y)
    'UAA': '*', 'UAG': '*', 'UGA': '*',  # Stop codons (*)
    'CAU': 'H', 'CAC': 'H',  # Histidine (H)
    'CAA': 'Q', 'CAG': 'Q',  # Glutamine (Q)
    'AAU': 'N', 'AAC': 'N',  # Asparagine (N)
    'AAA': 'K', 'AAG': 'K',  # Lysine (K)
    'GAU': 'D', 'GAC': 'D',  # Aspartic acid (D)
    'GAA': 'E', 'GAG': 'E',  # Glutamic acid (E)
    'UGU': 'C', 'UGC': 'C',  # Cysteine (C)
    'UGG': 'W',  # Tryptophan (W)
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',  # Arginine (R)
    'AGU': 'S', 'AGC': 'S',  # Serine (S)
    'AGA': 'R', 'AGG': 'R',  # Arginine (R)
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'   # Glycine (G)
}

nucleotides = ['A', 'C', 'G', 'U']
codons = [a+b+c for a in nucleotides for b in nucleotides for c in nucleotides]
STOP_codons = ['UAA', 'UAG', 'UGA']

# Function to translate a codon
def translate_codon(codon: str) -> str:
    return genetic_code.get(codon.upper(), "Invalid codon")

# Function to answer wether a mutation was a transition or a transversion
def is_transition(original_codon: str, mutated_codon: str) -> bool:
    transitions = [('A', 'G'), ('G', 'A'), ('C', 'U'), ('U', 'C')]
    return (original_codon[0], mutated_codon[0]) in transitions or (original_codon[1], mutated_codon[1]) in transitions or (original_codon[2], mutated_codon[2]) in transitions

# Creating a hash directory with all the possible point mutations from each codon
def create_hash_mutations() -> dict:
    possible_mutations = {}

    for codon in codons:
        original_aa = translate_codon(codon)
        curr_mutations = []
        for position in range(3):
            for nt in nucleotides:
                if nt == codon[position]:
                    continue
                mutated_codon = codon[:position] + nt + codon[position+1:]
                if (mutated_codon in STOP_codons): 
                    continue

                mutated_aa = translate_codon(mutated_codon)
                if (mutated_aa != original_aa) and (mutated_aa != 'Invalid codon'): # and (mutated_aa != '*'):
                    curr_mutations.append(mutated_aa)
                
        possible_mutations[codon] = curr_mutations

    return possible_mutations

def get_mutation_ratios (mutation_set : dict, codon : str) -> dict:
    mutation_ratios = {}
    for mutation in mutation_set[codon]:
        mutation_ratios[mutation] = mutation_set[codon].count(mutation)
    return mutation_ratios

# Reading input data into a pandas data frame
def data_transformation(file_path : str, mutation_hash : dict) -> pd.DataFrame:
    data = pd.read_csv(file_path)

    # Changing DNA codons to RNA codons
    data['Codon'] = data['Codon'].str.replace('T', 'U')

    # Adding a new column for number of point mutations that can lead to this change of amino acid
    data['Mutation_Pathways'] = 0

    # Filling the new column with the mutation ratios from possible_mutations that i prepared above
    i = 0
    for codon in data['Codon']:
        ratios = get_mutation_ratios(mutation_hash, codon)
        new_aa = data.loc[i, 'Alt_aa']
        data.loc[i, 'Mutation_Pathways'] = ratios.get(new_aa, 0)

        i += 1
    
    return data

# Function to check if mutation counts are different
def has_different_counts(ratios: dict) -> bool:
    counts = list(set(ratios.values()))
    return len(set(counts)) > 1

def seperate_interesting_codons(mutation_hash: dict) -> list:
    # Creating a subset of interesting codons - those with distinct mutation counts and not stop codons
    codons_of_interest = []
    for codon in codons:
        if codon in STOP_codons:
            continue
        if has_different_counts(get_mutation_ratios(mutation_hash, codon)):
            codons_of_interest.append(codon)
        #else:
        #    print(get_mutation_ratios(mutation_hash, codon))
    return codons_of_interest

# Sanity checks for the data
def sanity_checks(data, mutation_hash : dict):
    # 1. Check if the translation of the codon is correct
    for _, row in data.iterrows():
        codon = row['Codon']
        expected_aa = row['Initiall_aa']
        translated_aa = translate_codon(codon)
        
        if translated_aa != expected_aa:
            print(f"Sanity check failed for Codon {codon}: Expected {expected_aa}, but got {translated_aa}")
    
    # 2. Check for missing mutations (check that every codon has all possible alt_aa combinations except the original one)    
    for codon in codons:
        alt_aa_values = set(data[data['Codon'] == codon]['Alt_aa'].unique())
        expected_aa_set = set(mutation_hash[codon])
        
        if not alt_aa_values.issubset(expected_aa_set):
            print(f"Sanity check failed for Codon {codon}: Missing or incorrect Alt_aa values. Expected one of {expected_aa_set}, but found {alt_aa_values}")
    
    # 3. Check if 'Mutation_Pathways' is zero (it shouldn't be)
    for _, row in data.iterrows():
        if row['Mutation_Pathways'] == 0:
            print(f"Sanity check failed for Codon {row['Codon']}: 'Mutation_Pathways' is 0, but mutation exists.")

"""
def is_it_transition(init_base: str, alt_base: str) -> bool:
    purines = ['A', 'G']
    pyrimidines = ['C', 'T']
    
    if init_base in purines and alt_base in purines:
        print(f"Transition: {init_base} -> {alt_base}")
        return True
    elif init_base in pyrimidines and alt_base in pyrimidines:
        print(f"Transition: {init_base} -> {alt_base}")
        return True
    else:
        print(f"Not a transition: {init_base} -> {alt_base}")
        return False

def transition_tryout(file_path: str):
    # Read the data
    data = pd.read_csv(file_path)
    data["Transition"] = 0
    for idx, row in data.iterrows():
        if is_it_transition(row['Ref.base'], row['Alt.base']):
            data.at[idx, 'Transition'] = 1
        else:
            data.at[idx, 'Transition'] = 0
    
    import matplotlib.pyplot as plt
    import seaborn as sns

    plt.figure(figsize=(7,5))
    sns.boxplot(x = data['Transition'],
                y = data['Events'],
                linewidth=1)
    counts, bins, patches = plt.hist(data['Transition'], bins=2, edgecolor='white')
    summed_events = data.groupby('Transition')['Events'].sum().sort_index()
    print(summed_events)
    print(summed_events[1] / summed_events[0])
    sns.barplot(x=summed_events.index, y=summed_events.values, palette='pastel')
    plt.xticks([0, 1], ['Transversion (0)', 'Transition (1)'])
    plt.ylabel('Number of mutations')

    plt.show()

transition_tryout("../data/data_nuc_changes.csv")
"""