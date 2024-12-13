import os
import math
import csv

"""
PROGRAM DETAILS:
- Performs Jukes Cantor Distance algorithm to determine distance between DNA sequences (how different/similar they are)
- Performs the Neighbor Joining algorithm to construct a phylogenetic tree from a distance matrix.
- Tree is represented by a list, where each tuple contains two nodes and the distance between them
"""

"""
DNA SEQUENCES: 
'German_Neanderthal.fasta' = 0 
'Russian_Neanderthal.fasta' = 1
'European_Human.fasta' = 2
'Mountain_Gorilla_Rwanda.fasta' = 3
'Chimp_Troglodytes.fasta' = 4 
'Puti_Orangutan.fasta' = 5
'Jari_Orangutan.fasta' = 6 
'Western_Lowland_Gorilla.fasta' = 7
'Eastern_Lowland_Gorilla.fasta' = 8
'Chimp_Schweinfurthii.fasta' = 9 
'Chimp_Vellerosus.fasta' = 10
'Chimp_Verus.fasta' = 11
"""

# Define the valid characters for the DNA sequences
valid_chars = {'A', 'T', 'C', 'G'}

# Define a function to read in a FASTA file and extract the DNA sequence
def read_fasta(filename):
    """Read a FASTA file and return the DNA sequence."""
    with open(filename) as file:
        lines = file.readlines()

    seq = ''
    for line in lines:
        if line.startswith('>'):
            continue
        line = line.strip()
        if any(c not in 'ATCG' for c in line):
            continue
        seq += line

    return seq

# Define a function to read in multiple FASTA files and create an array of sequences
def create_seq_array(filenames):
    seq_array = []
    for filename in filenames:
        # Validate the file extension
        if os.path.splitext(filename)[1] != '.fasta':
            raise ValueError('Invalid file extension')
        seq = read_fasta(filename)
        seq_array.append(seq)
    return seq_array

# Define a function to calculate the Jukes-Cantor distance between two sequences
def jukes_cantor_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        return None
    num_diff = sum(1 for i in range(len(seq1)) if seq1[i] != seq2[i])
    p = abs(num_diff / len(seq1))
    if p >= 3/4:
        return None
    dist = (-3/4) * math.log(1 - (4/3 * p))
    return dist

# Define a function to create a distance matrix using the Jukes-Cantor distances
def create_distance_matrix(seq_array):
    n = len(seq_array)
    dist_matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            dist = jukes_cantor_distance(seq_array[i], seq_array[j])
            if dist is not None:
                dist_matrix[i][j] = dist
                dist_matrix[j][i] = dist
    return dist_matrix

def neighbor_joining(dist_matrix):
    n = len(dist_matrix)

    # Base case: If only two sequences remain, return the direct connection
    if n == 2:
        return [(0, 1, dist_matrix[0][1])]

    # Step 1: Compute the Q-matrix
    q_matrix = []
    for i in range(n):
        for j in range(i + 1, n):
            q_value = (dist_matrix[i][j] * (n - 2)) - sum(dist_matrix[i]) - sum(dist_matrix[j])
            q_matrix.append((q_value, i, j))

    # Step 2: Find the pair of nodes with the minimum Q-value
    min_q_value, min_i, min_j = min(q_matrix)

    # Step 3: Calculate the distance of the new node to the chosen pair
    total_dist_i = sum(dist_matrix[min_i])
    total_dist_j = sum(dist_matrix[min_j])
    delta = (total_dist_i - total_dist_j) / (n - 2)
    dist_to_new_node_i = (dist_matrix[min_i][min_j] + delta) / 2
    dist_to_new_node_j = dist_matrix[min_i][min_j] - dist_to_new_node_i

    # Step 4: Create the new distance matrix
    new_dist_matrix = []
    for k in range(n):
        if k != min_i and k != min_j:
            new_dist = (dist_matrix[min_i][k] + dist_matrix[min_j][k] - dist_matrix[min_i][min_j]) / 2
            new_dist_matrix.append(new_dist)

    # Add the new node distances
    new_dist_matrix.append([dist_to_new_node_i, dist_to_new_node_j])

    # Recursively perform Neighbor Joining on the reduced matrix
    reduced_matrix = [
        [dist_matrix[x][y] for y in range(n) if y != min_i and y != min_j]
        for x in range(n) if x != min_i and x != min_j
    ]
    tree = neighbor_joining(reduced_matrix)

    # Step 5: Add the new branches to the tree
    new_node_index = n - 2
    tree.append((min_i, new_node_index, dist_to_new_node_i))
    tree.append((min_j, new_node_index, dist_to_new_node_j))

    return tree

def write_csv(filename, tree):
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sequence 1', 'Sequence 2', 'Distance'])
        for branch in tree:
            writer.writerow(branch)


# Define the list of filenames
filenames = ['German_Neanderthal.fasta', 'Russian_Neanderthal.fasta', 'European_Human.fasta',
             'Mountain_Gorilla_Rwanda.fasta', 'Chimp_Troglodytes.fasta', 'Puti_Orangutan.fasta',
             'Jari_Orangutan.fasta', 'Western_Lowland_Gorilla.fasta', 'Eastern_Lowland_Gorilla.fasta',
             'Chimp_Schweinfurthii.fasta', 'Chimp_Vellerosus.fasta', 'Chimp_Verus.fasta']

# Create the array of sequences
seq_array = create_seq_array(filenames)

# Create the distance matrix
dist_matrix = create_distance_matrix(seq_array)

# Perform the neighbor joining algorithm
tree = neighbor_joining(dist_matrix)
write_csv('neighbor_joining_result.csv', tree)
