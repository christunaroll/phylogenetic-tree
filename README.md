# DNA Sequence Phylogenetic Tree Construction

## Program Overview:
This program performs two key tasks:
1. **Jukes-Cantor Distance Algorithm**: Calculates the evolutionary distance between pairs of DNA sequences to determine how similar or different they are.
2. **Neighbor Joining Algorithm**: Constructs a phylogenetic tree based on the Jukes-Cantor distance matrix. The phylogenetic tree is built step by step by successively joining the most similar nodes, resulting in a tree representation of the genetic relationships between the DNA sequences.

The tree is represented as a list of tuples, where each tuple contains two nodes and the distance between them. 

## Features:
- **Jukes-Cantor Distance**: Calculates the genetic distance between two sequences based on nucleotide differences.
- **Neighbor Joining**: Uses the distance matrix to iteratively join the closest sequences, constructing a phylogenetic tree.
- **CSV Output**: The resulting tree is saved in a CSV format for easy interpretation and visualization.

## DNA Sequences:
The program processes these DNA sequences from their corresponding FASTA files and constructs the phylogenetic tree:
German Neanderthal, Russian Neanderthal, European Human, Mountain Gorilla Rwanda, Chimp Troglodytes, Puti Orangutan, Jari Orangutan, Western Lowland Gorilla, Eastern Lowland Gorilla, Chimp Schweinfurthii, Chimp Vellerosus, Chimp Verus

After running, you will get the phylogenetic tree saved in `neighbor_joining_result.csv`.

## Notes:
- Ensure all FASTA files are formatted correctly (each file should have a `.fasta` extension).
- The program uses the Jukes-Cantor model for distance calculation, which assumes equal mutation rates across all positions in the sequences.
- The program includes a safeguard for files that do not have the `.fasta` extension.
