from Bio import SeqIO
import pandas as pd
import os
import uuid

# Function to count amino acids and calculate percentages
def count_amino_acids(seq):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    count_dict = {aa: 0 for aa in amino_acids}
    total_count = 0

    # Count occurrences of amino acids
    for aa in amino_acids:
        count = seq.count(aa)
        count_dict[aa] += count
        total_count += count

    # Calculate percentages
    percentages = {aa: (count / total_count) * 100 for aa, count in count_dict.items()}
    return percentages

# List of filenames
file_list = ['K01190_modified_fasta.fasta', 'K01220_modified_fasta.fasta', 'K01224_modified_fasta.fasta', 'K12111_modified_fasta.fasta', 'K12308_modified_fasta.fasta', 'K12309_modified_fasta.fasta']  # Replace with your file names

# Process each file and store the amino acid percentages for each sequence
results = {}
for file_index, file in enumerate(file_list, start=1):
    sequences = SeqIO.parse(file, 'fasta')
    for seq_index, seq_record in enumerate(sequences, start=1):
        amino_acid_percentages = count_amino_acids(str(seq_record.seq))
        header_words = seq_record.description.split()[:2]  # Extract first two words from the header
        base_filename = os.path.splitext(os.path.basename(file))[0]
        seq_uuid = uuid.uuid4().hex  # Generate a random UUID
        modified_identifier = f"{base_filename}_{'_'.join(header_words)}_{seq_uuid}"
        results[modified_identifier] = amino_acid_percentages

# Convert the results to a DataFrame and save as CSV
df = pd.DataFrame(results).T
output_file = 'individual_sequence_amino_acid_percentages.csv'
df.to_csv(output_file)
print(f"File '{output_file}' has been generated successfully.")






