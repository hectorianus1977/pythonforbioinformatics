import re

# Read the list of KEGG codes and corresponding names
kegg_data = {}
with open('lista_KEGG_completa.txt', 'r') as kegg_file:
    for line in kegg_file:
        parts = line.strip().split('\t:')
        if len(parts) >= 2:
            kegg_code = parts[0].strip()
            species_name = parts[1].strip()
            kegg_data[kegg_code] = species_name  # Store KEGG code as key

# Regular expression pattern to match KEGG codes in the header
kegg_pattern = re.compile(r'>([a-zA-Z0-9_]+)')

# Process the FASTA file
new_fasta_lines = []
with open('K18579.fasta', 'r') as fasta_file:
    for line in fasta_file:
        match = kegg_pattern.match(line)
        if match:  # Header line
            kegg_code = match.group(1)  # Extract KEGG code from the header

            # Check if KEGG code exists in the kegg_data
            if kegg_code in kegg_data:
                species_name = kegg_data[kegg_code]
                new_header = f">{species_name} {' '.join(line.strip().split()[1:])}"  # Update the header
                new_fasta_lines.append(new_header + '\n')
                continue  # Move to the next line
        new_fasta_lines.append(line)  # Keep the original line if KEGG code not found or line is sequence

# Write the modified FASTA content to a new file
with open('K18579_modified_fasta.fasta', 'w') as output_file:
    output_file.writelines(new_fasta_lines)
