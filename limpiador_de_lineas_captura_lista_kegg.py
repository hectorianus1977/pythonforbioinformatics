import re

# Define the regular expression pattern to match and keep the content after the KEGG code format
pattern = re.compile(r'^.*?([a-zA-Z]{2,5}\s+\t:\s+.*)$')

# Read the content of the file
with open('lista_KEGG.txt', 'r') as file:
    lines = file.readlines()

# Process the lines, keeping content after the KEGG code, discarding everything before it
processed_lines = []
for line in lines:
    match = pattern.sub(r'\1', line)
    processed_lines.append(match)

# Write the processed lines back to the file
with open('lista_KEGG_completa.txt', 'w') as file:
    file.writelines(processed_lines)
