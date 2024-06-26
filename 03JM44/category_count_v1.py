# Open the GFF3 file for reading
with open('npJM44.gff', 'r') as gff_file:
    # Initialize a dictionary to store category counts
    category_counts = {}

    # Iterate through each line in the file
    for line in gff_file:
        # Skip comment lines starting with #
        if line.startswith('#'):
            continue

        # Split the line into fields
        fields = line.strip().split('\t')

        # Check if the line has enough fields
        if len(fields) >= 9:
            # Extract the attributes field
            attributes = fields[8]

            # Split the attributes field into key-value pairs
            attribute_pairs = [pair.split('=') for pair in attributes.split(';')]
            
            # skip poor identity
            skip_line = False;
            for key, value in attribute_pairs:
                if key == 'Identity':
                    if value == 'NA' or float(value) < 0.8:
                        skip_line = True
            if skip_line:
                continue

            # Extract the category/type from the attribute pairs
            for key, value in attribute_pairs:              
                if key == 'Classification':
                    category = value 
                    break
                else:
                    # If 'Type' is not found, use the last key-value pair
                    category = attribute_pairs[-1][0]
            
            # Update the category count in the dictionary
            category_counts[category] = category_counts.get(category, 0) + 1

# Print the category counts
sum_count = 0
for category, count in category_counts.items():
    print(f'{category}: {count}')
    sum_count = sum_count + count
print(sum_count)
