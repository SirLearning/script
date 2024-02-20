# Open the GFF3 file for reading
with open('npJM44.gff') as gff_file:
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

            skip_line = False

            # Extract the category/type from the attribute pairs
            for key, value in attribute_pairs:
                if key == 'Classification':
                    category = value
                else:
                    # If 'Type' is not found, use the last key-value pair
                    category = attribute_pairs[-1][0]
                if key == 'Identity' and float(str(value)) < 0.8:
                    skip_line = True
            # skip poor identity
            if skip_line:
                continue
            # Update the category count in the dictionary
            category_counts[category] = category_counts.get(category, 0) + 1

# Print the category counts
for category, count in category_counts.items():
    print(f'{category}: {count}')