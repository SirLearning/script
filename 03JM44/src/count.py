import re

# Initialize a dictionary to store the sums
sums = {}

# Open the file and read it line by line
with open('newplot.txt', 'r') as file:
    for line in file:
        # Split the line into parts by spaces or colons
        parts = line.split()
        last_part = ''
        max_percentage = 0
        max_word = ''
        count_number = 0

        # For each part, check if it contains a percentage
        for part in parts:
            match = re.match(r'(\w+)_\w+', last_part)
            if match:
                word = match.groups()
                partc = part.split(":")
                percentage = float(partc[0])
                # Keep track of the word with the highest percentage
                if percentage > max_percentage:
                    max_percentage = percentage
                    max_word = word
            last_part = part

        # Extract the count number from the last part
        count_number = int(parts[-1])

        print(f'parts: {parts}')
        print(f'max_percentage: {max_percentage}')
        print(f'max_word: {max_word}')
        print(f'count_number: {count_number}')

        # Add the count number to the total for the first three letters of the max word
        if max_word:
            sums[max_word[:3]] = sums.get(max_word[:3], 0) + count_number

# Print the sums
for key, value in sums.items():
    print(f'{key}: {value}')