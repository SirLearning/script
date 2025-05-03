import matplotlib.pyplot as plt
import numpy as np
import biopython as bp

# Generate random DNA sequence
sequence = ''.join(np.random.choice(['A', 'C', 'G', 'T'], size=100))

# Compute position-specific scoring matrix
pwm = np.random.rand(4, 10)

# Convert the string array to a float64 array
sequence = np.array(sequence, dtype=np.float64)
sequence[np.isnan(sequence)] = np.nan

# Compute logo
logo = np.zeros((100, 4))
for i in range(4):
    logo[:, i] = pwm[i, :] * sequence

# Plot logo
plt.imshow(logo, cmap='Greys_r')
plt.xlabel('Position')
plt.ylabel('Base')
plt.show()