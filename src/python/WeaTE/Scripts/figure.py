# Define the range and step size for x
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

x_min = 0
x_max = 20
step_size = 0.01

# Generate x values
x = np.arange(x_min, x_max, step_size)

sigma1 = 0.6
sigma2 = 0.2
sigma3 = 0.05


mu1 = 0
mu2 = 1.6
mu3 = 2.8 # Shift the third distribution further to the right

# Recalculate the distributions with the new parameters

# y1 = (1 / (sigma1 * np.sqrt(2 * np.pi))) * np.exp(-(x - mu1)**2 / (2 * sigma1**2))
# y2 = (1 / (sigma2 * np.sqrt(2 * np.pi))) * np.exp(-(x - mu2)**2 / (2 * sigma2**2))
# y3 = (1 / (sigma3 * np.sqrt(2 * np.pi))) * np.exp(-(x - mu3)**2 / (2 * sigma3**2))
y1 = (1 / (x * sigma1 * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu1)**2 / (2 * sigma1**2))
y2 = (1 / (x * sigma2 * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu2)**2 / (2 * sigma2**2))
y3 = (1 / (x * sigma3 * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu3)**2 / (2 * sigma3**2))

# Plot with filled areas and transparency

plt.figure(figsize=(6, 6))
plt.fill_between(x, y1, color='#f38121', alpha=0.5)
plt.fill_between(x, y2, color='#1b78b3', alpha=0.5)
plt.fill_between(x, y3, color='#2da340', alpha=0.5)
# plt.fill_between(x, y4, color='purple', alpha=0.5)

# Set x-axis limits
plt.xlim(-1, 22)
plt.ylim(0, 1)
# Remove grid
plt.grid(False)

# Remove axis ticks and labels
plt.xticks([])
plt.yticks([])
# Add title and labels
plt.title('TE age distribution')
plt.xlabel('time')
plt.ylabel('TE activity')

# Show plot
plt.legend().set_visible(False)
plt.show()
# plt.savefig('plot.pdf')