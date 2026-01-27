import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Settings
images = [
    # "/data/home/tusr1/git/script/out/king_ibs_dist_full_py.png",
    "/data/home/tusr1/git/script/out/raw_missing_filtered_ibs_distribution_linear.png",
    "/data/home/tusr1/git/script/out/raw_missing_filtered_ibs_distribution_log.png"
]
output_file = "/data/home/tusr1/git/script/out/raw_missing_filtered_ibs_combined_dist.png"

# Create a figure with subplots
fig, axes = plt.subplots(1, len(images), figsize=(6 * len(images), 4))

for ax, img_path in zip(axes, images):
    try:
        img = mpimg.imread(img_path)
        ax.imshow(img)
        ax.axis('off')  # Hide axis
    except FileNotFoundError:
        print(f"Error: File not found {img_path}")
        ax.text(0.5, 0.5, 'Image Not Found', ha='center', va='center')
        ax.axis('off')

plt.tight_layout()
plt.savefig(output_file, dpi=300)
print(f"Combined image saved to {output_file}")
