import matplotlib.pyplot as plt
import matplotlib.image as mpimg

def combine_plots(
    image1 = "/data/home/tusr1/git/script/out/site_depth_variant_reg_log_mean_vs_log_var.png",
    image2 = "/data/home/tusr1/git/script/out/site_depth_variant_reg_log_mean_vs_log_cv.png",
    output_file = "/data/home/tusr1/git/script/out/site_depth_variant_reg_log_mean_vs_log_combined.png"
):
    # Settings
    images = [
        image1,
        image2
    ]

    # Create a figure with subplots
    fig, axes = plt.subplots(1, len(images), figsize=(6 * len(images), 5))

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

if __name__ == "__main__":
    combine_plots()
