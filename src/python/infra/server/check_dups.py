import os
import hashlib

import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Check for duplicate BAM files across USB drives.")
    parser.add_argument("--server", type=str, default="s115", help="Server directory name (e.g., s115)")
    parser.add_argument("--usb_dirs", nargs="+", default=["/mnt/usb", "/mnt/usb-2", "/mnt/usb-3"], help="List of USB mount points")
    return parser.parse_args()

def build_file_map(mounts, server_subfolder):
    print(f"Scanning USB drives {mounts} for files in folder '{server_subfolder}'...")
    file_map = {} # {filename: [path1, path2, ...]}
    
    for mount in mounts:
        # Construct the full path where nextflow copies files
        search_path = os.path.join(mount, server_subfolder)
        
        if os.path.exists(search_path):
            print(f"Scanning {search_path}...")
            # Simple listdir since nextflow structure is flat under server_subfolder
            try:
                files = os.listdir(search_path)
                for file in files:
                    if file.endswith(".bam"): 
                        full_path = os.path.join(search_path, file)
                        if file in file_map:
                            file_map[file].append(full_path)
                        else:
                            file_map[file] = [full_path]
            except OSError as e:
                print(f"Error accessing {search_path}: {e}")
                
    return file_map

def find_duplicates():
    args = parse_args()
    file_map = build_file_map(args.usb_dirs, args.server)
    
    total_dupes = 0
    dupe_size = 0
    
    print("\n--- Duplicate Report ---")
    duplicate_found = False
    for filename, paths in file_map.items():
        if len(paths) > 1:
            duplicate_found = True
            total_dupes += 1
            # Get size of one copy
            try:
                size = os.path.getsize(paths[0])
                dupe_size += size * (len(paths) - 1) # Count extra copies
                print(f"DUPLICATE: {filename}")
                for p in paths:
                    print(f"  - {p}")
            except OSError:
                pass
    
    if not duplicate_found:
        print("No duplicates found.")

    print(f"\nTotal Duplicates found: {total_dupes}")
    print(f"Wasted Space: {dupe_size / (1024**3):.2f} GB ({dupe_size / (1024**4):.2f} TB)")

if __name__ == "__main__":
    find_duplicates()
