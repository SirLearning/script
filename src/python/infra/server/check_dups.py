import os
import hashlib

import shutil
import argparse

def get_available_space(path):
    """Returns available space in bytes for the given path."""
    try:
        total, used, free = shutil.disk_usage(path)
        return free
    except OSError:
        return 0

def parse_args():
    parser = argparse.ArgumentParser(description="Check for duplicate BAM files across USB drives.")
    parser.add_argument("--server", type=str, default="s115", help="Server directory name (e.g., s115)")
    parser.add_argument("--usb_dirs", nargs="+", default=["/mnt/usb", "/mnt/usb-2", "/mnt/usb-3"], help="List of USB mount points")
    parser.add_argument("--delete", action="store_true", help="Delete duplicate files (keep one copy)")
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
    
    if args.delete and total_dupes > 0:
        print("\n--- Deleting Duplicates ---")
        print("Strategy: Keep the copy on the drive with the MOST free space.")
        
        for filename, paths in file_map.items():
            if len(paths) > 1:
                # Determine which copy to keep
                # We want to keep the one on the drive with the most space to balance usage?
                # Or keep the one on the drive with the LEAST space to fill it up?
                # Wait, usually we want to FREE space on full drives.
                # So we KEEP the file on the empty drive, DELETE from full drive.
                
                path_spaces = []
                for p in paths:
                    # simplistic assumption about mount point being 2 levels up or just check disk usage of folder
                    mount_point = os.path.dirname(os.path.dirname(p)) # /mnt/usb/s115 -> /mnt/usb
                    free = get_available_space(mount_point)
                    path_spaces.append((p, free))
                
                # Sort by free space descending. Keep the first one (most free space).
                # This moves data to the emptier drive effectively by deleting from fuller ones.
                path_spaces.sort(key=lambda x: x[1], reverse=True)
                
                keep_tuple = path_spaces[0]
                keep_path = keep_tuple[0]
                keep_space = keep_tuple[1]
                
                delete_tuples = path_spaces[1:]
                
                print(f"Keeping: {keep_path} (Drive has {keep_space/1024**3:.2f} GB free)")
                
                for p_tup in delete_tuples:
                    p = p_tup[0]
                    try:
                        print(f"  Deleting: {p}")
                        os.remove(p)
                        
                        # Also delete associated .bai and .md5
                        # 1. foo.bam.bai / foo.bam.md5
                        for ext in [".bai", ".md5"]:
                            assoc_file = p + ext
                            if os.path.exists(assoc_file):
                                os.remove(assoc_file)
                                print(f"    Deleted associated: {assoc_file}")
                            
                            # 2. foo.bai / foo.md5 (if different)
                            assoc_file_alt = os.path.splitext(p)[0] + ext
                            if assoc_file_alt != assoc_file and os.path.exists(assoc_file_alt):
                                os.remove(assoc_file_alt)
                                print(f"    Deleted associated: {assoc_file_alt}")

                    except OSError as e:
                        print(f"    Error deleting {p}: {e}")
            
    elif total_dupes > 0:
        print("\nTo delete duplicates and free up space, run with --delete")


if __name__ == "__main__":
    find_duplicates()
