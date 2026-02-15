import os
import shutil
import subprocess
import concurrent.futures

def get_available_space(path):
    """Returns available space in bytes for the given path."""
    try:
        total, used, free = shutil.disk_usage(path)
        return free
    except OSError:
        return 0

def get_file_size(path):
    """Returns file size in bytes."""
    try:
        path = os.path.realpath(path)
        return os.path.getsize(path)
    except OSError:
        return 0

def copy_file_task(task):
    """
    Task to copy a single file. (Helper for ThreadPoolExecutor)
    task: tuple of (src, dst, size)
    """
    src, dst, size = task
    try:
        # Check if file already exists at destination with same size
        if os.path.exists(dst):
            dst_size = os.path.getsize(dst)
            # Allow small variance or exact match. For strictness use exact match.
            if dst_size == size:
                print(f"Skipping {src} -> {dst} (Already exists)")
                return (src, dst, size, "SKIPPED")
            else:
                 print(f"Overwriting {src} -> {dst} (Size mismatch: {dst_size} vs {size})")

        # Create destination directory if it doesn't exist
        dst_dir = os.path.dirname(dst)
        if not os.path.exists(dst_dir):
            os.makedirs(dst_dir, exist_ok=True)
            
        print(f"Copying {src} ({size/1024/1024:.2f} MB) -> {dst}")
        # Use -L to dereference symlinks (copy content, not link) because USB drives (exFAT/FAT32) usually don't support symlinks.
        subprocess.run(["cp", "-L", src, dst], check=True)
        return (src, dst, size, "SUCCESS")
    except Exception as e:
        print(f"Failed to copy {src} to {dst}: {e}")
        return (src, dst, size, f"FAILED: {str(e)}")

def get_group_size(item):
    """Returns size of a file or total size of a list of files."""
    if isinstance(item, str):
        return get_file_size(item)
    elif isinstance(item, (list, tuple)):
        return sum(get_end_file_size(f) for f in item)
    return 0

def get_end_file_size(f):
    """Helper to get size of individual file"""
    if os.path.exists(f):
        return get_file_size(f)
    print(f"Warning: File not found {f}")
    return 0

def distribute_files_to_usbs(files, usb_dirs, safety_margin_mb=100):
    """
    Distribute files to USB drives based on available space.
    'files' can be a list of file paths (strs) or a list of groups (lists of strs).
    Returns: (tasks_list, failed_list)
    """
    # 1. Gather sizes
    items_with_size = []
    for item in files:
        size = get_group_size(item)
        if size > 0:
            items_with_size.append((item, size))

    # 2. Sort by size descending
    items_with_size.sort(key=lambda x: x[1], reverse=True)

    # 3. Allocate files to drives (Logically)
    # We maintain a 'virtual' free space tracker to avoid race conditions
    drive_free_space = {}
    for d in usb_dirs:
        if os.path.exists(d):
            drive_free_space[d] = get_available_space(d)
        else:
            if not os.path.exists(d):
                try:
                    os.makedirs(d, exist_ok=True)
                    drive_free_space[d] = get_available_space(d)
                except:
                    print(f"Warning: USB dir not found/creatable {d}")
            else:
                drive_free_space[d] = get_available_space(d)

    copy_tasks = [] #(src, dst, size)
    failed_assignments = []
    
    # 3b. Check for existing files first to prevent duplicates across drives
    # If a file already exists on ANY drive, we must assign it there (to be skipped) 
    # or remove it from the list.
    
    already_done_items = []
    
    # Map to track where files are currently located
    # This is a bit expensive but necessary for resume without duplicates
    print("Checking for existing files on USB drives to prevent duplicates...")
    
    # Create a set of existing filenames on each drive for fast lookup
    existing_files_map = {} # { drive_path: set(filenames) }
    for d in usb_dirs:
        existing_files_map[d] = set()
        if os.path.exists(d):
            try:
                # Only list the top level of the target dir (no recursion for now as we copy flat or single subdir)
                # If we use subdirs like 's115', 'd' should point to that subdir based on nextflow logic
                for content in os.listdir(d):
                    existing_files_map[d].add(content)
            except OSError:
                pass

    items_to_process = []
    
    for item, size in items_with_size:
        # Item can be a single file path or a list of file paths
        sub_files = [item] if isinstance(item, str) else item
        
        # Check if ALL files in this group exist on a specific drive
        found_drive = None
        
        # We need to find a drive that contains the MAIN file (usually the bam)
        # or all files. Let's look for the first file in the group.
        first_file_name = os.path.basename(sub_files[0])
        
        for drive in usb_dirs:
            if first_file_name in existing_files_map.get(drive, set()):
                # Found the file on this drive. 
                # We assume if the first file is here, the group is intended to be here.
                found_drive = drive
                break
        
        if found_drive:
            # The file exists on 'found_drive'. We assign it there.
            # verify sizes in the actual task, but here we just assign it 
            # so the copy logic can say "SKIPPED".
            # We do NOT deduct availability from 'drive_free_space' because it's already consuming that space.
            
            for f in sub_files:
                dest_path = os.path.join(found_drive, os.path.basename(f))
                copy_tasks.append((f, dest_path, get_file_size(f)))
                
        else:
            # Not found anywhere. Needs allocation.
            items_to_process.append((item, size))

    # Now allocate ONLY the new/missing items
    # Note: items_to_process contains (item, size) for items NOT found on any drive
    for item, size in items_to_process:
        assigned = False
        for drive in usb_dirs:
            if drive not in drive_free_space: continue
            
            # Use safety margin (MB to bytes)
            if drive_free_space[drive] > (size + safety_margin_mb * 1024 * 1024):
                # Update virtual free space
                drive_free_space[drive] -= size
                assigned = True

                # Generate tasks
                sub_files = [item] if isinstance(item, str) else item
                for f in sub_files:
                    dest_path = os.path.join(drive, os.path.basename(f))
                    copy_tasks.append((f, dest_path, get_file_size(f)))
                break
        
        if not assigned:
            failed_assignments.append((str(item), size))
            print(f"ERROR: No space for item ({size/1024/1024:.2f} MB)")
            
    return copy_tasks, failed_assignments

def run_copy_process(files, usb_dirs, log_file, threads=1):
    """
    Main function to run the copy process.
    """
    copy_tasks, failed_assignments = distribute_files_to_usbs(
        files, usb_dirs
    )

    # 4. Execute copies in parallel
    results = []
    
    print(f"Starting transfer of {len(copy_tasks)} files using {threads} threads...")
    
    # Pre-check for existing files to adjust free space calculation if we were resuming?
    # No, simple resume: iterate tasks, if exists sync size, skip.
    
    with open(log_file, "a") as log: # Append mode for partial logs if we implement checkpointing
        # But here we just write at the end.
        pass

    # Real-time logging to file would be better for resume, 
    # but concurrent writes to single file might be messy without a lock.
    # We will rely on stdout (which is redirected to .log in nextflow) and checks in copy_file_task.

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(copy_file_task, task) for task in copy_tasks]
        for future in concurrent.futures.as_completed(futures):
            res = future.result()
            results.append(res)
            # Optional: print progress

    # 5. Write Log (Summary)
    # Note: If process killed, this part is not reached, hence missing log file.
    # But because we added "SKIPPED" check in copy_file_task, re-running will be fast.
    with open(log_file, "w") as log:
        for src, dst, size, status in results:
            log.write(f"{src}\t{dst}\t{size}\t{status}\n")
        for src, size in failed_assignments:
            log.write(f"{src}\tFAILED_NO_SPACE\t{size}\tFAILED\n")
