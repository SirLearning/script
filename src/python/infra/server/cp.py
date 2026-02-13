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
        # Create destination directory if it doesn't exist
        dst_dir = os.path.dirname(dst)
        if not os.path.exists(dst_dir):
            os.makedirs(dst_dir, exist_ok=True)
            
        print(f"Copying {src} ({size/1024/1024:.2f} MB) -> {dst}")
        subprocess.run(["cp", "-ax", src, dst], check=True)
        return (src, dst, size, "SUCCESS")
    except Exception as e:
        print(f"Failed to copy {src} to {dst}: {e}")
        return (src, dst, size, f"FAILED: {str(e)}")

def distribute_files_to_usbs(files, usb_dirs, safety_margin_mb=100):
    """
    Distribute files to USB drives based on available space.
    Returns: (tasks_list, failed_list)
    """
    # 1. Gather file sizes
    file_sizes = []
    for f in files:
        if os.path.exists(f):
            file_sizes.append((f, get_file_size(f)))
        else:
            print(f"Warning: File not found {f}")

    # 2. Sort by size descending
    file_sizes.sort(key=lambda x: x[1], reverse=True)

    # 3. Allocate files to drives (Logically)
    # We maintain a 'virtual' free space tracker to avoid race conditions
    drive_free_space = {}
    for d in usb_dirs:
        if os.path.exists(d):
            drive_free_space[d] = get_available_space(d)
        else:
             print(f"Warning: USB dir not found {d}")

    copy_tasks = [] #(src, dst, size)
    failed_assignments = []

    for file_path, size in file_sizes:
        assigned = False
        for drive in usb_dirs:
            if drive not in drive_free_space: continue
            
            # Use safety margin
            if drive_free_space[drive] > (size + safety_margin_mb * 1024 * 1024):
                dest_path = os.path.join(drive, os.path.basename(file_path))
                copy_tasks.append((file_path, dest_path, size))
                
                # Update virtual free space
                drive_free_space[drive] -= size
                assigned = True
                break
        
        if not assigned:
            failed_assignments.append((file_path, size))
            print(f"ERROR: No space for {file_path} ({size/1024/1024:.2f} MB)")
            
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
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(copy_file_task, task) for task in copy_tasks]
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())

    # 5. Write Log
    with open(log_file, "w") as log:
        for src, dst, size, status in results:
            log.write(f"{src}\t{dst}\t{size}\t{status}\n")
        for src, size in failed_assignments:
            log.write(f"{src}\tFAILED_NO_SPACE\t{size}\tFAILED\n")
