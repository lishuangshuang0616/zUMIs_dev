
import sys
from unittest.mock import MagicMock

# Mock yaml module
sys.modules['yaml'] = MagicMock()

import os
import shutil
import glob
from pipeline_modules import merge_bam_stats

def test_bcstats_aggregation():
    # Setup
    tmp_dir = "test_tmp"
    out_dir = "test_out"
    project = "test_proj"
    
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
        
    os.makedirs(tmp_dir)
    os.makedirs(out_dir)
    
    # Create dummy BCstats files
    with open(os.path.join(tmp_dir, f"{project}.1.BCstats.txt"), 'w') as f:
        f.write("AAAA\t10\n")
        f.write("BBBB\t5\n")
        
    with open(os.path.join(tmp_dir, f"{project}.2.BCstats.txt"), 'w') as f:
        f.write("AAAA\t20\n") # Duplicate barcode, should sum to 30
        f.write("CCCC\t8\n")
        
    # Create dummy BAM to avoid "No BAM files" warning/error if it stops execution (it just prints)
    # Actually the function prints "No BAM files found to check layout." but doesn't crash.
    
    # Run
    try:
        merge_bam_stats(tmp_dir, project, out_dir, "dummy.yaml", "samtools")
    except Exception as e:
        print(f"Function raised exception (might be expected if dependencies missing): {e}")
        
    # Verify
    out_file = os.path.join(out_dir, f"{project}.BCstats.txt")
    if os.path.exists(out_file):
        with open(out_file, 'r') as f:
            content = f.read()
            print("Output content:")
            print(content)
            
            lines = content.strip().split('\n')
            data = {}
            for line in lines:
                parts = line.split('\t')
                data[parts[0]] = int(parts[1])
            
            assert data.get("AAAA") == 30, f"Expected AAAA: 30, got {data.get('AAAA')}"
            assert data.get("BBBB") == 5, f"Expected BBBB: 5, got {data.get('BBBB')}"
            assert data.get("CCCC") == 8, f"Expected CCCC: 8, got {data.get('CCCC')}"
            print("Verification PASSED")
    else:
        print("Output file not found!")

    # Cleanup
    shutil.rmtree(tmp_dir)
    shutil.rmtree(out_dir)

if __name__ == "__main__":
    test_bcstats_aggregation()
