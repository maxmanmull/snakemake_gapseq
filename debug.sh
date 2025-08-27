#!/bin/bash

# Correct script to fix thread allocations in Snakefile
# This maximizes parallelization across samples

echo "Backing up original Snakefile..."
cp snakefile snakefile.backup.$(date +%Y%m%d_%H%M%S)

echo "Fixing thread allocations..."

# Create a Python script to do the replacement correctly
cat > fix_threads.py << 'EOF'
import re

with open('snakefile', 'r') as f:
    content = f.read()

# Replace thread counts in rules (not in comments)
# Pattern matches "threads: NUMBER" not preceded by #
content = re.sub(r'^(\s+threads:\s+)16\s*$', r'\g<1>1', content, flags=re.MULTILINE)
content = re.sub(r'^(\s+threads:\s+)8\s*$', r'\g<1>1', content, flags=re.MULTILINE)
content = re.sub(r'^(\s+threads:\s+)4\s*$', r'\g<1>1', content, flags=re.MULTILINE)
content = re.sub(r'^(\s+threads:\s+)10\s*$', r'\g<1>1', content, flags=re.MULTILINE)

# Also fix thread references in shell commands
content = re.sub(r'-t 8\b', '-t 1', content)
content = re.sub(r'-t 16\b', '-t 1', content)
content = re.sub(r'-t 4\b', '-t 1', content)
content = re.sub(r'--threads 8\b', '--threads 1', content)
content = re.sub(r'--threads 4\b', '--threads 1', content)
content = re.sub(r'--cpus 8\b', '--cpus 1', content)
content = re.sub(r'-@ 8\b', '-@ 1', content)

# Special handling for GTDB-Tk batch - keep it at 32
# Find the gtdbtk_classify_batch rule and preserve its thread count
pattern = r'(rule gtdbtk_classify_batch:.*?threads:\s+)1(\s)'
content = re.sub(pattern, r'\g<1>32\g<2>', content, flags=re.DOTALL)

with open('snakefile.fixed', 'w') as f:
    f.write(content)

print("Fixed thread allocations")
EOF

python fix_threads.py

echo ""
echo "Thread allocations in original:"
grep -E '^\s+threads:' snakefile | sort | uniq -c

echo ""
echo "Thread allocations in fixed version:"
grep -E '^\s+threads:' snakefile.fixed | sort | uniq -c

echo ""
echo "Sample of changes made:"
diff -u snakefile snakefile.fixed | grep -E "threads:|^[\+\-]\s+-t |^[\+\-]\s+--threads" | head -20

echo ""
echo "Do you want to apply these changes? (y/n)"
read -r response

if [[ "$response" == "y" ]]; then
    mv snakefile.fixed snakefile
    rm fix_threads.py
    echo "✓ Snakefile updated!"
    echo ""
    echo "Now you can run:"
    echo "snakemake --use-singularity --cores 40 --jobs 40 --rerun-incomplete"
    echo ""
    echo "This will run up to 40 samples in parallel with 1 core each!"
else
    rm snakefile.fixed fix_threads.py
    echo "Changes discarded."
fi
