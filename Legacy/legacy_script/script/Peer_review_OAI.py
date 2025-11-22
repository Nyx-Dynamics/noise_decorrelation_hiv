# Create a single ZIP archive with all the key files and also copy the manuscript to a simpler filename
import os, shutil, zipfile, pathlib

base = "/mnt/data"
files = [
    "Noise_Mediated_Neuroprotection-3.pdf",
    "IMPLEMENTATION_SUMMARY.md",
    "README_v2.md",
    "chronic_NAA_underprediction_analysis.md",
    "EVOLUTIONARY_FRAMEWORK.md",
]

# Verify files exist and collect absolute paths
abs_files = []
for f in files:
    p = os.path.join(base, f)
    if os.path.exists(p):
        abs_files.append(p)

# Copy manuscript to a simpler filename for convenience
manuscript_src = os.path.join(base, "Noise_Mediated_Neuroprotection-3.pdf")
manuscript_dst = os.path.join(base, "Noise_Mediated_Neuroprotection_Manuscript.pdf")
if os.path.exists(manuscript_src):
    shutil.copyfile(manuscript_src, manuscript_dst)

# Create ZIP
zip_path = os.path.join(base, "Neuroprotection_Manuscript_and_Notes.zip")
with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zf:
    for f in abs_files:
        zf.write(f, arcname=os.path.basename(f))
    # Also include the simplified manuscript copy if created
    if os.path.exists(manuscript_dst):
        zf.write(manuscript_dst, arcname=os.path.basename(manuscript_dst))

zip_path, os.path.exists(zip_path), manuscript_dst, os.path.exists(manuscript_dst)
