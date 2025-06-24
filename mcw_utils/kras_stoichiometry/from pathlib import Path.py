from pathlib import Path
import shutil

p = Path("/scratch/g/dseo/tempus_rna_seq/")

for f in p.iterdir():
    if f.is_dir():
        star = f / "star_index"
        if star.is_dir():
            shutil.rmtree(star)
        rsem = f / "rsem_index"
        if rsem.is_dir():
            shutil.rmtree(rsem)
