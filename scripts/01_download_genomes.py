'''
01_download_genomes.py
- Fetch candidate genomes from NCBI
- Flatten them to a single FASTA per accession for easy looping over later processes
- Pull the metadata into a neat table for later lookup

Run from project root directory with: conda run -n kp-tools python scripts/01_download_genomes.py
'''

# --- Imports ---
from pathlib import Path
import subprocess
import zipfile
import shutil
import sys

# --- Configuration settings ---
ACCESSION_LIST = Path('data/accession.txt')
RAW_DIR = Path('data/raw')
GENOME_DIR = Path('data/genomes')
METADATA_OUT = Path('data/downloaded_metadata.tsv')
ZIP_PATH = RAW_DIR/'ncbi_download.zip'


def main():
    # --- 1. Setup ---

    # Check accession list exists
    if not ACCESSION_LIST.exists():
        sys.exit(f'ERROR: {ACCESSION_LIST} not found')

    # Read file and kill any whitespace
    accessions = ACCESSION_LIST.read_text().split()
    print(f'Read {len(accessions)} accessions from {ACCESSION_LIST}')

    # Create target directories if not already present
    RAW_DIR.mkdir(parents = True, exist_ok = True)
    GENOME_DIR.mkdir(parents = True, exist_ok = True)


    # --- 2. Download from NCBI ---

    # Skip this step if zipped downloads are already present to allow quick re-runs during debugging
    if ZIP_PATH.exists() and ZIP_PATH.stat().st_size > 0:
        print(f'File archive already present ({ZIP_PATH}) - skipping download')
    else:
        print('Downloading from NCBI...')
        subprocess.run(['datasets', 'download', 'genome', 'accession', '--inputfile', str(ACCESSION_LIST),
                        '--include', 'genome', '--filename', str(ZIP_PATH)], check = True)


    # --- 3. Unpack downloads ---
    print('Unpacking downloads...')
    with zipfile.ZipFile(ZIP_PATH) as archive:
        archive.extractall(RAW_DIR)

    data_dir = RAW_DIR/'ncbi_dataset'/'data'
    if not data_dir.is_dir():
        sys.exit(f'ERROR: Unexpected NCBI record layout; expected {data_dir} after unzip.')


    # --- 4. Flatten & rename FASTA files ---

    # Each genome will be stored in data/<ACCESSION>/<ACCESSION>_<ASM>_genomic.fna
    # Each should end up being stored in GENOME_DIR as <ACCESSION>.fna
    
    print(f"Flattening into {GENOME_DIR}/")
    n_copied = 0
    for genome_dir in sorted(data_dir.glob("GCF*")):
        if not genome_dir.is_dir():
            continue
        accession = genome_dir.name                    # e.g. "GCF_009025895.1"

        # list(...glob("*.fna")) collects the FASTA(s) in this folder.
        fasta_files = list(genome_dir.glob("*.fna"))
        if not fasta_files:
            print(f"  WARNING: no .fna for {accession}, skipping")
            continue

        # Copy the first (there should only be one) to <accession>.fna.
        shutil.copy(fasta_files[0], GENOME_DIR / f"{accession}.fna")
        n_copied += 1
    print(f"Copied {n_copied} genomes")


    # --- 5. Extract metadata ---
    report = data_dir/'assembly_data_report.jsonl'
    if report.exists():
        print('Extracting metadata...')
        report_tsv = subprocess.run(['dataformat', 'tsv', 'genome', '--inputfile', str(report)],
                                    check = True, capture_output = True, text = True)
        METADATA_OUT.write_text(report_tsv.stdout)
        print(f'Wrote {METADATA_OUT}')
    else:
        print(f'WARNING: {report} not found; metadata skipped.')

    print('\nDone. Next: scripts/02_mlst.sh')


if __name__ == '__main__':
    main()