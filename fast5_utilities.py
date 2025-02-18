from pathlib import Path
import shutil
from ont_fast5_api.fast5_interface import get_fast5_file
import logging
import subprocess
import h5py
import logging
from pathlib import Path

# Define the location and name of the log file
log_file_path = Path("logs/old_brain.log")

# Create the directory if it doesn't exist
log_file_path.parent.mkdir(parents=True, exist_ok=True)

# Configure the logging module
logging.basicConfig(
    level=logging.DEBUG,  # Set the logging level (DEBUG, INFO, WARNING, etc.)
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",  # Time and log message format
    datefmt="%Y-%m-%d %H:%M:%S",  # Time format (year-month-day hour:minute:second)
    handlers=[
        logging.FileHandler(log_file_path, mode='w'),  # Log to a file. with 'w' the file is overwritten.
        logging.StreamHandler()  # Also log to the console
    ]
)

# Example of logging
logger = logging.getLogger(__name__)

def multi_to_single(in_dir: str, out_dir: str):
    subprocess.run(f'multi_to_single_fast5 --recursive -i {in_dir} -s {out_dir} -t 10', shell=True)


def remove_basecalling_from_single_fast5s(single_f5_dir: str):
    for fast5_filepath in Path(single_f5_dir).glob('**/*fast5'):
        read_id = fast5_filepath.stem
        with h5py.File(fast5_filepath, mode='r+') as f5_file:  # Open in read-write mode
            try:
                del f5_file['Analyses']
                logger.info(f"Basecalling data removed for read ID {read_id} in {fast5_filepath}")
            except KeyError as e:
                logger.info(f"Basecalling group not found for read {read_id} in file {fast5_filepath}")

def make_single_and_remove_basecalling(nanopore_dir: str):
    final_summary = next(Path(nanopore_dir).glob('*/*/final_summary*'), None)
    out_dir = Path(final_summary).parent / 'single_fast5s'
    if out_dir.exists():
        return None
    logger.info(f'{out_dir} to single and remove basecalling.')
    multi_to_single(nanopore_dir, out_dir)
    remove_basecalling_from_single_fast5s(out_dir)


def basecalling_guppyv6(input_dir: str, output_dir: str):
    input_dir_PSCA0047 = Path("/media/owner/809f403f-5b66-4d70-be53-a585528402c5/Minion_220607/201020_cas/47/20201020_0453_MN32877_FAO42372_086073be")
    if not input_dir:
        input_dir = input_dir_PSCA0047
    if not output_dir:
        output_dir = Path("test_files/guppyv6_output_PSCA0047")
    output_dir.mkdir(exist_ok=True)
    cmd = f'~/Softwares/ont-guppy_6.4.6_linux64/ont-guppy/bin/guppy_basecaller -i {input_dir} -s {output_dir} -c dna_r9.4.1_450bps_modbases_5mc_cg_sup.cfg --device cuda:0 --recursive'
    subprocess.run(cmd, shell=True)


def guppy_calling_v5(input_dir: str = "", output_dir: str = ""):
    input_dir_PSCA0047 = Path("/media/owner/809f403f-5b66-4d70-be53-a585528402c5/Minion_220607/201020_cas/47/20201020_0453_MN32877_FAO42372_086073be")
    if not input_dir:
        input_dir = input_dir_PSCA0047
    if not output_dir:
        output_dir = Path("test_files/guppyv5_fast5output_prom_PSCA0047")
    output_dir.mkdir(exist_ok=True)
    cmd = f'~/Softwares/ont-guppy_5.0.16_linux64/ont-guppy/bin/guppy_basecaller -i {input_dir} -s {output_dir} -c dna_r9.4.1_450bps_modbases_5mc_hac.cfg --device cuda:0 --recursive --fast5_out'
    subprocess.run(cmd, shell=True)


if __name__ == '__main__':
    for nanopore_dir in Path('/media/owner/809f403f-5b66-4d70-be53-a585528402c5/old_brain_data/').glob('*'):
        input_dir = Path(next(nanopore_dir.glob('*/*/single_fast5s'), None))
        basecalling_guppyv6(input_dir, input_dir.parent / 'single_fast5s_bc')