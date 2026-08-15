import hashlib
from pathlib import Path

import requests
from rich.console import Console

from linkapy.logger import setup_logger


class Linkapy_Example:
    '''
    Linkapy_Example class for downloading the example data, and generate an example command to run Linkapy.
    Note that this command requires internet access.

    :param str output: The output directory to download the data into.
    '''

    def __init__(self, output: str = 'linkapy_example') -> None:
        self.output = Path(output)
        self.output.mkdir(parents=True, exist_ok=True)
        self.url = 'https://zenodo.org/record/17117295/files/linkapy_example.tar.gz'
        self.md5 = '204bc39ea2f2442f1f6dc7a993f6e44c'

        # Set up log
        self.logfile = self.output / 'linkapy_example.log'
        self.logger = setup_logger(self.logfile, verbose=False)
        
        console = Console()
        console.rule("[bold green]Linkapy Parser[/bold green]")
        self.logger.info(f"Logging under {self.logfile}")

        # Download file
        self._download_and_verify()
        console.print("\n[bold]Generate mudata/matrices with: [/bold]\n")
        console.print("linkapy parsing \\")
        console.print(f"    -m {self.output / 'test_data'} \\")
        console.print(f"    -t {self.output / 'test_data'} \\")
        console.print(f"    -o {self.output / 'linkapy_out'} \\")
        console.print("    -j 4 \\")
        console.print(f"    -r {self.output / 'test_data' / 'genes.bed'} \\")
        console.print(f"    -r {self.output / 'test_data' / 'promoters.bed'} \\")
        console.print("    --methylation_pattern *GCHN*tsv.gz \\")
        console.print("    --methylation_pattern_names ATAC \\")
        console.print("    --methylation_pattern *WCGN*tsv.gz \\")
        console.print("    --methylation_pattern_names METH \\")
        console.print("    --transcriptome_pattern *rna.tsv")

    def _download_and_verify(self):
        '''
        Download the actual data.
        '''
        import tarfile
        self.logger.info(f"Downloading example data from zenodo to {self.output.resolve()}")
        response = requests.get(self.url)
        tarball_path = self.output / 'linkapy_example.tar.gz'
        with open(tarball_path, 'wb') as f:
            f.write(response.content)
        
        self.logger.info("Download finished, verifying md5")
        md5_hash = hashlib.md5()
        with open(tarball_path, "rb") as f:
            # Read and update hash string value in blocks of 4K
            for byte_block in iter(lambda: f.read(4096), b""):
                md5_hash.update(byte_block)
        if md5_hash.hexdigest() != self.md5:
            self.logger.error("MD5 checksum does not match, Exiting.")
            raise ValueError("MD5 checksum does not match, something went wrong during download.")
        self.logger.info(f"Unpacking tarball into {self.output / 'test_data'}")
        with tarfile.open(tarball_path, "r:gz") as tar:
            tar.extractall(path=self.output)
        Path(tarball_path).unlink()
        self.logger.info(f"Example data is ready under {self.output}")        
