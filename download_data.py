"""Download the real data by accession number into data/ and check it.

Usage:
    uv run python download_data.py
"""

import hashlib
import urllib.request
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"

# accession -> MD5 of the FASTA file as served by ENA. ENA does not
# publish an MD5 for this FASTA, so it was computed once after the
# first download; it guards against truncated or changed downloads.
GENOMES = {
    "U00096.3": "805d1558950b76c15883a718da46d3e1",  # E. coli K-12 MG1655
}


def md5(path):
    """Return the MD5 checksum of a file.

    Parameters
    ----------
    path : pathlib.Path
        The file; it is read in blocks of 1 MiB.

    Returns
    -------
    str
        The checksum as 32 hexadecimal digits.
    """
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(1 << 20)  # 1 MiB
            if not block:  # b"" = end of file
                break
            digest.update(block)
    return digest.hexdigest()


def download(accession, expected_md5):
    """Download one genome from ENA unless a correct copy exists.

    Parameters
    ----------
    accession : str
        ENA accession of the genome, e.g. "U00096.3".
    expected_md5 : str
        MD5 checksum of the correct FASTA file.

    Returns
    -------
    pathlib.Path
        The FASTA file in data/.

    Raises
    ------
    RuntimeError
        If the downloaded file does not have the expected checksum.
    """
    path = DATA_DIR / f"{accession}.fasta"
    if path.exists() and md5(path) == expected_md5:
        print(f"{path.name}: already downloaded")
        return path
    url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}"
    print(f"{path.name}: downloading {url}")
    urllib.request.urlretrieve(url, path)
    if md5(path) != expected_md5:
        raise RuntimeError(
            f"{path.name}: MD5 does not match, the download is broken"
        )
    return path


def main():
    """Download every genome in GENOMES.

    Returns
    -------
    list of pathlib.Path
        The FASTA files, in the order of GENOMES.
    """
    DATA_DIR.mkdir(exist_ok=True)
    return [
        download(accession, checksum)
        for accession, checksum in GENOMES.items()
    ]


if __name__ == "__main__":
    main()
