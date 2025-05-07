
from pathlib import Path
import subprocess

from micall.utils.dir_path import DirPath
from .find_file import find_file
from .logger import logger


def stitch_contigs(info_file: Path, output: Path) -> None:
    assert info_file.name.endswith(".json")
    directory = DirPath(info_file.with_suffix(""))

    try:
        the_unstitched_contigs_path = find_file(directory, ".*unstitched.*contig.*[.]csv$")
    except ValueError:
        try:
            the_unstitched_contigs_path = find_file(directory, ".*contigs.*[.]csv$")
        except ValueError as ex:
            the_unstitched_contigs_path = None
            logger.warning("%s", ex)

    if the_unstitched_contigs_path:
        plot = directory / "stitcher_plot.svg"
        log = directory / "stitcher.log"
        if log.exists() and plot.exists():
            logger.debug("Run %r already stitched.", info_file.name)

        else:
            with open(log, "w") as log_writer:
                subprocess.check_call(
                    ["micall", "contig_stitcher",
                     "with-references",
                     "--debug",
                     "--plot", str(plot),
                     str(the_unstitched_contigs_path),
                     "/dev/null",
                     ],
                    stderr=log_writer,
                )

    output.touch()
