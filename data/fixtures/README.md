# Small input fixtures

`drive-probe.json` is the 153-byte synthetic payload written and read in Colab by `notebooks/00_drive_setup.ipynb`. It was downloaded from the mounted project Drive and checked against the printed SHA-256. GitHub Actions verifies that SHA-256 on the checked-out file. The fixture contains no electrophysiology or credentials.

`mouse_population_40ms.mat` contains derived 40 ms population counts grouped by session and state from the 13 datasets specified by Sooter et al. The authors' [Figshare deposit](https://doi.org/10.6084/m9.figshare.32012181.v1) provides the source files under CC BY 4.0. The preprocessing script and [manifest](../../manifests/mouse-preparation.json) record the bin convention, durations, provenance, and fixture digest. This fixture contains aggregated counts.

`zebra_lman_population.mat` contains derived LMAN population counts and sensitivity inputs from the eligible ZF1 and ZF2 recordings. The [Lorenz et al. data release](https://doi.org/10.12751/g-node.vsxyi8) licenses source data under CC BY 4.0. Session selection, provenance and checksums are recorded in the [preparation manifest](../../manifests/zebra-finch-preparation.json).
