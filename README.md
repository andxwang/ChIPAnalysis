## Setup

First, install (keeping simple as of now): `pip install -r requirements.txt`. Needed for analysis, but not UI.

## How to run analysis: set file paths and constants (see below) by using command line arguments OR filling in config.

### Option 1: Command line:

```sh
python3 main.py --help
```
will show 

```sh
usage: main.py [-h] [-c CONFIG] [-g GENES] [-p PEAKS] [-o OUTPUT] [--proximity PROXIMITY] [--location-edge-cutoff LOCATION_EDGE_CUTOFF]

Main script to annotate ChIP-seq peaks.

options:
  -h, --help            show this help message and exit
  -c CONFIG, --config CONFIG
                        Path to configuration JSON file (default: config.json)
  -g GENES, --genes GENES
                        Path to gene annotation GFF file
  -p PEAKS, --peaks PEAKS
                        Path to peaks GFF file
  -o OUTPUT, --output OUTPUT
                        Output CSV filename
  --proximity PROXIMITY
                        maximum threshold of distance between a peak and gene to consider that peak to regulate that gene (overrides config)
  --location-edge-cutoff LOCATION_EDGE_CUTOFF
                        bp distance from boundary of gene to consider a peak at start/end of that gene (overrides config)
```

Example command:

```bash
python3 main.py 
  -g ~/Documents/ChIPAnalysis/MabATCC19977_gff.gff 
  -p ~/Documents/ChIPAnalysis/SigHP1/SigHP1_FDR0.01combo.gff 
  -o ~/Documents/ChIPAnalysis/regulated_peaks_out.csv
  --proximity 600 
  --location-edge-cutoff 500
```


### Option 2: use config.json

Fill in the correct file paths in `config.json`. Also set **constants**. Then run `python main.py`. E.g.:

```json
{
  "paths": {
    "genes": "~/Documents/ChIPAnalysis/MabATCC19977_gff.gff",
    "peaks": "~/Documents/ChIPAnalysis/SigHP1/SigHP1_FDR0.01combo.gff",
    "output_path": "~/Documents/ChIPAnalysis/regulated_peaks_out.csv"
  },
  "analysis": {
    "proximity": 600,
    "location_edge_cutoff": 500
  }
}
```

- `proximity`: the maximum threshold of distance between a peak and gene to consider that peak to regulate that gene
- `location_edge_cutoff`: bp distance from boundary of gene to consider a peak at start/end of that gene


## To run the UI page:
```sh
cd ChIPAnalysis/
python -m http.server 8000
```
Go to your browser: http://127.0.0.1:8000/UI/
