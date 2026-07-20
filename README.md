## Setup

First, install (keeping simple as of now): `pip install -r requirements.txt`

Fill in the correct file paths in `config.json`. E.g.:

```json
{
  "paths": {
    "genes": "~/Documents/ChIPAnalysis/MabATCC19977_gff.gff",
    "peaks": "~/Documents/ChIPAnalysis/SigHP1/SigHP1_FDR0.01combo.gff",
    "output_path": "~/Documents/ChIPAnalysis/regulated_peaks_out.csv"
  },
  "analysis": {
    "proximity": 600,
    "operon_gap": 30
  }
}
```

`proximity`: the maximum threshold of distance between a peak and gene to consider that peak to regulate that gene
`operon_gap`: maximum threshold of distance between two genes to consider them an operon

EXAMPLE folder structure: in folder like this:

```
ChIPAnalysis/
│   config.json  <-- fill in paths here if different
│   MabATCC19977_gff.gff
├───SigHP1
│       SigHP1_FDR0.01combo.gff
├───SigHP2
│       SigHP2_FDR0.01combo.gff
```

## To run the analysis:
```bash
cd ChIPAnalysis/
python main.py
```

Which will save the output to the output path you specified in `config.json`.

## To run the UI page:
```sh
cd ChIPAnalysis/
python -m http.server 8000
```
Go to your browser: http://127.0.0.1:8000/UI/
