# GdcAnalysis
Tools for analysing NIH GDC data

# To get started

## On the commandline, clone the repo from GitHub:
```
git clone https://github.com/awr01/GdcAnalysis.git
```

## On the commandline, set up python:
```
pip install tqdm requests numpy scipy openpyxl
```

## Download the raw NIH GDC dataset:
```
python GdcDownload.py --dest NIH_GDC_dataset.tgz
```
where the `NIH_GDC_dataset.tgz` can be any valid filename. This will take ~3 hours, download ~12GB of raw data and produce a 1.5GB output file.

The downloader keeps a cache in the `.cache` folder in case it is interrupted. If an error occurs on a cached operation, you can try deleting the cached file being used, or failing that, the entire `.cache` folder.

Once the file is successfully downloaded, the `.cache` folder can be deleted to free 12GB of disk space.

## Run the co-mutation analysis
```
python GdcAnalysis.py --src NIH_GDC_dataset.tgz --dest Analysis.tsv --mutations ATRX
```
or
```
python GdcAnalysis.py --src NIH_GDC_dataset.tgz --dest Analysis.xlsx --mutations ATRX
```
This produces a tab-delimited text file or Microsoft Excel spreadsheet with the co-mutation statistics for different cancer-types.

## Run the volcano plot analysis, highlighting DRG2
```
python GdcVolcano.py --src NIH_GDC_dataset.tgz --dest DUMMY.tsv --mutations ATRX
```
This produces volcano plots of the statistical significance of star-counts vs. the log1.5(fold-ratio) of ATRX mutant over ATRX wild-type in different cancer-types, specifically highlighting DRG2, saved to the file `DRG2-volcano.pdf`.

## Run the DRG2 box-plot analysis
```
python GdcBoxPlots.py --src NIH_GDC_dataset.tgz --dest DUMMY.tsv --mutations ATRX
```
This produces combined scatter and box+whisker plots of the DRG2 star-counts for ATRX mutant vs ATRX wild-type in different cancer-types saved to the file `DRG2-boxplot.pdf`.
