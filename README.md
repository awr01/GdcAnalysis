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
where the `NIH_GDC_dataset.tgz` can be any valid filename. This will take ~3 hours and produce a 1.5GB file.

The downloader keeps a cache in the `.cache` folder in case it is interrupted. If an error occurs on a cached operation, you can try deleting the cached file being used, or failing that, the entire `.cache` folder.

## Run the DRG2 box-plot analysis
```
python GdcBoxPlots.py --src NIH_GDC_dataset.tgz --dest DUMMY.tsv --mutations ATRX
```
This produces combined scatter and box+whisker plots of the DRG2 star-counts for ATRX mutant vs ATRX wild-type in different cancer-types saved to the file DRG2.pdf.
