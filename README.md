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
python GdcDownload.py --dest NIH_GDC_dataset.tar
```
where the `NIH_GDC_dataset.tar` can be any valid filename. This will take ~4 hours, download ~18GB of raw data and produce a 9GB output file.

The downloader keeps a cache in the `.cache` folder in case it is interrupted. If an error occurs on a cached operation, you can try deleting the cached file being used, or failing that, the entire `.cache` folder.

Once the file is successfully downloaded, the `.cache` folder can be deleted to free 18GB of disk space.

## Run the co-mutation analysis
```
python GdcAnalysis.py --src NIH_GDC_dataset.tgz --dest Analysis.xlsx --mutation ATRX
```
This produces a Microsoft Excel spreadsheet with one sheet per different cancer-types, each giving the star-count statistics for mutant and wild-type cases of the specified gene, here ATRX.

## Run the various plot analyses, highlighting DRG2
```
python GdcVolcano.py          --src NIH_GDC_dataset.tgz --dest [Target].pdf --mutation ATRX
python GdcBoxPlots.py         --src NIH_GDC_dataset.tgz --dest [Target].pdf --mutation ATRX
```
For each cancer type, these produce respectively
 - volcano plots of the statistical significance of star-counts vs. the log1.5(fold-ratio) of ATRX mutant over ATRX wild-type in different cancer-types, specifically highlighting DRG2, saved to the specified file.
 - combined scatter and box+whisker plots of the DRG2 star-counts for ATRX mutant vs ATRX wild-type in different cancer-types saved to the specified file.

```
python GdcVolcanoCombined.py  --src NIH_GDC_dataset.tgz --dest [Target].pdf --mutation ATRX
python GdcBoxPlotsCombined.py --src NIH_GDC_dataset.tgz --dest [Target].pdf --mutation ATRX
```

As above, but with broader categories, rather than specific cancer types.
