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

## Run the various analysis tools
```
python Gdc.py --src [SOURCE FILE] --dest [OUTPUT FILE] --output [ANALYSIS TYPE] --classification [DATA CLASSIFICATION]
```

* ```[SOURCE FILE]``` - The tar file containing the dataset
* ```[OUTPUT FILE]``` - Optional. Explicitly specify output file name, otherwise a default file name is defined from the other flags.
* ```[ANALYSIS TYPE]``` - Specify the analysis output type. Three options:
    * **Volcano** - Volcano plots (p-value vs. fold-ratio) of star-counts for all available genes for ATRX Mutant vs. Wildtype, for each classification
    * **BoxPlot** - Box plots of DRG2 star-counts for ATRX Mutant vs. Wildtype, for each classification
    * **Excel**   - Microsoft Excel spreadsheet, with one sheet per classification, each giving the star-count statistics for all available genes for ATRX Mutant vs. Wildtype
* ```[DATA CLASSIFICATION]``` - How the data is classified for analysis
    * **PerDisease** - Classify the data according to NIH GDC diagnosis label
    * **BroadClasses** - Classify the NIH GDC diagnosis label more broadly as "Carcinoma", "Melanoma", "Leukemia", "None", "Unspecified" and "Other"

Examples:
```
python Gdc.py --src NIH_GDC_dataset.tar --output BoxPlot --classification PerDisease
python Gdc.py --src NIH_GDC_dataset.tar --output BoxPlot --classification BroadClasses
python Gdc.py --src NIH_GDC_dataset.tar --output Volcano --classification PerDisease
python Gdc.py --src NIH_GDC_dataset.tar --output Volcano --classification BroadClasses
python Gdc.py --src NIH_GDC_dataset.tar --output Excel   --classification PerDisease
python Gdc.py --src NIH_GDC_dataset.tar --output Excel   --classification BroadClasses
```
