# NBLAST_Scripts [![LinkToJanelia](../image/HHMI_Janelia_Color_Alternate_180x40.png)](https://www.janelia.org)
These scripts allow users to create NBLAST databases and to perform NBLAST clustering without any coding.

## Precompiled Binaries
https://github.com/JaneliaSciComp/NBLAST_Scripts/releases

## Prerequisite
 1. R installation from here (https://cran.cnr.berkeley.edu/) 

## Installation
 1. Download the latest precompiled binaries from https://github.com/JaneliaSciComp/NBLAST_Scripts/releases.
 2. Decompress and copy all files into `Fiji.app/plugins` directory. The .R scrips are needed to place directly under the Plugins folder.

## Create NBLAST Database (Also adding neurons into existing NBLAST database)
1. Launch Fiji.
1. Run `Plugins>build nblastdb`.
1. Choose an input directory. This plugin supports swc and nrrd file formats.
1. Choose a location to save the database.
1. Set parameters.
   1. **RScript**: Set the file path to the Rscript. The default path to the RScript is `C:\Program Files\R\R-x.x.x\bin\RScript.exe` (Windows) or `/usr/local/bin/RScript` (Mac and Linux).
   1. **Resample**: A new length to which all segmented edges will be resampled.
   1. **K**: Number of nearest neighbours to use for tangent vector calculation.
1. You can use the generated NBLAST database on VVDViewer.

## NBLAST Search
1. Launch Fiji.
1. Run `Plugins>nblast search`.
1. Choose an input directory. This plugin supports swc and nrrd file formats.
1. Choose a location of the NBLAST database.
1. Choose an output directory.
1. Set parameters.
   1. **RScript**: Set the file path to the Rscript. The default path to the RScript is `C:\Program Files\R\R-x.x.x\bin\RScript.exe` (Windows) or `/usr/local/bin/RScript` (Mac and Linux).
   1. **Resample**: A new length to which all segmented edges will be resampled.
   1. **K**: Number of nearest neighbours to use for tangent vector calculation.
   1. **Scoring Method**: Scoring method for NBLAST. The "mean" score is an avarage of normalized forward and reverse scores.
   1. **Number of Results**: Maximum number of search results.
   1. **Threads**: Maximum thread number for NBLAST search.
1. You can import the results on VVDViewer.

## Perform NBLAST Clustering
### Generate Score Matrix
1. Launch Fiji.
1. Run `Plugins>nblast scoremat nlfh`.
1. Choose a NBLAST database.
1. Set parameters.
   1. **RScript**: Set the file path to the Rscript. The default path to the RScript is `C:\Program Files\R\R-x.x.x\bin\RScript.exe` (Windows) or `/usr/local/bin/RScript` (Mac and Linux).
   1. **Threads**: Maximum thread number for NBLAST score matrix calculation.
1. The score matrix will be saved in the selected NBLAST database directory.
### Run NBLAST Clustering
1. Launch Fiji.
1. Run `Plugins>nblast clustering`.
1. Choose an input score matrix.
1. Choose an output directory.
1. Set parameters.
   1. **RScript**: Set the file path to the Rscript. The default path to the RScript is `C:\Program Files\R\R-x.x.x\bin\RScript.exe` (Windows) or `/usr/local/bin/RScript` (Mac and Linux).
   1. **Method**: Clustering method.
   1. **N Clusters**: Number of clusters.
   1. **Height cutoff value**: Height at which to cut tree.
   1. **"Set number of clusters"**: If you choose this, `N Clusters` will be used for clustering. (`Height cutoff value` will be ignored.)
   1. **"Use hight cutoff value"**: If you choose this, `Height cutoff value` will be used for clustering. (`N Clusters` will be ignored.)
