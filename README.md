# NBLAST_Scripts

These scripts allow users to create NBLAST databases and to perform NBLAST clustering without any coding.

## Precompiled Binaries
https://github.com/JaneliaSciComp/NBLAST_Scripts/releases

## Installation
 1. Download the latest precompiled binaries from https://github.com/JaneliaSciComp/NBLAST_Scripts/releases.
 2. Decompress and copy all files into `Fiji.app/plugins` directory.

## Create NBLAST Database
1. Launch Fiji.
1. Run `Plugins>build nblastdb`.
1. Choose an input directory. This plugin support swc and nrrd file formats.
1. Choose a location to save the database.
1. Set parameters.
   1. **RScript**: Set the file path to the Rscript. The default path to the RScript is `C:\Program Files\R\R-x.x.x\bin\RScript.exe` (Windows) or `/usr/local/bin/RScript` (Mac and Linux).
   1. **Resample**: A new length to which all segmented edges will be resampled.
   1. **K**: Number of nearest neighbours to use for tangent vector calculation.
1. You can use the generated NBLAST database on VVDViewer.

## Perform NBLAST Clustering
### Generate Score Matrix
1. Launch Fiji.
1. Run `Plugins>nblast scoremat`.
1. Choose an input directory. This plugin support swc and nrrd file formats.
1. Choose a location to save the score matrix.
1. Set parameters.
   1. **RScript**: Set the file path to the Rscript. The default path to the RScript is `C:\Program Files\R\R-x.x.x\bin\RScript.exe` (Windows) or `/usr/local/bin/RScript` (Mac and Linux).
   1. **Resample**: A new length to which all segmented edges will be resampled.
   1. **K**: Number of nearest neighbours to use for tangent vector calculation.
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
