#NeuroSynth ROI Extraction

A Python workflow for converting thresholded NeuroSynth or other statistical NIfTI maps into atlas-defined regions of interest (ROIs).

The script identifies connected statistical clusters, intersects them with anatomical atlases, separates disconnected ROI components, and exports individual binary ROI masks and CSV summary tables.

The resulting masks can be used in downstream neuroimaging analyses with tools such as FSL, Nilearn, SPM, and CONN.

The GitHub version makes the workflow easier to configure and reuse while preserving the extraction methodology and findings of the original working script.

## Overview

The workflow:

1. Loads a 3D thresholded statistical NIfTI image.
2. Identifies positive-valued voxels.
3. Detects connected clusters using configurable voxel connectivity.
4. Removes clusters below a minimum voxel-count threshold.
5. Resamples anatomical atlas label images to the input image grid.
6. Intersects retained clusters with atlas-defined regions.
7. Splits disconnected ROI components.
8. Filters ROI components by minimum size and cluster fraction.
9. Assigns anatomical labels and hemisphere classifications.
10. Ranks ROI components by voxel count.
11. Calculates ROI centroid coordinates.
12. Saves binary NIfTI masks and atlas-specific CSV summaries.

**Important:** The script assumes that the input statistical map has already been thresholded. It does not perform statistical thresholding, FDR correction, or generation of the original statistical map.

## Atlases

The workflow uses three FSL maximum-probability atlas images:

- Harvard–Oxford Cortical Atlas
- Harvard–Oxford Subcortical Atlas
- Cerebellar Atlas (MNI FNIRT)

The default atlas images use the 25% threshold variants where applicable.

Anatomical labels are read from the corresponding FSL XML files.

The script accounts for the relationship between XML label indices and integer label values in the maximum-probability atlas images.

## Requirements

- Python 3.9 or later recommended
- NumPy
- SciPy
- NiBabel
- Nilearn
- FSL atlas files

FSL itself is not called as a command-line program by the script, but the required atlas images and XML files must be available.

### Installation

Create and activate a Python environment if necessary.

Install the dependencies:

```bash
pip install -r requirements.txt
```

The `requirements.txt` file should contain:

```text
numpy
scipy
nibabel
nilearn
```

For reproducible analyses, consider recording the package versions used in your validated environment.

## Atlas configuration

By default, the script looks for atlas files under:

```text
$FSLDIR/data/atlases/
```

The expected files are:

```text
atlases/
├── HarvardOxford-Cortical.xml
├── HarvardOxford-Subcortical.xml
├── HarvardOxford/
│   ├── HarvardOxford-cort-maxprob-thr25-2mm.nii.gz
│   └── HarvardOxford-sub-maxprob-thr25-2mm.nii.gz
├── Cerebellum_MNIfnirt.xml
└── Cerebellum/
    └── Cerebellum-MNIfnirt-maxprob-thr25-2mm.nii.gz
```

If FSL is installed in a nonstandard location, specify the atlas directory using `--atlas_dir`.

For example:

```bash
python extract_rois.py \
    --input map.nii.gz \
    --outdir results \
    --atlas_dir /path/to/fsl/data/atlases
```

The specified directory must contain the expected atlas files.

## Usage

### Basic example

```bash
python extract_rois.py \
    --input motor_uniformity-test_z_FDR_0.01.nii.gz
```

The default output directory is `outputs/`.

### Specify an output directory

```bash
python extract_rois.py \
    --input map.nii.gz \
    --outdir results
```

### Adjust cluster and ROI thresholds

```bash
python extract_rois.py \
    --input map.nii.gz \
    --outdir results \
    --min_cluster_size 100 \
    --min_roi_size 50 \
    --min_cluster_fraction 0.02
```

### Change cluster connectivity

```bash
python extract_rois.py \
    --input map.nii.gz \
    --outdir results \
    --connectivity 26
```

Supported connectivity values are:

- `6`: Voxels sharing a face.
- `18`: Voxels sharing a face or edge.
- `26`: Voxels sharing a face, edge, or corner.

The same connectivity setting is used for initial cluster detection and separation of disconnected ROI components.

The default is 26-neighbour connectivity, matching the original extraction workflow.

### Disable the cluster-fraction filter

```bash
python extract_rois.py \
    --input map.nii.gz \
    --outdir results \
    --min_cluster_fraction 0
```

### Command-line parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | Required | Input 3D thresholded NIfTI image |
| `--outdir` | `outputs` | Output directory |
| `--atlas_dir` | `$FSLDIR/data/atlases/` | FSL atlas directory |
| `--connectivity` | `26` | Cluster connectivity: 6, 18, or 26 |
| `--min_cluster_size` | `100` | Minimum cluster size in voxels |
| `--min_roi_size` | `50` | Minimum ROI component size in voxels |
| `--min_cluster_fraction` | `0.02` | Minimum ROI component fraction of its parent cluster |
| `--bilateral_ratio` | `0.30` | Minimum hemisphere voxel ratio for bilateral classification |
| `--midline_tolerance` | `2.0` | Central-region threshold in mm based on mean MNI X |
| `--clean_output` | Disabled | Delete the existing output directory before processing |

All voxel-count thresholds are expressed in voxels, not mm³.

The `--clean_output` option permanently deletes the specified output directory before processing. Use it only when you are certain that the directory contains no files you need to preserve.

Without `--clean_output`, files from previous runs may remain in an existing output directory if they are not overwritten. Using a separate output directory for each analysis is recommended.

## Output structure

The script produces the following directory structure:

```text
outputs/
├── clusters_combined.nii.gz
├── cortical/
│   ├── cortical_summary.csv
│   └── masks/
├── subcortical/
│   ├── subcortical_summary.csv
│   └── masks/
└── cerebellum/
    ├── cerebellum_summary.csv
    └── masks/
```

### Combined cluster map

`clusters_combined.nii.gz` contains the retained statistical clusters.

Each cluster is assigned a unique positive integer identifier.

- `0`: Background
- `1`, `2`, `3`, etc.: Retained clusters

Clusters smaller than `--min_cluster_size` are excluded.

### ROI masks

Individual binary NIfTI masks are saved under each atlas directory.

Each mask contains:

- `1`: Voxels belonging to the ROI component.
- `0`: Voxels outside the component.

Filenames identify the parent cluster, component rank, atlas region ID, and anatomical region name.

For example:

```text
cluster_001_rank_01_region_007_Precentral_Gyrus_L.nii.gz
```

The atlas region ID included in the filename helps identify the corresponding anatomical label.

Output filenames may differ from those generated by the original script, even when the extracted ROI masks and scientific findings are identical.

### Summary tables

Each atlas produces a CSV summary:

- `cortical_summary.csv`
- `subcortical_summary.csv`
- `cerebellum_summary.csv`

The tables contain the following fields:

| Column | Description |
|---|---|
| `ClusterID` | Parent statistical cluster identifier |
| `Atlas` | Atlas used for extraction |
| `RegionRank` | Component rank within the parent cluster and atlas |
| `PrimaryRegion` | Whether the component is ranked first |
| `RegionID` | Integer label value in the atlas image |
| `RegionName` | Anatomical region name |
| `VoxelCount` | Number of voxels in the ROI component |
| `PercentOfCluster` | Percentage of the parent cluster represented by the component |
| `ClusterSize` | Number of voxels in the parent cluster |
| `CentroidX` | Mean world-space X coordinate in mm |
| `CentroidY` | Mean world-space Y coordinate in mm |
| `CentroidZ` | Mean world-space Z coordinate in mm |
| `Filename` | Corresponding ROI mask filename |

## ROI extraction and ranking

Each retained statistical cluster is intersected with the atlas labels.

For each atlas region, disconnected components are identified separately using the selected connectivity setting.

Components below the minimum ROI size or minimum cluster fraction are excluded.

Remaining components are ranked by voxel count in descending order within each parent cluster and atlas.

Rank 1 is designated the primary region for organizational purposes. This designation does not imply statistical significance or biological importance.

Voxels outside labelled atlas regions are not assigned to an atlas-defined ROI. Consequently, the sum of ROI voxel counts may be smaller than the parent cluster size.

## Hemisphere assignment

For atlas labels that do not already specify a hemisphere, laterality is estimated from the component's world-space X coordinates.

Possible suffixes include:

- `_L`: Predominantly left-sided component.
- `_R`: Predominantly right-sided component.
- `_B`: Component with substantial representation in both hemispheres.
- `_C`: Central component according to the mean-X criterion.

Atlas labels that already specify a hemisphere retain their original names.

The central-region threshold and bilateral ratio can be adjusted using `--midline_tolerance` and `--bilateral_ratio`.

Laterality is a geometric classification and should not be interpreted as a separate anatomical segmentation procedure.

## Centroid coordinates

Centroids are calculated by averaging the world-space coordinates of the voxels in each ROI component.

Coordinates are reported in millimeters according to the input image's NIfTI affine transformation.

They represent MNI coordinates only when the input image is appropriately aligned to MNI space.

## Reproducibility and validation

The GitHub version was developed from an existing extraction script that had already been used to generate the original results.

The purpose of the updates was to improve usability, configuration, input validation, and documentation without changing the intended extraction methodology.

The updated script retains the principal processing steps:

1. Selecting positive-valued voxels from an already-thresholded statistical map.
2. Identifying connected clusters using 26-neighbour connectivity by default.
3. Excluding clusters smaller than the specified minimum size.
4. Resampling atlas label images to the input image grid using nearest-neighbour interpolation.
5. Intersecting clusters with atlas-defined regions.
6. Separating disconnected ROI components.
7. Applying the minimum ROI size and cluster-fraction filters.
8. Assigning anatomical labels and laterality.
9. Calculating voxel counts and centroid coordinates.
10. Exporting binary ROI masks and CSV summary tables.

**The updated workflow is intended to preserve the findings of the original extraction.** Output filenames, directory organization, and some implementation details may differ.

For reproducible results, use the same:

- Input statistical NIfTI image.
- FSL atlas images and XML label files.
- Cluster connectivity.
- Minimum cluster size.
- Minimum ROI size.
- Minimum cluster fraction.
- Laterality parameters.

When comparing different versions of the script, the most informative checks are:

- Voxel-by-voxel agreement between the combined cluster maps.
- Agreement between corresponding binary ROI masks.
- Agreement in ROI voxel counts.
- Agreement in anatomical labels and cluster assignments.
- Agreement in centroid coordinates, allowing for any rounding differences.

Differences in filenames alone do not indicate differences in the extracted data.

For formal validation, compare the outputs generated by both versions using the same input image, atlas files, and parameter values.

## Methods description

Thresholded statistical maps were represented as binary masks by selecting positive-valued voxels. Connected clusters were identified using 26-neighbour voxel connectivity, and clusters smaller than the specified minimum size were excluded.

The remaining clusters were intersected with Harvard–Oxford cortical, Harvard–Oxford subcortical, and cerebellar maximum-probability atlases. Atlas images were resampled to the input image grid using nearest-neighbour interpolation.

Each cluster was subdivided into atlas-defined regions, and disconnected components were extracted separately. Components below the minimum ROI size or minimum parent-cluster fraction were discarded.

Retained components were ranked by voxel count, assigned anatomical labels and laterality where applicable, and saved as individual binary masks.

Atlas-specific CSV files contain ROI sizes, parent-cluster information, centroid coordinates, and corresponding filenames.

When reporting an analysis, specify the actual parameter values used, including connectivity and the minimum cluster and ROI sizes.

## Methodological considerations

### Input thresholding

The script selects voxels with values greater than zero.

It assumes that the input map has already undergone the intended statistical thresholding procedure. The script does not independently perform FDR correction or verify the threshold indicated by a filename.

Negative-valued effects are excluded by the positive-voxel selection rule.

### Spatial alignment

Atlas images are resampled to the input image grid using nearest-neighbour interpolation, which preserves discrete label values.

Resampling does not perform anatomical registration. The input image and atlases must already be in compatible anatomical coordinate spaces.

### Atlas coverage

Voxels outside labelled atlas regions are not assigned an atlas-defined ROI. Therefore, the sum of ROI voxel counts may be smaller than the parent cluster size.

### Threshold dependence

Minimum cluster size, minimum ROI size, and minimum cluster fraction influence which regions are retained.

These parameters should be documented and kept consistent when comparing analyses.

### Interpretation

The workflow provides anatomical characterization of statistical clusters. It does not independently establish statistical significance or demonstrate that a region is functionally specific to a task.

## Citation

Please cite the original sources of the statistical maps, the atlases, and the software used in the analysis.

For NeuroSynth, see:

Yarkoni, T., Poldrack, R. A., Nichols, T. E., Van Essen, D. C., & Wager, T. D. (2011). Large-scale automated synthesis of human functional neuroimaging data. *Nature Methods*, 8, 665–670.

https://doi.org/10.1038/nmeth.1635

Also consult the relevant publications and documentation for the Harvard–Oxford atlas, the cerebellar atlas, FSL, NiBabel, Nilearn, NumPy, and SciPy.

Check the licensing terms of the atlases separately before redistributing atlas files or derived materials.
