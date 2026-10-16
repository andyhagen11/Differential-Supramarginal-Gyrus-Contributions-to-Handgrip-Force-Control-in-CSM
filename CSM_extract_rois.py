#!/usr/bin/env python3
"""
Extract atlas-defined ROIs from a thresholded statistical NIfTI map.

Workflow
--------
1. Load a 3D statistical map.
2. Select positive-valued voxels.
3. Identify connected clusters.
4. Remove clusters below the minimum cluster size.
5. Resample atlas label images to the input image grid.
6. Intersect clusters with atlas-defined regions.
7. Split disconnected ROI components.
8. Filter components by voxel count and cluster fraction.
9. Assign region names and laterality.
10. Save individual ROI masks and CSV summaries.

Important
---------
The input is assumed to be an already-thresholded statistical map.
This script does not perform statistical thresholding or FDR correction.

Atlas XML indices are converted to label values by adding one, as
appropriate for FSL 3D maximum-probability summary images.

All voxel-count thresholds are expressed in voxels, not mm^3.

This version retains the extraction logic of the original script
while adding command-line configuration and input validation.
"""

import os
import csv
import argparse
import shutil
import xml.etree.ElementTree as ET

import nibabel as nib
import numpy as np

from scipy.ndimage import label, generate_binary_structure
from nilearn.image import resample_to_img


# ==========================================================
# DEFAULT ATLAS CONFIGURATION
# ==========================================================

DEFAULT_FSLDIR = os.environ.get(
    "FSLDIR",
    "/usr/local/fsl"
)


def get_atlases(atlas_dir):
    """Return atlas image and XML paths for the supported atlases."""

    return {
        "cortical": {
            "img": os.path.join(
                atlas_dir,
                "HarvardOxford",
                "HarvardOxford-cort-maxprob-thr25-2mm.nii.gz"
            ),
            "xml": os.path.join(
                atlas_dir,
                "HarvardOxford-Cortical.xml"
            ),
        },
        "subcortical": {
            "img": os.path.join(
                atlas_dir,
                "HarvardOxford",
                "HarvardOxford-sub-maxprob-thr25-2mm.nii.gz"
            ),
            "xml": os.path.join(
                atlas_dir,
                "HarvardOxford-Subcortical.xml"
            ),
        },
        "cerebellum": {
            "img": os.path.join(
                atlas_dir,
                "Cerebellum",
                "Cerebellum-MNIfnirt-maxprob-thr25-2mm.nii.gz"
            ),
            "xml": os.path.join(
                atlas_dir,
                "Cerebellum_MNIfnirt.xml"
            ),
        },
    }


# ==========================================================
# FUNCTIONS
# ==========================================================

def get_structure(connectivity):
    """Return a 3D connectivity structure: 6, 18, or 26."""

    if connectivity == 6:
        return generate_binary_structure(3, 1)

    if connectivity == 18:
        return generate_binary_structure(3, 2)

    if connectivity == 26:
        return np.ones((3, 3, 3), dtype=np.uint8)

    raise ValueError(
        "Connectivity must be 6, 18, or 26."
    )


def parse_fsl_xml(xml_file):
    """
    Parse FSL XML labels for a 3D maximum-probability atlas.

    FSL XML indices correspond to probabilistic atlas volume indices.
    The corresponding label values in a maximum-probability summary
    image are index + 1.
    """

    tree = ET.parse(xml_file)
    root = tree.getroot()

    labels = {}

    for lab in root.findall(".//label"):
        text = (lab.text or "").strip()

        if not text:
            continue

        label_value = int(lab.attrib["index"]) + 1
        labels[label_value] = text

    return labels


def label_has_hemisphere(text):
    """Return True when an atlas label already specifies laterality."""

    txt = text.lower()

    return any(
        term in txt
        for term in ("left", "right", "vermis")
    )


def assign_laterality(
    mask,
    affine,
    bilateral_ratio=0.30,
    midline_tolerance=2.0
):
    """
    Assign a component a laterality label.

    Returns
    -------
    L : predominantly left-hemisphere component
    R : predominantly right-hemisphere component
    B : component with substantial representation in both hemispheres
    C : central component or empty mask

    Coordinates are calculated in world space using the image affine.
    """

    coords = np.array(np.where(mask)).T

    if coords.size == 0:
        return "C"

    mni = nib.affines.apply_affine(
        affine,
        coords
    )

    x = mni[:, 0]

    left = int(np.sum(x < 0))
    right = int(np.sum(x > 0))

    if left == 0 and right == 0:
        return "C"

    if abs(float(np.mean(x))) < midline_tolerance:
        return "C"

    if left > 0 and right > 0:
        ratio = min(left, right) / max(left, right)

        if ratio >= bilateral_ratio:
            return "B"

    return "L" if left > right else "R"


def compute_centroid(mask, affine):
    """Return the mean voxel coordinate in world space (X, Y, Z)."""

    coords = np.array(np.where(mask)).T

    if coords.size == 0:
        return (np.nan, np.nan, np.nan)

    world_coords = nib.affines.apply_affine(
        affine,
        coords
    )

    centroid = np.mean(world_coords, axis=0)

    return tuple(np.round(centroid, 1))


def safe_filename(text):
    """Convert a region name into a filename-friendly string."""

    text = text.replace(" ", "_")
    text = text.replace(",", "")
    text = text.replace("/", "_")
    text = text.replace("\\", "_")

    return text


def save_mask(mask, reference_img, output_path):
    """Save a binary ROI mask using the reference image geometry."""

    header = reference_img.header.copy()
    header.set_data_dtype(np.uint8)

    out_img = nib.Nifti1Image(
        mask.astype(np.uint8),
        reference_img.affine,
        header
    )

    nib.save(out_img, output_path)


def validate_files(input_file, atlases):
    """Check that the input and required atlas files exist."""

    if not os.path.isfile(input_file):
        raise FileNotFoundError(
            f"Input NIfTI file not found: {input_file}"
        )

    for atlas_name, atlas_info in atlases.items():
        for key in ("img", "xml"):
            path = atlas_info[key]

            if not os.path.isfile(path):
                raise FileNotFoundError(
                    f"Missing {atlas_name} atlas {key}: {path}"
                )


# ==========================================================
# ARGUMENTS
# ==========================================================

def parse_arguments():
    """Configure command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Extract atlas-defined ROIs from a thresholded "
            "statistical NIfTI image."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input 3D thresholded statistical NIfTI image."
    )

    parser.add_argument(
        "--outdir",
        default="outputs",
        help="Output directory. Default: outputs"
    )

    parser.add_argument(
        "--atlas_dir",
        default=os.path.join(
            DEFAULT_FSLDIR,
            "data",
            "atlases"
        ),
        help="FSL atlas directory."
    )

    parser.add_argument(
        "--connectivity",
        type=int,
        choices=[6, 18, 26],
        default=26,
        help="Cluster connectivity. Default: 26"
    )

    parser.add_argument(
        "--min_cluster_size",
        type=int,
        default=100,
        help="Minimum cluster size in voxels. Default: 100"
    )

    parser.add_argument(
        "--min_roi_size",
        type=int,
        default=50,
        help="Minimum ROI component size in voxels. Default: 50"
    )

    parser.add_argument(
        "--min_cluster_fraction",
        type=float,
        default=0.02,
        help=(
            "Minimum ROI component / parent cluster voxel ratio. "
            "Set to 0 to disable. Default: 0.02"
        )
    )

    parser.add_argument(
        "--bilateral_ratio",
        type=float,
        default=0.30,
        help=(
            "Minimum smaller/larger hemisphere voxel ratio for "
            "a bilateral label. Default: 0.30"
        )
    )

    parser.add_argument(
        "--midline_tolerance",
        type=float,
        default=2.0,
        help=(
            "Central-region threshold in mm, based on mean world X. "
            "Default: 2.0"
        )
    )

    parser.add_argument(
        "--clean_output",
        action="store_true",
        help=(
            "Delete the existing output directory before processing. "
            "Use with caution."
        )
    )

    return parser.parse_args()


# ==========================================================
# MAIN WORKFLOW
# ==========================================================

def main():

    args = parse_arguments()

    # ------------------------------------------------------
    # VALIDATE PARAMETERS
    # ------------------------------------------------------

    if args.min_cluster_size < 1:
        raise ValueError(
            "--min_cluster_size must be at least 1."
        )

    if args.min_roi_size < 1:
        raise ValueError(
            "--min_roi_size must be at least 1."
        )

    if not 0 <= args.min_cluster_fraction <= 1:
        raise ValueError(
            "--min_cluster_fraction must be between 0 and 1."
        )

    if not 0 <= args.bilateral_ratio <= 1:
        raise ValueError(
            "--bilateral_ratio must be between 0 and 1."
        )

    if args.midline_tolerance < 0:
        raise ValueError(
            "--midline_tolerance cannot be negative."
        )

    out_base = os.path.abspath(args.outdir)

    # Prevent accidental deletion of the filesystem root.
    if (
        args.clean_output
        and out_base == os.path.abspath(os.path.sep)
    ):
        raise ValueError(
            "Refusing to clean the filesystem root."
        )

    # ------------------------------------------------------
    # PREPARE OUTPUT DIRECTORY
    # ------------------------------------------------------

    if args.clean_output and os.path.isdir(out_base):
        print(f"Cleaning output directory: {out_base}")
        shutil.rmtree(out_base)

    os.makedirs(out_base, exist_ok=True)

    # ------------------------------------------------------
    # LOAD ATLAS CONFIGURATION
    # ------------------------------------------------------

    atlases = get_atlases(args.atlas_dir)

    validate_files(args.input, atlases)

    structure = get_structure(args.connectivity)

    # ------------------------------------------------------
    # LOAD INPUT MAP
    # ------------------------------------------------------

    print("Loading input map...")

    zimg = nib.load(args.input)

    if len(zimg.shape) != 3:
        raise ValueError(
            "The input image must be 3D. "
            f"Received shape: {zimg.shape}"
        )

    zdata = zimg.get_fdata()

    if not np.all(np.isfinite(zdata)):
        raise ValueError(
            "Input image contains NaN or infinite values."
        )

    # Preserve the original method:
    # select positive voxels from the thresholded statistical map.
    thresholded_mask = zdata > 0

    if not np.any(thresholded_mask):
        raise ValueError(
            "No positive-valued voxels found in the input image."
        )

    # ------------------------------------------------------
    # CLUSTER DETECTION
    # ------------------------------------------------------

    print("Finding clusters...")

    labeled, n_clusters = label(
        thresholded_mask,
        structure=structure
    )

    cluster_map = np.zeros(
        labeled.shape,
        dtype=np.int32
    )

    new_id = 1

    for cid in range(1, n_clusters + 1):

        cluster = labeled == cid

        if int(np.sum(cluster)) >= args.min_cluster_size:
            cluster_map[cluster] = new_id
            new_id += 1

    cluster_ids = np.unique(cluster_map)
    cluster_ids = cluster_ids[cluster_ids > 0]

    print(f"Detected clusters: {n_clusters}")
    print(f"Retained clusters: {len(cluster_ids)}")

    cluster_img = nib.Nifti1Image(
        cluster_map,
        zimg.affine,
        zimg.header
    )

    nib.save(
        cluster_img,
        os.path.join(
            out_base,
            "clusters_combined.nii.gz"
        )
    )

    # ------------------------------------------------------
    # LOAD AND RESAMPLE ATLASES
    # ------------------------------------------------------

    atlas_data = {}
    atlas_labels = {}

    for atlas_name, atlas_info in atlases.items():

        print(f"Loading {atlas_name} atlas...")

        atlas_img = nib.load(atlas_info["img"])

        if len(atlas_img.shape) != 3:
            raise ValueError(
                f"{atlas_name} atlas must be a 3D label image."
            )

        # Preserve original nearest-neighbour resampling.
        resampled_img = resample_to_img(
            atlas_img,
            zimg,
            interpolation="nearest",
            force_resample=True,
            copy_header=True
        )

        atlas_values = resampled_img.get_fdata()

        if not np.all(np.isfinite(atlas_values)):
            raise ValueError(
                f"{atlas_name} atlas contains non-finite values."
            )

        # Preserve the original rounded integer label values.
        atlas_data[atlas_name] = np.rint(
            atlas_values
        ).astype(np.int32)

        atlas_labels[atlas_name] = parse_fsl_xml(
            atlas_info["xml"]
        )

        os.makedirs(
            os.path.join(
                out_base,
                atlas_name,
                "masks"
            ),
            exist_ok=True
        )

    # ------------------------------------------------------
    # SUBCORTICAL EXCLUSIONS
    # ------------------------------------------------------

    exclude_subcortical_labels = {
        "Left Cerebral White Matter",
        "Right Cerebral White Matter",
        "Left Cerebral Cortex",
        "Right Cerebral Cortex",
    }

    # ------------------------------------------------------
    # ROI EXTRACTION
    # ------------------------------------------------------

    summary_rows = []

    for atlas_name in atlases:

        print(f"\nProcessing {atlas_name}...")

        atlas_vol = atlas_data[atlas_name]
        labels_dict = atlas_labels[atlas_name]

        for cid in cluster_ids:

            cluster_mask = cluster_map == cid

            cluster_size = int(
                np.sum(cluster_mask)
            )

            region_ids = np.unique(
                atlas_vol[
                    cluster_mask & (atlas_vol > 0)
                ]
            )

            roi_entries = []

            for rid in region_ids:

                rid = int(rid)

                label_text = labels_dict.get(
                    rid,
                    f"Unknown_{rid}"
                )

                if (
                    atlas_name == "subcortical"
                    and label_text in exclude_subcortical_labels
                ):
                    continue

                atlas_region_mask = (
                    cluster_mask & (atlas_vol == rid)
                )

                # Preserve the original behaviour:
                # split disconnected pieces within each atlas region.
                sub_labels, n_sub = label(
                    atlas_region_mask,
                    structure=structure
                )

                for sub_id in range(1, n_sub + 1):

                    component = sub_labels == sub_id

                    roi_size = int(
                        np.sum(component)
                    )

                    if roi_size < args.min_roi_size:
                        continue

                    fraction = roi_size / cluster_size

                    if (
                        args.min_cluster_fraction > 0
                        and fraction < args.min_cluster_fraction
                    ):
                        continue

                    roi_entries.append({
                        "RegionID": rid,
                        "LabelText": label_text,
                        "VoxelCount": roi_size,
                        "Fraction": fraction,
                        "Mask": component,
                    })

            # IMPORTANT:
            # Restore the original sorting rule.
            # Sort by ROI size in descending order.
            # Python's sort is stable, so equal-sized entries retain
            # their original order from the extraction loop.
            roi_entries.sort(
                key=lambda item: item["VoxelCount"],
                reverse=True
            )

            for rank, entry in enumerate(
                roi_entries,
                start=1
            ):

                rid = entry["RegionID"]
                label_text = entry["LabelText"]
                roi_mask = entry["Mask"]

                if label_has_hemisphere(label_text):

                    region_name = label_text

                else:

                    side = assign_laterality(
                        roi_mask,
                        zimg.affine,
                        bilateral_ratio=args.bilateral_ratio,
                        midline_tolerance=args.midline_tolerance
                    )

                    region_name = f"{label_text}_{side}"

                cx, cy, cz = compute_centroid(
                    roi_mask,
                    zimg.affine
                )

                safe_name = safe_filename(region_name)

                # The region ID is included to make filenames clearer.
                # This does not alter the ROI mask itself.
                out_name = (
                    f"cluster_{cid:03d}_"
                    f"rank_{rank:02d}_"
                    f"region_{rid:03d}_"
                    f"{safe_name}.nii.gz"
                )

                out_path = os.path.join(
                    out_base,
                    atlas_name,
                    "masks",
                    out_name
                )

                save_mask(
                    roi_mask,
                    zimg,
                    out_path
                )

                summary_rows.append({
                    "ClusterID": int(cid),
                    "Atlas": atlas_name,
                    "RegionRank": rank,
                    "PrimaryRegion": rank == 1,
                    "RegionID": rid,
                    "RegionName": region_name,
                    "VoxelCount": entry["VoxelCount"],
                    "PercentOfCluster": round(
                        100 * entry["Fraction"],
                        2
                    ),
                    "ClusterSize": cluster_size,
                    "CentroidX": cx,
                    "CentroidY": cy,
                    "CentroidZ": cz,
                    "Filename": out_name,
                })

    # ------------------------------------------------------
    # SAVE SUMMARY TABLES
    # ------------------------------------------------------

    fieldnames = [
        "ClusterID",
        "Atlas",
        "RegionRank",
        "PrimaryRegion",
        "RegionID",
        "RegionName",
        "VoxelCount",
        "PercentOfCluster",
        "ClusterSize",
        "CentroidX",
        "CentroidY",
        "CentroidZ",
        "Filename",
    ]

    for atlas_name in atlases:

        rows = [
            row for row in summary_rows
            if row["Atlas"] == atlas_name
        ]

        csv_file = os.path.join(
            out_base,
            atlas_name,
            f"{atlas_name}_summary.csv"
        )

        with open(
            csv_file,
            "w",
            newline="",
            encoding="utf-8"
        ) as f:

            writer = csv.DictWriter(
                f,
                fieldnames=fieldnames
            )

            writer.writeheader()
            writer.writerows(rows)

        print(
            f"{atlas_name}: {len(rows)} ROI components saved."
        )

    print("\nProcessing complete.")
    print(f"Output directory: {out_base}")


if __name__ == "__main__":
    main()

