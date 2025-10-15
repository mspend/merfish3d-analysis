import napari
from tifffile import imread
import pandas as pd
import numpy as np
from pathlib import Path
from numpy.typing import ArrayLike
from scipy.spatial import cKDTree

def greedy_match_indices_within_radius(
    qi2lab_coords: ArrayLike,
    qi2lab_gene_ids: ArrayLike,
    gt_coords: ArrayLike,
    gt_gene_ids: ArrayLike,
    radius: float,
):
    """
    Greedy closest-first, within-radius, same-gene, one-to-one matching.

    Parameters
    ----------
    qi2lab_coords : ArrayLike
        (Nq, 3) z, y, x coordinates for found spots (microns; world coords).
    qi2lab_gene_ids : ArrayLike
        (Nq,) gene IDs for found spots.
    gt_coords : ArrayLike
        (Ng, 3) z, y, x coordinates for ground-truth spots (microns; world coords).
    gt_gene_ids : ArrayLike
        (Ng,) gene IDs for ground-truth spots.
    radius : float
        Maximum 3D distance (same units as coordinates) allowed for matching.

    Returns
    -------
    matched_q: ArrayLike
        matched qi2lab spot
    matched_g: ArrayLike
        match GT spot
    q_fp_idx: ArrayLike
        qi2lab false positives
    g_fn_idx: ArrayLike
        qi2lab false negatives    
    """

    qi2lab_coords = np.asarray(qi2lab_coords)
    qi2lab_gene_ids = np.asarray(qi2lab_gene_ids)
    gt_coords = np.asarray(gt_coords)
    gt_gene_ids = np.asarray(gt_gene_ids)

    Nq, Ng = len(qi2lab_coords), len(gt_coords)
    if Nq == 0 or Ng == 0:
        return (
            np.empty((0,), int),  # q_matched_idx
            np.empty((0,), int),  # g_matched_idx
            np.arange(Nq, dtype=int),  # q_fp_idx (all qi2lab if Ng==0)
            np.arange(Ng, dtype=int),  # g_fn_idx (all GT if Nq==0)
        )

    # Collect candidate pairs (global indices) within radius for shared genes
    pair_q, pair_g, pair_d = [], [], []
    for gene in np.intersect1d(np.unique(qi2lab_gene_ids), np.unique(gt_gene_ids)):
        q_idx = np.flatnonzero(qi2lab_gene_ids == gene)
        g_idx = np.flatnonzero(gt_gene_ids == gene)
        if q_idx.size == 0 or g_idx.size == 0:
            continue
        q_pts = qi2lab_coords[q_idx]
        g_pts = gt_coords[g_idx]
        q_tree = cKDTree(q_pts)
        g_tree = cKDTree(g_pts)
        try:
            coo = q_tree.sparse_distance_matrix(g_tree, max_distance=radius, output_type="coo_matrix")
        except TypeError:
            coo = q_tree.sparse_distance_matrix(g_tree, max_distance=radius).tocoo()
        if coo.nnz == 0:
            continue
        pair_q.append(q_idx[coo.row])
        pair_g.append(g_idx[coo.col])
        pair_d.append(coo.data)

    if not pair_q:
        return np.empty((0,), int), np.empty((0,), int), np.arange(Nq, int), np.arange(Ng, int)

    pair_q = np.concatenate(pair_q)
    pair_g = np.concatenate(pair_g)
    pair_d = np.concatenate(pair_d)

    # Greedy: accept closest first
    order = np.argsort(pair_d, kind="stable")
    pair_q, pair_g = pair_q[order], pair_g[order]

    q_used = np.zeros(Nq, bool)
    g_used = np.zeros(Ng, bool)
    matched_q, matched_g = [], []
    for qi, gi in zip(pair_q, pair_g):
        if q_used[qi] or g_used[gi]:
            continue
        q_used[qi] = True
        g_used[gi] = True
        matched_q.append(qi)
        matched_g.append(gi)

    matched_q = np.array(matched_q, int)
    matched_g = np.array(matched_g, int)

    q_fp_idx = np.flatnonzero(~q_used)   # False Positives (qi2lab unmatched)
    g_fn_idx = np.flatnonzero(~g_used)   # False Negatives (GT unmatched)

    return matched_q, matched_g, q_fp_idx, g_fn_idx

def greedy_F1_within_radius(
    qi2lab_coords: ArrayLike,
    qi2lab_gene_ids: ArrayLike,
    gt_coords: ArrayLike,
    gt_gene_ids: ArrayLike,
    radius: float,
) -> dict:
    """
    Compute F1 using greedy closest-first matching within a max radius, with same-gene and
    one-to-one constraints.

    Algorithm
    ---------
    1) For each gene present in both sets, find all qi2lab↔GT pairs within `radius` via KD-trees.
    2) Concatenate all candidates across genes; sort by distance (ascending).
    3) Greedily accept pairs if both endpoints are unused; remove both from consideration.
    4) TP = #accepted pairs; FP = #qi2lab unmatched; FN = #GT unmatched.

    Parameters
    ----------
    qi2lab_coords : ArrayLike
        (Nq, 3) z, y, x coordinates for found spots (microns; world coords).
    qi2lab_gene_ids : ArrayLike
        (Nq,) gene IDs for found spots.
    gt_coords : ArrayLike
        (Ng, 3) z, y, x coordinates for ground-truth spots (microns; world coords).
    gt_gene_ids : ArrayLike
        (Ng,) gene IDs for ground-truth spots.
    radius : float
        Maximum 3D distance (same units as coordinates) allowed for matching.

    Returns
    -------
    results : dict
        {
          "F1 Score": float,
          "Precision": float,
          "Recall": float,
          "True Positives": int,
          "False Positives": int,
          "False Negatives": int,
        }
    """
    qi2lab_coords = np.asarray(qi2lab_coords)
    qi2lab_gene_ids = np.asarray(qi2lab_gene_ids)
    gt_coords = np.asarray(gt_coords)
    gt_gene_ids = np.asarray(gt_gene_ids)

    Nq = qi2lab_coords.shape[0]
    Ng = gt_coords.shape[0]

    # Short-circuit trivial cases
    if Nq == 0 and Ng == 0:
        return {
            "F1 Score": 1.0,
            "Precision": 1.0,
            "Recall": 1.0,
            "True Positives": 0,
            "False Positives": 0,
            "False Negatives": 0,
        }
    if Nq == 0:
        return {
            "F1 Score": 0.0,
            "Precision": 0.0,
            "Recall": 0.0 if Ng > 0 else 1.0,
            "True Positives": 0,
            "False Positives": 0,
            "False Negatives": int(Ng),
        }
    if Ng == 0:
        return {
            "F1 Score": 0.0,
            "Precision": 0.0,
            "Recall": 0.0,
            "True Positives": 0,
            "False Positives": int(Nq),
            "False Negatives": 0,
        }

    # Build candidate pairs within radius, per gene; merge globally
    pair_q_idx_all: list[np.ndarray] = []
    pair_g_idx_all: list[np.ndarray] = []
    pair_dist_all:  list[np.ndarray] = []

    common_genes = np.intersect1d(np.unique(qi2lab_gene_ids), np.unique(gt_gene_ids))
    for gene in common_genes:
        q_idx = np.flatnonzero(qi2lab_gene_ids == gene)
        g_idx = np.flatnonzero(gt_gene_ids == gene)
        if q_idx.size == 0 or g_idx.size == 0:
            continue

        q_pts = qi2lab_coords[q_idx]
        g_pts = gt_coords[g_idx]

        q_tree = cKDTree(q_pts)
        g_tree = cKDTree(g_pts)
        try:
            dist_coo = q_tree.sparse_distance_matrix(
                g_tree, max_distance=radius, output_type="coo_matrix"
            )
        except TypeError:
            dist_coo = q_tree.sparse_distance_matrix(g_tree, max_distance=radius).tocoo()

        if dist_coo.nnz == 0:
            continue

        # Local -> global index mapping
        q_local = dist_coo.row
        g_local = dist_coo.col
        dists   = dist_coo.data

        pair_q_idx_all.append(q_idx[q_local])
        pair_g_idx_all.append(g_idx[g_local])
        pair_dist_all.append(dists)

    if not pair_q_idx_all:
        # No candidates within radius at all
        tp = 0
        fp = int(Nq)
        fn = int(Ng)
    else:
        pair_q_idx = np.concatenate(pair_q_idx_all)
        pair_g_idx = np.concatenate(pair_g_idx_all)
        pair_dist  = np.concatenate(pair_dist_all)

        # Sort globally by distance (closest first)
        order = np.argsort(pair_dist, kind="stable")
        pair_q_idx = pair_q_idx[order]
        pair_g_idx = pair_g_idx[order]

        # Greedy selection with one-to-one constraint
        q_used = np.zeros(Nq, dtype=bool)
        g_used = np.zeros(Ng, dtype=bool)

        tp = 0
        for qi, gi in zip(pair_q_idx, pair_g_idx):
            if q_used[qi] or g_used[gi]:
                continue
            q_used[qi] = True
            g_used[gi] = True
            tp += 1

        fp = int(Nq - tp)
        fn = int(Ng - g_used.sum())

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall    = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1        = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    return {
        "F1 Score": f1,
        "Precision": precision,
        "Recall": recall,
        "True Positives": int(tp),
        "False Positives": int(fp),
        "False Negatives": int(fn),
    }


def main(root_path: Path, sim_path: Path = "cells"):

    image_0315_voxel = np.array([0.315, 0.1625, .1625])
    image_10_voxel = np.array([1.0, 0.1625, .1625])
    image_15_voxel = np.array([1.5, 0.1625, .1625])

    z_spacings = ["0.315", "1.0", "1.5"]

    top_path = root_path / Path(sim_path)
    images = []
    status = []
    points = []
    f1 = []
    for z in z_spacings:
        temp = imread(top_path / Path(str(z)) / Path("data_r0001_tile0000_1") / Path("data_r0001_tile0000.tif"))
        if not(z == "0.315"):
            temp_filled = np.zeros_like(images[0])
            if z == "1.0":
                for z_idx in range(temp.shape[1]):
                    temp_filled[:, int(np.floor(z_idx*image_10_voxel[0]/image_0315_voxel[0])),:,:] = temp[:,z_idx,:,:]
            elif z == "1.5":
                for z_idx in range(temp.shape[1]):
                    temp_filled[:, int(np.floor(z_idx*image_15_voxel[0]/image_0315_voxel[0])),:,:] = temp[:,z_idx,:,:]
            images.append(temp_filled)
        else:
            images.append(temp)

        gt_spots = pd.read_csv(top_path / Path(str(z)) / Path("GT_spots.csv"))
        codebook = pd.read_csv(top_path / Path(str(z)) / Path("codebook.csv"))
        codebook_genes = codebook['gene_id'].to_numpy()
        
        decoded_spots = pd.read_parquet(top_path / Path(str(z)) / Path("decoded_features.parquet"))
        qi2lab_coords = decoded_spots[['global_z', 'global_y', 'global_x']].to_numpy()
        qi2lab_gene_ids = decoded_spots['gene_id'].to_numpy()
        gt_coords = gt_spots[['Z','X','Y']].to_numpy()
        gt_gene_ids = codebook_genes[(gt_spots['Gene_label'].to_numpy(dtype=int)-1)]
        gt_offset = [
            0, 0, 0 
            # (1*images[0].shape[-2]/2)*image_0315_voxel[1]-image_0315_voxel[1]/2,
            # (1*images[0].shape[-1]/2)*image_0315_voxel[2]-image_0315_voxel[2]/2
        ]

        gt_coords = gt_coords + gt_offset

        if z == "0.315":
            mq, mg, q_fp, g_fn = greedy_match_indices_within_radius(
                qi2lab_coords,  
                qi2lab_gene_ids,
                gt_coords, 
                gt_gene_ids,
                1.5
            )
            z_f1 = greedy_F1_within_radius(
                qi2lab_coords,  
                qi2lab_gene_ids,
                gt_coords, 
                gt_gene_ids,
                1.5
            )

        elif z == "1.0":
            mq, mg, q_fp, g_fn = greedy_match_indices_within_radius(
                qi2lab_coords,  
                qi2lab_gene_ids,
                gt_coords, 
                gt_gene_ids,
                1.5
            )
            z_f1 = greedy_F1_within_radius(
                qi2lab_coords,  
                qi2lab_gene_ids,
                gt_coords, 
                gt_gene_ids,
                1.5
            )
        elif z == "1.5":
            mq, mg, q_fp, g_fn = greedy_match_indices_within_radius(
                qi2lab_coords,  
                qi2lab_gene_ids,
                gt_coords, 
                gt_gene_ids,
                1.5
            )
            z_f1 = greedy_F1_within_radius(
                qi2lab_coords,  
                qi2lab_gene_ids,
                gt_coords, 
                gt_gene_ids,
                1.5
            )

        # Build one points array + one categorical status property
        z_parts = []
        z_labels = []

        if mq.size:
            z_parts.append(qi2lab_coords[mq])          # True Positives → plot qi2lab side
            z_labels.append(np.full(mq.size, "TP", dtype=object))
        if q_fp.size:
            z_parts.append(qi2lab_coords[q_fp])        # False Positives → qi2lab unmatched
            z_labels.append(np.full(q_fp.size, "FP", dtype=object))
        if g_fn.size:
            z_parts.append(gt_coords[g_fn])            # False Negatives → GT unmatched
            z_labels.append(np.full(g_fn.size, "FN", dtype=object))

        if z_parts:
            z_points = np.vstack(z_parts)
            z_status = np.concatenate(z_labels)
        else:
            z_points = np.empty((0, 3), dtype=float)
            z_status = np.empty((0,), dtype=object)

        status.append(z_status)
        points.append(z_points)
        f1.append(z_f1)

    print(f"F1 results z=0.315: {f1[0]}")
    print(f"F1 results z=1.0: {f1[1]}")
    print(f"F1 results z=1.5: {f1[2]}")

    viewer = napari.Viewer()

    layer6 = viewer.add_points(points[2][:, [0, 2]], name="1.5 RNA XZ", scale=[1,1], size=.75, properties={"status": status[2]}, face_color="status")
    layer6.face_color_cycle = {"TP": "gray", "FP": "cyan", "FN": "orange"}
    layer5 = viewer.add_points(points[1][:, [0, 2]], name="1.0 RNA XZ", scale=[1,1], size=.75, properties={"status": status[1]}, face_color="status")
    layer5.face_color_cycle = {"TP": "gray", "FP": "cyan", "FN": "orange"}
    layer4 = viewer.add_points(points[0][:, [0, 2]], name=".315 RNA XZ", scale=[1,1], size=.75, properties={"status": status[0]}, face_color="status")
    layer4.face_color_cycle = {"TP": "gray", "FP": "cyan", "FN": "orange"}

    viewer.add_image(np.squeeze(np.max(np.swapaxes(images[2],1,2),axis=1)), name="1.5 bit 1 XZ", scale=[image_0315_voxel[0],image_0315_voxel[1]], colormap='gray', blending='additive', contrast_limits=[0,4000])
    viewer.add_image(np.squeeze(np.max(np.swapaxes(images[1],1,2),axis=1)), name="1.0 bit 1 XZ", scale=[image_0315_voxel[0],image_0315_voxel[1]], colormap='gray', blending='additive', contrast_limits=[0,4000])
    viewer.add_image(np.squeeze(np.max(np.swapaxes(images[0],1,2),axis=1)), name="0.315 bit 1 XZ", scale=[image_0315_voxel[0],image_0315_voxel[1]], colormap='gray', blending='additive', contrast_limits=[0,4000])

    layer3 = viewer.add_points(points[2][:, [1, 2]], name="1.5 RNA XY", scale=[1,1], size=.75, properties={"status": status[2]}, face_color="status")
    layer3.face_color_cycle = {"TP": "gray", "FP": "cyan", "FN": "orange"}
    layer2 = viewer.add_points(points[1][:, [1, 2]], name="1.0 RNA XY", scale=[1,1], size=.75, properties={"status": status[1]}, face_color="status")
    layer2.face_color_cycle = {"TP": "gray", "FP": "cyan", "FN": "orange"}
    layer1 = viewer.add_points(points[0][:, [1, 2]], name=".315 RNA XY", scale=[1,1], size=.75, properties={"status": status[0]}, face_color="status")
    layer1.face_color_cycle = {"TP": "gray", "FP": "cyan", "FN": "orange"}

    viewer.add_image(np.max(images[2],axis=1), name="1.5 bit 1 XY", scale=image_0315_voxel[1:], colormap='gray', blending='additive', contrast_limits=[0,4000])
    viewer.add_image(np.max(images[1],axis=1), name="1.0 bit 1 XY", scale=image_0315_voxel[1:], colormap='gray', blending='additive', contrast_limits=[0,4000])
    viewer.add_image(np.max(images[0],axis=1), name="0.315 bit 1 XY", scale=image_0315_voxel[1:], colormap='gray', blending='additive', contrast_limits=[0,4000])


    viewer.scale_bar.visible = True
    viewer.scale_bar.unit = 'um'


    napari.run()

if __name__ == "__main__":
    root_path = Path(r"/data/smFISH/simulated_data/092225/paper_figure")
    sim_path = "cells" # "cells" or "flat" or ""
    main(root_path,sim_path)