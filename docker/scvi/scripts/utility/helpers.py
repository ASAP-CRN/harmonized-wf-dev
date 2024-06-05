import anndata as ad
from scipy.sparse import csr_matrix
import tables
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata
from typing import Dict, Optional

from scipy.stats import multinomial
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import (
    accuracy_score,
    roc_auc_score,
    average_precision_score,
    roc_curve,
    precision_recall_curve,
)
from scipy.special import softmax
from scvi.external import SOLO
from scvi.model import SCVI
import umap

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

matplotlib.use("Agg")


def update_validation_metrics(adata: ad.AnnData, step: str, val_metrics: pd.DataFrame):
    """Update the validation metrics dataframe with new metrics for adata object."""
    new_metrics = get_validation_metrics(adata, step)
    val_metrics = pd.concat([val_metrics, new_metrics], ignore_index=True)
    return val_metrics


def get_validation_metrics(adata: ad.AnnData, step: str):
    """Log validation metrics for adata object."""
    n_samples = adata.obs["sample"].nunique()
    n_cells = adata.n_obs
    n_genes = adata.n_vars

    # fractional metrics
    n_mito_contaminated = (adata.obs["pct_counts_mt"] > 10).sum()
    n_predicted_doublet = (adata.obs["doublet_score"] >= 0.2).sum()
    n_low_counts = (adata.obs["total_counts"] < 500).sum()
    n_high_counts = (adata.obs["total_counts"] > 100000).sum()
    n_low_genes = (adata.obs["n_genes_by_counts"] < 300).sum()
    n_high_genes = (adata.obs["n_genes_by_counts"] > 10000).sum()

    # create a dataframe with the metrics plus the column step= "concatenation"
    val_metrics = pd.DataFrame(
        {
            "step": [step],
            "n_samples": [n_samples],
            "n_cells": [n_cells],
            "n_genes": [n_genes],
            "n_mito_contaminated": [n_mito_contaminated],
            "n_predicted_doublet": [n_predicted_doublet],
            "n_low_counts": [n_low_counts],
            "n_high_counts": [n_high_counts],
            "n_low_genes": [n_low_genes],
            "n_high_genes": [n_high_genes],
        }
    )

    return val_metrics


def minify_adata(adata: ad.AnnData):
    """Minify the adata object to save space."""

    all_zeros = csr_matrix(adata.X.shape)
    layers = {layer: all_zeros.copy() for layer in adata.layers}
    bdata = ad.AnnData(
        X=all_zeros,
        layers=layers,
        uns=adata.uns.copy(),
        obs=adata.obs,
        var=adata.var,
        varm=adata.varm,
        obsm=adata.obsm,
        obsp=adata.obsp,
    )
    return bdata


# copied from cellbender.remove_background.downstream which only works with python 3.7
"""Functions for downstream work with outputs of remove-background."""


def dict_from_h5(file: str) -> Dict[str, np.ndarray]:
    """Read in everything from an h5 file and put into a dictionary.

    Args:
        file: The h5 file

    Returns:
        Dictionary containing all the information from the h5 file
    """
    d = {}
    with tables.open_file(file) as f:
        # read in everything
        for array in f.walk_nodes("/", "Array"):
            d[array.name] = array.read()
    return d


def anndata_from_h5(file: str, analyzed_barcodes_only: bool = True) -> anndata.AnnData:
    """Load an output h5 file into an AnnData object for downstream work.

    Args:
        file: The h5 file
        analyzed_barcodes_only: False to load all barcodes, so that the size of
            the AnnData object will match the size of the input raw count matrix.
            True to load a limited set of barcodes: only those analyzed by the
            algorithm. This allows relevant latent variables to be loaded
            properly into adata.obs and adata.obsm, rather than adata.uns.

    Returns:
        anndata.AnnData: The anndata object, populated with inferred latent variables
            and metadata.

    """

    d = dict_from_h5(file)
    X = (
        sp.csc_matrix(
            (d.pop("data"), d.pop("indices"), d.pop("indptr")), shape=d.pop("shape")
        )
        .transpose()
        .tocsr()
    )

    # check and see if we have barcode index annotations, and if the file is filtered
    barcode_key = [k for k in d.keys() if (("barcode" in k) and ("ind" in k))]
    if len(barcode_key) > 0:
        max_barcode_ind = d[barcode_key[0]].max()
        filtered_file = max_barcode_ind >= X.shape[0]
    else:
        filtered_file = True

    if analyzed_barcodes_only:
        if filtered_file:
            # filtered file being read, so we don't need to subset
            print('Assuming we are loading a "filtered" file that contains only cells.')
            pass
        elif "barcode_indices_for_latents" in d.keys():
            X = X[d["barcode_indices_for_latents"], :]
            d["barcodes"] = d["barcodes"][d["barcode_indices_for_latents"]]
        elif "barcodes_analyzed_inds" in d.keys():
            X = X[d["barcodes_analyzed_inds"], :]
            d["barcodes"] = d["barcodes"][d["barcodes_analyzed_inds"]]
        else:
            print(
                "Warning: analyzed_barcodes_only=True, but the key "
                '"barcodes_analyzed_inds" or "barcode_indices_for_latents" '
                "is missing from the h5 file. "
                "Will output all barcodes, and proceed as if "
                "analyzed_barcodes_only=False"
            )

    # Construct the anndata object.
    adata = anndata.AnnData(
        X=X,
        obs={"barcode": d.pop("barcodes").astype(str)},
        var={
            "gene_name": (
                d.pop("gene_names") if "gene_names" in d.keys() else d.pop("name")
            ).astype(str)
        },
        dtype=X.dtype,
    )
    adata.obs.set_index("barcode", inplace=True)
    adata.var.set_index("gene_name", inplace=True)

    # For CellRanger v2 legacy format, "gene_ids" was called "genes"... rename this
    if "genes" in d.keys():
        d["id"] = d.pop("genes")

    # For purely aesthetic purposes, rename "id" to "gene_id"
    if "id" in d.keys():
        d["gene_id"] = d.pop("id")

    # If genomes are empty, try to guess them based on gene_id
    if "genome" in d.keys():
        if np.array([s.decode() == "" for s in d["genome"]]).all():
            if "_" in d["gene_id"][0].decode():
                print(
                    "Genome field blank, so attempting to guess genomes based on gene_id prefixes"
                )
                d["genome"] = np.array(
                    [s.decode().split("_")[0] for s in d["gene_id"]], dtype=str
                )

    # Add other information to the anndata object in the appropriate slot.
    _fill_adata_slots_automatically(adata, d)

    # Add a special additional field to .var if it exists.
    if "features_analyzed_inds" in adata.uns.keys():
        adata.var["cellbender_analyzed"] = [
            True if (i in adata.uns["features_analyzed_inds"]) else False
            for i in range(adata.shape[1])
        ]
    elif "features_analyzed_inds" in adata.var.keys():
        adata.var["cellbender_analyzed"] = [
            True if (i in adata.var["features_analyzed_inds"].values) else False
            for i in range(adata.shape[1])
        ]

    if analyzed_barcodes_only:
        for col in adata.obs.columns[
            adata.obs.columns.str.startswith("barcodes_analyzed")
            | adata.obs.columns.str.startswith("barcode_indices")
        ]:
            try:
                del adata.obs[col]
            except Exception:
                pass
    else:
        # Add a special additional field to .obs if all barcodes are included.
        if "barcodes_analyzed_inds" in adata.uns.keys():
            adata.obs["cellbender_analyzed"] = [
                True if (i in adata.uns["barcodes_analyzed_inds"]) else False
                for i in range(adata.shape[0])
            ]
        elif "barcodes_analyzed_inds" in adata.obs.keys():
            adata.obs["cellbender_analyzed"] = [
                True if (i in adata.obs["barcodes_analyzed_inds"].values) else False
                for i in range(adata.shape[0])
            ]

    return adata


def _fill_adata_slots_automatically(adata, d):
    """Add other information to the adata object in the appropriate slot."""

    # TODO: what about "features_analyzed_inds"?  If not all features are analyzed, does this work?

    for key, value in d.items():
        try:
            if value is None:
                continue
            value = np.asarray(value)
            if len(value.shape) == 0:
                adata.uns[key] = value
            elif value.shape[0] == adata.shape[0]:
                if (len(value.shape) < 2) or (value.shape[1] < 2):
                    adata.obs[key] = value
                else:
                    adata.obsm[key] = value
            elif value.shape[0] == adata.shape[1]:
                if value.dtype.name.startswith("bytes"):
                    adata.var[key] = value.astype(str)
                else:
                    adata.var[key] = value
            else:
                adata.uns[key] = value
        except Exception:
            print("Unable to load data into AnnData: ", key, value, type(value))


# from calico/solo repo
def get_solo_results(
    solo: SOLO,
    adata: ad.AnnData,
    vae: SCVI,
    doublet_ratio: float | None = 0.05,
    expected_number_of_doublets: Optional[int] = None,
    gen_report: bool = False,
):
    """ """

    num_cells, num_genes = adata.shape

    latent = vae.get_latent_representation()

    logit_predictions = solo.predict(include_simulated_doublets=True)

    is_doublet_known = solo.adata.obs._solo_doub_sim == "doublet"
    is_doublet_pred = logit_predictions.idxmin(axis=1) == "singlet"

    validation_is_doublet_known = is_doublet_known[solo.validation_indices]
    validation_is_doublet_pred = is_doublet_pred[solo.validation_indices]
    training_is_doublet_known = is_doublet_known[solo.train_indices]
    training_is_doublet_pred = is_doublet_pred[solo.train_indices]

    valid_as = accuracy_score(validation_is_doublet_known, validation_is_doublet_pred)
    valid_roc = roc_auc_score(validation_is_doublet_known, validation_is_doublet_pred)
    valid_ap = average_precision_score(
        validation_is_doublet_known, validation_is_doublet_pred
    )

    train_as = accuracy_score(training_is_doublet_known, training_is_doublet_pred)
    train_roc = roc_auc_score(training_is_doublet_known, training_is_doublet_pred)
    train_ap = average_precision_score(
        training_is_doublet_known, training_is_doublet_pred
    )

    print(f"Training results")
    print(f"AUROC: {train_roc}, Accuracy: {train_as}, Average precision: {train_ap}")

    print(f"Validation results")
    print(f"AUROC: {valid_roc}, Accuracy: {valid_as}, Average precision: {valid_ap}")

    # write predictions
    # softmax predictions
    softmax_predictions = pd.DataFrame(
        softmax(logit_predictions, axis=1),
        columns=logit_predictions.columns,
        index=logit_predictions.index,
    )

    doublet_score = softmax_predictions.loc[:, "doublet"]

    # logit predictions
    logit_doublet_score = logit_predictions.loc[:, "doublet"]

    # update threshold as a function of Solo's estimate of the number of
    # doublets
    # essentially a log odds update
    # TODO put in a function
    # currently overshrinking softmaxes

    if doublet_ratio is None:
        expected_number_of_doublets = None
    else:
        expected_number_of_doublets = int(doublet_ratio * num_cells)

    solo_scores = doublet_score[:num_cells]

    if expected_number_of_doublets is not None:
        k = len(solo_scores) - expected_number_of_doublets
        if expected_number_of_doublets / len(solo_scores) > 0.5:
            print(
                """Make sure you actually expect more than half your cells
                   to be doublets. If not change your
                   -e parameter value"""
            )
        assert k > 0
        idx = np.argpartition(solo_scores, k)
        threshold = np.max(solo_scores[idx[:k]])
        is_solo_doublet = solo_scores > threshold
    else:
        is_solo_doublet = solo_scores > 0.5

    smoothed_preds = knn_smooth_pred_class(
        X=latent, pred_class=is_doublet_pred[:num_cells]
    )

    ret_obs = adata.obs.copy()
    ret_obs["is_doublet"] = is_solo_doublet[:num_cells].values.astype(bool)
    ret_obs["logit_scores"] = logit_doublet_score[:num_cells].values.astype(float)
    ret_obs["softmax_scores"] = solo_scores[:num_cells].values.astype(float)

    ret_figs = []
    if gen_report:

        train_solo_scores = doublet_score[solo.train_indices]
        validation_solo_scores = doublet_score[solo.validation_indices]

        train_fpr, train_tpr, _ = roc_curve(
            training_is_doublet_known, train_solo_scores
        )
        val_fpr, val_tpr, _ = roc_curve(
            validation_is_doublet_known, validation_solo_scores
        )
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        # plot ROC
        ax.plot(train_fpr, train_tpr, label="Train")
        ax.plot(val_fpr, val_tpr, label="Validation")
        ax.set_xlabel("False positive rate")
        ax.set_ylabel("True positive rate")
        ax.legend()
        # plt.savefig(os.path.join(args.out_dir, "roc.pdf"))
        ret_figs.append(fig)

        train_precision, train_recall, _ = precision_recall_curve(
            training_is_doublet_known, train_solo_scores
        )
        val_precision, val_recall, _ = precision_recall_curve(
            validation_is_doublet_known, validation_solo_scores
        )

        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        # plot accuracy
        ax.plot(train_recall, train_precision, label="Train")
        ax.plot(val_recall, val_precision, label="Validation")
        ax.set_xlabel("Recall")
        ax.set_ylabel("pytPrecision")
        ax.legend()
        # plt.savefig(os.path.join(args.out_dir, "precision_recall.pdf"))
        # plt.show()
        ret_figs.append(fig)

        # plot distributions
        obs_indices = solo.validation_indices[solo.validation_indices < num_cells]
        sim_indices = solo.validation_indices[solo.validation_indices > num_cells]

        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        sns.displot(doublet_score[sim_indices], label="Simulated", ax=ax)
        sns.displot(doublet_score[obs_indices], label="Observed", ax=ax)
        ax.legend()
        # plt.savefig(os.path.join(args.out_dir, "sim_vs_obs_dist.pdf"))
        # plt.show()
        ret_figs.append(fig)

        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        sns.distplot(solo_scores[:num_cells], label="Observed (transformed)", ax=ax)
        ax.legend()
        # plt.savefig(os.path.join(args.out_dir, "real_cells_dist.pdf"))
        # plt.show()
        ret_figs.append(fig)

        scvi_umap = umap.UMAP(n_neighbors=16).fit_transform(latent)
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        cax = ax.scatter(
            scvi_umap[:, 0],
            scvi_umap[:, 1],
            c=doublet_score[:num_cells],
            s=8,
            cmap="GnBu",
        )
        plt.colorbar(cax, ax=ax, pad=0.01, fraction=0.08, aspect=30, location="right")
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title("SOLO doublet_score")
        # fig.savefig(os.path.join(args.out_dir, "umap_solo_scores.pdf"))
        # plt.show()
        ret_figs.append(fig)
        return ret_obs, ret_figs

    return ret_obs


# https://github.com/calico/solo/blob/master/solo/utils.py
def knn_smooth_pred_class(
    X: np.ndarray,
    pred_class: np.ndarray,
    grouping: np.ndarray = None,
    k: int = 15,
) -> np.ndarray:
    """
    Smooths class predictions by taking the modal class from each cell's
    nearest neighbors.
    Parameters
    ----------
    X : np.ndarray
        [N, Features] embedding space for calculation of nearest neighbors.
    pred_class : np.ndarray
        [N,] array of unique class labels.
    groupings : np.ndarray
        [N,] unique grouping labels for i.e. clusters.
        if provided, only considers nearest neighbors *within the cluster*.
    k : int
        number of nearest neighbors to use for smoothing.
    Returns
    -------
    smooth_pred_class : np.ndarray
        [N,] unique class labels, smoothed by kNN.
    Examples
    --------
    >>> smooth_pred_class = knn_smooth_pred_class(
    ...     X = X,
    ...     pred_class = raw_predicted_classes,
    ...     grouping = louvain_cluster_groups,
    ...     k = 15,)
    Notes
    -----
    scNym classifiers do not incorporate neighborhood information.
    By using a simple kNN smoothing heuristic, we can leverage neighborhood
    information to improve classification performance, smoothing out cells
    that have an outlier prediction relative to their local neighborhood.
    """
    if grouping is None:
        # do not use a grouping to restrict local neighborhood
        # associations, create a universal pseudogroup `0`.
        grouping = np.zeros(X.shape[0])

    smooth_pred_class = np.zeros_like(pred_class)
    for group in np.unique(grouping):
        # identify only cells in the relevant group
        group_idx = np.where(grouping == group)[0].astype("int")
        X_group = X[grouping == group, :]
        # if there are < k cells in the group, change `k` to the
        # group size
        if X_group.shape[0] < k:
            k_use = X_group.shape[0]
        else:
            k_use = k
        # compute a nearest neighbor graph and identify kNN
        nns = NearestNeighbors(
            n_neighbors=k_use,
        ).fit(X_group)
        dist, idx = nns.kneighbors(X_group)

        # for each cell in the group, assign a class as
        # the majority class of the kNN
        for i in range(X_group.shape[0]):
            classes = pred_class[group_idx[idx[i, :]]]
            uniq_classes, counts = np.unique(classes, return_counts=True)
            maj_class = uniq_classes[int(np.argmax(counts))]
            smooth_pred_class[group_idx[i]] = maj_class
    return smooth_pred_class


# copied from scib: https://github.com/theislab/scib/blob/main/scib/preprocessing.py
# Cell Cycle
def score_cell_cycle(adata, organism="mouse"):
    """Score cell cycle score given an organism

    Wrapper function for `scanpy.tl.score_genes_cell_cycle`_

    .. _scanpy.tl.score_genes_cell_cycle: https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes_cell_cycle.html

    Tirosh et al. cell cycle marker genes downloaded from
    https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt

    s_genes and g2m_genes extracted from cell_cycle_genes like this:
    ``` python
    cell_cycle_genes = [x.strip() for x in open('./data/regev_lab_cell_cycle_genes.txt')]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    ```

    For human, mouse genes are capitalised and used directly. This is under the assumption that cell cycle genes are
    well conserved across species.

    :param adata: anndata object containing
    :param organism: organism of gene names to match cell cycle genes
    :return: tuple of ``(s_genes, g2m_genes)`` of S-phase genes and G2- and M-phase genes scores
    """
    import pathlib

    root = pathlib.Path(__file__).parent

    cc_files = {
        "mouse": [
            "/opt/resources/s_genes_tirosh.txt",
            "/opt/resources/g2m_genes_tirosh.txt",
        ],
        "human": [
            "/opt/resources/s_genes_tirosh_hm.txt",
            "/opt/resources/g2m_genes_tirosh_hm.txt",
        ],
    }

    with open(cc_files[organism][0]) as f:
        s_genes = [x.strip() for x in f.readlines() if x.strip() in adata.var.index]
    with open(cc_files[organism][1]) as f:
        g2m_genes = [x.strip() for x in f.readlines() if x.strip() in adata.var.index]

    if (len(s_genes) == 0) or (len(g2m_genes) == 0):
        rand_choice = np.random.randint(1, adata.n_vars, 10)
        rand_genes = adata.var_names[rand_choice].tolist()
        raise ValueError(
            f"cell cycle genes not in adata\n organism: {organism}\n varnames: {rand_genes}"
        )

    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
