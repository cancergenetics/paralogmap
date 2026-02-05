#!/usr/bin/env python3
"""Build static JSON products for the paralog interaction browser.

This script:
- Computes per-pair mean/std across all cell lines
- Generates top-100 lists for each cell line (by score and z)
- Generates top-100 lists for each gene pair (by score and z), bucketed
- Builds gene aggregates by mean score, bucketed
- Writes search/index JSON files and alias maps

Run once per dataset version; outputs to dist/.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from collections import defaultdict
import heapq
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


def normalize_alias(value: str) -> str:
    return "".join(ch for ch in value.upper() if ch.isalnum())


def parse_pair_column(col: str) -> Tuple[str, str]:
    parts = col.split("_", 1)
    if len(parts) != 2:
        raise ValueError(f"Unexpected pair column format: {col}")
    return parts[0], parts[1]


def pair_id_for(col: str) -> Tuple[str, str, str]:
    gene1, gene2 = parse_pair_column(col)
    pair_id = f"{gene1}__{gene2}"
    return pair_id, gene1, gene2


def pair_bucket(pair_id: str) -> str:
    digest = hashlib.md5(pair_id.encode("utf-8")).digest()
    return f"{digest[0]:02x}"


def disease_bucket(name: str) -> str:
    norm = normalize_alias(name)
    digest = hashlib.md5(norm.encode("utf-8")).digest()
    return f"{digest[0]:02x}"


def gene_bucket(symbol: str) -> str:
    norm = normalize_alias(symbol)
    if not norm:
        return "00"
    bucket = ord(norm[0]) % 64
    return f"{bucket:02d}"


def load_model_aliases(model_csv: Path) -> Dict[str, Dict[str, List[str]]]:
    if not model_csv.exists():
        return {}

    model_map: Dict[str, Dict[str, List[str]]] = {}
    with model_csv.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            depmap_id = row.get("ModelID") or row.get("DepMap_ID")
            if not depmap_id:
                continue
            aliases = []
            for key in ("ModelID", "DepMap_ID", "CellLineName", "StrippedCellLineName"):
                val = row.get(key)
                if val:
                    aliases.append(val)
            display = row.get("CellLineName") or depmap_id
            lineage = row.get("OncotreeLineage") or ""
            primary_disease = row.get("OncotreePrimaryDisease") or ""
            model_map[depmap_id] = {
                "display": display,
                "aliases": list(dict.fromkeys(aliases)),
                "lineage": lineage,
                "primary_disease": primary_disease,
            }
    return model_map


def load_gene_metadata(gene_json: Path) -> Dict[str, Dict[str, str]]:
    if not gene_json.exists():
        return {}
    with gene_json.open() as f:
        data = json.load(f)
    docs = data.get("response", {}).get("docs", [])
    meta = {}
    for doc in docs:
        symbol = doc.get("symbol")
        name = doc.get("name")
        hgnc_id = doc.get("hgnc_id")
        if symbol:
            meta[symbol] = {
                "name": name or "",
                "hgnc_id": hgnc_id or "",
            }
    return meta


def top_k_indices(arr: np.ndarray, k: int) -> np.ndarray:
    if k >= arr.size:
        return np.argsort(arr)[::-1]
    idx = np.argpartition(arr, -k)[-k:]
    return idx[np.argsort(arr[idx])[::-1]]


def top_k_indices_masked(arr: np.ndarray, k: int) -> np.ndarray:
    valid = np.isfinite(arr)
    if not valid.any():
        return np.array([], dtype=np.int64)
    if valid.sum() <= k:
        valid_idx = np.where(valid)[0]
        order = np.argsort(arr[valid])[::-1]
        return valid_idx[order]
    arr2 = arr.copy()
    arr2[~valid] = -np.inf
    idx = np.argpartition(arr2, -k)[-k:]
    return idx[np.argsort(arr2[idx])[::-1]]


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(payload, f, separators=(",", ":"))


def rankdata(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty_like(order, dtype=np.float64)
    n = len(values)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg_rank = (i + j + 2) / 2.0
        ranks[order[i : j + 1]] = avg_rank
        i = j + 1
    return ranks


def mann_whitney_pvalue(x: np.ndarray, y: np.ndarray) -> float:
    try:
        from scipy.stats import mannwhitneyu

        result = mannwhitneyu(x, y, alternative="two-sided")
        return float(result.pvalue)
    except Exception:
        pass

    n1 = x.size
    n2 = y.size
    if n1 == 0 or n2 == 0:
        return 1.0
    data = np.concatenate([x, y])
    ranks = rankdata(data)
    r1 = ranks[:n1].sum()
    u1 = r1 - n1 * (n1 + 1) / 2.0
    u2 = n1 * n2 - u1
    u = min(u1, u2)
    mean_u = n1 * n2 / 2.0
    _, counts = np.unique(data, return_counts=True)
    tie_term = np.sum(counts**3 - counts)
    denom = (n1 + n2) * (n1 + n2 - 1)
    if denom == 0:
        return 1.0
    var_u = n1 * n2 / 12.0 * ((n1 + n2 + 1) - (tie_term / denom))
    if var_u <= 0:
        return 1.0
    z = (u - mean_u) / math.sqrt(var_u)
    p = math.erfc(abs(z) / math.sqrt(2))
    return float(p)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build static JSON files for paralog browser")
    parser.add_argument("--matrix", default="full_SLprediction_matrix.csv")
    parser.add_argument("--model", default="Model.csv")
    parser.add_argument("--gene-info", default="gene_with_protein_product.json")
    parser.add_argument("--out", default="dist")
    parser.add_argument("--tmp", default="dist/_tmp")
    parser.add_argument(
        "--limit-rows",
        type=int,
        default=0,
        help="If >0, only process the first N cell lines (for test builds).",
    )
    parser.add_argument(
        "--limit-pairs",
        type=int,
        default=0,
        help="If >0, only process the first N gene pair columns (for test builds).",
    )
    args = parser.parse_args()

    matrix_path = Path(args.matrix)
    model_path = Path(args.model)
    gene_info_path = Path(args.gene_info)
    out_dir = Path(args.out)
    tmp_dir = Path(args.tmp)
    limit_rows = int(args.limit_rows) if args.limit_rows else 0
    limit_pairs = int(args.limit_pairs) if args.limit_pairs else 0

    out_index = out_dir / "index"
    out_top = out_dir / "top"
    out_agg = out_dir / "agg"
    out_disease = out_dir / "disease"

    out_index.mkdir(parents=True, exist_ok=True)
    out_top.mkdir(parents=True, exist_ok=True)
    out_agg.mkdir(parents=True, exist_ok=True)
    out_disease.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    model_map = load_model_aliases(model_path)
    gene_meta = load_gene_metadata(gene_info_path)

    # Read header to get pair list
    with matrix_path.open() as f:
        reader = csv.reader(f)
        header = next(reader)

    if not header or header[0] != "DepMap_ID":
        raise ValueError("Expected first column to be DepMap_ID")

    pair_columns = header[1:]
    if limit_pairs and limit_pairs < len(pair_columns):
        pair_columns = pair_columns[:limit_pairs]
    pair_ids: List[str] = []
    pair_genes: List[Tuple[str, str]] = []
    for col in pair_columns:
        pair_id, gene1, gene2 = pair_id_for(col)
        pair_ids.append(pair_id)
        pair_genes.append((gene1, gene2))

    n_pairs = len(pair_ids)

    # Precompute pair bucket list
    pair_bucket_list = [pair_bucket(pid) for pid in pair_ids]
    pair_buckets_map = {pid: bucket for pid, bucket in zip(pair_ids, pair_bucket_list)}

    # Pass 1: compute mean and std per pair (ignore missing values)
    mean = np.zeros(n_pairs, dtype=np.float64)
    m2 = np.zeros(n_pairs, dtype=np.float64)
    counts = np.zeros(n_pairs, dtype=np.int64)
    depmap_ids: List[str] = []

    with matrix_path.open() as f:
        reader = csv.reader(f)
        next(reader)  # header
        processed = 0
        for row in reader:
            if not row:
                continue
            depmap_id = row[0]
            depmap_ids.append(depmap_id)
            processed += 1
            row_values = row[1:1 + n_pairs]
            arr = np.array([float(x) if x != "" else np.nan for x in row_values], dtype=np.float64)
            mask = ~np.isnan(arr)
            if not mask.any():
                continue
            counts[mask] += 1
            delta = arr - mean
            mean[mask] += delta[mask] / counts[mask]
            delta2 = arr - mean
            m2[mask] += delta[mask] * delta2[mask]
            if limit_rows and processed >= limit_rows:
                break

    std = np.zeros(n_pairs, dtype=np.float64)
    valid_var = counts > 1
    var = np.zeros(n_pairs, dtype=np.float64)
    var[valid_var] = m2[valid_var] / (counts[valid_var] - 1)
    std[valid_var] = np.sqrt(var[valid_var])

    inv_std = np.zeros_like(std)
    nonzero = std > 0
    inv_std[nonzero] = 1.0 / std[nonzero]

    # Pass 2: top-100 per cell line + spill to bucket temp CSVs
    top_pairs_by_cell: Dict[str, List[Dict[str, float]]] = {}
    top_zpairs_by_cell: Dict[str, List[Dict[str, float]]] = {}

    tmp_files = {}
    for i in range(256):
        bucket = f"{i:02x}"
        tmp_path = tmp_dir / f"pair_bucket_{bucket}.csv"
        tmp_files[bucket] = tmp_path.open("w")

    with matrix_path.open() as f:
        reader = csv.reader(f)
        next(reader)
        processed = 0
        for row in reader:
            if not row:
                continue
            depmap_id = row[0]
            processed += 1
            row_values = row[1:1 + n_pairs]
            arr = np.array([float(x) if x != "" else np.nan for x in row_values], dtype=np.float64)
            z = (arr - mean) * inv_std

            # Top 100 by score and by z-score
            top_k = 100
            idx_score = top_k_indices_masked(arr, top_k)
            idx_z = top_k_indices_masked(z, top_k)

            top_pairs_by_cell[depmap_id] = [
                {"pair_id": pair_ids[i], "score": float(arr[i])}
                for i in idx_score
                if math.isfinite(arr[i])
            ]
            top_zpairs_by_cell[depmap_id] = [
                {"pair_id": pair_ids[i], "score": float(arr[i]), "z": float(z[i])}
                for i in idx_z
                if math.isfinite(arr[i]) and math.isfinite(z[i])
            ]

            # Spill per-pair values into bucket temp files
            for i, pair_id in enumerate(pair_ids):
                score_val = arr[i]
                if not math.isfinite(score_val):
                    continue
                bucket = pair_bucket_list[i]
                z_val = z[i]
                tmp_files[bucket].write(
                    f"{pair_id},{depmap_id},{score_val},{z_val}\n"
                )
            if limit_rows and processed >= limit_rows:
                break

    for f in tmp_files.values():
        f.close()

    # Write top-by-cell-line outputs
    write_json(out_top / "top_pairs_by_cell_line.json", top_pairs_by_cell)
    write_json(out_top / "top_zpairs_by_cell_line.json", top_zpairs_by_cell)

    # Reduce each bucket to top 100 per pair
    for i in range(256):
        bucket = f"{i:02x}"
        tmp_path = tmp_dir / f"pair_bucket_{bucket}.csv"
        if not tmp_path.exists():
            continue

        score_heaps: Dict[str, List[Tuple[float, str]]] = defaultdict(list)
        z_heaps: Dict[str, List[Tuple[float, float, str]]] = defaultdict(list)

        with tmp_path.open() as f:
            reader = csv.reader(f)
            for pair_id, depmap_id, score_s, z_s in reader:
                score = float(score_s)
                z_val = float(z_s)
                if not math.isfinite(score):
                    continue

                # score heap
                heap = score_heaps[pair_id]
                item = (score, depmap_id)
                if len(heap) < 100:
                    heapq.heappush(heap, item)
                else:
                    if score > heap[0][0]:
                        heapq.heapreplace(heap, item)

                # z heap
                if math.isfinite(z_val):
                    zheap = z_heaps[pair_id]
                    zitem = (z_val, score, depmap_id)
                    if len(zheap) < 100:
                        heapq.heappush(zheap, zitem)
                    else:
                        if z_val > zheap[0][0]:
                            heapq.heapreplace(zheap, zitem)

        # Convert to sorted outputs
        top_by_score = {}
        for pair_id, heap in score_heaps.items():
            heap.sort(key=lambda x: x[0], reverse=True)
            top_by_score[pair_id] = [
                {"depmap_id": depmap_id, "score": float(score)}
                for score, depmap_id in heap
            ]

        write_json(out_top / f"pair_topcells_bucket_{bucket}.json", top_by_score)

        top_by_z = {}
        for pair_id, heap in z_heaps.items():
            heap.sort(key=lambda x: x[0], reverse=True)
            top_by_z[pair_id] = [
                {"depmap_id": depmap_id, "score": float(score), "z": float(z_val)}
                for z_val, score, depmap_id in heap
            ]

        write_json(out_top / f"pair_topcells_z_bucket_{bucket}.json", top_by_z)

    # Gene aggregates by mean score
    gene_pairs: Dict[str, List[Dict[str, float]]] = defaultdict(list)
    for i, (gene1, gene2) in enumerate(pair_genes):
        mean_val = float(mean[i])
        pair_id = pair_ids[i]
        gene_pairs[gene1].append(
            {"pair_id": pair_id, "partner": gene2, "mean": mean_val}
        )
        gene_pairs[gene2].append(
            {"pair_id": pair_id, "partner": gene1, "mean": mean_val}
        )

    gene_bucket_map: Dict[str, Dict[str, List[Dict[str, float]]]] = defaultdict(dict)
    for gene, entries in gene_pairs.items():
        entries.sort(key=lambda x: x["mean"], reverse=True)
        bucket = gene_bucket(gene)
        gene_bucket_map[bucket][gene] = entries

    for bucket, payload in gene_bucket_map.items():
        write_json(out_agg / f"gene_pairs_by_mean_bucket_{bucket}.json", payload)

    # Disease-specific dependencies (OncotreePrimaryDisease)
    disease_labels = []
    for depmap_id in depmap_ids:
        info = model_map.get(depmap_id)
        disease = ""
        if info:
            disease = (info.get("primary_disease") or "").strip()
        disease_labels.append(disease)

    from collections import Counter

    disease_counts = Counter([d for d in disease_labels if d and d != "Non-Cancerous"])
    diseases = sorted([d for d, c in disease_counts.items() if c > 9])
    if diseases:
        # Load full matrix into memory for disease stats
        data = np.empty((len(depmap_ids), n_pairs), dtype=np.float32)
        with matrix_path.open() as f:
            reader = csv.reader(f)
            next(reader)
            for i, row in enumerate(reader):
                if i >= len(depmap_ids):
                    break
                row_values = row[1:1 + n_pairs]
                arr = np.array([float(x) if x != "" else np.nan for x in row_values], dtype=np.float32)
                data[i] = arr

        labels_arr = np.array(disease_labels)
        known_mask = labels_arr != ""
        pair_active_mask = np.sum(data > 0.1, axis=0) >= 3

        disease_bucket_map: Dict[str, Dict[str, List[Dict[str, float]]]] = defaultdict(dict)
        for disease in diseases:
            idx_a = np.where(labels_arr == disease)[0]
            idx_b = np.where((labels_arr != disease) & known_mask)[0]
            bucket = disease_bucket(disease)
            if idx_a.size < 2 or idx_b.size < 2:
                disease_bucket_map[bucket][disease] = []
                continue

            group_a = data[idx_a]
            group_b = data[idx_b]
            valid_a = np.isfinite(group_a).any(axis=0)
            valid_b = np.isfinite(group_b).any(axis=0)
            valid_cols = valid_a & valid_b

            median_a = np.nanmedian(group_a, axis=0)
            median_b = np.nanmedian(group_b, axis=0)
            diff = median_a - median_b

            candidates = np.where((diff > 0.1) & pair_active_mask & valid_cols)[0]
            results = []
            for i in candidates:
                x = data[idx_a, i]
                y = data[idx_b, i]
                x = x[np.isfinite(x)]
                y = y[np.isfinite(y)]
                if x.size < 2 or y.size < 2:
                    continue
                p_val = mann_whitney_pvalue(x, y)
                if p_val < 0.05:
                    results.append(
                        {
                            "pair_id": pair_ids[i],
                            "median_disease": float(median_a[i]),
                            "median_other": float(median_b[i]),
                            "diff": float(diff[i]),
                            "p_value": float(p_val),
                        }
                    )

            results.sort(key=lambda x: x["diff"], reverse=True)
            disease_bucket_map[bucket][disease] = results

        for bucket, payload in disease_bucket_map.items():
            write_json(out_disease / f"disease_pairs_bucket_{bucket}.json", payload)

    # Index files
    cell_lines_index = []
    alias_cell_lines = {}
    for depmap_id in depmap_ids:
        info = model_map.get(depmap_id)
        if info:
            display = info["display"]
            aliases = info["aliases"]
            lineage = info.get("lineage", "")
            primary_disease = info.get("primary_disease", "")
        else:
            display = depmap_id
            aliases = [depmap_id]
            lineage = ""
            primary_disease = ""

        unique_aliases = list(dict.fromkeys([a for a in aliases if a]))
        cell_lines_index.append(
            {
                "depmap_id": depmap_id,
                "display": display,
                "aliases": unique_aliases,
                "lineage": lineage,
                "primary_disease": primary_disease,
            }
        )

        for alias in unique_aliases:
            norm = normalize_alias(alias)
            if norm and norm not in alias_cell_lines:
                alias_cell_lines[norm] = depmap_id

    # Genes
    genes_index = []
    alias_genes = {}
    gene_symbols = sorted(set([g for pair in pair_genes for g in pair]))
    for symbol in gene_symbols:
        aliases = [symbol]
        meta = gene_meta.get(symbol, {})
        name = meta.get("name")
        if name:
            aliases.append(name)
        unique_aliases = list(dict.fromkeys([a for a in aliases if a]))
        genes_index.append({
            "symbol": symbol,
            "aliases": unique_aliases,
            "name": name or "",
            "hgnc_id": meta.get("hgnc_id", ""),
        })
        for alias in unique_aliases:
            norm = normalize_alias(alias)
            if norm and norm not in alias_genes:
                alias_genes[norm] = symbol

    # Pairs
    pairs_index = []
    alias_pairs = {}
    for pair_id, (gene1, gene2) in zip(pair_ids, pair_genes):
        pairs_index.append({"pair_id": pair_id, "gene1": gene1, "gene2": gene2})
        alias_variants = [
            f"{gene1}_{gene2}",
            f"{gene2}_{gene1}",
            pair_id,
            f"{gene2}__{gene1}",
        ]
        for alias in alias_variants:
            norm = normalize_alias(alias)
            if norm and norm not in alias_pairs:
                alias_pairs[norm] = pair_id

    alias_maps = {
        "cell_lines": alias_cell_lines,
        "genes": alias_genes,
        "pairs": alias_pairs,
        "diseases": {},
    }

    write_json(out_index / "cell_lines.json", cell_lines_index)
    write_json(out_index / "genes.json", genes_index)
    write_json(out_index / "pairs.json", pairs_index)
    write_json(out_index / "pair_buckets.json", pair_buckets_map)

    # Disease index + alias maps
    diseases_index = []
    disease_bucket_map = {}
    alias_diseases = {}
    for disease in diseases:
        diseases_index.append({"disease": disease, "aliases": [disease]})
        norm = normalize_alias(disease)
        if norm and norm not in alias_diseases:
            alias_diseases[norm] = disease
        disease_bucket_map[disease] = disease_bucket(disease)

    alias_maps["diseases"] = alias_diseases
    write_json(out_index / "alias_maps.json", alias_maps)
    write_json(out_index / "diseases.json", diseases_index)
    write_json(out_index / "disease_buckets.json", disease_bucket_map)


if __name__ == "__main__":
    main()
