#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
from pathlib import Path
from glob import glob
from statistics import median
from collections import Counter, defaultdict

import pandas as pd
import matplotlib.pyplot as plt

def read_fasta_lengths(path):
    """Вернёт список длин последовательностей из FASTA (быстро, без biopython)."""
    if not os.path.exists(path):
        return []
    lens = []
    with open(path, 'r') as f:
        length = 0
        for line in f:
            if line.startswith('>'):
                if length:
                    lens.append(length)
                    length = 0
            else:
                s = line.strip()
                if s:
                    length += len(s)
        if length:
            lens.append(length)
    return lens

def gc_percent_from_fasta(path):
    """GC% по нуклеотидным FASTA (genes.fna)."""
    if not os.path.exists(path):
        return None
    gc = 0
    total = 0
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            s = line.strip().upper()
            if not s: 
                continue
            total += len(s)
            gc += s.count('G') + s.count('C')
    return (100.0 * gc / total) if total else None

def n50(lengths):
    """N50 по списку длин (bp/aa)."""
    if not lengths:
        return None
    arr = sorted(lengths, reverse=True)
    half = sum(arr) / 2
    cum = 0
    for L in arr:
        cum += L
        if cum >= half:
            return L
    return arr[-1]

def parse_gff_counts(path):
    """Подсчёт фич из GFF: количество CDS и распределение по цепям."""
    if not os.path.exists(path):
        return 0, 0, 0
    cds = 0
    plus = 0
    minus = 0
    with open(path, 'r') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            typ = parts[2]
            strand = parts[6] if len(parts) > 6 else '.'
            if typ.lower() == 'cds':
                cds += 1
                if strand == '+':
                    plus += 1
                elif strand == '-':
                    minus += 1
    return cds, plus, minus

def summarize_one(sample, assembler, caller, base_dir):
    """Сбор метрик для одной тройки (sample/assembler/caller)."""
    root = Path(base_dir) / assembler / sample / caller
    proteins = root / "proteins.faa"
    genes = root / "genes.fna"
    gff = root / "genes.gff"

    aa_lens = read_fasta_lengths(proteins)
    nt_lens = read_fasta_lengths(genes)

    cds_count, cds_plus, cds_minus = parse_gff_counts(gff)
    gc = gc_percent_from_fasta(genes)

    row = {
        "sample": sample,
        "assembler": assembler,
        "caller": caller,
        # proteins
        "n_proteins": len(aa_lens),
        "aa_total_len": sum(aa_lens) if aa_lens else 0,
        "aa_mean_len": (sum(aa_lens)/len(aa_lens)) if aa_lens else None,
        "aa_median_len": median(aa_lens) if aa_lens else None,
        "aa_min_len": min(aa_lens) if aa_lens else None,
        "aa_max_len": max(aa_lens) if aa_lens else None,
        "aa_N50": n50(aa_lens),
        # genes (nt)
        "n_genes": len(nt_lens),
        "nt_total_len": sum(nt_lens) if nt_lens else 0,
        "nt_mean_len": (sum(nt_lens)/len(nt_lens)) if nt_lens else None,
        "nt_median_len": median(nt_lens) if nt_lens else None,
        "nt_min_len": min(nt_lens) if nt_lens else None,
        "nt_max_len": max(nt_lens) if nt_lens else None,
        "nt_N50": n50(nt_lens),
        "nt_GC_percent": gc,
        # gff
        "gff_cds": cds_count,
        "gff_cds_plus": cds_plus,
        "gff_cds_minus": cds_minus,
        "gff_plus_frac": (cds_plus/cds_count) if cds_count else None,
    }
    return row

def discover_triplets(orfs_root):
    """Найдёт все (assembler/sample/caller), где есть хотя бы proteins.faa."""
    triplets = set()
    for proteins in glob(os.path.join(orfs_root, "*", "*", "*", "proteins.faa")):
        p = Path(proteins)
        caller = p.parent.name
        sample = p.parent.parent.name
        assembler = p.parent.parent.parent.name
        triplets.add((sample, assembler, caller))
    return sorted(triplets)

def make_plots(df, outdir):
    os.makedirs(outdir, exist_ok=True)
    # 1) количество белков по сборщикам/предикторам
    if {"assembler","caller","n_proteins"}.issubset(df.columns):
        agg = df.groupby(["assembler","caller"], as_index=False)["n_proteins"].median()
        plt.figure()
        x = range(len(agg))
        plt.bar(x, agg["n_proteins"])
        plt.xticks(x, [f"{a}\n{c}" for a,c in agg[["assembler","caller"]].itertuples(index=False, name=None)], rotation=0)
        plt.ylabel("Медиана # белков")
        plt.title("Медианное число предсказанных белков по методам")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "median_n_proteins_by_method.png"), dpi=150)
        plt.close()

    # 2) boxplot длин белков по методам (агрегировано по семплам медианами)
    if {"assembler","caller","aa_median_len"}.issubset(df.columns):
        groups = []
        labels = []
        for (a,c), g in df.groupby(["assembler","caller"]):
            labels.append(f"{a}\n{c}")
            groups.append(g["aa_median_len"].dropna())
        if groups:
            plt.figure()
            plt.boxplot(groups, labels=labels)
            plt.ylabel("Медианная длина белков (aa)")
            plt.title("Распределение медианной длины белков")
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, "box_aa_median_len_by_method.png"), dpi=150)
            plt.close()

    # 3) scatter: nt_total_len vs n_genes
    if {"nt_total_len","n_genes"}.issubset(df.columns):
        plt.figure()
        for (a,c), g in df.groupby(["assembler","caller"]):
            plt.scatter(g["nt_total_len"], g["n_genes"], alpha=0.7, label=f"{a}/{c}")
        plt.xlabel("Суммарная длина генов (nt)")
        plt.ylabel("# генов")
        plt.title("nt_total_len vs n_genes")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "scatter_nt_total_vs_n_genes.png"), dpi=150)
        plt.close()

def main():
    ap = argparse.ArgumentParser(description="Сбор статистики по ORF-предсказаниям (proteins.faa / genes.fna / genes.gff).")
    ap.add_argument("--orfs-root", default="orfs", help="Корень с результатами ORF: orfs/{assembler}/{sample}/{caller}/")
    ap.add_argument("--outdir", default="qc/orf_stats", help="Куда писать сводку и графики")
    args = ap.parse_args()

    orfs_root = os.path.abspath(os.path.expanduser(args.orfs_root))
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    os.makedirs(outdir, exist_ok=True)

    triplets = discover_triplets(orfs_root)
    if not triplets:
        raise SystemExit(f"Не найдено triplets в {orfs_root}")

    rows = []
    for sample, assembler, caller in triplets:
        rows.append(summarize_one(sample, assembler, caller, orfs_root))
    df = pd.DataFrame(rows)

    # Сохраняем сводку
    csv_path = os.path.join(outdir, "orf_stats_summary.csv")
    df.to_csv(csv_path, index=False)

    # Пара быстрых агрегаций для удобства
    by_method = df.groupby(["assembler","caller"], as_index=False).agg(
        n_proteins_median=("n_proteins","median"),
        aa_median_len_median=("aa_median_len","median"),
        nt_GC_percent_median=("nt_GC_percent","median"),
        n_genes_median=("n_genes","median"),
    )
    by_method.to_csv(os.path.join(outdir, "by_method_summary.csv"), index=False)

    # Графики
    make_plots(df, os.path.join(outdir, "plots"))

    print(f"[OK] Samples×Methods: {len(df)} записей")
    print(f"[OK] CSV: {csv_path}")
    print(f"[OK] By-method CSV: {os.path.join(outdir, 'by_method_summary.csv')}")
    print(f"[OK] Plots dir: {os.path.join(outdir, 'plots')}")

if __name__ == "__main__":
    main()
