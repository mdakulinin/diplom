#!/usr/bin/env python3
import re
from pathlib import Path
from collections import defaultdict

# Регулярка для извлечения организма из заголовка UniProt:
# берём всё после " OS=" до следующего тега (OX=|GN=|PE=|SV=) или конца строки
OS_RE = re.compile(r"\sOS=([^=]+?)\s(?:OX=|GN=|PE=|SV=|$)")

def extract_protein_name(p: Path) -> str:
    # из имени файла uniprot_{PROT}.fasta получить {PROT}
    name = p.stem  # "uniprot_{PROT}"
    if name.startswith("uniprot_"):
        return name[len("uniprot_"):]
    return name

def extract_os(header: str) -> str:
    m = OS_RE.search(header)
    return m.group(1).strip() if m else None

def main():
    files = sorted(Path(".").glob("uniprot_*.fasta"))
    protein2taxa = defaultdict(set)

    for f in files:
        prot = extract_protein_name(f)
        with f.open("r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith(">"):
                    os_val = extract_os(line)
                    if os_val:
                        protein2taxa[prot].add(os_val)

    # Печать TSV-таблицы: protein \t taxa1, taxa2, ...
    with open('protein_taxa.tsv', 'w') as f:
        f.write("protein\ttaxa\n")
        for prot in sorted(protein2taxa):
            taxa_list = ", ".join(sorted(protein2taxa[prot]))
            f.write(f"{prot}\t{taxa_list}\n")

if __name__ == "__main__":
    main()
