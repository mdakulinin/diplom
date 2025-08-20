#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Пути
ROOT="/home/galanova/students/misha/diplom/ORF_pred/orfs"
DB="/home/galanova/students/misha/diplom/genebase/final_base.fasta"
OUTROOT="/home/galanova/students/misha/diplom/blastp"
THREADS=32

# Проверяем, есть ли готовый индекс для базы
if [ ! -f "${DB}.pin" ]; then
    echo "[INFO] Создаю BLAST базу..."
    makeblastdb -in "$DB" -dbtype prot -out "$DB"
fi

# Обходим все proteins.faa
find "$ROOT" -type f -name "proteins.faa" | while read -r faa; do
    # Относительный путь относительно ROOT
    relpath="${faa#$ROOT/}"
    # Убираем файл, оставляем только подкаталоги
    subdir="$(dirname "$relpath")"
    # Папка для результата
    outdir="$OUTROOT/$subdir"
    mkdir -p "$outdir"
    outfile="$outdir/results.tsv"

    echo "[INFO] BLASTP: $faa -> $outfile"

    blastp -query "$faa" -db "$DB" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle" \
        -num_threads $THREADS \
    | awk '($3 > 60) && ($4 >= 0.9 * (($13 < $14) ? $13 : $14))' > "$outfile"
done

echo "ГОТОВО. Результаты в $OUTROOT"
