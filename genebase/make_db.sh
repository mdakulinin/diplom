#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import time
import hashlib
import requests
from pathlib import Path
from typing import Optional, Tuple, Iterator, Dict, List, Iterable, DefaultDict
from collections import defaultdict
import urllib.parse

INPUT_GENES_FILE = "headers_unique.txt"                # по одному названию в строке
OUT_DIR = Path("uniprot_downloads")                    # сюда сложим результаты
COMBINED_FASTA_RAW = OUT_DIR / "uniprot_all.raw.fasta" # промежуточный общий файл до дедупликации
COMBINED_FASTA = OUT_DIR / "uniprot_all.dedup.fasta"   # общий итоговый файл после дедупликации
DUP_REPORT = OUT_DIR / "uniprot_all.dedup_report.tsv"  # отчёт о дублях (по хешу и числу записей)
PAGE_SIZE = 500                                        # сколько записей за страницу (макс ~500)
REVIEWED = True                                        # True -> сначала тянем только Swiss-Prot (reviewed:true)

# --- Fallback TrEMBL ---
TREMBL_PER_TAXON_LIMIT = 5         # максимум unreviewed-акцессий на один taxon
SEARCH_PAGE_SIZE = 500             # размер страницы в JSON-поиске (для TrEMBL отбора)
STREAM_BATCH_SIZE = 200            # сколько acc за раз в stream-запрос (FASTA)

TAX_QUERY = (
    "taxonomy_id:1654 OR taxonomy_id:1678 OR taxonomy_id:103892 OR taxonomy_id:84111 OR "
    "taxonomy_id:644652 OR taxonomy_id:239759 OR taxonomy_id:816 OR taxonomy_id:40544 OR "
    "taxonomy_id:375288 OR taxonomy_id:838 OR taxonomy_id:904 OR taxonomy_id:909 OR "
    "taxonomy_id:39148 OR taxonomy_id:207244 OR taxonomy_id:39491 OR taxonomy_id:572511 OR "
    "taxonomy_id:1988 OR taxonomy_id:830 OR taxonomy_id:85028 OR taxonomy_id:1485 OR "
    "taxonomy_id:44258 OR taxonomy_id:39948 OR taxonomy_id:39491 OR taxonomy_id:1350 OR "
    "taxonomy_id:1730 OR taxonomy_id:853 OR taxonomy_id:84108 OR taxonomy_id:1578 OR "
    "taxonomy_id:1243 OR taxonomy_id:1637 OR taxonomy_id:836 OR taxonomy_id:12954 OR "
    "taxonomy_id:906 OR taxonomy_id:543311 OR taxonomy_id:1253 OR taxonomy_id:84127 OR "
    "taxonomy_id:1263 OR taxonomy_id:1301 OR taxonomy_id:409438 OR taxonomy_id:191303 OR "
    "taxonomy_id:29465 OR taxonomy_id:46255 OR taxonomy_id:544 OR taxonomy_id:872 OR "
    "taxonomy_id:547 OR taxonomy_id:561 OR taxonomy_id:209 OR taxonomy_id:570 OR "
    "taxonomy_id:846 OR taxonomy_id:583 OR taxonomy_id:586"
)

# UniProt REST v3
BASE_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
BASE_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"  # для пакетной выгрузки FASTA по запросу
BASE_ENTRY_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"  # запасной вариант

# -------------------- Утилиты запросов --------------------

def sanitize_filename(s: str) -> str:
    s = s.strip().replace(" ", "_")
    return re.sub(r"[^0-9A-Za-z._\\-]+", "_", s)

def build_query(term: str, reviewed: Optional[bool]) -> str:
    """
    Формируем поисковый запрос:
      - фильтр по таксонам (TAX_QUERY)
      - reviewed:true/false (если задано)
      - ищем по точной фразе: название белка, имя гена либо общий текстовый поиск
    """
    term = term.strip().strip('"').strip("'")
    quoted = f'"{term}"'
    parts = [f"({TAX_QUERY})"]
    if reviewed is True:
        parts.append("reviewed:true")
    elif reviewed is False:
        parts.append("reviewed:false")
    name_part = f'(protein_name:{quoted} OR gene_exact:{quoted} OR {quoted})'
    parts.append(name_part)
    return " AND ".join(parts)

def get_next_link(r: requests.Response) -> Optional[str]:
    link = r.headers.get("Link")
    if not link:
        return None
    for chunk in link.split(","):
        if 'rel="next"' in chunk:
            m = re.search(r"<([^>]+)>", chunk)
            if m:
                return m.group(1)
    return None

# -------------------- Прямая выгрузка FASTA (reviewed pass) --------------------

def fetch_fasta(query: str, out_path: Path) -> int:
    """
    Вытягиваем все страницы FASTA по запросу (как сейчас у вас) и пишем в out_path (append mode).
    Возвращаем число записей (приблизительно по количеству заголовков '>').
    """
    params = {"query": query, "format": "fasta", "size": str(PAGE_SIZE)}
    url = BASE_SEARCH_URL
    headers = {"Accept": "text/x-fasta"}
    n_headers = 0

    with requests.Session() as sess, open(out_path, "ab") as out_f:
        while True:
            for attempt in range(6):  # backoff
                try:
                    r = sess.get(url, params=params, headers=headers, timeout=60)
                except requests.RequestException:
                    time.sleep(2 ** attempt)
                    continue

                if r.status_code in (429, 500, 502, 503, 504):
                    time.sleep(2 ** attempt)
                    continue

                if r.status_code == 200:
                    chunk = r.content
                    out_f.write(chunk)
                    n_headers += chunk.count(b"\n>")
                    nxt = get_next_link(r)
                    if nxt:
                        url = nxt
                        params = {}
                        break
                    else:
                        return n_headers + (1 if chunk.startswith(b">") else 0)
                else:
                    time.sleep(2 ** attempt)
            else:
                raise RuntimeError(f"UniProt FASTA request failed after retries. Last URL: {url}")

# -------------------- Fallback TrEMBL: отбор до 5 на таксон --------------------

def search_unreviewed_accessions_grouped(term: str) -> Dict[int, List[str]]:
    """
    Делает JSON-поиск по unreviewed:false? Нет — именно reviewed:false.
    Группирует найденные accession по organismId (taxon id) и оставляет не более TREMBL_PER_TAXON_LIMIT на таксон.
    Возвращает dict: {taxon_id -> [accession, ...]}.
    """
    query = build_query(term, reviewed=False)
    headers = {"Accept": "application/json"}
    url = BASE_SEARCH_URL
    params = {
        "query": query,
        "format": "json",
        "size": str(SEARCH_PAGE_SIZE),
        "fields": "accession,organism_id,reviewed"  # экономим трафик
    }

    grouped: DefaultDict[int, List[str]] = defaultdict(list)

    with requests.Session() as sess:
        while True:
            for attempt in range(6):
                try:
                    r = sess.get(url, params=params, headers=headers, timeout=60)
                except requests.RequestException:
                    time.sleep(2 ** attempt)
                    continue

                if r.status_code in (429, 500, 502, 503, 504):
                    time.sleep(2 ** attempt)
                    continue

                if r.status_code == 200:
                    data = r.json()
                    results = data.get("results", [])
                    for rec in results:
                        # подстраховка: берём только reviewed == False
                        if rec.get("reviewed") is True:
                            continue
                        acc = rec.get("primaryAccession")
                        tax = rec.get("organism", {}).get("taxonId")
                        if not acc or tax is None:
                            continue
                        bucket = grouped[tax]
                        if len(bucket) < TREMBL_PER_TAXON_LIMIT:
                            bucket.append(acc)

                    nxt = get_next_link(r)
                    if nxt:
                        # при переходе по Link-пагинации убираем params и ходим по абсолютному URL
                        url = nxt
                        params = {}
                        # если по всем таксонам уже «добрали» по лимиту — можно прекратить раньше
                        # (но мы не знаем список всех таксонов заранее, поэтому продолжаем до конца или пока ответ не иссякнет)
                        break
                    else:
                        return dict(grouped)
                else:
                    time.sleep(2 ** attempt)
            else:
                # не получилось достать страницу
                return dict(grouped)

def fetch_fasta_for_accessions(accessions: List[str], out_path: Path) -> int:
    """
    Скачивает FASTA по списку accession’ов.
    Пытаемся использовать stream endpoint батчами; если что-то пойдёт не так — отпадаем на поштучную загрузку.
    Возвращает число записей по количеству заголовков.
    """
    if not accessions:
        return 0

    n_headers = 0
    with requests.Session() as sess, open(out_path, "ab") as out_f:
        for i in range(0, len(accessions), STREAM_BATCH_SIZE):
            batch = accessions[i:i + STREAM_BATCH_SIZE]
            # query вида: (accession:ACC1) OR (accession:ACC2) ...
            q = " OR ".join(f"accession:{urllib.parse.quote_plus(acc)}" for acc in batch)
            params = {"query": q, "format": "fasta", "size": str(STREAM_BATCH_SIZE)}
            for attempt in range(6):
                try:
                    r = sess.get(BASE_STREAM_URL, params=params, timeout=60)
                except requests.RequestException:
                    time.sleep(2 ** attempt)
                    continue
                if r.status_code in (429, 500, 502, 503, 504):
                    time.sleep(2 ** attempt)
                    continue
                if r.status_code == 200:
                    chunk = r.content
                    out_f.write(chunk)
                    n_headers += chunk.count(b"\n>") + (1 if chunk.startswith(b">") else 0)
                    break
                else:
                    time.sleep(2 ** attempt)
            else:
                # fallback поштучно
                for acc in batch:
                    url = BASE_ENTRY_FASTA.format(acc=acc)
                    for attempt2 in range(4):
                        try:
                            rr = sess.get(url, timeout=60)
                        except requests.RequestException:
                            time.sleep(2 ** attempt2)
                            continue
                        if rr.status_code == 200:
                            out_f.write(rr.content)
                            n_headers += rr.content.count(b"\n>") + (1 if rr.content.startswith(b">") else 0)
                            break
                        elif rr.status_code in (429, 500, 502, 503, 504):
                            time.sleep(2 ** attempt2)
                            continue
                        else:
                            break
    return n_headers

# ---------------------- FASTA утилиты ----------------------

def fasta_iter(path: Path) -> Iterator[Tuple[str, str]]:
    header = None
    seq_chunks: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks).replace(" ", "").replace("\r", "").replace("\n", "")
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks).replace(" ", "").replace("\r", "").replace("\n", "")

def seq_hash(seq: str) -> str:
    return hashlib.sha256(seq.upper().encode("ascii", errors="ignore")).hexdigest()

def deduplicate_fasta(in_path: Path, out_path: Path, report_path: Path) -> Tuple[int, int]:
    seen: Dict[str, Tuple[str, int]] = {}
    kept = 0
    total = 0

    with open(out_path, "w", encoding="utf-8") as fout:
        for hdr, seq in fasta_iter(in_path):
            total += 1
            h = seq_hash(seq)
            if h not in seen:
                fout.write(f">{hdr}\n")
                for i in range(0, len(seq), 60):
                    fout.write(seq[i:i+60] + "\n")
                seen[h] = (hdr, 1)
                kept += 1
            else:
                kept_hdr, cnt = seen[h]
                seen[h] = (kept_hdr, cnt + 1)

    with open(report_path, "w", encoding="utf-8") as rep:
        rep.write("hash\tkept_header\tn_total\n")
        for h, (hdr, cnt) in seen.items():
            rep.write(f"{h}\t{hdr}\t{cnt}\n")

    return total, kept

# ---------------------- ОСНОВНАЯ ЛОГИКА ----------------------

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    for p in (COMBINED_FASTA_RAW, COMBINED_FASTA, DUP_REPORT, OUT_DIR / "not_found.txt"):
        if p.exists():
            p.unlink()

    not_found_file = OUT_DIR / "not_found.txt"

    with open(INPUT_GENES_FILE, "r", encoding="utf-8") as f:
        terms = [ln.strip() for ln in f if ln.strip()]

    total_records = 0

    for i, term in enumerate(terms, 1):
        base_name = sanitize_filename(term)
        out_path = OUT_DIR / f"uniprot_{base_name}.fasta"

        print(f"[{i}/{len(terms)}] '{term}' -> запрашиваем UniProt (reviewed:true)...")
        q_reviewed = build_query(term, reviewed=True)
        try:
            n = fetch_fasta(q_reviewed, out_path)
        except Exception as e:
            print(f"   ! Ошибка для '{term}' (reviewed): {e}")
            n = 0

        if n == 0:
            # --- Fallback TrEMBL ---
            print(f"   · Reviewed не найдено. Делаю fallback на TrEMBL (reviewed:false) c лимитом {TREMBL_PER_TAXON_LIMIT}/taxon...")
            grouped = search_unreviewed_accessions_grouped(term)
            # собираем все отобранные acc
            accs: List[str] = []
            for tax, lst in grouped.items():
                accs.extend(lst)
            if accs:
                # перед скачиванием убедимся, что файл пуст (мы писали в append в reviewed-проходе)
                if out_path.exists() and out_path.stat().st_size > 0:
                    out_path.unlink()
                n = fetch_fasta_for_accessions(accs, out_path)
                print(f"   · TrEMBL добавлен: {n} последовательностей.")
            else:
                print(f"   · TrEMBL тоже ничего не дал.")
                # фиксация "не найдено"
                try:
                    if out_path.exists():
                        out_path.unlink()
                except Exception:
                    pass
                with open(not_found_file, "a", encoding="utf-8") as nf:
                    nf.write(term + "\n")
                continue
        else:
            print(f"   · Найдено (reviewed) записей: {n}")

        # добавим в общий сырой файл
        with open(COMBINED_FASTA_RAW, "ab") as combined, open(out_path, "rb") as part:
            combined.write(part.read())

        total_records += n
        time.sleep(0.2)

    if not COMBINED_FASTA_RAW.exists() or COMBINED_FASTA_RAW.stat().st_size == 0:
        print("\nГотово, но записей не найдено.")
        if not_found_file.exists():
            print(f"Список генов/белков без результатов: {not_found_file}")
        return

    print("\nДедупликация по аминокислотной последовательности...")
    n_in, n_out = deduplicate_fasta(COMBINED_FASTA_RAW, COMBINED_FASTA, DUP_REPORT)

    print("\nГотово!")
    print(f"Суммарно скачано записей (оценка по заголовкам): ~{total_records}")
    print(f"FASTA до дедупликации: {COMBINED_FASTA_RAW} (записей прочитано: {n_in})")
    print(f"Итоговый файл без дублей: {COMBINED_FASTA} (уникальных последовательностей: {n_out})")
    print(f"Отчёт по дублям: {DUP_REPORT}")
    if not_found_file.exists():
        print(f"Список генов/белков без результатов: {not_found_file}")
    print(f"Отдельные файлы по запросам лежат в: {OUT_DIR.resolve()}")

if __name__ == "__main__":
    main()
