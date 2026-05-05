#!/usr/bin/env bash
# =============================================================================
# SCRIPT 1/2 — Import QIIME2, Cutadapt, DADA2, courbes de raréfaction
# Projet  : vague_project | 16S V4-V5 | Pointe de l'Artillerie
# Données : déjà nettoyées en amont, démarrage direct à QIIME2
# Chemin  : /nvme/bio/data_fungi/vague_project
# Usage   : bash 01_qiime2_preprocess.sh
# =============================================================================

set -euo pipefail

export ROOTDIR="/nvme/bio/data_fungi/vague_project"
export NTHREADS=16
export QIIME2_ENV="qiime2-amplicon-2025.7"
export TMPDIR="${ROOTDIR}/tmp"
export RAWDATA="${ROOTDIR}/00_raw_data"
export DBDIR="${ROOTDIR}/98_databasefiles"
export QDIR="${ROOTDIR}/05_QIIME2"

PRIMER_F="GTGYCAGCMGCCGCGGTAA"
PRIMER_R="CCGYCAATTYMTTTRAGTTT"

log() { echo -e "\n[$(date +'%F %T')] === $* ===\n"; }

mkdir -p "$TMPDIR" "$DBDIR" "${QDIR}/core" "${QDIR}/visual" "${QDIR}/subtables"
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
set -u

log "SCRIPT 1 DÉMARRÉ"

log "Génération manifest et métadonnées"
python3 << 'PYEOF'
import os, re, csv
from collections import Counter

ROOTDIR = "/nvme/bio/data_fungi/vague_project"
RAWDATA = os.path.join(ROOTDIR, "00_raw_data")
DBDIR   = os.path.join(ROOTDIR, "98_databasefiles")
manifest_path = os.path.join(DBDIR, "manifest")
metadata_path = os.path.join(DBDIR, "sample-metadata.tsv")

samples = []
for f in sorted(os.listdir(RAWDATA)):
    if not f.endswith("_R1_001.fastq.gz"):
        continue
    base = f.replace("_R1_001.fastq.gz", "")
    sample_id = re.sub(r'_S\d+$', '', base)
    r1 = os.path.join(RAWDATA, f)
    r2 = os.path.join(RAWDATA, f.replace("_R1_001", "_R2_001"))
    if os.path.exists(r2):
        samples.append((sample_id, r1, r2))

with open(manifest_path, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(["sample-id", "forward-absolute-filepath", "reverse-absolute-filepath"])
    for sid, r1, r2 in samples:
        w.writerow([sid, r1, r2])

headers = [
    "sample-id", "sample_type", "depth_cm", "condition", "location",
    "is_control", "control_type", "description"
]
types = [
    "categorical", "numeric", "categorical", "categorical",
    "categorical", "categorical", "categorical"
]

rows = []
for sid, _, _ in samples:
    sl = sid.lower()
    if "sed" in sl and not sl.startswith("t_"):
        m = re.match(r'^(\d+)', sid)
        depth = m.group(1) if m else ""
        rows.append([sid, "sediment", depth, "NA", "PA", "no", "none", f"Carotte sédiment tranche {depth} cm - Pointe de l'Artillerie"])
    elif "calm" in sl:
        m = re.match(r'^(\d+)', sid)
        rep = m.group(1) if m else "1"
        rows.append([sid, "seawater", "", "calm", "PA", "no", "none", f"Eau de mer temps calme réplicat {rep} - Pointe de l'Artillerie"])
    elif "storm" in sl or "strom" in sl:
        m = re.match(r'^(\d+)', sid)
        rep = m.group(1) if m else "1"
        rows.append([sid, "seawater", "", "storm", "PA", "no", "none", f"Eau de mer tempête réplicat {rep} - Pointe de l'Artillerie"])
    elif "blanc_colonne" in sl:
        rows.append([sid, "negative_control", "", "blank_seawater", "PA", "yes", "blank_seawater", "Blanc colonne d'eau"])
    elif sl == "t_sed_pa":
        rows.append([sid, "negative_control", "", "sediment_extraction", "PA", "yes", "sediment_extraction", "Contrôle négatif extraction sédiment"])
    elif "t_1_filter" in sl:
        rows.append([sid, "negative_control", "", "water_extraction", "PA", "yes", "water_extraction", "Contrôle négatif extraction filtre eau"])
    else:
        rows.append([sid, "unknown", "", "NA", "PA", "no", "none", f"Type inconnu : {sid}"])

with open(metadata_path, 'w') as fh:
    fh.write('\t'.join(headers) + '\n')
    fh.write('#q2:types\t' + '\t'.join(types) + '\n')
    for row in rows:
        fh.write('\t'.join(map(str, row)) + '\n')

print(f"Manifest écrit : {manifest_path}")
print(f"Métadonnées écrites : {metadata_path}")
print(f"Total échantillons : {len(rows)}")
counts = Counter(r[1] for r in rows)
for k, v in counts.items():
    print(f"  {k}: {v}")
PYEOF

cd "${QDIR}/core"
log "Import QIIME2"
conda run -n "$QIIME2_ENV" qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "${DBDIR}/manifest" \
  --output-path "demux_paired.qza" \
  --input-format PairedEndFastqManifestPhred33V2

conda run -n "$QIIME2_ENV" qiime demux summarize \
  --i-data demux_paired.qza \
  --o-visualization ../visual/demux-summary.qzv

log "Cutadapt - suppression amorces"
conda run -n "$QIIME2_ENV" qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux_paired.qza \
  --p-front-f "$PRIMER_F" \
  --p-front-r "$PRIMER_R" \
  --p-discard-untrimmed \
  --p-no-indels \
  --p-overlap 10 \
  --p-cores "$NTHREADS" \
  --o-trimmed-sequences demux_trimmed.qza \
  --verbose 2> "${TMPDIR}/cutadapt.log"

conda run -n "$QIIME2_ENV" qiime demux summarize \
  --i-data demux_trimmed.qza \
  --o-visualization ../visual/demux-trimmed-summary.qzv

log "DADA2"
TRUNC_F=0
TRUNC_R=0
conda run -n "$QIIME2_ENV" qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux_trimmed.qza \
  --p-trunc-len-f "$TRUNC_F" \
  --p-trunc-len-r "$TRUNC_R" \
  --p-n-threads "$NTHREADS" \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

conda run -n "$QIIME2_ENV" qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization ../visual/denoising-stats.qzv

conda run -n "$QIIME2_ENV" qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization ../visual/rep-seqs.qzv

conda run -n "$QIIME2_ENV" qiime feature-table summarize \
  --i-table table.qza \
  --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
  --o-visualization ../visual/table-summary.qzv

log "Arbre phylogénétique pour rarefaction Faith PD"
conda run -n "$QIIME2_ENV" qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --p-n-threads "$NTHREADS" \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree tree.qza

log "Courbes de raréfaction"
MAX_DEPTH=50000
conda run -n "$QIIME2_ENV" qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny tree.qza \
  --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
  --p-max-depth "$MAX_DEPTH" \
  --p-steps 20 \
  --o-visualization ../visual/alpha-rarefaction-curves.qzv

log "SCRIPT 1 TERMINÉ"
echo "Consulter ${QDIR}/visual/table-summary.qzv et alpha-rarefaction-curves.qzv puis choisir RAREFACTION_DEPTH pour le script 2."
