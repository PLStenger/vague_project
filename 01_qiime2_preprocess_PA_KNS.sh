#!/usr/bin/env bash
# =============================================================================
# SCRIPT 1/2 — Import QIIME2, Cutadapt, DADA2, décontamination par témoins
# Projet : vague_project | 16S V4-V5 | Pointe d'Artillerie + Koniambo (KNS)
# Données : déjà nettoyées en amont, démarrage direct à QIIME2
# Chemin : /nvme/bio/data_fungi/vague_project
# Usage : bash 01_qiime2_preprocess.sh
# =============================================================================

#set -euo pipefail

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

log "Génération manifest et métadonnées (PA + KNS)"
python3 << 'PYEOF'
import os
import re
import csv
from collections import Counter

ROOTDIR = "/nvme/bio/data_fungi/vague_project"
RAWDATA = os.path.join(ROOTDIR, "00_raw_data")
DBDIR = os.path.join(ROOTDIR, "98_databasefiles")
manifest_path = os.path.join(DBDIR, "manifest")
metadata_path = os.path.join(DBDIR, "sample-metadata.tsv")

SEARCH_DIRS = [
    RAWDATA,
    os.path.join(RAWDATA, "02_Carotte_KNS"),
]

def infer_location(sample_id: str, rel_dir: str) -> str:
    sl = sample_id.lower()
    rdl = rel_dir.lower()
    if "_kns" in sl or sl.startswith("kns_") or "carotte_kns" in rdl:
        return "KNS"
    if "_pa" in sl:
        return "PA"
    return "unknown"

def classify_sample(sample_id: str, rel_dir: str):
    sl = sample_id.lower()
    location = infer_location(sample_id, rel_dir)

    # Témoins négatifs
    if "blanc_colonne" in sl:
        return {
            "sample_type": "negative_control",
            "depth_cm": "",
            "condition": "blank_seawater",
            "location": location if location != "unknown" else "PA",
            "is_control": "yes",
            "control_type": "blank_seawater",
            "description": "Blanc colonne d'eau"
        }

    if sl == "t_sed_pa":
        return {
            "sample_type": "negative_control",
            "depth_cm": "",
            "condition": "sediment_extraction",
            "location": "PA",
            "is_control": "yes",
            "control_type": "sediment_extraction",
            "description": "Contrôle négatif extraction sédiment - Pointe d'Artillerie"
        }

    if "t_1_filter" in sl:
        return {
            "sample_type": "negative_control",
            "depth_cm": "",
            "condition": "water_extraction",
            "location": "PA",
            "is_control": "yes",
            "control_type": "water_extraction",
            "description": "Contrôle négatif extraction filtre eau - Pointe d'Artillerie"
        }

    # Nouveaux témoins KNS notés T-
    if sl.startswith("t-_sed"):
        return {
            "sample_type": "negative_control",
            "depth_cm": "",
            "condition": "sediment_extraction",
            "location": "KNS",
            "is_control": "yes",
            "control_type": "sediment_extraction",
            "description": "Contrôle négatif extraction sédiment - Koniambo"
        }

    if sl.startswith("t-_filtre") or sl.startswith("t-_filter"):
        return {
            "sample_type": "negative_control",
            "depth_cm": "",
            "condition": "water_extraction",
            "location": "KNS",
            "is_control": "yes",
            "control_type": "water_extraction",
            "description": "Contrôle négatif extraction filtre eau - Koniambo"
        }

    # Eau de mer KNS
    if sl.startswith("kns_sw"):
        m = re.search(r'kns_sw(\d+)', sl)
        rep = m.group(1) if m else ""
        return {
            "sample_type": "seawater",
            "depth_cm": "",
            "condition": "NA",
            "location": "KNS",
            "is_control": "no",
            "control_type": "none",
            "description": f"Eau de mer Koniambo réplicat {rep}".strip()
        }

    # Eau de mer PA calme/tempête
    if "calm" in sl:
        m = re.match(r'^(\d+)', sample_id)
        rep = m.group(1) if m else "1"
        return {
            "sample_type": "seawater",
            "depth_cm": "",
            "condition": "calm",
            "location": "PA",
            "is_control": "no",
            "control_type": "none",
            "description": f"Eau de mer temps calme réplicat {rep} - Pointe d'Artillerie"
        }

    if "storm" in sl or "strom" in sl:
        m = re.match(r'^(\d+)', sample_id)
        rep = m.group(1) if m else "1"
        return {
            "sample_type": "seawater",
            "depth_cm": "",
            "condition": "storm",
            "location": "PA",
            "is_control": "no",
            "control_type": "none",
            "description": f"Eau de mer tempête réplicat {rep} - Pointe d'Artillerie"
        }

    # Sédiments PA et KNS
    if "sed" in sl:
        m = re.match(r'^(\d+)', sample_id)
        depth = m.group(1) if m else ""
        if location == "PA":
            desc = f"Carotte sédiment tranche {depth} cm - Pointe d'Artillerie".strip()
        elif location == "KNS":
            desc = f"Carotte sédiment tranche {depth} cm - Koniambo".strip()
        else:
            desc = f"Carotte sédiment tranche {depth} cm".strip()

        return {
            "sample_type": "sediment",
            "depth_cm": depth,
            "condition": "NA",
            "location": location,
            "is_control": "no",
            "control_type": "none",
            "description": desc
        }

    return {
        "sample_type": "unknown",
        "depth_cm": "",
        "condition": "NA",
        "location": location,
        "is_control": "no",
        "control_type": "none",
        "description": f"Type inconnu : {sample_id}"
    }

samples = []
seen = set()

for base_dir in SEARCH_DIRS:
    if not os.path.isdir(base_dir):
        continue

    for f in sorted(os.listdir(base_dir)):
        if not f.endswith("_R1_001.fastq.gz"):
            continue

        r1 = os.path.join(base_dir, f)
        r2 = os.path.join(base_dir, f.replace("_R1_001", "_R2_001"))
        if not os.path.exists(r2):
            continue

        base = f.replace("_R1_001.fastq.gz", "")
        sample_id = re.sub(r'_S\d+$', '', base)

        if sample_id in seen:
            raise RuntimeError(f"Sample ID dupliqué détecté : {sample_id}")

        seen.add(sample_id)
        rel_dir = os.path.relpath(base_dir, RAWDATA)
        samples.append((sample_id, r1, r2, rel_dir))

samples.sort(key=lambda x: x[0])

with open(manifest_path, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(["sample-id", "forward-absolute-filepath", "reverse-absolute-filepath"])
    for sid, r1, r2, _ in samples:
        w.writerow([sid, r1, r2])

headers = [
    "sample-id",
    "sample_type",
    "depth_cm",
    "condition",
    "location",
    "is_control",
    "control_type",
    "description"
]

types = [
    "categorical",
    "numeric",
    "categorical",
    "categorical",
    "categorical",
    "categorical",
    "categorical"
]

rows = []
for sid, _, _, rel_dir in samples:
    ann = classify_sample(sid, rel_dir)
    rows.append([
        sid,
        ann["sample_type"],
        ann["depth_cm"],
        ann["condition"],
        ann["location"],
        ann["is_control"],
        ann["control_type"],
        ann["description"]
    ])

with open(metadata_path, 'w') as fh:
    fh.write('\t'.join(headers) + '\n')
    fh.write('#q2:types\t' + '\t'.join(types) + '\n')
    for row in rows:
        fh.write('\t'.join(map(str, row)) + '\n')

print(f"Manifest écrit : {manifest_path}")
print(f"Métadonnées écrites : {metadata_path}")
print(f"Total échantillons : {len(rows)}")

counts_type = Counter(r[1] for r in rows)
counts_loc = Counter(r[4] for r in rows)
counts_ctrl = Counter(r[6] for r in rows if r[5] == "yes")

print("\nComptage par type :")
for k, v in sorted(counts_type.items()):
    print(f"  {k}: {v}")

print("\nComptage par site :")
for k, v in sorted(counts_loc.items()):
    print(f"  {k}: {v}")

print("\nContrôles négatifs :")
for k, v in sorted(counts_ctrl.items()):
    print(f"  {k}: {v}")
PYEOF

cd "${QDIR}/core"

log "Import QIIME2"
conda run -n "$QIIME2_ENV" qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "${DBDIR}/manifest" \
  --output-path demux_paired.qza \
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
  --o-table table_raw.qza \
  --o-representative-sequences rep-seqs_raw.qza \
  --o-denoising-stats denoising-stats.qza

conda run -n "$QIIME2_ENV" qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization ../visual/denoising-stats.qzv

conda run -n "$QIIME2_ENV" qiime feature-table tabulate-seqs \
  --i-data rep-seqs_raw.qza \
  --o-visualization ../visual/rep-seqs-raw.qzv

conda run -n "$QIIME2_ENV" qiime feature-table summarize \
  --i-table table_raw.qza \
  --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
  --o-visualization ../visual/table-summary-raw.qzv

log "Extraction des ASV présentes dans les témoins négatifs"
NEG_IDS="${TMPDIR}/negative-control-ids.txt"
NEG_FEATURES_TSV="${TMPDIR}/negative-control-features.tsv"
NEG_FEATURES_QZA="${TMPDIR}/negative-control-features.qza"

awk -F '\t' '
BEGIN {OFS="\t"}
NR==1 {
  for (i=1; i<=NF; i++) col[$i]=i
  next
}
NR==2 && $1=="#q2:types" {next}
$col["is_control"]=="yes" {print $col["sample-id"]}
' "${DBDIR}/sample-metadata.tsv" | sort -u > "${NEG_IDS}"

if [[ ! -s "${NEG_IDS}" ]]; then
  echo "Aucun témoin négatif détecté dans les métadonnées." >&2
  exit 1
fi

conda run -n "$QIIME2_ENV" qiime feature-table filter-samples \
  --i-table table_raw.qza \
  --m-metadata-file "${NEG_IDS}" \
  --o-filtered-table table_negative_controls_only.qza

conda run -n "$QIIME2_ENV" qiime feature-table summarize \
  --i-table table_negative_controls_only.qza \
  --o-visualization ../visual/table-negative-controls-summary.qzv

conda run -n "$QIIME2_ENV" qiime tools export \
  --input-path table_negative_controls_only.qza \
  --output-path "${TMPDIR}/table_negative_controls_export"

BIOM_FILE=$(find "${TMPDIR}/table_negative_controls_export" -name 'feature-table.biom' | head -n 1)
if [[ -z "${BIOM_FILE}" ]]; then
  echo "Impossible de trouver le fichier BIOM exporté des témoins négatifs." >&2
  exit 1
fi

biom convert \
  -i "${BIOM_FILE}" \
  -o "${TMPDIR}/table_negative_controls.tsv" \
  --to-tsv

python3 << 'PYEOF'
import os

tmpdir = "/nvme/bio/data_fungi/vague_project/tmp"
infile = os.path.join(tmpdir, "table_negative_controls.tsv")
outfile = os.path.join(tmpdir, "negative-control-features.tsv")

features = []
with open(infile) as fh:
    for line in fh:
        if not line.strip():
            continue
        if line.startswith("#"):
            continue
        if line.lower().startswith("feature id"):
            continue
        parts = line.rstrip("\n").split("\t")
        fid = parts[0]
        vals = parts[1:]
        total = 0
        for x in vals:
            try:
                total += float(x)
            except ValueError:
                pass
        if total > 0:
            features.append(fid)

with open(outfile, "w") as out:
    out.write("Feature ID\n")
    for fid in features:
        out.write(f"{fid}\n")

print(f"Nombre d'ASV détectées dans les témoins négatifs : {len(features)}")
PYEOF

if [[ ! -s "${NEG_FEATURES_TSV}" ]]; then
  echo "La liste des ASV de témoins négatifs est vide ou absente." >&2
  exit 1
fi

conda run -n "$QIIME2_ENV" qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path "${NEG_FEATURES_TSV}" \
  --output-path "${NEG_FEATURES_QZA}"

log "Suppression des ASV détectées dans les témoins négatifs"
conda run -n "$QIIME2_ENV" qiime feature-table filter-features \
  --i-table table_raw.qza \
  --m-metadata-file "${NEG_FEATURES_QZA}" \
  --p-exclude-ids \
  --o-filtered-table table.qza

conda run -n "$QIIME2_ENV" qiime feature-table filter-seqs \
  --i-data rep-seqs_raw.qza \
  --m-metadata-file "${NEG_FEATURES_QZA}" \
  --p-exclude-ids \
  --o-filtered-data rep-seqs.qza

conda run -n "$QIIME2_ENV" qiime feature-table summarize \
  --i-table table.qza \
  --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
  --o-visualization ../visual/table-summary-decontaminated.qzv

conda run -n "$QIIME2_ENV" qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization ../visual/rep-seqs-decontaminated.qzv

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
echo "Consulter ${QDIR}/visual/table-summary-decontaminated.qzv et alpha-rarefaction-curves.qzv puis choisir RAREFACTION_DEPTH pour le script 2."
