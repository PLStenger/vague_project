#!/usr/bin/env bash
# =============================================================================
# SCRIPT 1/2 — Import QIIME2, Cutadapt, DADA2, décontamination par témoins
# Projet : vague_project | 16S V4-V5 | Pointe d'Artillerie + Koniambo (KNS)
# Usage : bash 01_qiime2_preprocess_PA_KNS.sh
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

# ===========================================================================
# BLOC 1 — Génération manifest et métadonnées (PA + KNS)
# DÉJÀ EXÉCUTÉ AVEC SUCCÈS : manifest et sample-metadata.tsv existent
# ===========================================================================
# log "Génération manifest et métadonnées (PA + KNS)"
# python3 << 'PYEOF'
# ... (code Python inchangé)
# PYEOF

# ===========================================================================
# BLOC 2 — Import QIIME2
# DÉJÀ EXÉCUTÉ AVEC SUCCÈS : demux_paired.qza et demux-summary.qzv existent
# ===========================================================================
# log "Import QIIME2"
# conda run -n "$QIIME2_ENV" qiime tools import \
#   --type 'SampleData[PairedEndSequencesWithQuality]' \
#   --input-path "${DBDIR}/manifest" \
#   --output-path "${QDIR}/core/demux_paired.qza" \
#   --input-format PairedEndFastqManifestPhred33V2
#
# conda run -n "$QIIME2_ENV" qiime demux summarize \
#   --i-data "${QDIR}/core/demux_paired.qza" \
#   --o-visualization "${QDIR}/visual/demux-summary.qzv"

# ===========================================================================
# BLOC 3 — Cutadapt
# DÉJÀ EXÉCUTÉ AVEC SUCCÈS : demux_trimmed.qza et demux-trimmed-summary.qzv existent
# ===========================================================================
# log "Cutadapt - suppression amorces"
# conda run -n "$QIIME2_ENV" qiime cutadapt trim-paired \
#   --i-demultiplexed-sequences "${QDIR}/core/demux_paired.qza" \
#   --p-front-f "$PRIMER_F" \
#   --p-front-r "$PRIMER_R" \
#   --p-discard-untrimmed \
#   --p-no-indels \
#   --p-overlap 10 \
#   --p-cores "$NTHREADS" \
#   --o-trimmed-sequences "${QDIR}/core/demux_trimmed.qza" \
#   --verbose 2> "${TMPDIR}/cutadapt.log"
#
# conda run -n "$QIIME2_ENV" qiime demux summarize \
#   --i-data "${QDIR}/core/demux_trimmed.qza" \
#   --o-visualization "${QDIR}/visual/demux-trimmed-summary.qzv"

# ===========================================================================
# BLOC 4 — DADA2
# DÉJÀ EXÉCUTÉ AVEC SUCCÈS : table_raw.qza, rep-seqs_raw.qza, denoising-stats.qza existent
# ===========================================================================
# log "DADA2"
# conda run -n "$QIIME2_ENV" qiime dada2 denoise-paired \
#   --i-demultiplexed-seqs "${QDIR}/core/demux_trimmed.qza" \
#   --p-trunc-len-f 0 \
#   --p-trunc-len-r 0 \
#   --p-n-threads "$NTHREADS" \
#   --o-table "${QDIR}/core/table_raw.qza" \
#   --o-representative-sequences "${QDIR}/core/rep-seqs_raw.qza" \
#   --o-denoising-stats "${QDIR}/core/denoising-stats.qza"

# ===========================================================================
# BLOC 5 — Visualisations post-DADA2 (raw)
# DÉJÀ EXÉCUTÉ AVEC SUCCÈS : denoising-stats.qzv, rep-seqs-raw.qzv, table-summary-raw.qzv existent
# ===========================================================================
# conda run -n "$QIIME2_ENV" qiime metadata tabulate \
#   --m-input-file "${QDIR}/core/denoising-stats.qza" \
#   --o-visualization "${QDIR}/visual/denoising-stats.qzv"
#
# conda run -n "$QIIME2_ENV" qiime feature-table tabulate-seqs \
#   --i-data "${QDIR}/core/rep-seqs_raw.qza" \
#   --o-visualization "${QDIR}/visual/rep-seqs-raw.qzv"
#
# conda run -n "$QIIME2_ENV" qiime feature-table summarize \
#   --i-table "${QDIR}/core/table_raw.qza" \
#   --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
#   --o-visualization "${QDIR}/visual/table-summary-raw.qzv"

# ===========================================================================
# BLOC 6 — Extraction des ASV des témoins négatifs (CORRIGÉ)
# ERREUR PRÉCÉDENTE : negative-control-ids.txt sans header TSV valide
# CORRECTION : --p-where "[is_control]='yes'" directement sur sample-metadata.tsv
# À LANCER
# ===========================================================================

cd "${QDIR}/core"

log "Extraction des ASV présentes dans les témoins négatifs (méthode corrigée)"

# Étape 6a : filtrer la table pour ne garder que les témoins négatifs
# On utilise directement le sample-metadata.tsv existant avec --p-where
conda run -n "$QIIME2_ENV" qiime feature-table filter-samples \
  --i-table table_raw.qza \
  --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
  --p-where "[is_control]='yes'" \
  --o-filtered-table table_negative_controls_only.qza

conda run -n "$QIIME2_ENV" qiime feature-table summarize \
  --i-table table_negative_controls_only.qza \
  --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
  --o-visualization ../visual/table-negative-controls-summary.qzv

# Étape 6b : exporter la table des témoins en BIOM puis TSV
conda run -n "$QIIME2_ENV" qiime tools export \
  --input-path table_negative_controls_only.qza \
  --output-path "${TMPDIR}/table_negative_controls_export"

biom convert \
  -i "${TMPDIR}/table_negative_controls_export/feature-table.biom" \
  -o "${TMPDIR}/table_negative_controls.tsv" \
  --to-tsv

# Étape 6c : extraire la liste des feature-id présentes dans les témoins
# On génère un TSV avec header "feature-id" valide pour QIIME2
python3 << 'PYEOF'
import os

tmpdir = "/nvme/bio/data_fungi/vague_project/tmp"
infile  = os.path.join(tmpdir, "table_negative_controls.tsv")
outfile = os.path.join(tmpdir, "neg_ctrl_features_metadata.tsv")

features = []
with open(infile) as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        # Sauter les lignes de commentaire BIOM (#Constructed from ...)
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        fid = parts[0]
        # Sauter la ligne de header "OTU ID  sample1  sample2 ..."
        if fid.lower() in ("otu id", "feature id", "feature-id", "id", "sample-id"):
            continue
        vals = parts[1:]
        total = sum(float(x) for x in vals if x not in ("", "NA"))
        if total > 0:
            features.append(fid)

# Écrire un TSV avec header "feature-id" + colonne factice pour QIIME2 metadata
with open(outfile, "w") as out:
    out.write("feature-id\tin_negative_control\n")
    for fid in features:
        out.write(f"{fid}\tyes\n")

print(f"Nombre d'ASV détectées dans les témoins négatifs : {len(features)}")
print(f"Fichier écrit : {outfile}")
PYEOF

# Vérifier que le fichier est non vide
NEG_FEATURES_META="${TMPDIR}/neg_ctrl_features_metadata.tsv"
N_FEATURES=$(tail -n +2 "${NEG_FEATURES_META}" | wc -l)
echo "Nombre d'ASV contaminantes à exclure : ${N_FEATURES}"

if [[ "${N_FEATURES}" -eq 0 ]]; then
  echo "AVERTISSEMENT : aucune ASV détectée dans les témoins négatifs."
  echo "Vérifier table_negative_controls.tsv et sample-metadata.tsv"
  echo "La table brute sera copiée sans filtrage."
  cp table_raw.qza table.qza
  cp rep-seqs_raw.qza rep-seqs.qza
else
  # ===========================================================================
  # BLOC 7 — Suppression des ASV contaminantes de la table principale
  # À LANCER
  # ===========================================================================
  log "Suppression des ASV détectées dans les témoins négatifs"

  # Filtrer les features (ASV) de la table
  conda run -n "$QIIME2_ENV" qiime feature-table filter-features \
    --i-table table_raw.qza \
    --m-metadata-file "${NEG_FEATURES_META}" \
    --p-exclude-ids \
    --o-filtered-table table.qza

  # Filtrer les séquences représentatives en conséquence
  conda run -n "$QIIME2_ENV" qiime feature-table filter-seqs \
    --i-data rep-seqs_raw.qza \
    --m-metadata-file "${NEG_FEATURES_META}" \
    --p-exclude-ids \
    --o-filtered-data rep-seqs.qza
fi

# ===========================================================================
# BLOC 8 — Visualisations post-décontamination
# À LANCER
# ===========================================================================
log "Visualisations post-décontamination"

conda run -n "$QIIME2_ENV" qiime feature-table summarize \
  --i-table table.qza \
  --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
  --o-visualization ../visual/table-summary-decontaminated.qzv

conda run -n "$QIIME2_ENV" qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization ../visual/rep-seqs-decontaminated.qzv

# ===========================================================================
# BLOC 9 — Arbre phylogénétique
# À LANCER
# ===========================================================================
log "Arbre phylogénétique (mafft + fasttree)"

conda run -n "$QIIME2_ENV" qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --p-n-threads "$NTHREADS" \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree tree.qza

# ===========================================================================
# BLOC 10 — Courbes de raréfaction
# Adapter MAX_DEPTH après avoir consulté table-summary-decontaminated.qzv
# À LANCER
# ===========================================================================
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
echo ""
echo "Fichiers produits dans ${QDIR}/visual/ :"
echo "  - table-negative-controls-summary.qzv  → vérifier les témoins"
echo "  - table-summary-decontaminated.qzv     → choisir RAREFACTION_DEPTH pour script 2"
echo "  - alpha-rarefaction-curves.qzv         → confirmer le seuil de raréfaction"
