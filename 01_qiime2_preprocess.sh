#!/usr/bin/env bash
# =============================================================================
# SCRIPT 1/2 — Import QIIME2, Cutadapt, DADA2, Courbes de raréfaction
# Projet  : vague_project | 16S V4-V5 | Pointe de l'Artillerie
# Usage   : bash 01_qiime2_preprocess.sh
# Résultat: Visualisations dans 05_QIIME2/visual/ à consulter avant le script 2
# =============================================================================

set -euo pipefail

# ============================================================
# VARIABLES — modifier si besoin
# ============================================================
export ROOTDIR="/nvme/bio/data_fungi/vague_project"
export NTHREADS=16
export QIIME2_ENV="qiime2-amplicon-2025.7"
export TMPDIR="${ROOTDIR}/tmp"

export RAWDATA="${ROOTDIR}/00_raw_data"
export DBDIR="${ROOTDIR}/98_databasefiles"
export QDIR="${ROOTDIR}/05_QIIME2"

# Amorces V4-V5 (515F / 926R)
PRIMER_F="GTGYCAGCMGCCGCGGTAA"
PRIMER_R="CCGYCAATTYMTTTRAGTTT"

# ============================================================
# INITIALISATION
# ============================================================
log() { echo -e "\n[$(date +'%F %T')] ===  $* ===\n"; }

mkdir -p "$TMPDIR" "$DBDIR" \
    "${QDIR}/core" \
    "${QDIR}/visual" \
    "${QDIR}/subtables"

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
set -u

log "SCRIPT 1 DÉMARRÉ"

# ============================================================
# ÉTAPE 00 — GÉNÉRATION MANIFEST & MÉTADONNÉES QIIME2
# ============================================================
log "Génération manifest et métadonnées"

python3 << 'PYEOF'
import os, re, csv

ROOTDIR = "/nvme/bio/data_fungi/vague_project"
RAWDATA = os.path.join(ROOTDIR, "00_raw_data")
DBDIR   = os.path.join(ROOTDIR, "98_databasefiles")

manifest_path = os.path.join(DBDIR, "manifest")
metadata_path = os.path.join(DBDIR, "sample-metadata.tsv")

# ---- Scan des paires R1/R2 dans 00_raw_data ----
samples = []
for f in sorted(os.listdir(RAWDATA)):
    if not f.endswith("_R1_001.fastq.gz"):
        continue
    # Extraire le sample-id : retirer _S{num}_R1_001.fastq.gz
    base = f.replace("_R1_001.fastq.gz", "")
    sample_id = re.sub(r'_S\d+$', '', base)
    r1 = os.path.join(RAWDATA, f)
    r2 = os.path.join(RAWDATA, f.replace("_R1_001", "_R2_001"))
    if os.path.exists(r2):
        samples.append((sample_id, r1, r2))
    else:
        print(f"[WARN] R2 manquant pour {f}, ignoré")

print(f"  {len(samples)} paires R1/R2 trouvées")

# ---- Écriture du manifest QIIME2 ----
with open(manifest_path, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(["sample-id", "forward-absolute-filepath", "reverse-absolute-filepath"])
    for sid, r1, r2 in samples:
        w.writerow([sid, r1, r2])
print(f"  Manifest écrit : {manifest_path}")

# ---- Construction des métadonnées ----
# Colonnes : sample-id | sample_type | depth_cm | condition | location | description
headers = ["sample-id", "sample_type", "depth_cm", "condition", "location", "description"]
# Types QIIME2 (hors sample-id)
types   = ["categorical",  "numeric",   "categorical", "categorical", "categorical"]

rows = []
for sid, _, _ in samples:
    sl = sid.lower()

    # --- Sédiments (ex: 11_sed_PA, 50_sed_PA …) ---
    if "sed" in sl and not sl.startswith("t_"):
        m = re.match(r'^(\d+)', sid)
        depth = m.group(1) if m else "NA"
        rows.append({
            "sample-id":   sid,
            "sample_type": "sediment",
            "depth_cm":    depth,
            "condition":   "NA",
            "location":    "PA",
            "description": f"Carotte sédiment tranche {depth} cm - Pointe de l'Artillerie"
        })

    # --- Colonne d'eau calme ---
    elif "calm" in sl:
        m = re.match(r'^(\d+)', sid)
        rep = m.group(1) if m else "1"
        rows.append({
            "sample-id":   sid,
            "sample_type": "seawater",
            "depth_cm":    "NA",
            "condition":   "calm",
            "location":    "PA",
            "description": f"Eau de mer temps calme réplicat {rep} - Pointe de l'Artillerie"
        })

    # --- Colonne d'eau tempête (storm ou strom typo) ---
    elif "storm" in sl or "strom" in sl:
        m = re.match(r'^(\d+)', sid)
        rep = m.group(1) if m else "1"
        rows.append({
            "sample-id":   sid,
            "sample_type": "seawater",
            "depth_cm":    "NA",
            "condition":   "storm",
            "location":    "PA",
            "description": f"Eau de mer tempête réplicat {rep} - Pointe de l'Artillerie"
        })

    # --- Blanc colonne d'eau ---
    elif "blanc_colonne" in sl or "blanc" in sl:
        rows.append({
            "sample-id":   sid,
            "sample_type": "seawater",
            "depth_cm":    "NA",
            "condition":   "blank_seawater",
            "location":    "PA",
            "description": "Blanc extraction colonne d'eau - Pointe de l'Artillerie"
        })

    # --- Contrôle négatif sédiment ---
    elif sl == "t_sed_pa":
        rows.append({
            "sample-id":   sid,
            "sample_type": "negative_control",
            "depth_cm":    "NA",
            "condition":   "sediment_extraction",
            "location":    "PA",
            "description": "Contrôle négatif extraction sédiment"
        })

    # --- Contrôle négatif filtre eau ---
    elif "t_1_filter" in sl:
        rows.append({
            "sample-id":   sid,
            "sample_type": "negative_control",
            "depth_cm":    "NA",
            "condition":   "water_extraction",
            "location":    "PA",
            "description": "Contrôle négatif extraction filtre eau"
        })

    else:
        rows.append({
            "sample-id":   sid,
            "sample_type": "unknown",
            "depth_cm":    "NA",
            "condition":   "NA",
            "location":    "PA",
            "description": f"Type inconnu : {sid}"
        })

# Écriture du fichier métadonnées (format QIIME2 avec #q2:types)
with open(metadata_path, 'w') as fh:
    fh.write('\t'.join(headers) + '\n')
    fh.write('#q2:types\t' + '\t'.join(types) + '\n')
    for row in rows:
        fh.write('\t'.join([str(row[h]) for h in headers]) + '\n')

print(f"  Métadonnées écrites : {metadata_path}")
print(f"  Total : {len(rows)} échantillons")

# Résumé
from collections import Counter
counts = Counter(r['sample_type'] for r in rows)
for k, v in counts.items():
    print(f"    {k}: {v}")
PYEOF

log "Manifest et métadonnées générés"
head -5 "${DBDIR}/manifest"

# ============================================================
# ÉTAPE 01 — IMPORT QIIME2 (données déjà nettoyées)
# ============================================================
log "Import QIIME2 - données pré-nettoyées (format PairedEnd Manifest)"

cd "${QDIR}/core"

if [[ ! -f "demux_paired.qza" ]]; then
    conda run -n "$QIIME2_ENV" qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "${DBDIR}/manifest" \
        --output-path "demux_paired.qza" \
        --input-format PairedEndFastqManifestPhred33V2 || {
            log "ERREUR import QIIME2"
            echo "Vérifiez le manifest :"
            cat "${DBDIR}/manifest"
            exit 1
        }
    log "Import QIIME2 réussi"
else
    log "demux_paired.qza déjà présent — import ignoré"
fi

# Visualisation du fichier importé
conda run -n "$QIIME2_ENV" qiime demux summarize \
    --i-data "demux_paired.qza" \
    --o-visualization "../visual/demux-summary.qzv"

log "→ Visualisez ../visual/demux-summary.qzv pour vérifier la qualité des reads"

# ============================================================
# ÉTAPE 02 — CUTADAPT : suppression des amorces V4-V5
# ============================================================
# Les amorces PCR doivent être retirées AVANT DADA2
# même si trimmomatic a déjà retiré les adaptateurs Illumina.
# ============================================================
log "Cutadapt — suppression amorces V4-V5 (${PRIMER_F} / ${PRIMER_R})"

cd "${QDIR}/core"

if [[ ! -f "demux_trimmed.qza" ]]; then
    conda run -n "$QIIME2_ENV" qiime cutadapt trim-paired \
        --i-demultiplexed-sequences "demux_paired.qza" \
        --p-front-f "$PRIMER_F" \
        --p-front-r "$PRIMER_R" \
        --p-discard-untrimmed \
        --p-no-indels \
        --p-overlap 10 \
        --p-cores "$NTHREADS" \
        --o-trimmed-sequences "demux_trimmed.qza" \
        --verbose 2> "${ROOTDIR}/tmp/cutadapt.log" || {
            log "ERREUR cutadapt"
            tail -20 "${ROOTDIR}/tmp/cutadapt.log"
            exit 1
        }
    log "Cutadapt terminé — log : ${ROOTDIR}/tmp/cutadapt.log"
else
    log "demux_trimmed.qza déjà présent — cutadapt ignoré"
fi

# Résumé post-cutadapt
conda run -n "$QIIME2_ENV" qiime demux summarize \
    --i-data "demux_trimmed.qza" \
    --o-visualization "../visual/demux-trimmed-summary.qzv"

log "→ Visualisez ../visual/demux-trimmed-summary.qzv pour vérifier post-cutadapt"

# ============================================================
# ÉTAPE 03 — DADA2 : débruitage et génération ASV
# ============================================================
# Paramètres de troncature :
#   --p-trunc-len-f et --p-trunc-len-r à 0 = pas de troncature positionnelle.
#   Ajustez ces valeurs selon le profil qualité de demux-trimmed-summary.qzv
#   (règle générale : tronquer quand Q < 25-30 de façon continue)
# ============================================================
log "DADA2 — débruitage paired-end"

TRUNC_F=0   # ← modifier si besoin après inspection demux-trimmed-summary.qzv
TRUNC_R=0   # ← modifier si besoin après inspection demux-trimmed-summary.qzv

cd "${QDIR}/core"

if [[ ! -f "table.qza" ]]; then
    conda run -n "$QIIME2_ENV" qiime dada2 denoise-paired \
        --i-demultiplexed-seqs "demux_trimmed.qza" \
        --p-trunc-len-f "$TRUNC_F" \
        --p-trunc-len-r "$TRUNC_R" \
        --p-n-threads "$NTHREADS" \
        --o-table "table.qza" \
        --o-representative-sequences "rep-seqs.qza" \
        --o-denoising-stats "denoising-stats.qza" || {
            log "ERREUR DADA2"
            exit 1
        }
    log "DADA2 terminé"
else
    log "table.qza déjà présent — DADA2 ignoré"
fi

# Visualisation des statistiques de débruitage
conda run -n "$QIIME2_ENV" qiime metadata tabulate \
    --m-input-file "denoising-stats.qza" \
    --o-visualization "../visual/denoising-stats.qzv"

# Visualisation des séquences représentatives
conda run -n "$QIIME2_ENV" qiime feature-table tabulate-seqs \
    --i-data "rep-seqs.qza" \
    --o-visualization "../visual/rep-seqs.qzv"

log "→ Visualisez ../visual/denoising-stats.qzv pour les stats DADA2"

# ============================================================
# ÉTAPE 04 — RÉSUMÉ DE LA TABLE ASV
# ============================================================
log "Résumé feature-table (fréquences par échantillon)"

cd "${QDIR}/core"

conda run -n "$QIIME2_ENV" qiime feature-table summarize \
    --i-table "table.qza" \
    --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
    --o-visualization "../visual/table-summary.qzv"

log "→ Visualisez ../visual/table-summary.qzv pour voir les effectifs par échantillon"

# ============================================================
# ÉTAPE 05 — COURBES DE RARÉFACTION ALPHA
# ============================================================
# Ces courbes permettent de choisir la profondeur de raréfaction.
# max_depth : fixé à 50000 ; ajuster si vos échantillons ont plus de reads.
# ============================================================
log "Calcul des courbes de raréfaction alpha (patience…)"

MAX_DEPTH=50000    # ← augmenter si vos échantillons sont plus profonds

cd "${QDIR}/core"

conda run -n "$QIIME2_ENV" qiime diversity alpha-rarefaction \
    --i-table "table.qza" \
    --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
    --p-max-depth "$MAX_DEPTH" \
    --p-steps 20 \
    --p-metrics "shannon" \
    --p-metrics "observed_features" \
    --p-metrics "faith_pd" \
    --o-visualization "../visual/alpha-rarefaction-curves.qzv" || {
        log "WARN: faith_pd nécessite un arbre — nouvelle tentative sans faith_pd"
        conda run -n "$QIIME2_ENV" qiime diversity alpha-rarefaction \
            --i-table "table.qza" \
            --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
            --p-max-depth "$MAX_DEPTH" \
            --p-steps 20 \
            --p-metrics "shannon" \
            --p-metrics "observed_features" \
            --o-visualization "../visual/alpha-rarefaction-curves.qzv"
    }

log "Courbes de raréfaction générées"

# ============================================================
# RÉSUMÉ FINAL — SCRIPT 1 TERMINÉ
# ============================================================
log "SCRIPT 1 TERMINÉ"

echo ""
echo "======================================================================="
echo "  SCRIPT 1 TERMINÉ — ACTION REQUISE AVANT LANCER LE SCRIPT 2"
echo "======================================================================="
echo ""
echo "  Visualisations à consulter dans : ${QDIR}/visual/"
echo ""
echo "  1. table-summary.qzv         → effectifs reads par échantillon"
echo "  2. alpha-rarefaction-curves.qzv → courbes de raréfaction Shannon/Observed"
echo "  3. denoising-stats.qzv       → stats DADA2 (reads filtrés/mergés/chimères)"
echo "  4. demux-trimmed-summary.qzv → qualité après cutadapt"
echo ""
echo "  → Ouvrir avec : https://view.qiime2.org"
echo ""
echo "  Choisissez une profondeur de raréfaction en observant :"
echo "   - Le plateau des courbes alpha (Shannon / observed_features)"
echo "   - Le compromis nombre de reads / nombre d'échantillons conservés"
echo "   - La valeur minimale pour l'échantillon le moins séquencé que"
echo "     vous souhaitez conserver (exclure les contrôles négatifs)"
echo ""
echo "  → Editez ensuite RAREFACTION_DEPTH dans 02_diversity_analysis.sh"
echo "     puis lancez : bash 02_diversity_analysis.sh"
echo "======================================================================="
