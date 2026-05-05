#!/usr/bin/env bash
# =============================================================================
# SCRIPT 2/2 — Classification SILVA, Diversité & Exports
# Projet  : vague_project | 16S V4-V5 | Pointe de l'Artillerie
# Usage   : bash 02_diversity_analysis.sh
#           ou : bash 02_diversity_analysis.sh --depth 10000
# Prérequis : avoir terminé le Script 1 et choisi la profondeur de raréfaction
# =============================================================================

set -euo pipefail

# ============================================================
# ► PARAMÈTRE PRINCIPAL — À MODIFIER PAR L'UTILISATEUR ◄
# ============================================================
RAREFACTION_DEPTH=10000   # ← REMPLACER par la valeur choisie sur les courbes

# Lecture en argument optionnel : bash 02_diversity_analysis.sh --depth XXXX
while [[ $# -gt 0 ]]; do
    case "$1" in
        --depth) RAREFACTION_DEPTH="$2"; shift 2 ;;
        *) echo "Argument inconnu : $1"; exit 1 ;;
    esac
done

# ============================================================
# VARIABLES GLOBALES
# ============================================================
export ROOTDIR="/nvme/bio/data_fungi/vague_project"
export NTHREADS=16
export QIIME2_ENV="qiime2-amplicon-2025.7"
export TMPDIR="${ROOTDIR}/tmp"

export DBDIR="${ROOTDIR}/98_databasefiles"
export QDIR="${ROOTDIR}/05_QIIME2"

# Classifier SILVA 138.2 — région V4-V5 (515F/926R)
# S'il n'existe pas, il sera créé automatiquement (long ~4-6h)
CLASSIFIER="${DBDIR}/silva-138.2-ssu-nr99-515f-926r-classifier.qza"

PRIMER_F="GTGYCAGCMGCCGCGGTAA"
PRIMER_R="CCGYCAATTYMTTTRAGTTT"

# ============================================================
# INITIALISATION
# ============================================================
log() { echo -e "\n[$(date +'%F %T')] ===  $* ===\n"; }

mkdir -p "$TMPDIR" \
    "${QDIR}/core/diversity" \
    "${QDIR}/core/pcoa" \
    "${QDIR}/subtables" \
    "${QDIR}/visual" \
    "${QDIR}/export"

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
set -u

log "SCRIPT 2 DÉMARRÉ — Profondeur raréfaction : ${RAREFACTION_DEPTH} reads"

# Vérification prérequis
for f in "${QDIR}/core/table.qza" "${QDIR}/core/rep-seqs.qza"; do
    if [[ ! -f "$f" ]]; then
        echo "[ERREUR] Fichier manquant : $f"
        echo "  → Lancez d'abord 01_qiime2_preprocess.sh"
        exit 1
    fi
done

cd "${QDIR}/core"

# ============================================================
# ÉTAPE 06 — ARBRE PHYLOGÉNÉTIQUE (MAFFT + FastTree)
# ============================================================
log "Génération arbre phylogénétique (MAFFT + FastTree)"

if [[ ! -f "tree.qza" ]]; then
    conda run -n "$QIIME2_ENV" qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences "rep-seqs.qza" \
        --p-n-threads "$NTHREADS" \
        --o-alignment "aligned-rep-seqs.qza" \
        --o-masked-alignment "masked-aligned-rep-seqs.qza" \
        --o-tree "unrooted-tree.qza" \
        --o-rooted-tree "tree.qza" || {
            log "ERREUR génération arbre phylogénétique"
            exit 1
        }
    log "Arbre phylogénétique généré"
else
    log "tree.qza déjà présent — arbre ignoré"
fi

# ============================================================
# ÉTAPE 07 — CLASSIFIER SILVA 138.2 (V4-V5)
# ============================================================
log "Vérification / création classifier SILVA 138.2 (V4-V5)"

if ! conda run -n "$QIIME2_ENV" qiime tools validate "$CLASSIFIER" 2>/dev/null; then
    log "Classifier absent — création en cours (peut prendre 4-6h)…"

    cd "$DBDIR"

    # Téléchargement SILVA 138.2
    conda run -n "$QIIME2_ENV" qiime rescript get-silva-data \
        --p-version '138.2' \
        --p-target 'SSURef_NR99' \
        --o-silva-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
        --o-silva-taxonomy silva-138.2-ssu-nr99-tax.qza

    # Conversion RNA → DNA
    conda run -n "$QIIME2_ENV" qiime rescript reverse-transcribe \
        --i-rna-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
        --o-dna-sequences silva-138.2-ssu-nr99-seqs.qza

    # Nettoyage séquences ambiguës
    conda run -n "$QIIME2_ENV" qiime rescript cull-seqs \
        --i-sequences silva-138.2-ssu-nr99-seqs.qza \
        --o-clean-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza

    # Filtre par longueur et taxonomie
    conda run -n "$QIIME2_ENV" qiime rescript filter-seqs-length-by-taxon \
        --i-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza \
        --i-taxonomy silva-138.2-ssu-nr99-tax.qza \
        --p-labels Archaea Bacteria Eukaryota \
        --p-min-lens 900 1200 1400 \
        --o-filtered-seqs silva-138.2-ssu-nr99-seqs-filt.qza \
        --o-discarded-seqs silva-138.2-ssu-nr99-seqs-discard.qza

    # Déréplication uniq
    conda run -n "$QIIME2_ENV" qiime rescript dereplicate \
        --i-sequences silva-138.2-ssu-nr99-seqs-filt.qza \
        --i-taxa silva-138.2-ssu-nr99-tax.qza \
        --p-mode 'uniq' \
        --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
        --o-dereplicated-taxa silva-138.2-ssu-nr99-tax-derep-uniq.qza

    # Extraction région V4-V5
    conda run -n "$QIIME2_ENV" qiime feature-classifier extract-reads \
        --i-sequences silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
        --p-f-primer "$PRIMER_F" \
        --p-r-primer "$PRIMER_R" \
        --p-n-jobs "$NTHREADS" \
        --p-read-orientation 'forward' \
        --o-reads silva-138.2-ssu-nr99-seqs-515f-926r.qza

    # Déréplication région extraite
    conda run -n "$QIIME2_ENV" qiime rescript dereplicate \
        --i-sequences silva-138.2-ssu-nr99-seqs-515f-926r.qza \
        --i-taxa silva-138.2-ssu-nr99-tax-derep-uniq.qza \
        --p-mode uniq \
        --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-515f-926r-uniq.qza \
        --o-dereplicated-taxa silva-138.2-ssu-nr99-tax-515f-926r-derep-uniq.qza

    # Entraînement Naive Bayes
    conda run -n "$QIIME2_ENV" qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads silva-138.2-ssu-nr99-seqs-515f-926r-uniq.qza \
        --i-reference-taxonomy silva-138.2-ssu-nr99-tax-515f-926r-derep-uniq.qza \
        --o-classifier "$CLASSIFIER"

    log "Classifier SILVA 138.2 V4-V5 créé : $CLASSIFIER"
else
    log "Classifier SILVA 138.2 V4-V5 déjà présent — réutilisé"
fi

# ============================================================
# ÉTAPE 08 — CLASSIFICATION TAXONOMIQUE
# ============================================================
log "Classification taxonomique SILVA 138.2"

cd "${QDIR}/core"

if [[ ! -f "taxonomy.qza" ]]; then
    conda run -n "$QIIME2_ENV" qiime feature-classifier classify-sklearn \
        --i-classifier "$CLASSIFIER" \
        --i-reads "rep-seqs.qza" \
        --p-n-jobs "$NTHREADS" \
        --o-classification "taxonomy.qza" || {
            log "ERREUR classification taxonomique"
            exit 1
        }
    log "Classification terminée"
else
    log "taxonomy.qza déjà présent — classification ignorée"
fi

# Visualisation taxonomie
conda run -n "$QIIME2_ENV" qiime metadata tabulate \
    --m-input-file "taxonomy.qza" \
    --o-visualization "../visual/taxonomy.qzv"

# Barplot taxa
conda run -n "$QIIME2_ENV" qiime taxa barplot \
    --i-table "table.qza" \
    --i-taxonomy "taxonomy.qza" \
    --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
    --o-visualization "../visual/taxa-bar-plots.qzv"

log "→ Visualisez ../visual/taxa-bar-plots.qzv"

# ============================================================
# ÉTAPE 09 — RARÉFACTION DE LA TABLE
# ============================================================
log "Raréfaction de la table à ${RAREFACTION_DEPTH} reads"

conda run -n "$QIIME2_ENV" qiime feature-table rarefy \
    --i-table "table.qza" \
    --p-sampling-depth "$RAREFACTION_DEPTH" \
    --o-rarefied-table "../subtables/RarTable-depth${RAREFACTION_DEPTH}.qza" || {
        log "ERREUR raréfaction"
        exit 1
    }

# Résumé table raréfiée
conda run -n "$QIIME2_ENV" qiime feature-table summarize \
    --i-table "../subtables/RarTable-depth${RAREFACTION_DEPTH}.qza" \
    --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
    --o-visualization "../visual/table-rarefied-summary.qzv"

log "Table raréfiée générée : subtables/RarTable-depth${RAREFACTION_DEPTH}.qza"

# ============================================================
# ÉTAPE 10 — MÉTRIQUES DE DIVERSITÉ (phylogénétiques)
# ============================================================
log "Calcul des métriques de diversité (core-metrics-phylogenetic)"

cd "${QDIR}/core"

# Nettoyage résultats précédents si existants
rm -f \
    diversity/Vector-faith_pd.qza \
    diversity/Vector-observed_asv.qza \
    diversity/Vector-shannon.qza \
    diversity/Vector-evenness.qza \
    diversity/Matrix-unweighted_unifrac.qza \
    diversity/Matrix-weighted_unifrac.qza \
    diversity/Matrix-jaccard.qza \
    diversity/Matrix-braycurtis.qza \
    pcoa/PCoA-unweighted_unifrac.qza \
    pcoa/PCoA-weighted_unifrac.qza \
    pcoa/PCoA-jaccard.qza \
    pcoa/PCoA-braycurtis.qza \
    ../visual/Emperor-*.qzv 2>/dev/null || true

conda run -n "$QIIME2_ENV" qiime diversity core-metrics-phylogenetic \
    --i-table "table.qza" \
    --i-phylogeny "tree.qza" \
    --p-sampling-depth "$RAREFACTION_DEPTH" \
    --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
    --o-rarefied-table "rarefied_table.qza" \
    --o-faith-pd-vector "diversity/Vector-faith_pd.qza" \
    --o-observed-features-vector "diversity/Vector-observed_asv.qza" \
    --o-shannon-vector "diversity/Vector-shannon.qza" \
    --o-evenness-vector "diversity/Vector-evenness.qza" \
    --o-unweighted-unifrac-distance-matrix "diversity/Matrix-unweighted_unifrac.qza" \
    --o-weighted-unifrac-distance-matrix "diversity/Matrix-weighted_unifrac.qza" \
    --o-jaccard-distance-matrix "diversity/Matrix-jaccard.qza" \
    --o-bray-curtis-distance-matrix "diversity/Matrix-braycurtis.qza" \
    --o-unweighted-unifrac-pcoa-results "pcoa/PCoA-unweighted_unifrac.qza" \
    --o-weighted-unifrac-pcoa-results "pcoa/PCoA-weighted_unifrac.qza" \
    --o-jaccard-pcoa-results "pcoa/PCoA-jaccard.qza" \
    --o-bray-curtis-pcoa-results "pcoa/PCoA-braycurtis.qza" \
    --o-unweighted-unifrac-emperor "../visual/Emperor-unweighted_unifrac.qzv" \
    --o-weighted-unifrac-emperor "../visual/Emperor-weighted_unifrac.qzv" \
    --o-jaccard-emperor "../visual/Emperor-jaccard.qzv" \
    --o-bray-curtis-emperor "../visual/Emperor-braycurtis.qzv" || {

    log "WARN: core-metrics-phylogenetic échoué — tentative sans phylogénie"
    conda run -n "$QIIME2_ENV" qiime diversity core-metrics \
        --i-table "table.qza" \
        --p-sampling-depth "$RAREFACTION_DEPTH" \
        --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
        --o-rarefied-table "rarefied_table.qza" \
        --o-observed-features-vector "diversity/Vector-observed_asv.qza" \
        --o-shannon-vector "diversity/Vector-shannon.qza" \
        --o-evenness-vector "diversity/Vector-evenness.qza" \
        --o-jaccard-distance-matrix "diversity/Matrix-jaccard.qza" \
        --o-bray-curtis-distance-matrix "diversity/Matrix-braycurtis.qza" \
        --o-jaccard-pcoa-results "pcoa/PCoA-jaccard.qza" \
        --o-bray-curtis-pcoa-results "pcoa/PCoA-braycurtis.qza" \
        --o-jaccard-emperor "../visual/Emperor-jaccard.qzv" \
        --o-bray-curtis-emperor "../visual/Emperor-braycurtis.qzv"
}

log "Métriques de diversité calculées"

# ---- Diversité alpha — visualisations ----
log "Visualisations diversité alpha"

for metric in shannon evenness observed_features faith_pd; do
    vec_file="diversity/Vector-${metric/observed_features/observed_asv}.qza"
    [[ "$metric" == "observed_features" ]] && vec_file="diversity/Vector-observed_asv.qza"
    [[ "$metric" == "faith_pd" ]] && vec_file="diversity/Vector-faith_pd.qza"
    if [[ -f "$vec_file" ]]; then
        conda run -n "$QIIME2_ENV" qiime diversity alpha-group-significance \
            --i-alpha-diversity "$vec_file" \
            --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
            --o-visualization "../visual/alpha-${metric}-significance.qzv" 2>/dev/null || true
    fi
done

# ---- Diversité beta — PERMANOVA sur sample_type ----
log "Tests PERMANOVA beta diversité (sample_type)"

for matrix in braycurtis jaccard; do
    mat_file="diversity/Matrix-${matrix}.qza"
    if [[ -f "$mat_file" ]]; then
        conda run -n "$QIIME2_ENV" qiime diversity beta-group-significance \
            --i-distance-matrix "$mat_file" \
            --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
            --m-metadata-column "sample_type" \
            --p-method "permanova" \
            --p-permutations 999 \
            --o-visualization "../visual/beta-${matrix}-permanova-sampletype.qzv" 2>/dev/null || true
    fi
done

# ---- Core features ----
log "Core features analysis"

conda run -n "$QIIME2_ENV" qiime feature-table core-features \
    --i-table "../subtables/RarTable-depth${RAREFACTION_DEPTH}.qza" \
    --p-min-fraction 0.1 \
    --p-max-fraction 1.0 \
    --p-steps 10 \
    --o-visualization "../visual/CoreBiom-all.qzv" 2>/dev/null || true

# ============================================================
# ÉTAPE 11 — EXPORTS
# ============================================================
log "Export des fichiers QIIME2 vers formats tabulaires"

mkdir -p \
    "${QDIR}/export/core/table" \
    "${QDIR}/export/core/rep-seqs" \
    "${QDIR}/export/core/taxonomy" \
    "${QDIR}/export/subtables/RarTable" \
    "${QDIR}/export/diversity_tsv"

cd "${QDIR}"

# Table ASV principale (BIOM)
conda run -n "$QIIME2_ENV" qiime tools export \
    --input-path "core/table.qza" \
    --output-path "export/core/table"

# Séquences représentatives
conda run -n "$QIIME2_ENV" qiime tools export \
    --input-path "core/rep-seqs.qza" \
    --output-path "export/core/rep-seqs"

# Taxonomie
conda run -n "$QIIME2_ENV" qiime tools export \
    --input-path "core/taxonomy.qza" \
    --output-path "export/core/taxonomy"

# Table raréfiée
conda run -n "$QIIME2_ENV" qiime tools export \
    --input-path "subtables/RarTable-depth${RAREFACTION_DEPTH}.qza" \
    --output-path "export/subtables/RarTable"

# Stats DADA2
conda run -n "$QIIME2_ENV" qiime tools export \
    --input-path "core/denoising-stats.qza" \
    --output-path "export/diversity_tsv" || true

# Métriques de diversité
export_tsv() {
    local qza="$1"
    local name="$2"
    if [[ -f "$qza" ]]; then
        local tmp="${QDIR}/export/diversity_tsv/${name}_tmp"
        rm -rf "$tmp"
        conda run -n "$QIIME2_ENV" qiime tools export \
            --input-path "$qza" \
            --output-path "$tmp" && {
            find "$tmp" -name "*.tsv" -o -name "*.txt" | while read -r f; do
                cp "$f" "${QDIR}/export/diversity_tsv/${name}_$(basename $f)"
            done
            rm -rf "$tmp"
        } || log "WARN: export échoué pour $name"
    fi
}

export_tsv "core/diversity/Vector-faith_pd.qza"         "faith_pd"
export_tsv "core/diversity/Vector-shannon.qza"          "shannon"
export_tsv "core/diversity/Vector-observed_asv.qza"     "observed_features"
export_tsv "core/diversity/Vector-evenness.qza"         "evenness"
export_tsv "core/diversity/Matrix-braycurtis.qza"       "bray_curtis"
export_tsv "core/diversity/Matrix-jaccard.qza"          "jaccard"
export_tsv "core/diversity/Matrix-unweighted_unifrac.qza" "unweighted_unifrac"
export_tsv "core/diversity/Matrix-weighted_unifrac.qza"   "weighted_unifrac"
export_tsv "core/pcoa/PCoA-braycurtis.qza"              "pcoa_braycurtis"
export_tsv "core/pcoa/PCoA-unweighted_unifrac.qza"      "pcoa_unweighted_unifrac"

# ---- Conversion BIOM → TSV + fusion taxonomie ----
log "Conversion BIOM table raréfiée + fusion taxonomie SILVA"

BIOM_FILE="${QDIR}/export/subtables/RarTable/feature-table.biom"
TAX_FILE="${QDIR}/export/core/taxonomy/taxonomy.tsv"

if [[ -f "$BIOM_FILE" ]]; then
    conda run -n "$QIIME2_ENV" biom convert \
        -i "$BIOM_FILE" \
        -o "${QDIR}/export/subtables/RarTable/table-from-biom.tsv" \
        --to-tsv

    # Nettoyer header BIOM
    sed '1d; s/#OTU ID/ASV_ID/' \
        "${QDIR}/export/subtables/RarTable/table-from-biom.tsv" \
        > "${QDIR}/export/subtables/RarTable/ASV.tsv"

    log "ASV.tsv créé ($(wc -l < "${QDIR}/export/subtables/RarTable/ASV.tsv") lignes)"

    # Fusion avec taxonomie SILVA → ASV_taxonomy.tsv
    if [[ -f "$TAX_FILE" ]]; then
        python3 << PYEOF2
import csv, re, os

asv_path = "${QDIR}/export/subtables/RarTable/ASV.tsv"
tax_path = "${QDIR}/export/core/taxonomy/taxonomy.tsv"
out_path = "${QDIR}/export/subtables/RarTable/ASV_taxonomy.tsv"

# Charger taxonomie
taxonomy = {}
with open(tax_path) as fh:
    reader = csv.reader(fh, delimiter='\t')
    next(reader)  # skip header
    for row in reader:
        if len(row) >= 2:
            taxonomy[row[0]] = row[1]

def parse_silva(tax_str):
    levels = ["Kingdom","Phylum","Class","Order","Family","Genus","Species"]
    patterns = [r'D_0__([^;]+)',r'D_1__([^;]+)',r'D_2__([^;]+)',
                r'D_3__([^;]+)',r'D_4__([^;]+)',r'D_5__([^;]+)',r'D_6__([^;]+)']
    result = []
    for p in patterns:
        m = re.search(p, tax_str)
        result.append(m.group(1).strip() if m else "Unassigned")
    return result

# Lire ASV table
with open(asv_path) as fi, open(out_path, 'w') as fo:
    reader = csv.reader(fi, delimiter='\t')
    header = next(reader)
    tax_levels = ["Kingdom","Phylum","Class","Order","Family","Genus","Species"]
    fo.write('\t'.join(["ASV_ID"] + tax_levels + header[1:]) + '\n')
    for row in reader:
        asv_id = row[0]
        tax_string = taxonomy.get(asv_id, "")
        tax_parsed = parse_silva(tax_string)
        fo.write('\t'.join([asv_id] + tax_parsed + row[1:]) + '\n')

print(f"ASV_taxonomy.tsv créé : {out_path}")
PYEOF2
    fi
fi

# ============================================================
# RAPPORT FINAL
# ============================================================
log "SCRIPT 2 TERMINÉ"

echo ""
echo "======================================================================="
echo "  SCRIPT 2 TERMINÉ — ANALYSES TERMINÉES"
echo "  Profondeur raréfaction utilisée : ${RAREFACTION_DEPTH} reads"
echo "======================================================================="
echo ""
echo "  VISUALISATIONS (https://view.qiime2.org) :"
echo "    ${QDIR}/visual/"
echo "    ├── taxa-bar-plots.qzv            → composition taxonomique"
echo "    ├── alpha-rarefaction-curves.qzv  → courbes raréfaction"
echo "    ├── alpha-shannon-significance.qzv → test alpha diversité"
echo "    ├── Emperor-*.qzv                 → PCoA interactif"
echo "    ├── beta-*-permanova-*.qzv        → PERMANOVA"
echo "    └── CoreBiom-all.qzv              → core features"
echo ""
echo "  FICHIERS TABULAIRES :"
echo "    ${QDIR}/export/"
echo "    ├── subtables/RarTable/ASV_taxonomy.tsv   → table ASV + taxonomie SILVA"
echo "    ├── core/taxonomy/taxonomy.tsv             → classifications"
echo "    ├── core/rep-seqs/dna-sequences.fasta      → séquences ASV"
echo "    └── diversity_tsv/                         → métriques alpha/beta"
echo ""
echo "======================================================================="
