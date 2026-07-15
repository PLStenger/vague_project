#!/usr/bin/env bash
# =============================================================================
# SCRIPT 2/2 — Taxonomie + Diversité | PA + KNS
# Double classification 16S (SILVA 138.2) + 18S (PR2 v5) pour eucaryotes
# Usage : bash 02_diversity_analysis.sh [--depth 42105]
# =============================================================================

set -euo pipefail

RAREFACTION_DEPTH=42105
while [[ $# -gt 0 ]]; do
    case "$1" in
        --depth) RAREFACTION_DEPTH="$2"; shift 2 ;;
        *) echo "Argument inconnu: $1"; exit 1 ;;
    esac
done

export ROOTDIR="/nvme/bio/data_fungi/vague_project"
export NTHREADS=16
export QIIME2_ENV="qiime2-amplicon-2025.7"
export TMPDIR="${ROOTDIR}/tmp"
export DBDIR="${ROOTDIR}/98_databasefiles"
export QDIR="${ROOTDIR}/05_QIIME2"
export PRIMER_F="GTGYCAGCMGCCGCGGTAA"
export PRIMER_R="CCGYCAATTYMTTTRAGTTT"

# Chemins PR2 (déjà générés manuellement)
PR2_CLASSIFIER="${DBDIR}/pr2-v5-classifier-515F-926R.qza"
PR2_SEQS_RAW="${DBDIR}/pr2_v5_SSU_seqs.qza"
PR2_TAX_RAW="${DBDIR}/pr2_v5_SSU_tax.qza"
PR2_SEQS_TRIMMED="${DBDIR}/pr2_v5_515F_926R_seqs.qza"

RAREFACTION_DEPTH_EUK=1000   # À ajuster après inspection de table-euk-candidates-summary.qzv

log() { echo -e "\n[$(date +'%F %T')] === $* ===\n"; }

mkdir -p "$TMPDIR" "${QDIR}/core/diversity" "${QDIR}/core/pcoa" \
         "${QDIR}/visual" "${QDIR}/subtables" "${QDIR}/export"

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
set -u

for f in "${QDIR}/core/table-decontam.qza" \
          "${QDIR}/core/rep-seqs-decontam.qza" \
          "${DBDIR}/sample-metadata.tsv" \
          "$PR2_CLASSIFIER"; do
    [[ -f "$f" ]] || { echo "Fichier manquant: $f"; exit 1; }
done

cd "${QDIR}/core"

# =============================================================================
# ÉTAPES DÉJÀ EFFECTUÉES — commentées
# =============================================================================

# log "Préparation listes contrôles / non-contrôles"
# awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} NR==2 && $1=="#q2:types"{next} $c["is_control"]=="yes"{print $1}' \
#     "${DBDIR}/sample-metadata.tsv" > "${DBDIR}/control-samples.txt"
# awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} NR==2 && $1=="#q2:types"{next} $c["is_control"]!="yes"{print $1}' \
#     "${DBDIR}/sample-metadata.tsv" > "${DBDIR}/non-control-samples.txt"

# log "Extraction table des contrôles"
# conda run -n "$QIIME2_ENV" qiime feature-table filter-samples \
#     --i-table table.qza \
#     --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
#     --p-where "[is_control]='yes'" \
#     --o-filtered-table table-controls-only.qza
# conda run -n "$QIIME2_ENV" qiime feature-table summarize \
#     --i-table table-controls-only.qza \
#     --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
#     --o-visualization ../visual/table-controls-only-summary.qzv

# log "Extraction liste des ASV présents dans les contrôles (359 ASV)"
# mkdir -p "${TMPDIR}/controls_export"
# rm -rf "${TMPDIR}/controls_export"/* 2>/dev/null || true
# conda run -n "$QIIME2_ENV" qiime tools export \
#     --input-path table-controls-only.qza \
#     --output-path "${TMPDIR}/controls_export"
# conda run -n "$QIIME2_ENV" biom convert \
#     -i "${TMPDIR}/controls_export/feature-table.biom" \
#     -o "${TMPDIR}/controls_export/controls_table.tsv" \
#     --to-tsv
# { echo -e "feature-id"; \
#   awk 'BEGIN{FS="\t"} NR>2 {sum=0; for(i=2;i<=NF;i++) sum+=$i; if(sum>0) print $1}' \
#   "${TMPDIR}/controls_export/controls_table.tsv"; \
# } > "${DBDIR}/features-in-controls.txt"

# log "Extraction table échantillons biologiques uniquement"
# conda run -n "$QIIME2_ENV" qiime feature-table filter-samples \
#     --i-table table.qza \
#     --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
#     --p-where "[is_control]!='yes'" \
#     --o-filtered-table table-non-controls.qza
# conda run -n "$QIIME2_ENV" qiime feature-table summarize \
#     --i-table table-non-controls.qza \
#     --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
#     --o-visualization ../visual/table-non-controls-summary.qzv

# log "Retrait des ASV présents dans les contrôles"
# conda run -n "$QIIME2_ENV" qiime feature-table filter-features \
#     --i-table table-non-controls.qza \
#     --m-metadata-file "${DBDIR}/features-in-controls.txt" \
#     --p-exclude-ids \
#     --o-filtered-table table-decontam.qza
# conda run -n "$QIIME2_ENV" qiime feature-table summarize \
#     --i-table table-decontam.qza \
#     --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
#     --o-visualization ../visual/table-decontam-summary.qzv

# log "Filtrage rep-seqs sur ASV décontaminés"
# conda run -n "$QIIME2_ENV" qiime feature-table filter-seqs \
#     --i-data rep-seqs.qza \
#     --i-table table-decontam.qza \
#     --o-filtered-data rep-seqs-decontam.qza

# =============================================================================
# ÉTAPE 07a — CLASSIFICATION SILVA 138.2 (16S — bactéries + archées)
# À LANCER si taxonomy-silva-raw.qza absent
# =============================================================================
log "Classification taxonomique SILVA 138.2 (16S)"

CLASSIFIER_SILVA="${DBDIR}/silva-138.2-ssu-nr99-515f-926r-classifier.qza"
CLASSIFIER_SOURCE="/nvme/bio/data_fungi/valormicro_nc/98_databasefiles/silva-138.2-ssu-nr99-515f-926r-classifier.qza"

if [[ ! -f "$CLASSIFIER_SILVA" ]]; then
    log "Copie du classifier SILVA depuis valormicro_nc"
    cp "$CLASSIFIER_SOURCE" "$CLASSIFIER_SILVA" || \
        { log "ERREUR: Classifier SILVA source introuvable"; exit 1; }
fi

conda run -n "$QIIME2_ENV" qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER_SILVA" \
    --i-reads rep-seqs-decontam.qza \
    --p-n-jobs "$NTHREADS" \
    --o-classification taxonomy-silva-raw.qza

conda run -n "$QIIME2_ENV" qiime metadata tabulate \
    --m-input-file taxonomy-silva-raw.qza \
    --o-visualization ../visual/taxonomy-silva-raw.qzv

# =============================================================================
# ÉTAPE 07b — SÉPARATION procaryotes / candidats eucaryotes post-SILVA
# =============================================================================
log "Séparation ASVs 16S (procaryotes) / candidats 18S (eucaryotes)"

# Table + séquences procaryotes strictes (exclut Eukaryota, Mito, Chloro, Unassigned)
conda run -n "$QIIME2_ENV" qiime taxa filter-table \
    --i-table table-decontam.qza \
    --i-taxonomy taxonomy-silva-raw.qza \
    --p-exclude "Eukaryota,Mitochondria,Chloroplast,Unassigned" \
    --o-filtered-table table-prokaryotes.qza

conda run -n "$QIIME2_ENV" qiime taxa filter-seqs \
    --i-sequences rep-seqs-decontam.qza \
    --i-taxonomy taxonomy-silva-raw.qza \
    --p-exclude "Eukaryota,Mitochondria,Chloroplast,Unassigned" \
    --o-filtered-sequences rep-seqs-prokaryotes.qza

conda run -n "$QIIME2_ENV" qiime feature-table summarize \
    --i-table table-prokaryotes.qza \
    --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
    --o-visualization ../visual/table-prokaryotes-summary.qzv

# Table + séquences candidats eucaryotes (Eukaryota SILVA + Unassigned)
conda run -n "$QIIME2_ENV" qiime taxa filter-table \
    --i-table table-decontam.qza \
    --i-taxonomy taxonomy-silva-raw.qza \
    --p-include "Eukaryota,Unassigned" \
    --o-filtered-table table-euk-candidates.qza

conda run -n "$QIIME2_ENV" qiime taxa filter-seqs \
    --i-sequences rep-seqs-decontam.qza \
    --i-taxonomy taxonomy-silva-raw.qza \
    --p-include "Eukaryota,Unassigned" \
    --o-filtered-sequences rep-seqs-euk-candidates.qza

conda run -n "$QIIME2_ENV" qiime feature-table summarize \
    --i-table table-euk-candidates.qza \
    --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
    --o-visualization ../visual/table-euk-candidates-summary.qzv

# =============================================================================
# ÉTAPE 07c — Préparation PR2 DÉJÀ EFFECTUÉE (manuellement)
# Les fichiers suivants existent déjà :
#   ${DBDIR}/pr2_v5_SSU_seqs.qza          (import mothur fasta)
#   ${DBDIR}/pr2_v5_SSU_tax.qza           (import taxonomie)
#   ${DBDIR}/pr2_v5_515F_926R_seqs.qza    (extract-reads V4-V5)
#   ${DBDIR}/pr2-v5-classifier-515F-926R.qza  (fit-classifier)
# =============================================================================

# log "Import séquences PR2 dans QIIME2"
# conda run -n "$QIIME2_ENV" qiime tools import \
#     --type 'FeatureData[Sequence]' \
#     --input-path "${DBDIR}/pr2_version_5.0.0_SSU_mothur.fasta" \
#     --output-path "$PR2_SEQS_RAW"

# log "Import taxonomie PR2 dans QIIME2"
# conda run -n "$QIIME2_ENV" qiime tools import \
#     --type 'FeatureData[Taxonomy]' \
#     --input-format HeaderlessTSVTaxonomyFormat \
#     --input-path "${DBDIR}/pr2_v5_tax_qiime2.tsv" \
#     --output-path "$PR2_TAX_RAW"

# log "Extraction in-silico région V4-V5 sur PR2 (515F/926R)"
# conda run -n "$QIIME2_ENV" qiime feature-classifier extract-reads \
#     --i-sequences "$PR2_SEQS_RAW" \
#     --p-f-primer "$PRIMER_F" \
#     --p-r-primer "$PRIMER_R" \
#     --p-min-length 200 \
#     --p-max-length 600 \
#     --p-n-jobs "$NTHREADS" \
#     --o-reads "$PR2_SEQS_TRIMMED"

# log "Entraînement classifier PR2 V4-V5"
# conda run -n "$QIIME2_ENV" qiime feature-classifier fit-classifier-naive-bayes \
#     --i-reference-reads "$PR2_SEQS_TRIMMED" \
#     --i-reference-taxonomy "$PR2_TAX_RAW" \
#     --o-classifier "$PR2_CLASSIFIER"

# =============================================================================
# ÉTAPE 07d — CLASSIFICATION PR2 (18S — eucaryotes)
# =============================================================================
log "Classification taxonomique PR2 v5 (18S eucaryotes)"

conda run -n "$QIIME2_ENV" qiime feature-classifier classify-sklearn \
    --i-classifier "$PR2_CLASSIFIER" \
    --i-reads rep-seqs-euk-candidates.qza \
    --p-n-jobs "$NTHREADS" \
    --p-confidence 0.7 \
    --o-classification taxonomy-pr2-eukaryotes.qza

conda run -n "$QIIME2_ENV" qiime metadata tabulate \
    --m-input-file taxonomy-pr2-eukaryotes.qza \
    --o-visualization ../visual/taxonomy-pr2-eukaryotes.qzv

conda run -n "$QIIME2_ENV" qiime taxa barplot \
    --i-table table-euk-candidates.qza \
    --i-taxonomy taxonomy-pr2-eukaryotes.qza \
    --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
    --o-visualization ../visual/taxa-bar-plots-eukaryotes-pr2.qzv

# =============================================================================
# ÉTAPE 07e — Barplot procaryotes avec taxonomie SILVA (sans conta eucaryotes)
# =============================================================================
log "Barplot procaryotes (SILVA, Eukaryota/Mito/Chloro exclus)"

conda run -n "$QIIME2_ENV" qiime taxa barplot \
    --i-table table-prokaryotes.qza \
    --i-taxonomy taxonomy-silva-raw.qza \
    --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
    --o-visualization ../visual/taxa-bar-plots-prokaryotes-silva.qzv

# =============================================================================
# ÉTAPE 08 — RARÉFACTION
# Consulter table-prokaryotes-summary.qzv et table-euk-candidates-summary.qzv
# pour ajuster les seuils avant de lancer
# =============================================================================
log "Raréfaction procaryotes à ${RAREFACTION_DEPTH} reads"

conda run -n "$QIIME2_ENV" qiime feature-table rarefy \
    --i-table table-prokaryotes.qza \
    --p-sampling-depth "$RAREFACTION_DEPTH" \
    --o-rarefied-table ../subtables/RarTable-prokaryotes-depth${RAREFACTION_DEPTH}.qza

conda run -n "$QIIME2_ENV" qiime feature-table summarize \
    --i-table ../subtables/RarTable-prokaryotes-depth${RAREFACTION_DEPTH}.qza \
    --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
    --o-visualization ../visual/table-rarefied-prokaryotes-summary.qzv

log "Raréfaction eucaryotes à ${RAREFACTION_DEPTH_EUK} reads (ajuster si nécessaire)"

conda run -n "$QIIME2_ENV" qiime feature-table rarefy \
    --i-table table-euk-candidates.qza \
    --p-sampling-depth "$RAREFACTION_DEPTH_EUK" \
    --o-rarefied-table ../subtables/RarTable-eukaryotes-depth${RAREFACTION_DEPTH_EUK}.qza || \
    echo "AVERTISSEMENT: raréfaction eucaryotes échouée — vérifier seuil dans table-euk-candidates-summary.qzv"

# =============================================================================
# ÉTAPE 09 — ARBRE PHYLOGÉNÉTIQUE (procaryotes)
# =============================================================================
log "Arbre phylogénétique — procaryotes"

conda run -n "$QIIME2_ENV" qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences rep-seqs-prokaryotes.qza \
    --p-n-threads "$NTHREADS" \
    --o-alignment aligned-rep-seqs-prokaryotes.qza \
    --o-masked-alignment masked-aligned-rep-seqs-prokaryotes.qza \
    --o-tree unrooted-tree-prokaryotes.qza \
    --o-rooted-tree tree-prokaryotes.qza

# =============================================================================
# ÉTAPE 10 — CORE METRICS PHYLOGENETIC (procaryotes)
# =============================================================================
log "Core metrics phylogenetic — procaryotes"

conda run -n "$QIIME2_ENV" qiime diversity core-metrics-phylogenetic \
    --i-table table-prokaryotes.qza \
    --i-phylogeny tree-prokaryotes.qza \
    --p-sampling-depth "$RAREFACTION_DEPTH" \
    --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
    --o-rarefied-table rarefied_table_prokaryotes.qza \
    --o-faith-pd-vector diversity/Vector-faith_pd.qza \
    --o-observed-features-vector diversity/Vector-observed_asv.qza \
    --o-shannon-vector diversity/Vector-shannon.qza \
    --o-evenness-vector diversity/Vector-evenness.qza \
    --o-unweighted-unifrac-distance-matrix diversity/Matrix-unweighted_unifrac.qza \
    --o-weighted-unifrac-distance-matrix diversity/Matrix-weighted_unifrac.qza \
    --o-jaccard-distance-matrix diversity/Matrix-jaccard.qza \
    --o-bray-curtis-distance-matrix diversity/Matrix-braycurtis.qza \
    --o-unweighted-unifrac-pcoa-results pcoa/PCoA-unweighted_unifrac.qza \
    --o-weighted-unifrac-pcoa-results pcoa/PCoA-weighted_unifrac.qza \
    --o-jaccard-pcoa-results pcoa/PCoA-jaccard.qza \
    --o-bray-curtis-pcoa-results pcoa/PCoA-braycurtis.qza \
    --o-unweighted-unifrac-emperor ../visual/Emperor-unweighted_unifrac.qzv \
    --o-weighted-unifrac-emperor ../visual/Emperor-weighted_unifrac.qzv \
    --o-jaccard-emperor ../visual/Emperor-jaccard.qzv \
    --o-bray-curtis-emperor ../visual/Emperor-braycurtis.qzv

# Diversité alpha — tests de significativité
for metric in shannon evenness observed_features faith_pd; do
    vec="diversity/Vector-${metric}.qza"
    [[ "$metric" == "observed_features" ]] && vec="diversity/Vector-observed_asv.qza"
    [[ -f "$vec" ]] && conda run -n "$QIIME2_ENV" qiime diversity alpha-group-significance \
        --i-alpha-diversity "$vec" \
        --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
        --o-visualization "../visual/alpha-${metric}-significance.qzv" || true
done

# Diversité beta — PERMANOVA par sample_type et par site (PA vs KNS)
for matrix in braycurtis jaccard unweighted_unifrac weighted_unifrac; do
    mf="diversity/Matrix-${matrix}.qza"
    for col in sample_type location; do
        [[ -f "$mf" ]] && conda run -n "$QIIME2_ENV" qiime diversity beta-group-significance \
            --i-distance-matrix "$mf" \
            --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
            --m-metadata-column "$col" \
            --p-method permanova \
            --p-permutations 999 \
            --o-visualization "../visual/beta-${matrix}-permanova-${col}.qzv" || true
    done
done

conda run -n "$QIIME2_ENV" qiime feature-table core-features \
    --i-table ../subtables/RarTable-prokaryotes-depth${RAREFACTION_DEPTH}.qza \
    --p-min-fraction 0.1 \
    --p-max-fraction 1.0 \
    --p-steps 10 \
    --o-visualization ../visual/CoreBiom-prokaryotes.qzv || true

# =============================================================================
# ÉTAPE 11 — EXPORTS
# =============================================================================
log "Exports"

mkdir -p "${QDIR}/export/core/table-prokaryotes" \
         "${QDIR}/export/core/table-eukaryotes" \
         "${QDIR}/export/core/rep-seqs-prokaryotes" \
         "${QDIR}/export/core/rep-seqs-eukaryotes" \
         "${QDIR}/export/core/taxonomy-silva" \
         "${QDIR}/export/core/taxonomy-pr2" \
         "${QDIR}/export/subtables/RarTable-prokaryotes" \
         "${QDIR}/export/subtables/RarTable-eukaryotes" \
         "${QDIR}/export/diversity_tsv"

cd "${QDIR}"

conda run -n "$QIIME2_ENV" qiime tools export \
    --input-path core/table-prokaryotes.qza \
    --output-path export/core/table-prokaryotes

conda run -n "$QIIME2_ENV" qiime tools export \
    --input-path core/table-euk-candidates.qza \
    --output-path export/core/table-eukaryotes

conda run -n "$QIIME2_ENV" qiime tools export \
    --input-path core/rep-seqs-prokaryotes.qza \
    --output-path export/core/rep-seqs-prokaryotes

conda run -n "$QIIME2_ENV" qiime tools export \
    --input-path core/rep-seqs-euk-candidates.qza \
    --output-path export/core/rep-seqs-eukaryotes

conda run -n "$QIIME2_ENV" qiime tools export \
    --input-path core/taxonomy-silva-raw.qza \
    --output-path export/core/taxonomy-silva

conda run -n "$QIIME2_ENV" qiime tools export \
    --input-path core/taxonomy-pr2-eukaryotes.qza \
    --output-path export/core/taxonomy-pr2

conda run -n "$QIIME2_ENV" qiime tools export \
    --input-path subtables/RarTable-prokaryotes-depth${RAREFACTION_DEPTH}.qza \
    --output-path export/subtables/RarTable-prokaryotes

[[ -f "subtables/RarTable-eukaryotes-depth${RAREFACTION_DEPTH_EUK}.qza" ]] && \
    conda run -n "$QIIME2_ENV" qiime tools export \
        --input-path subtables/RarTable-eukaryotes-depth${RAREFACTION_DEPTH_EUK}.qza \
        --output-path export/subtables/RarTable-eukaryotes || true

# Fonction export diversité (accolade fermante corrigée)
export_tsv() {
    local qza="$1"
    local name="$2"
    local tmp="${QDIR}/export/diversity_tsv/${name}_tmp"
    [[ -f "$qza" ]] || return 0
    rm -rf "$tmp"
    conda run -n "$QIIME2_ENV" qiime tools export --input-path "$qza" --output-path "$tmp"
    find "$tmp" \( -name "*.tsv" -o -name "*.txt" -o -name "*.csv" \) -type f | \
        while read -r f; do
            cp "$f" "${QDIR}/export/diversity_tsv/${name}_$(basename "$f")"
        done
    rm -rf "$tmp"
}  # <- accolade manquante dans le script original

export_tsv core/diversity/Vector-faith_pd.qza         faith_pd
export_tsv core/diversity/Vector-shannon.qza           shannon
export_tsv core/diversity/Vector-observed_asv.qza      observed_features
export_tsv core/diversity/Vector-evenness.qza          evenness
export_tsv core/diversity/Matrix-braycurtis.qza        bray_curtis
export_tsv core/diversity/Matrix-jaccard.qza           jaccard
export_tsv core/diversity/Matrix-unweighted_unifrac.qza unweighted_unifrac
export_tsv core/diversity/Matrix-weighted_unifrac.qza   weighted_unifrac
export_tsv core/pcoa/PCoA-braycurtis.qza               pcoa_braycurtis
export_tsv core/pcoa/PCoA-jaccard.qza                  pcoa_jaccard
export_tsv core/pcoa/PCoA-unweighted_unifrac.qza       pcoa_unweighted_unifrac
export_tsv core/pcoa/PCoA-weighted_unifrac.qza         pcoa_weighted_unifrac

# Conversion BIOM → TSV + fusion taxonomies (SILVA pour prokaryotes, PR2 pour eucaryotes)
for SUBSET in prokaryotes eukaryotes; do
    if [[ "$SUBSET" == "prokaryotes" ]]; then
        BIOM_FILE="${QDIR}/export/subtables/RarTable-prokaryotes/feature-table.biom"
        TAX_FILE="${QDIR}/export/core/taxonomy-silva/taxonomy.tsv"
        OUT_PREFIX="${QDIR}/export/subtables/RarTable-prokaryotes"
    else
        BIOM_FILE="${QDIR}/export/subtables/RarTable-eukaryotes/feature-table.biom"
        TAX_FILE="${QDIR}/export/core/taxonomy-pr2/taxonomy.tsv"
        OUT_PREFIX="${QDIR}/export/subtables/RarTable-eukaryotes"
    fi

    [[ -f "$BIOM_FILE" ]] || { echo "BIOM absent pour $SUBSET, skip export TSV"; continue; }

    conda run -n "$QIIME2_ENV" biom convert \
        -i "$BIOM_FILE" \
        -o "${OUT_PREFIX}/table-from-biom.tsv" \
        --to-tsv

    sed '1d; s/#OTU ID/ASV_ID/' "${OUT_PREFIX}/table-from-biom.tsv" \
        > "${OUT_PREFIX}/ASV.tsv"

    [[ -f "$TAX_FILE" ]] || { echo "Taxonomie absente pour $SUBSET, skip fusion"; continue; }

    python3 << PYEOF
import csv, re

asv_path   = "${OUT_PREFIX}/ASV.tsv"
tax_path   = "${TAX_FILE}"
out_path   = "${OUT_PREFIX}/ASV_taxonomy.tsv"
subset     = "${SUBSET}"

taxonomy = {}
with open(tax_path) as fh:
    r = csv.reader(fh, delimiter='\t')
    next(r)
    for row in r:
        if len(row) >= 2:
            taxonomy[row[0]] = row[1]

def parse_silva(s):
    pats = [r'D_0__([^;]+)', r'D_1__([^;]+)', r'D_2__([^;]+)',
            r'D_3__([^;]+)', r'D_4__([^;]+)', r'D_5__([^;]+)', r'D_6__([^;]+)']
    out = []
    for p in pats:
        m = re.search(p, s or "")
        out.append(m.group(1).strip() if m else "Unassigned")
    return out

def parse_pr2(s):
    # PR2 mothur format: Kingdom;Supergroup;Division;Class;Order;Family;Genus;Species
    levels = [l.strip() for l in (s or "").split(";")]
    while len(levels) < 8:
        levels.append("Unassigned")
    return levels[:8]

if subset == "prokaryotes":
    headers  = ["ASV_ID","Kingdom","Phylum","Class","Order","Family","Genus","Species"]
    parse_fn = parse_silva
else:
    headers  = ["ASV_ID","Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"]
    parse_fn = parse_pr2

with open(asv_path) as fi, open(out_path, 'w') as fo:
    r = csv.reader(fi, delimiter='\t')
    h = next(r)
    fo.write('\t'.join(headers + h[1:]) + '\n')
    for row in r:
        taxon_fields = parse_fn(taxonomy.get(row[0], ""))
        fo.write('\t'.join([row[0]] + taxon_fields + row[1:]) + '\n')

print(f"Table {subset} exportée : {out_path}")
PYEOF
done

log "SCRIPT 2 TERMINÉ"
echo ""
echo "======================================================================="
echo "  RÉSULTATS DANS : ${QDIR}/"
echo ""
echo "  PROCARYOTES (SILVA 138.2):"
echo "  → visual/taxa-bar-plots-prokaryotes-silva.qzv"
echo "  → export/subtables/RarTable-prokaryotes/ASV_taxonomy.tsv"
echo ""
echo "  EUCARYOTES (PR2 v5):"
echo "  → visual/taxa-bar-plots-eukaryotes-pr2.qzv"
echo "  → visual/taxonomy-pr2-eukaryotes.qzv"
echo "  → export/subtables/RarTable-eukaryotes/ASV_taxonomy.tsv"
echo ""
echo "  DIVERSITÉ:"
echo "  → visual/Emperor-*.qzv                   PCoA interactifs"
echo "  → visual/alpha-*-significance.qzv         tests alpha"
echo "  → visual/beta-*-permanova-*.qzv           PA vs KNS + sample_type"
echo "  → export/diversity_tsv/                   matrices TSV brutes"
echo "======================================================================="
