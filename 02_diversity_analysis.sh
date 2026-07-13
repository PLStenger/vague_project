#!/usr/bin/env bash
# =============================================================================
# SCRIPT 2/2 — Taxonomie + Diversité | PA + KNS
# MODIFIÉ : double classification 16S (SILVA) + 18S (PR2) pour eucaryotes
# Usage   : bash 02_diversity_analysis.sh --depth 42105
# =============================================================================

#set -euo pipefail

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

log() { echo -e "\n[$(date +'%F %T')] === $* ===\n"; }
mkdir -p "$TMPDIR" "${QDIR}/core/diversity" "${QDIR}/core/pcoa" \
         "${QDIR}/visual" "${QDIR}/subtables" "${QDIR}/export"
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
set -u

for f in "${QDIR}/core/table-decontam.qza" \
          "${QDIR}/core/rep-seqs-decontam.qza" \
          "${DBDIR}/sample-metadata.tsv"; do
    [[ -f "$f" ]] || { echo "Fichier manquant: $f"; exit 1; }
done

cd "${QDIR}/core"

# =============================================================================
# ÉTAPES DÉJÀ EFFECTUÉES — commentées
# =============================================================================

### log "Préparation listes contrôles / non-contrôles"
### awk -F'\t' '...' "${DBDIR}/sample-metadata.tsv" > "${DBDIR}/control-samples.txt"
### awk -F'\t' '...' "${DBDIR}/sample-metadata.tsv" > "${DBDIR}/non-control-samples.txt"
##
### log "Extraction table des contrôles"
### conda run -n "$QIIME2_ENV" qiime feature-table filter-samples \
###     --i-table table.qza \
###     --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
###     --p-where "[is_control]='yes'" \
###     --o-filtered-table table-controls-only.qza
### conda run -n "$QIIME2_ENV" qiime feature-table summarize \
###     --i-table table-controls-only.qza \
###     --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
###     --o-visualization ../visual/table-controls-only-summary.qzv
##
### log "Extraction liste des ASV présents dans les contrôles"
### mkdir -p "${TMPDIR}/controls_export"
### rm -rf "${TMPDIR}/controls_export"/* 2>/dev/null || true
### conda run -n "$QIIME2_ENV" qiime tools export \
###     --input-path table-controls-only.qza \
###     --output-path "${TMPDIR}/controls_export"
### conda run -n "$QIIME2_ENV" biom convert \
###     -i "${TMPDIR}/controls_export/feature-table.biom" \
###     -o "${TMPDIR}/controls_export/controls_table.tsv" \
###     --to-tsv
### { echo -e "feature-id"; \
###   awk 'BEGIN{FS="\t"} NR>2 {sum=0; for(i=2;i<=NF;i++) sum+=$i; if(sum>0) print $1}' \
###   "${TMPDIR}/controls_export/controls_table.tsv"; \
### } > "${DBDIR}/features-in-controls.txt"
##
### log "Extraction table échantillons biologiques"
### conda run -n "$QIIME2_ENV" qiime feature-table filter-samples \
###     --i-table table.qza \
###     --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
###     --p-where "[is_control]!='yes'" \
###     --o-filtered-table table-non-controls.qza
### conda run -n "$QIIME2_ENV" qiime feature-table summarize \
###     --i-table table-non-controls.qza \
###     --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
###     --o-visualization ../visual/table-non-controls-summary.qzv
##
### log "Retrait des ASV présents dans les contrôles"
### conda run -n "$QIIME2_ENV" qiime feature-table filter-features \
###     --i-table table-non-controls.qza \
###     --m-metadata-file "${DBDIR}/features-in-controls.txt" \
###     --p-exclude-ids \
###     --o-filtered-table table-decontam.qza
### conda run -n "$QIIME2_ENV" qiime feature-table summarize \
###     --i-table table-decontam.qza \
###     --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
###     --o-visualization ../visual/table-decontam-summary.qzv
##
### log "Filtrage rep-seqs sur ASV décontaminés"
### conda run -n "$QIIME2_ENV" qiime feature-table filter-seqs \
###     --i-data rep-seqs.qza \
###     --i-table table-decontam.qza \
###     --o-filtered-data rep-seqs-decontam.qza
##
### =============================================================================
### ÉTAPE 07a — CLASSIFICATION SILVA 138.2 (16S — bactéries + archées)
### À LANCER
### =============================================================================
##log "Classification taxonomique SILVA 138.2 (16S)"
##
##CLASSIFIER_SILVA="${DBDIR}/silva-138.2-ssu-nr99-515f-926r-classifier.qza"
##CLASSIFIER_SOURCE="/nvme/bio/data_fungi/valormicro_nc/98_databasefiles/silva-138.2-ssu-nr99-515f-926r-classifier.qza"
##
##if [[ ! -f "$CLASSIFIER_SILVA" ]]; then
##    log "Copie du classifier SILVA depuis valormicro_nc"
##    cp "$CLASSIFIER_SOURCE" "$CLASSIFIER_SILVA" || \
##        { log "ERREUR: Classifier SILVA source introuvable"; exit 1; }
##fi
##
##conda run -n "$QIIME2_ENV" qiime tools validate "$CLASSIFIER_SILVA" || \
##    { log "ERREUR: Classifier SILVA invalide"; exit 1; }
##
##conda run -n "$QIIME2_ENV" qiime feature-classifier classify-sklearn \
##    --i-classifier "$CLASSIFIER_SILVA" \
##    --i-reads rep-seqs-decontam.qza \
##    --p-n-jobs "$NTHREADS" \
##    --o-classification taxonomy-silva-raw.qza
##
##conda run -n "$QIIME2_ENV" qiime metadata tabulate \
##    --m-input-file taxonomy-silva-raw.qza \
##    --o-visualization ../visual/taxonomy-silva-raw.qzv
##
### =============================================================================
### ÉTAPE 07b — FILTRAGE POST-SILVA : séparer prokaryotes et candidats eucaryotes
### Les ASVs "Unassigned" ou classées Eukaryota par SILVA sont des candidats 18S
### À LANCER
### =============================================================================
##log "Séparation ASVs 16S (procaryotes) / candidats 18S (eucaryotes)"
##
### Table et séquences strictement procaryotes (bactéries + archées)
### Exclut Eukaryota, Mitochondria, Chloroplast, Unassigned
##conda run -n "$QIIME2_ENV" qiime taxa filter-table \
##    --i-table table-decontam.qza \
##    --i-taxonomy taxonomy-silva-raw.qza \
##    --p-exclude "Eukaryota,Mitochondria,Chloroplast,Unassigned" \
##    --o-filtered-table table-prokaryotes.qza
##
##conda run -n "$QIIME2_ENV" qiime taxa filter-seqs \
##    --i-sequences rep-seqs-decontam.qza \
##    --i-taxonomy taxonomy-silva-raw.qza \
##    --p-exclude "Eukaryota,Mitochondria,Chloroplast,Unassigned" \
##    --o-filtered-sequences rep-seqs-prokaryotes.qza
##
##conda run -n "$QIIME2_ENV" qiime feature-table summarize \
##    --i-table table-prokaryotes.qza \
##    --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
##    --o-visualization ../visual/table-prokaryotes-summary.qzv
##
### Table et séquences candidats eucaryotes :
### inclut Eukaryota assigné par SILVA + tout ce que SILVA n'a pas su assigner
##conda run -n "$QIIME2_ENV" qiime taxa filter-table \
##    --i-table table-decontam.qza \
##    --i-taxonomy taxonomy-silva-raw.qza \
##    --p-include "Eukaryota,Unassigned" \
##    --o-filtered-table table-euk-candidates.qza
##
##conda run -n "$QIIME2_ENV" qiime taxa filter-seqs \
##    --i-sequences rep-seqs-decontam.qza \
##    --i-taxonomy taxonomy-silva-raw.qza \
##    --p-include "Eukaryota,Unassigned" \
##    --o-filtered-sequences rep-seqs-euk-candidates.qza
##
##conda run -n "$QIIME2_ENV" qiime feature-table summarize \
##    --i-table table-euk-candidates.qza \
##    --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
##    --o-visualization ../visual/table-euk-candidates-summary.qzv
##
### =============================================================================
### ÉTAPE 07c — TÉLÉCHARGEMENT + PRÉPARATION CLASSIFIER PR2 V4-V5 POUR 18S
### PR2 v5 est la référence pour protistes / champignons / métazoaires marins
### À LANCER (téléchargement une seule fois, long ~15 min selon connexion)
### =============================================================================
##log "Préparation classifier PR2 v5 pour les eucaryotes (18S V4-V5)"
##
##PR2_CLASSIFIER="${DBDIR}/pr2-v5-classifier-515F-926R.qza"
##PR2_SEQS_RAW="${DBDIR}/pr2_v5_SSU_seqs.qza"
##PR2_TAX_RAW="${DBDIR}/pr2_v5_SSU_tax.qza"
##PR2_SEQS_TRIMMED="${DBDIR}/pr2_v5_515F_926R_seqs.qza"
##
##if [[ ! -f "$PR2_CLASSIFIER" ]]; then
##    log "Téléchargement PR2 v5.0.0 depuis GitHub releases"
##
##    # Séquences et taxonomie PR2 v5 (full-length SSU, format QIIME2)
##    wget -q -O "${DBDIR}/pr2_v5.0.0_SSU_seqs.fasta.gz" \
##        "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU.fasta.gz"
##    wget -q -O "${DBDIR}/pr2_v5.0.0_SSU_tax.tsv" \
##        "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU.tsv"
##
##    # Reformater la taxonomie PR2 au format QIIME2 (7 niveaux tab-séparés)
##    python3 << 'PYEOF'
##import csv, gzip, os
##
##tax_in  = "/nvme/bio/data_fungi/vague_project/98_databasefiles/pr2_v5.0.0_SSU_tax.tsv"
##tax_out = "/nvme/bio/data_fungi/vague_project/98_databasefiles/pr2_v5_tax_qiime2.tsv"
##
### Les colonnes PR2 : pr2_accession, kingdom, supergroup, division, class, order, family, genus, species, ...
### On construit une chaîne taxonomique à 8 niveaux compatible QIIME2
##with open(tax_in, newline='') as fi, open(tax_out, 'w') as fo:
##    fo.write("Feature ID\tTaxon\n")
##    r = csv.DictReader(fi, delimiter='\t')
##    for row in r:
##        acc = row.get('pr2_accession') or row.get('seqid') or list(row.values())[0]
##        levels = [
##            row.get('kingdom',''),
##            row.get('supergroup',''),
##            row.get('division',''),
##            row.get('class',''),
##            row.get('order',''),
##            row.get('family',''),
##            row.get('genus',''),
##            row.get('species',''),
##        ]
##        taxon = "; ".join(f"d__{l}" if i==0 else
##                          f"p__{l}" if i==1 else
##                          f"c__{l}" if i==3 else
##                          f"o__{l}" if i==4 else
##                          f"f__{l}" if i==5 else
##                          f"g__{l}" if i==6 else
##                          f"s__{l}" if i==7 else l
##                          for i,l in enumerate(levels) if l)
##        fo.write(f"{acc}\t{taxon}\n")
##print("Taxonomie PR2 reformatée")
##PYEOF
##
##    # Importer les séquences FASTA PR2 dans QIIME2
##    gunzip -f "${DBDIR}/pr2_v5.0.0_SSU_seqs.fasta.gz" -c \
##        > "${DBDIR}/pr2_v5.0.0_SSU_seqs.fasta"


# ============================================================
# IMPORT DANS QIIME2
# ============================================================
conda run -n qiime2-amplicon-2025.7 qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path /nvme/bio/data_fungi/vague_project/98_databasefiles/pr2_version_5.0.0_SSU_mothur.fasta \
    --output-path /nvme/bio/data_fungi/vague_project/98_databasefiles/pr2_v5_SSU_seqs.qza

conda run -n qiime2-amplicon-2025.7 qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path /nvme/bio/data_fungi/vague_project/98_databasefiles/pr2_v5_tax_qiime2.tsv \
    --output-path /nvme/bio/data_fungi/vague_project/98_databasefiles/pr2_v5_SSU_tax.qza

# ============================================================
# EXTRACTION IN-SILICO V4-V5 avec les amorces 515F/926R
# ============================================================
conda run -n qiime2-amplicon-2025.7 qiime feature-classifier extract-reads \
    --i-sequences /nvme/bio/data_fungi/vague_project/98_databasefiles/pr2_v5_SSU_seqs.qza \
    --p-f-primer GTGYCAGCMGCCGCGGTAA \
    --p-r-primer CCGYCAATTYMTTTRAGTTT \
    --p-min-length 200 \
    --p-max-length 600 \
    --p-n-jobs 16 \
    --o-reads /nvme/bio/data_fungi/vague_project/98_databasefiles/pr2_v5_515F_926R_seqs.qza

# ============================================================
# ENTRAÎNEMENT DU CLASSIFIER 
# ============================================================
conda run -n qiime2-amplicon-2025.7 qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads /nvme/bio/data_fungi/vague_project/98_databasefiles/pr2_v5_515F_926R_seqs.qza \
    --i-reference-taxonomy /nvme/bio/data_fungi/vague_project/98_databasefiles/pr2_v5_SSU_tax.qza \
    --o-classifier /nvme/bio/data_fungi/vague_project/98_databasefiles/pr2-v5-classifier-515F-926R.qza



echo "Suite ancien pipeline"


    conda run -n "$QIIME2_ENV" qiime tools import \
        --type 'FeatureData[Sequence]' \
        --input-path "${DBDIR}/pr2_v5.0.0_SSU_seqs.fasta" \
        --output-path "$PR2_SEQS_RAW"

    conda run -n "$QIIME2_ENV" qiime tools import \
        --type 'FeatureData[Taxonomy]' \
        --input-format HeaderlessTSVTaxonomyFormat \
        --input-path "${DBDIR}/pr2_v5_tax_qiime2.tsv" \
        --output-path "$PR2_TAX_RAW"

    # Extraction in-silico de la région V4-V5 avec nos amorces
    log "Extraction in-silico région V4-V5 sur PR2 (515F/926R)"
    conda run -n "$QIIME2_ENV" qiime feature-classifier extract-reads \
        --i-sequences "$PR2_SEQS_RAW" \
        --p-f-primer "${PRIMER_F}" \
        --p-r-primer "${PRIMER_R}" \
        --p-min-length 200 \
        --p-max-length 600 \
        --p-n-jobs "$NTHREADS" \
        --o-reads "$PR2_SEQS_TRIMMED"

    # Entraînement du classifier naïf-Bayes PR2 V4-V5
    log "Entraînement classifier PR2 V4-V5 (peut prendre 20-40 min)"
    conda run -n "$QIIME2_ENV" qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads "$PR2_SEQS_TRIMMED" \
        --i-reference-taxonomy "$PR2_TAX_RAW" \
        --o-classifier "$PR2_CLASSIFIER"

    log "Classifier PR2 V4-V5 prêt : $PR2_CLASSIFIER"
else
    log "Classifier PR2 déjà disponible, skip téléchargement"
fi

# =============================================================================
# ÉTAPE 07d — CLASSIFICATION PR2 (18S — eucaryotes)
# À LANCER
# =============================================================================
log "Classification taxonomique PR2 v5 (18S eucaryotes)"

cd "${QDIR}/core"

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
# ÉTAPE 07e — Barplot procaryotes avec taxonomie SILVA filtrée
# À LANCER
# =============================================================================
log "Barplot procaryotes (SILVA, sans Eukaryota/Mito/Chloro)"

conda run -n "$QIIME2_ENV" qiime taxa filter-table \
    --i-table table-decontam.qza \
    --i-taxonomy taxonomy-silva-raw.qza \
    --p-exclude "Eukaryota,Mitochondria,Chloroplast,Unassigned" \
    --o-filtered-table table-decontam-prokaryotes-only.qza

conda run -n "$QIIME2_ENV" qiime taxa barplot \
    --i-table table-decontam-prokaryotes-only.qza \
    --i-taxonomy taxonomy-silva-raw.qza \
    --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
    --o-visualization ../visual/taxa-bar-plots-prokaryotes-silva.qzv

# =============================================================================
# ÉTAPE 08 — RARÉFACTION (sur table procaryotes + eucaryotes séparément)
# À LANCER
# =============================================================================
log "Raréfaction à ${RAREFACTION_DEPTH} reads — procaryotes"

conda run -n "$QIIME2_ENV" qiime feature-table rarefy \
    --i-table table-prokaryotes.qza \
    --p-sampling-depth "$RAREFACTION_DEPTH" \
    --o-rarefied-table ../subtables/RarTable-prokaryotes-depth${RAREFACTION_DEPTH}.qza

conda run -n "$QIIME2_ENV" qiime feature-table summarize \
    --i-table ../subtables/RarTable-prokaryotes-depth${RAREFACTION_DEPTH}.qza \
    --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
    --o-visualization ../visual/table-rarefied-prokaryotes-summary.qzv

log "Raréfaction à ${RAREFACTION_DEPTH} reads — eucaryotes"

# Note : adapter le seuil pour les eucaryotes si nécessaire
# (souvent beaucoup moins de reads 18S que de reads 16S dans un jeu de données 16S)
RAREFACTION_DEPTH_EUK=1000
log "Profondeur de raréfaction eucaryotes : ${RAREFACTION_DEPTH_EUK} (à ajuster après avoir consulté table-euk-candidates-summary.qzv)"

conda run -n "$QIIME2_ENV" qiime feature-table rarefy \
    --i-table table-euk-candidates.qza \
    --p-sampling-depth "$RAREFACTION_DEPTH_EUK" \
    --o-rarefied-table ../subtables/RarTable-eukaryotes-depth${RAREFACTION_DEPTH_EUK}.qza || \
    echo "AVERTISSEMENT: raréfaction eucaryotes échouée — vérifier le seuil dans table-euk-candidates-summary.qzv"

# =============================================================================
# ÉTAPE 09 — ARBRE PHYLOGÉNÉTIQUE (procaryotes)
# À LANCER
# =============================================================================
log "Arbre phylogénétique sur séquences procaryotes décontaminées"

conda run -n "$QIIME2_ENV" qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences rep-seqs-prokaryotes.qza \
    --p-n-threads "$NTHREADS" \
    --o-alignment aligned-rep-seqs-prokaryotes.qza \
    --o-masked-alignment masked-aligned-rep-seqs-prokaryotes.qza \
    --o-tree unrooted-tree-prokaryotes.qza \
    --o-rooted-tree tree-prokaryotes.qza

# =============================================================================
# ÉTAPE 10 — CORE METRICS PHYLOGENETIC (procaryotes)
# À LANCER
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

# Diversité alpha — significativité
for metric in shannon evenness observed_features faith_pd; do
    vec="diversity/Vector-${metric}.qza"
    [[ "$metric" == "observed_features" ]] && vec="diversity/Vector-observed_asv.qza"
    [[ -f "$vec" ]] && conda run -n "$QIIME2_ENV" qiime diversity alpha-group-significance \
        --i-alpha-diversity "$vec" \
        --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
        --o-visualization "../visual/alpha-${metric}-significance.qzv" || true
done

# Diversité beta — PERMANOVA sample_type et location
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

# Export tables et séquences
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

# Fonction export diversité
export_tsv() {
    local qza="$1"; local name="$2"; local tmp="${QDIR}/export/diversity_tsv/${name}_tmp"
    [[ -f "$qza" ]] || return 0
    rm -rf "$tmp"
    conda run -n "$QIIME2_ENV" qiime tools export --input-path "$qza" --output-path "$tmp"
    find "$tmp" \( -name "*.tsv" -o -name "*.txt" -o -name "*.csv" \) -type f | \
        while read -r f; do cp "$f" "${QDIR}/export/diversity_tsv/${name}_$(basename "$f")"; done
    rm -rf "$tmp"
}

export_tsv core/diversity/Vector-faith_pd.qza       faith_pd
export_tsv core/diversity/Vector-shannon.qza         shannon
export_tsv core/diversity/Vector-observed_asv.qza    observed_features
export_tsv core/diversity/Vector-evenness.qza        evenness
export_tsv core/diversity/Matrix-braycurtis.qza      bray_curtis
export_tsv core/diversity/Matrix-jaccard.qza         jaccard
export_tsv core/diversity/Matrix-unweighted_unifrac.qza unweighted_unifrac
export_tsv core/diversity/Matrix-weighted_unifrac.qza   weighted_unifrac
export_tsv core/pcoa/PCoA-braycurtis.qza             pcoa_braycurtis
export_tsv core/pcoa/PCoA-jaccard.qza                pcoa_jaccard
export_tsv core/pcoa/PCoA-unweighted_unifrac.qza     pcoa_unweighted_unifrac
export_tsv core/pcoa/PCoA-weighted_unifrac.qza       pcoa_weighted_unifrac

# Conversion BIOM → TSV + fusion taxonomies (procaryotes SILVA + eucaryotes PR2)
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

    [[ -f "$BIOM_FILE" ]] || continue

    conda run -n "$QIIME2_ENV" biom convert \
        -i "$BIOM_FILE" \
        -o "${OUT_PREFIX}/table-from-biom.tsv" \
        --to-tsv

    sed '1d; s/#OTU ID/ASV_ID/' "${OUT_PREFIX}/table-from-biom.tsv" \
        > "${OUT_PREFIX}/ASV.tsv"

    [[ -f "$TAX_FILE" ]] && python3 << PYEOF
import csv, re
asv_path = "${OUT_PREFIX}/ASV.tsv"
tax_path = "${TAX_FILE}"
out_path = "${OUT_PREFIX}/ASV_taxonomy.tsv"
subset   = "${SUBSET}"

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
    # PR2 taxonomy : "d__Kingdom; p__Division; c__Class; o__Order; f__Family; g__Genus; s__Species"
    levels = [l.strip() for l in (s or "").split(";")]
    out = []
    for l in levels:
        m = re.match(r'[a-z]__(.+)', l)
        out.append(m.group(1).strip() if m else l)
    # Padder à 8 niveaux (Kingdom, Supergroup, Division, Class, Order, Family, Genus, Species)
    while len(out) < 8:
        out.append("Unassigned")
    return out[:8]

if subset == "prokaryotes":
    headers = ["ASV_ID","Kingdom","Phylum","Class","Order","Family","Genus","Species"]
    parse_fn = parse_silva
else:
    headers = ["ASV_ID","Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"]
    def parse_fn(s):
        return parse_pr2(s)

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
echo "  → visual/Emperor-*.qzv          PCoA interactifs"
echo "  → visual/alpha-*-significance.qzv"
echo "  → visual/beta-*-permanova-*.qzv  PA vs KNS + sample_type"
echo "======================================================================="
