#!/usr/bin/env bash
# =============================================================================
# SCRIPT 2/2 — Décontamination simple par contrôles + taxonomie + diversité
# Projet  : vague_project | 16S V4-V5 | Pointe de l'Artillerie
# Chemin  : /nvme/bio/data_fungi/vague_project
# Usage   : bash 02_diversity_analysis.sh --depth 42105
# Principe: les contrôles identifient les ASV contaminants, retirés des
#           échantillons biologiques via feature-table filter-features.
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

log() { echo -e "\n[$(date +'%F %T')] === $* ===\n"; }
mkdir -p "$TMPDIR" "${QDIR}/core/diversity" "${QDIR}/core/pcoa" "${QDIR}/visual" "${QDIR}/subtables" "${QDIR}/export"
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
set -u

for f in "${QDIR}/core/table.qza" "${QDIR}/core/rep-seqs.qza" "${DBDIR}/sample-metadata.tsv"; do
    [[ -f "$f" ]] || { echo "Fichier manquant: $f"; exit 1; }
done

cd "${QDIR}/core"

# =============================================================================
# ÉTAPES DÉJÀ EFFECTUÉES — conservées ici pour référence
# =============================================================================

# log "Préparation listes contrôles / non-contrôles"
# awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} NR==2 && $1=="#q2:types"{next} $c["is_control"]=="yes"{print $1}' "${DBDIR}/sample-metadata.tsv" > "${DBDIR}/control-samples.txt"
# awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} NR==2 && $1=="#q2:types"{next} $c["is_control"]!="yes"{print $1}' "${DBDIR}/sample-metadata.tsv" > "${DBDIR}/non-control-samples.txt"

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

# log "Extraction liste des ASV présents dans les contrôles (359 ASV trouvés)"
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

# log "Extraction table des échantillons biologiques uniquement"
# conda run -n "$QIIME2_ENV" qiime feature-table filter-samples \
#     --i-table table.qza \
#     --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
#     --p-where "[is_control]!='yes'" \
#     --o-filtered-table table-non-controls.qza
# conda run -n "$QIIME2_ENV" qiime feature-table summarize \
#     --i-table table-non-controls.qza \
#     --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
#     --o-visualization ../visual/table-non-controls-summary.qzv

# log "Retrait des ASV présents dans les contrôles (--p-exclude-ids)"
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
# ÉTAPE 07 — CLASSIFICATION TAXONOMIQUE (SILVA 138.2 V4-V5)
# Utilisation du classifier pré-entraîné du projet valormicro_nc
# =============================================================================
log "Classification taxonomique SILVA 138.2"

cd "${DBDIR}"

CLASSIFIER_SOURCE="/nvme/bio/data_fungi/valormicro_nc/98_databasefiles/silva-138.2-ssu-nr99-515f-926r-classifier.qza"
CLASSIFIER="${DBDIR}/silva-138.2-ssu-nr99-515f-926r-classifier.qza"

if [[ ! -f "$CLASSIFIER" ]]; then
    log "Copie du classifier depuis valormicro_nc"
    cp "$CLASSIFIER_SOURCE" "$CLASSIFIER" || { log "ERREUR: Classifier source introuvable ($CLASSIFIER_SOURCE)"; exit 1; }
fi

conda run -n "$QIIME2_ENV" qiime tools validate "$CLASSIFIER" || { log "ERREUR: Classifier invalide"; exit 1; }
log "Classifier SILVA 138.2 V4-V5 validé"

cd "${QDIR}/core"

conda run -n "$QIIME2_ENV" qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads rep-seqs-decontam.qza \
    --p-n-jobs "$NTHREADS" \
    --o-classification taxonomy-decontam.qza

log "Classification terminée"

conda run -n "$QIIME2_ENV" qiime metadata tabulate \
    --m-input-file taxonomy-decontam.qza \
    --o-visualization ../visual/taxonomy-decontam.qzv

conda run -n "$QIIME2_ENV" qiime taxa barplot \
    --i-table table-decontam.qza \
    --i-taxonomy taxonomy-decontam.qza \
    --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
    --o-visualization ../visual/taxa-bar-plots-decontam.qzv

# =============================================================================
# ÉTAPE 08 — RARÉFACTION
# =============================================================================
log "Raréfaction à ${RAREFACTION_DEPTH} reads"

conda run -n "$QIIME2_ENV" qiime feature-table rarefy \
    --i-table table-decontam.qza \
    --p-sampling-depth "$RAREFACTION_DEPTH" \
    --o-rarefied-table ../subtables/RarTable-decontam-depth${RAREFACTION_DEPTH}.qza

conda run -n "$QIIME2_ENV" qiime feature-table summarize \
    --i-table ../subtables/RarTable-decontam-depth${RAREFACTION_DEPTH}.qza \
    --m-sample-metadata-file "${DBDIR}/sample-metadata.tsv" \
    --o-visualization ../visual/table-rarefied-decontam-summary.qzv

# =============================================================================
# ÉTAPE 09 — ARBRE PHYLOGÉNÉTIQUE
# =============================================================================
log "Arbre phylogénétique sur séquences décontaminées"

conda run -n "$QIIME2_ENV" qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences rep-seqs-decontam.qza \
    --p-n-threads "$NTHREADS" \
    --o-alignment aligned-rep-seqs-decontam.qza \
    --o-masked-alignment masked-aligned-rep-seqs-decontam.qza \
    --o-tree unrooted-tree-decontam.qza \
    --o-rooted-tree tree-decontam.qza

# =============================================================================
# ÉTAPE 10 — CORE METRICS PHYLOGENETIC
# =============================================================================
log "Core metrics phylogenetic sur données décontaminées"

conda run -n "$QIIME2_ENV" qiime diversity core-metrics-phylogenetic \
    --i-table table-decontam.qza \
    --i-phylogeny tree-decontam.qza \
    --p-sampling-depth "$RAREFACTION_DEPTH" \
    --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
    --o-rarefied-table rarefied_table_decontam.qza \
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

# Diversité beta — PERMANOVA sample_type
for matrix in braycurtis jaccard unweighted_unifrac weighted_unifrac; do
    mf="diversity/Matrix-${matrix}.qza"
    [[ -f "$mf" ]] && conda run -n "$QIIME2_ENV" qiime diversity beta-group-significance \
        --i-distance-matrix "$mf" \
        --m-metadata-file "${DBDIR}/sample-metadata.tsv" \
        --m-metadata-column sample_type \
        --p-method permanova \
        --p-permutations 999 \
        --o-visualization "../visual/beta-${matrix}-permanova-sampletype.qzv" || true
done

conda run -n "$QIIME2_ENV" qiime feature-table core-features \
    --i-table ../subtables/RarTable-decontam-depth${RAREFACTION_DEPTH}.qza \
    --p-min-fraction 0.1 \
    --p-max-fraction 1.0 \
    --p-steps 10 \
    --o-visualization ../visual/CoreBiom-decontam.qzv || true

# =============================================================================
# ÉTAPE 11 — EXPORTS
# =============================================================================
log "Exports"

mkdir -p "${QDIR}/export/core/table" \
         "${QDIR}/export/core/rep-seqs" \
         "${QDIR}/export/core/taxonomy" \
         "${QDIR}/export/subtables/RarTable" \
         "${QDIR}/export/diversity_tsv"

cd "${QDIR}"

conda run -n "$QIIME2_ENV" qiime tools export \
    --input-path core/table-decontam.qza \
    --output-path export/core/table
conda run -n "$QIIME2_ENV" qiime tools export \
    --input-path core/rep-seqs-decontam.qza \
    --output-path export/core/rep-seqs
conda run -n "$QIIME2_ENV" qiime tools export \
    --input-path core/taxonomy-decontam.qza \
    --output-path export/core/taxonomy
conda run -n "$QIIME2_ENV" qiime tools export \
    --input-path subtables/RarTable-decontam-depth${RAREFACTION_DEPTH}.qza \
    --output-path export/subtables/RarTable

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

# Conversion BIOM → TSV + fusion taxonomie SILVA
BIOM_FILE="${QDIR}/export/subtables/RarTable/feature-table.biom"
TAX_FILE="${QDIR}/export/core/taxonomy/taxonomy.tsv"

if [[ -f "$BIOM_FILE" ]]; then
    conda run -n "$QIIME2_ENV" biom convert \
        -i "$BIOM_FILE" \
        -o "${QDIR}/export/subtables/RarTable/table-from-biom.tsv" \
        --to-tsv
    sed '1d; s/#OTU ID/ASV_ID/' \
        "${QDIR}/export/subtables/RarTable/table-from-biom.tsv" \
        > "${QDIR}/export/subtables/RarTable/ASV.tsv"

    if [[ -f "$TAX_FILE" ]]; then
python3 << PYEOF
import csv, re
asv_path = "${QDIR}/export/subtables/RarTable/ASV.tsv"
tax_path = "${QDIR}/export/core/taxonomy/taxonomy.tsv"
out_path = "${QDIR}/export/subtables/RarTable/ASV_taxonomy.tsv"
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
with open(asv_path) as fi, open(out_path, 'w') as fo:
    r = csv.reader(fi, delimiter='\t')
    h = next(r)
    fo.write('\t'.join(["ASV_ID","Kingdom","Phylum","Class","Order","Family","Genus","Species"] + h[1:]) + '\n')
    for row in r:
        fo.write('\t'.join([row[0]] + parse_silva(taxonomy.get(row[0], "")) + row[1:]) + '\n')
PYEOF
    fi
fi

log "SCRIPT 2 TERMINÉ"
echo ""
echo "======================================================================="
echo "  Résultats dans : ${QDIR}/"
echo "  → visual/taxa-bar-plots-decontam.qzv    barplots taxonomiques"
echo "  → visual/Emperor-*.qzv                  PCoA interactifs"
echo "  → visual/alpha-*-significance.qzv        tests diversité alpha"
echo "  → export/subtables/RarTable/ASV_taxonomy.tsv  table finale"
echo "======================================================================="
