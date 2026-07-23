from __future__ import annotations

import csv
import html
import json
import re
import subprocess
import time
import urllib.parse
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "manuscript" / "reference_audit"
NA = "NOT_APPLICABLE"
UNRESOLVED = "UNRESOLVED"


PUBLICATIONS = {
    "NovellaRausell2023MouseKidneyAtlas": {
        "title": "A comprehensive mouse kidney atlas enables rare cell population characterization and robust marker discovery",
        "first_author": "Novella-Rausell",
        "year": "2023",
        "journal": "iScience",
        "doi": "10.1016/j.isci.2023.106877",
    },
    "Kleshchevnikov2022Cell2location": {
        "title": "Cell2location maps fine-grained cell types in spatial transcriptomics",
        "first_author": "Kleshchevnikov",
        "year": "2022",
        "journal": "Nature Biotechnology",
        "doi": "10.1038/s41587-021-01139-4",
    },
    "TabulaSapiens2022Atlas": {
        "title": "The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans",
        "first_author": "Tabula Sapiens Consortium",
        "year": "2022",
        "journal": "Science",
        "doi": "10.1126/science.abl4896",
    },
    "Chen2022hECA": {
        "title": "hECA: The cell-centric assembly of a cell atlas",
        "first_author": "Chen",
        "year": "2022",
        "journal": "iScience",
        "doi": "10.1016/j.isci.2022.104318",
    },
    "GarciaAlonso2021Endometrium": {
        "title": "Mapping the temporal and spatial dynamics of the human endometrium in vivo and in vitro",
        "first_author": "Garcia-Alonso",
        "year": "2021",
        "journal": "Nature Genetics",
        "doi": "10.1038/s41588-021-00972-2",
    },
    "VentoTormo2018MaternalFetal": {
        "title": "Single-cell reconstruction of the early maternal-fetal interface in humans",
        "first_author": "Vento-Tormo",
        "year": "2018",
        "journal": "Nature",
        "doi": "10.1038/s41586-018-0698-6",
    },
    "Han2020HumanCellLandscape": {
        "title": "Construction of a human cell landscape at single-cell level",
        "first_author": "Han",
        "year": "2020",
        "journal": "Nature",
        "doi": "10.1038/s41586-020-2157-4",
    },
    "Reichart2022Cardiomyopathies": {
        "title": "Pathogenic variants damage cell composition and single cell transcription in cardiomyopathies",
        "first_author": "Reichart",
        "year": "2022",
        "journal": "Science",
        "doi": "10.1126/science.abo1984",
    },
    "Elmentaite2021IntestinalTract": {
        "title": "Cells of the human intestinal tract mapped across space and time",
        "first_author": "Elmentaite",
        "year": "2021",
        "journal": "Nature",
        "doi": "10.1038/s41586-021-03852-1",
    },
    "He2020AdultHumanAtlas": {
        "title": "Single-cell transcriptome profiling of an adult human cell atlas of 15 major organs",
        "first_author": "He",
        "year": "2020",
        "journal": "Genome Biology",
        "doi": "10.1186/s13059-020-02210-0",
    },
    "Suo2022DevelopingImmuneSystem": {
        "title": "Mapping the developing human immune system across organs",
        "first_author": "Suo",
        "year": "2022",
        "journal": "Science",
        "doi": "10.1126/science.abo0510",
    },
    "PijuanSala2019MouseGastrulation": {
        "title": "A single-cell molecular map of mouse gastrulation and early organogenesis",
        "first_author": "Pijuan-Sala",
        "year": "2019",
        "journal": "Nature",
        "doi": "10.1038/s41586-019-0933-9",
    },
    "Mittnenzweig2021MouseGastrulation": {
        "title": "A single-embryo, single-cell time-resolved model for mouse gastrulation",
        "first_author": "Mittnenzweig",
        "year": "2021",
        "journal": "Cell",
        "doi": "10.1016/j.cell.2021.04.004",
    },
    "Wu2021BreastCancerAtlas": {
        "title": "A single-cell and spatially resolved atlas of human breast cancers",
        "first_author": "Wu",
        "year": "2021",
        "journal": "Nature Genetics",
        "doi": "10.1038/s41588-021-00911-1",
    },
    "Lee2020ColorectalCancer": {
        "title": "Lineage-dependent gene expression programs influence the immune landscape of colorectal cancer",
        "first_author": "Lee",
        "year": "2020",
        "journal": "Nature Genetics",
        "doi": "10.1038/s41588-020-0636-z",
    },
    "Vahid2023CytoSPACE": {
        "title": "High-resolution alignment of single-cell and spatial transcriptomes with CytoSPACE",
        "first_author": "Vahid",
        "year": "2023",
        "journal": "Nature Biotechnology",
        "doi": "10.1038/s41587-023-01697-9",
    },
}


def joined_publication_fields(keys: list[str]) -> dict[str, str]:
    pubs = [PUBLICATIONS[key] for key in keys]
    return {
        "title": "; ".join(pub["title"] for pub in pubs),
        "first_author": "; ".join(pub["first_author"] for pub in pubs),
        "year": "; ".join(pub["year"] for pub in pubs),
        "journal": "; ".join(pub["journal"] for pub in pubs),
        "doi": "; ".join(pub["doi"] for pub in pubs),
    }


def vendor_spatial(
    *,
    tissue: str,
    disease: str,
    platform: str,
    assay: str,
    processing: str,
    name: str,
    sample: str,
    slide: str,
    area: str,
    url: str,
) -> dict[str, object]:
    return {
        "tissue": tissue,
        "disease": disease,
        "platform": platform,
        "assay": assay,
        "processing": processing,
        "name": name,
        "sample": sample,
        "slide": slide,
        "area": area,
        "database": "10x Genomics public datasets",
        "accession": name,
        "url": url,
        "paper_keys": [],
        "bibtex_key": "NOT_APPLICABLE_VENDOR_DATASET",
    }


SPATIAL = {
    "kidney": vendor_spatial(
        tissue="adult mouse kidney",
        disease="healthy",
        platform="10x Genomics Visium Spatial Gene Expression",
        assay="FFPE whole-transcriptome probe-based spatial expression",
        processing="FFPE direct placement",
        name="Visium_FFPE_Mouse_Kidney",
        sample="male C57BL/6 mouse, older than 8 weeks, sagittal kidney section",
        slide="V11A13-021",
        area="C1",
        url="https://www.10xgenomics.com/datasets/adult-mouse-kidney-ffpe-1-standard-1-3-0",
    ),
    "brain_ffpe": vendor_spatial(
        tissue="adult mouse sagittal brain",
        disease="healthy",
        platform="10x Genomics Visium CytAssist",
        assay="FFPE whole-transcriptome probe-based spatial expression",
        processing="FFPE CytAssist, 11 mm capture area",
        name="CytAssist_FFPE_Sagittal_Mouse_Brain",
        sample="adult mouse sagittal brain FFPE section",
        slide="V52B25-081",
        area="B",
        url="https://www.10xgenomics.com/datasets/preservation-method-comparison-on-visium-cytassist-ffpe-mouse-brain-sagittal-11-mm-capture-area-2-standard",
    ),
    "breast_ffpe": vendor_spatial(
        tissue="human breast",
        disease="ductal carcinoma in situ and invasive carcinoma, grade II",
        platform="10x Genomics Visium Spatial Gene Expression",
        assay="FFPE whole-transcriptome probe-based spatial expression",
        processing="FFPE direct placement",
        name="Visium_FFPE_Human_Breast_Cancer",
        sample="Block 738811QB, Section 1",
        slide="V11J26-008",
        area="B1",
        url="https://www.10xgenomics.com/datasets/human-breast-cancer-ductal-carcinoma-in-situ-invasive-carcinoma-ffpe-1-standard-1-3-0",
    ),
    "breast_ff": vendor_spatial(
        tissue="human breast",
        disease="invasive ductal carcinoma, T2N0M0",
        platform="10x Genomics Visium Spatial Gene Expression",
        assay="fresh-frozen whole-transcriptome spatial expression",
        processing="fresh frozen",
        name="Visium_Human_Breast_Cancer",
        sample="Block 1168993F",
        slide="V19B23-014",
        area="A1",
        url="https://www.10xgenomics.com/datasets/human-breast-cancer-visium-fresh-frozen-whole-transcriptome-1-standard",
    ),
    "breast_wta120": vendor_spatial(
        tissue="human breast",
        disease="invasive lobular carcinoma",
        platform="10x Genomics Visium Spatial Gene Expression",
        assay="fresh-frozen whole-transcriptome spatial expression",
        processing="fresh frozen; Space Ranger 1.2.0 release",
        name="Parent_Visium_Human_BreastCancer",
        sample="human invasive lobular carcinoma section",
        slide="V19L29-095",
        area="A1",
        url="https://www.10xgenomics.com/datasets/human-breast-cancer-whole-transcriptome-analysis-1-standard-1-2-0",
    ),
    "cervix": vendor_spatial(
        tissue="human cervix",
        disease="squamous cell carcinoma, T1bN0M0 stage IB",
        platform="10x Genomics Visium Spatial Gene Expression",
        assay="FFPE whole-transcriptome probe-based spatial expression",
        processing="FFPE direct placement",
        name="Visium_FFPE_Human_Cervical_Cancer",
        sample="Block C00084155.1a",
        slide="V10L13-019",
        area="A1",
        url="https://www.10xgenomics.com/datasets/human-cervical-cancer-1-standard",
    ),
    "heart": vendor_spatial(
        tissue="human heart ventricle",
        disease="healthy",
        platform="10x Genomics Visium Spatial Gene Expression",
        assay="fresh-frozen whole-transcriptome spatial expression",
        processing="fresh frozen",
        name="V1_Human_Heart",
        sample="human heart section",
        slide="V19S16-046",
        area="C1",
        url="https://www.10xgenomics.com/datasets/human-heart-1-standard-1-1-0",
    ),
    "intestine": vendor_spatial(
        tissue="human large intestine",
        disease="colorectal cancer",
        platform="10x Genomics Visium Spatial Gene Expression",
        assay="FFPE whole-transcriptome probe-based spatial expression",
        processing="FFPE direct placement",
        name="Visium_FFPE_Human_Intestinal_Cancer",
        sample="Block 1281585B",
        slide="V10L13-021",
        area="B1",
        url="https://www.10xgenomics.com/datasets/human-intestine-cancer-1-standard",
    ),
    "lymph": vendor_spatial(
        tissue="human lymph node",
        disease="healthy",
        platform="10x Genomics Visium Spatial Gene Expression",
        assay="fresh-frozen whole-transcriptome spatial expression",
        processing="fresh frozen",
        name="V1_Human_Lymph_Node",
        sample="human lymph node section",
        slide="V19L01-033",
        area="A1",
        url="https://www.10xgenomics.com/datasets/human-lymph-node-1-standard-1-0-0",
    ),
    "embryo": vendor_spatial(
        tissue="mouse embryo",
        disease="developmental tissue",
        platform="10x Genomics Visium CytAssist",
        assay="FFPE whole-transcriptome probe-based spatial expression",
        processing="5 micrometre FFPE section, 11 mm capture area",
        name="Visium_CytAssist_Mouse_Embryo_11mm_FFPE",
        sample="mouse embryo; developmental stage not stated by the public record",
        slide="V52Y09-019",
        area="B",
        url="https://www.10xgenomics.com/datasets/visium-cytassist-mouse-embryo-11-mm-capture-area-ffpe-2-standard",
    ),
    "tnbc": {
        "tissue": "human breast",
        "disease": "triple-negative breast cancer",
        "platform": "10x Genomics Visium Spatial Gene Expression",
        "assay": "fresh-frozen whole-transcriptome spatial expression",
        "processing": "fresh frozen",
        "name": "Wu breast-cancer spatial section CID4465",
        "sample": "CID4465",
        "slide": UNRESOLVED,
        "area": UNRESOLVED,
        "database": "Zenodo",
        "accession": "4739739",
        "url": "https://zenodo.org/records/4739739",
        "paper_keys": ["Wu2021BreastCancerAtlas"],
        "bibtex_key": "Wu2021BreastCancerAtlas",
    },
    "crc": vendor_spatial(
        tissue="human colorectum",
        disease="invasive adenocarcinoma, T4aN0M0 stage IIB",
        platform="10x Genomics Visium Spatial Gene Expression",
        assay="fresh-frozen whole-transcriptome spatial expression",
        processing="fresh frozen; Space Ranger 1.2.0 release",
        name="Parent_Visium_Human_ColorectalCancer",
        sample="human colorectal cancer section",
        slide=UNRESOLVED,
        area="C1",
        url="https://www.10xgenomics.com/datasets/human-colorectal-cancer-whole-transcriptome-analysis-1-standard-1-2-0",
    ),
}


REFERENCE = {
    "kidney_atlas": {
        "species": "Mus musculus",
        "tissue": "adult kidney",
        "assay": "integrated scRNA-seq atlas",
        "name": "Comprehensive Mouse Kidney Atlas",
        "cohort": "141,401 cells from eight locally retained study origins",
        "database": "CZ CELLxGENE plus source-study repositories",
        "accession": "CELLxGENE 42bb7f78-cef8-4b0d-9bba-50037d64d8c1; component accessions documented by Novella-Rausell et al.",
        "url": "https://cellxgene.cziscience.com/e/42bb7f78-cef8-4b0d-9bba-50037d64d8c1.cxg/",
        "paper_keys": ["NovellaRausell2023MouseKidneyAtlas"],
        "bibtex_key": "NovellaRausell2023MouseKidneyAtlas",
    },
    "brain_sn": {
        "species": "Mus musculus",
        "tissue": "adult brain",
        "assay": "single-nucleus RNA-seq",
        "name": "Cell2location paired adult mouse-brain snRNA-seq reference",
        "cohort": "40,572 nuclei; local sample identifiers 5705STDY8058280-5705STDY8058285",
        "database": "BioStudies / ArrayExpress",
        "accession": "E-MTAB-11115",
        "url": "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11115",
        "paper_keys": ["Kleshchevnikov2022Cell2location"],
        "bibtex_key": "Kleshchevnikov2022Cell2location",
    },
    "heca_breast": {
        "species": "Homo sapiens",
        "tissue": "healthy breast",
        "assay": "scRNA-seq and Smart-seq2 hECA export",
        "name": "hECA Tabula Sapiens breast project export",
        "cohort": "11,227 cells; donor TSP4; local study DOI 10.1101/2021.07.19.452956",
        "database": "hECA / Zenodo",
        "accession": "Zenodo 17008296",
        "url": "https://zenodo.org/records/17008296",
        "paper_keys": ["TabulaSapiens2022Atlas"],
        "bibtex_key": "TabulaSapiens2022Atlas;Chen2022hECA",
    },
    "heca_uterus": {
        "species": "Homo sapiens",
        "tissue": "healthy uterus and endometrium",
        "assay": "combined scRNA-seq hECA export",
        "name": "hECA healthy uterus composite reference",
        "cohort": "111,954 cells from Garcia-Alonso, Vento-Tormo, Han and Tabula Sapiens studies",
        "database": "BioStudies / hECA / Zenodo",
        "accession": "E-MTAB-10287; Zenodo 17008269; Zenodo 17008296",
        "url": "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10287; https://zenodo.org/records/17008269; https://zenodo.org/records/17008296",
        "paper_keys": ["GarciaAlonso2021Endometrium", "VentoTormo2018MaternalFetal", "Han2020HumanCellLandscape", "TabulaSapiens2022Atlas"],
        "bibtex_key": "GarciaAlonso2021Endometrium;VentoTormo2018MaternalFetal;Han2020HumanCellLandscape;TabulaSapiens2022Atlas;Chen2022hECA",
    },
    "heart_ref": {
        "species": "Homo sapiens",
        "tissue": "healthy ventricular heart",
        "assay": "single-cell and single-nucleus RNA-seq project export",
        "name": "Reichart et al. healthy ventricular reference",
        "cohort": "67,246 healthy ventricular cells/nuclei from six donors",
        "database": "European Genome-phenome Archive / hECA Zenodo",
        "accession": "EGAS00001006374; Zenodo 17010276",
        "url": "https://ega-archive.org/studies/EGAS00001006374; https://zenodo.org/records/17010276",
        "paper_keys": ["Reichart2022Cardiomyopathies"],
        "bibtex_key": "Reichart2022Cardiomyopathies;Chen2022hECA",
    },
    "intestine_ref": {
        "species": "Homo sapiens",
        "tissue": "healthy intestine",
        "assay": "combined scRNA-seq hECA export",
        "name": "hECA healthy intestine composite reference",
        "cohort": "465,681 cells from Elmentaite, Han and He studies",
        "database": "BioStudies / hECA / Zenodo",
        "accession": "E-MTAB-9543; E-MTAB-9536; Zenodo 17008296; Zenodo 17010276",
        "url": "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-9543; https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-9536; https://zenodo.org/records/17008296; https://zenodo.org/records/17010276",
        "paper_keys": ["Elmentaite2021IntestinalTract", "Han2020HumanCellLandscape", "He2020AdultHumanAtlas"],
        "bibtex_key": "Elmentaite2021IntestinalTract;Han2020HumanCellLandscape;He2020AdultHumanAtlas;Chen2022hECA",
    },
    "lymph_ref": {
        "species": "Homo sapiens",
        "tissue": "healthy lymph node",
        "assay": "scRNA-seq hECA project export",
        "name": "Suo et al. developing human immune-system lymph-node reference",
        "cohort": "5,872 cells; local donors F78 and F72",
        "database": "BioStudies / ArrayExpress and hECA Zenodo",
        "accession": "E-MTAB-11343; Zenodo 17010276",
        "url": "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11343; https://zenodo.org/records/17010276",
        "paper_keys": ["Suo2022DevelopingImmuneSystem"],
        "bibtex_key": "Suo2022DevelopingImmuneSystem;Chen2022hECA",
    },
    "embryo_ref": {
        "species": "Mus musculus",
        "tissue": "gastrulating and early-organogenesis embryo",
        "assay": "combined scRNA-seq atlas",
        "name": "Gottgens 2019 and Stelzer 2021 mouse-gastrulation reference",
        "cohort": "235,715 cells: 139,331 Pijuan-Sala/Gottgens plus 96,384 Mittnenzweig/Stelzer",
        "database": "BioStudies / ArrayExpress and NCBI GEO",
        "accession": "E-MTAB-6967; GSE169210",
        "url": "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6967; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169210",
        "paper_keys": ["PijuanSala2019MouseGastrulation", "Mittnenzweig2021MouseGastrulation"],
        "bibtex_key": "PijuanSala2019MouseGastrulation;Mittnenzweig2021MouseGastrulation",
    },
    "wu_breast": {
        "species": "Homo sapiens",
        "tissue": "breast cancer",
        "assay": "10x Genomics scRNA-seq",
        "name": "Wu et al. single-cell atlas of human breast cancers",
        "cohort": "26 breast tumours; project preparation selects relevant TNBC or HER2-positive cells",
        "database": "NCBI GEO",
        "accession": "GSE176078",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078",
        "paper_keys": ["Wu2021BreastCancerAtlas"],
        "bibtex_key": "Wu2021BreastCancerAtlas",
    },
    "lee_crc": {
        "species": "Homo sapiens",
        "tissue": "colorectal cancer",
        "assay": "10x Genomics scRNA-seq",
        "name": "Lee et al. colorectal-cancer single-cell atlas",
        "cohort": "colorectal cancer tumour and immune-cell atlas used by the CytoSPACE-derived preparation",
        "database": "NCBI GEO",
        "accession": "GSE132465",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465",
        "paper_keys": ["Lee2020ColorectalCancer"],
        "bibtex_key": "Lee2020ColorectalCancer",
    },
}


FIG2_SOURCE = (
    "result/real_profile_mask_foundation/foundation_manifest.csv; "
    "result/real_profile_mask_expression_recovery/expression_recovery_by_scenario.csv; "
    "visualizations/simulations/real_profile_mask_fig2d_only/fig2_panel_d_real_profile_mask_source_values.csv"
)
FIG2_PLOTS = (
    "scripts/build_masked_scenarios_real_stack.py; "
    "scripts/plot_real_profile_mask_fig2d_only.py; "
    "scripts/plot_real_profile_mask_expression_recovery.py"
)
FIG4_SOURCE = (
    "visualizations/stage3b_realdata_candidate_scan/spatial_9x2/"
    "stage3b_reference_dropout_spatial_stack_recommended_6x2_manifest.csv; "
    "visualizations/stage3b_realdata_candidate_scan/panel_a_blank_composition/"
    "stage3b_reference_dropout_panel_a_blank_region_composition_summary.csv"
)
FIG4_PLOTS = (
    "scripts/plot_stage3b_reference_dropout_spatial_stack.py; "
    "scripts/plot_stage3b_reference_dropout_panel_a_blank_composition.py"
)


SCENARIOS = [
    ("Fig. 2 low-resolution real profile masking", "Fig. 2", "A/B/E", 1, "adult_mouse_kidney_real_profile_mask_endo", "Endo", "profile masking", "kidney", "kidney_atlas", "HIGH"),
    ("Fig. 2 low-resolution real profile masking", "Fig. 2", "A/B/E", 2, "ffpe_mouse_brain_sagittal_real_profile_mask_microglia", "Microglia", "profile masking", "brain_ffpe", "brain_sn", "HIGH"),
    ("Fig. 2 low-resolution real profile masking", "Fig. 2", "A/B/E", 3, "human_breast_cancer_real_profile_mask_basal_cell", "Basal cell", "profile masking", "breast_ffpe", "heca_breast", "MEDIUM"),
    ("Fig. 2 low-resolution real profile masking", "Fig. 2", "A/B/E", 4, "human_breast_cancer_visium_ff_wta_real_profile_mask_macrophage", "Macrophage", "profile masking", "breast_ff", "heca_breast", "MEDIUM"),
    ("Fig. 2 low-resolution real profile masking", "Fig. 2", "A/B/E", 5, "human_breast_cancer_wta_120_real_profile_mask_endothelial_cell", "Endothelial cell", "profile masking", "breast_wta120", "heca_breast", "MEDIUM"),
    ("Fig. 2 low-resolution real profile masking", "Fig. 2", "A/B/E", 6, "human_cervical_cancer_real_profile_mask_epithelial_cell", "Epithelial cell", "profile masking", "cervix", "heca_uterus", "MEDIUM"),
    ("Fig. 2 low-resolution real profile masking", "Fig. 2", "A/B/E", 7, "human_heart_ff_real_profile_mask_endothelial_cell", "Endothelial cell", "profile masking", "heart", "heart_ref", "HIGH"),
    ("Fig. 2 low-resolution real profile masking", "Fig. 2", "A/B/E", 8, "human_intestine_cancer_real_profile_mask_endothelial_cell", "Endothelial cell", "profile masking", "intestine", "intestine_ref", "MEDIUM"),
    ("Fig. 2 low-resolution real profile masking", "Fig. 2", "A/B/E", 9, "human_lymph_node_real_profile_mask_b_cell", "B cell", "profile masking", "lymph", "lymph_ref", "HIGH"),
    ("Fig. 2 low-resolution real profile masking", "Fig. 2", "A/B/E", 10, "mouse_embryo_real_profile_mask_erythroid", "Erythroid", "profile masking", "embryo", "embryo_ref", "HIGH"),
    ("Fig. 4 real reference dropout", "Fig. 4", "F/G", 1, "mouse_embryo_real_sc_missing_endoderm_gut", "Endoderm/Gut", "single-cell-reference dropout", "embryo", "embryo_ref", "HIGH"),
    ("Fig. 4 real reference dropout", "Fig. 4", "F/G", 2, "mouse_embryo_real_sc_missing_erythroid", "Erythroid", "single-cell-reference dropout", "embryo", "embryo_ref", "HIGH"),
    ("Fig. 4 real reference dropout", "Fig. 4", "F/G", 3, "cytospace_fig2d_tme_brca_tnbc_fresh_frozen_sc_missing_plasma_cells", "Plasma cells", "single-cell-reference dropout", "tnbc", "wu_breast", "HIGH"),
    ("Fig. 4 real reference dropout", "Fig. 4", "F/G", 4, "cytospace_fig2d_tme_crc_fresh_frozen_sc_missing_b_cells", "B cells", "single-cell-reference dropout", "crc", "lee_crc", "HIGH"),
    ("Fig. 4 real reference dropout", "Fig. 4", "F/G", 5, "cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_plasma_cells", "Plasma cells", "single-cell-reference dropout", "breast_ffpe", "wu_breast", "HIGH"),
    ("Fig. 4 real reference dropout", "Fig. 4", "F/G", 6, "cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_epithelial_cells", "Epithelial cells", "single-cell-reference dropout", "breast_ffpe", "wu_breast", "HIGH"),
]


CSV_FIELDS = [
    "experiment_family", "figure", "figure_panel", "final_scenario_order", "internal_scenario_id",
    "config_path", "final_source_table_path", "final_plot_script_path", "target_population",
    "perturbation_type", "spatial_tissue", "spatial_disease", "spatial_platform", "spatial_assay",
    "spatial_processing_type", "spatial_formal_dataset_name", "spatial_sample_id", "spatial_slide_id",
    "spatial_capture_area", "spatial_accession_database", "spatial_accession", "spatial_public_url",
    "spatial_original_paper_title", "spatial_original_paper_first_author", "spatial_original_paper_year",
    "spatial_original_paper_journal", "spatial_original_paper_doi", "reference_species", "reference_tissue",
    "reference_assay", "reference_formal_dataset_name", "reference_sample_or_cohort",
    "reference_accession_database", "reference_accession", "reference_public_url",
    "reference_original_paper_title", "reference_original_paper_first_author", "reference_original_paper_year",
    "reference_original_paper_journal", "reference_original_paper_doi", "derived_from_cytospace_resource",
    "cytospace_source_note", "recommended_spatial_bibtex_key", "recommended_reference_bibtex_key",
    "confidence", "unresolved_fields", "repository_evidence", "public_source_evidence", "notes",
]


def publication_csv_fields(prefix: str, keys: list[str]) -> dict[str, str]:
    if not keys:
        return {
            f"{prefix}_original_paper_title": NA,
            f"{prefix}_original_paper_first_author": NA,
            f"{prefix}_original_paper_year": NA,
            f"{prefix}_original_paper_journal": NA,
            f"{prefix}_original_paper_doi": NA,
        }
    pubs = joined_publication_fields(keys)
    return {f"{prefix}_original_paper_{field}": value for field, value in pubs.items()}


def build_rows() -> tuple[list[dict[str, str]], list[dict[str, object]]]:
    rows: list[dict[str, str]] = []
    nested: list[dict[str, object]] = []
    for family, figure, panel, order, scenario_id, target, perturbation, spatial_key, reference_key, confidence in SCENARIOS:
        spatial = SPATIAL[spatial_key]
        reference = REFERENCE[reference_key]
        is_fig2 = figure == "Fig. 2"
        config = (
            "configs/datasets/mouse_embryo_real.yaml"
            if scenario_id == "mouse_embryo_real_profile_mask_erythroid"
            else f"configs/datasets/{scenario_id}.yaml"
        )
        source_table = FIG2_SOURCE if is_fig2 else FIG4_SOURCE
        plot_scripts = FIG2_PLOTS if is_fig2 else FIG4_PLOTS
        derived = "true" if scenario_id.startswith("cytospace_fig2d_tme_") else "false"
        cyto_note = (
            "The project preparation follows the CytoSPACE tumour-resource organization, but the cited data sources are the original Wu/Lee/10x records; CytoSPACE is a method/resource citation, not the original data publication."
            if derived == "true"
            else NA
        )
        unresolved = []
        for field_name, value in (("spatial_slide_id", spatial["slide"]), ("spatial_capture_area", spatial["area"])):
            if value == UNRESOLVED:
                unresolved.append(field_name)
        if confidence == "MEDIUM":
            unresolved.append("exact historical local hECA export checksum and merge manifest")
        repository = [config, source_table, plot_scripts]
        if is_fig2:
            base = scenario_id.removesuffix("_profile_mask_endo").removesuffix("_profile_mask_microglia").removesuffix("_profile_mask_basal_cell").removesuffix("_profile_mask_macrophage").removesuffix("_profile_mask_endothelial_cell").removesuffix("_profile_mask_epithelial_cell").removesuffix("_profile_mask_b_cell").removesuffix("_profile_mask_erythroid")
            repository.append(f"data/raw/low_resolution_experiments/{base}/real_input_info.json")
        else:
            repository.append("scripts/prepare_stage3b_realdata_reference_dropout_scenarios.py")
            if scenario_id.startswith("cytospace_fig2d_tme_"):
                repository.append("scripts/prepare_cytospace_fig2d_tme_stage1.py")
        public_urls = [part.strip() for part in str(spatial["url"]).split(";")]
        public_urls.extend(part.strip() for part in str(reference["url"]).split(";"))
        public_urls.extend(f"https://doi.org/{PUBLICATIONS[key]['doi']}" for key in spatial["paper_keys"])
        public_urls.extend(f"https://doi.org/{PUBLICATIONS[key]['doi']}" for key in reference["paper_keys"])
        resource_publication_keys: list[str] = []
        if "Chen2022hECA" in str(reference["bibtex_key"]):
            resource_publication_keys.append("Chen2022hECA")
        if derived == "true":
            resource_publication_keys.append("Vahid2023CytoSPACE")
        public_urls.extend(
            f"https://doi.org/{PUBLICATIONS[key]['doi']}" for key in resource_publication_keys
        )
        row = {
            "experiment_family": family,
            "figure": figure,
            "figure_panel": panel,
            "final_scenario_order": str(order),
            "internal_scenario_id": scenario_id,
            "config_path": config,
            "final_source_table_path": source_table,
            "final_plot_script_path": plot_scripts,
            "target_population": target,
            "perturbation_type": perturbation,
            "spatial_tissue": str(spatial["tissue"]),
            "spatial_disease": str(spatial["disease"]),
            "spatial_platform": str(spatial["platform"]),
            "spatial_assay": str(spatial["assay"]),
            "spatial_processing_type": str(spatial["processing"]),
            "spatial_formal_dataset_name": str(spatial["name"]),
            "spatial_sample_id": str(spatial["sample"]),
            "spatial_slide_id": str(spatial["slide"]),
            "spatial_capture_area": str(spatial["area"]),
            "spatial_accession_database": str(spatial["database"]),
            "spatial_accession": str(spatial["accession"]),
            "spatial_public_url": str(spatial["url"]),
            "reference_species": str(reference["species"]),
            "reference_tissue": str(reference["tissue"]),
            "reference_assay": str(reference["assay"]),
            "reference_formal_dataset_name": str(reference["name"]),
            "reference_sample_or_cohort": str(reference["cohort"]),
            "reference_accession_database": str(reference["database"]),
            "reference_accession": str(reference["accession"]),
            "reference_public_url": str(reference["url"]),
            "derived_from_cytospace_resource": derived,
            "cytospace_source_note": cyto_note,
            "recommended_spatial_bibtex_key": str(spatial["bibtex_key"]),
            "recommended_reference_bibtex_key": str(reference["bibtex_key"]),
            "confidence": confidence,
            "unresolved_fields": "; ".join(unresolved) if unresolved else "none",
            "repository_evidence": "; ".join(repository),
            "public_source_evidence": "; ".join(dict.fromkeys(public_urls)),
            "notes": (
                "The formal scenario is fixed by the final manifest/source table and plot order. No experiment was rerun. "
                + ("The reference is a documented multi-study composite; component citations must remain grouped." if len(reference["paper_keys"]) > 1 else "The repository identity and public record agree.")
            ),
        }
        row.update(publication_csv_fields("spatial", list(spatial["paper_keys"])))
        row.update(publication_csv_fields("reference", list(reference["paper_keys"])))
        missing_columns = [column for column in CSV_FIELDS if column not in row]
        if missing_columns:
            raise RuntimeError(f"Missing CSV fields for {scenario_id}: {missing_columns}")
        blank_columns = [column for column in CSV_FIELDS if not str(row[column]).strip()]
        if blank_columns:
            raise RuntimeError(f"Blank CSV fields for {scenario_id}: {blank_columns}")
        rows.append({column: row[column] for column in CSV_FIELDS})
        publication_keys = list(
            dict.fromkeys(
                list(spatial["paper_keys"])
                + list(reference["paper_keys"])
                + resource_publication_keys
            )
        )
        nested.append(
            {
                "scenario": {column: row[column] for column in CSV_FIELDS if column not in {"repository_evidence", "public_source_evidence"}},
                "evidence": {
                    "repository_evidence": repository,
                    "public_database_evidence": list(dict.fromkeys(public_urls)),
                    "publication_evidence": [dict(key=key, **PUBLICATIONS[key]) for key in publication_keys],
                },
                "confidence": confidence,
            }
        )
    return rows, nested


THALAMIC = {
    "role": "separate detailed case; excluded from the six Fig. 4 cross-dataset scenarios",
    "internal_scenario_id": "cell2loc_scan_st8059051_sc_missing_thalamic_excitatory",
    "figure_panels": "Fig. 4A-E and Fig. 4H",
    "target_population": "Ext_Thal_1 and Ext_Thal_2 (Thalamic excitatory)",
    "config_path": "configs/datasets/cell2loc_scan_st8059051_sc_missing_thalamic_excitatory.yaml",
    "preparation_script": "scripts/prepare_cell2location_mouse_brain_stage3b_case.py",
    "repository_case_record": "data/processed/cell2location_mouse_brain/cell2loc_scan_st8059051_sc_missing_thalamic_excitatory/stage1_preprocess/stage3b_case_info.json",
    "spatial_dataset": {
        "formal_name": "Cell2location adult mouse-brain 10x Visium spatial series",
        "sample_id": "ST8059051 = Visium-29B",
        "slide_id": "C05717-021",
        "capture_area": "B1",
        "accession": "E-MTAB-11114",
        "public_url": "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11114",
    },
    "reference_dataset": {
        "formal_name": "Cell2location paired adult mouse-brain snRNA-seq reference",
        "accession": "E-MTAB-11115",
        "public_url": "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11115",
    },
    "publication": dict(key="Kleshchevnikov2022Cell2location", **PUBLICATIONS["Kleshchevnikov2022Cell2location"]),
    "confidence": "HIGH",
    "interpretation_note": "Top-15% regions are reference-derived marker-region proxies, not independent biological ground truth.",
}


def scenario_table(rows: list[dict[str, str]], figure: str) -> str:
    subset = [row for row in rows if row["figure"] == figure]
    lines = [
        "| Order | Internal scenario ID | Target | Spatial source | Reference source | Confidence |",
        "|---:|---|---|---|---|---|",
    ]
    for row in subset:
        lines.append(
            f"| {row['final_scenario_order']} | `{row['internal_scenario_id']}` | {row['target_population']} | "
            f"{row['spatial_formal_dataset_name']} | {row['reference_formal_dataset_name']} | {row['confidence']} |"
        )
    return "\n".join(lines)


def detailed_scenarios(rows: list[dict[str, str]], prefix: str) -> str:
    lines: list[str] = []
    for row in rows:
        lines.extend(
            [
                f"### {prefix}{row['final_scenario_order']}. `{row['internal_scenario_id']}`",
                "",
                f"- **Perturbation and target:** {row['perturbation_type']}; `{row['target_population']}`.",
                f"- **Spatial evidence:** `{row['spatial_formal_dataset_name']}`; sample `{row['spatial_sample_id']}`; "
                f"slide `{row['spatial_slide_id']}`; capture area `{row['spatial_capture_area']}`; "
                f"public record: {row['spatial_public_url']}.",
                f"- **Reference evidence:** {row['reference_formal_dataset_name']}; {row['reference_sample_or_cohort']}; "
                f"accession(s): `{row['reference_accession']}`; public record(s): {row['reference_public_url']}.",
                f"- **Repository lock:** `{row['config_path']}`; `{row['final_source_table_path']}`; "
                f"`{row['final_plot_script_path']}`.",
                f"- **Publication mapping:** spatial DOI(s) `{row['spatial_original_paper_doi']}`; "
                f"reference DOI(s) `{row['reference_original_paper_doi']}`.",
                f"- **Confidence:** {row['confidence']}. Unresolved fields: {row['unresolved_fields']}.",
                "",
            ]
        )
    return "\n".join(lines)


def build_report(rows: list[dict[str, str]]) -> str:
    fig2 = [row for row in rows if row["figure"] == "Fig. 2"]
    fig4 = [row for row in rows if row["figure"] == "Fig. 4"]
    confidence = Counter(row["confidence"] for row in rows)
    paper_usage: Counter[str] = Counter()
    for _, _, _, _, _, _, _, spatial_key, reference_key, _ in SCENARIOS:
        for key in set(list(SPATIAL[spatial_key]["paper_keys"]) + list(REFERENCE[reference_key]["paper_keys"])):
            paper_usage[key] += 1
    shared = [(key, count) for key, count in paper_usage.items() if count > 1]
    publication_lines = [
        "| Key | Publication | DOI | Primary data role |",
        "|---|---|---|---|",
    ]
    source_keys = [key for key in PUBLICATIONS if key not in {"Chen2022hECA", "Vahid2023CytoSPACE"}]
    for key in source_keys:
        pub = PUBLICATIONS[key]
        publication_lines.append(f"| `{key}` | {pub['title']} ({pub['journal']}, {pub['year']}) | `{pub['doi']}` | Original source study |")
    shared_lines = [f"- `{key}` supports {count} of the 16 formal scenarios." for key, count in sorted(shared)]
    report = f"""# SVTuner Fig. 2 and Fig. 4 public-data source audit

## 1. Audit scope

This audit fixes and traces the 10 final low-resolution real profile-masking scenarios used by Fig. 2 and the 6 final real reference-dropout scenarios used by the Fig. 4 cross-dataset comparison. The thalamic-excitatory mouse-brain analysis is recorded separately and is not counted as a seventh cross-dataset scenario.

## 2. Frozen rules and non-actions

No Stage1, Stage3A, Stage3B, Stage4, mapping method, experiment, statistic or figure was rerun. Existing results, frozen configurations, source-value tables, formal manuscript files and the formal bibliography were not modified. Public identity was accepted only when repository paths, final source tables/plot order and authoritative public records were mutually consistent.

## 3. How the final 10 Fig. 2 scenarios were identified

The authoritative scenario lock is the 10-row `result/real_profile_mask_foundation/foundation_manifest.csv`, corroborated by `result/real_profile_mask_foundation/foundation_config.json`, `result/real_profile_mask_expression_recovery/expression_recovery_by_scenario.csv`, the final stack `visualizations/masked_scenarios/masked_scenarios_real_stack_10x4.png`, and its explicit ordering in `scripts/build_masked_scenarios_real_stack.py`. The separate 12-readout CytoSPACE-style state-enrichment benchmark is not this 10-scenario set.

## 4. Final Fig. 2 scenario list

{scenario_table(fig2, 'Fig. 2')}

## 5. How the final 6 Fig. 4 scenarios were identified

The authoritative lock is `visualizations/stage3b_realdata_candidate_scan/spatial_9x2/stage3b_reference_dropout_spatial_stack_recommended_6x2_manifest.csv`, including its six ordered rows and strict-validation fields. The same six feed the final spatial stack and blank-region composition summary. Other candidate scans and the thalamic case are not part of this six-scenario aggregate.

## 6. Final Fig. 4 scenario list

{scenario_table(fig4, 'Fig. 4')}

### Separate thalamic detailed case

`{THALAMIC['internal_scenario_id']}` uses ST8059051 (`Visium-29B`, slide `C05717-021`, capture area `B1`) from `E-MTAB-11114` and the paired mouse-brain snRNA-seq reference `E-MTAB-11115`. Both were published with Kleshchevnikov et al. (`10.1038/s41587-021-01139-4`). It supports Fig. 4A-E/H, remains outside the six cross-dataset scenarios, and its marker-region accuracy is reference-relative rather than independent biological ground truth.

## 7. Detailed spatial-data provenance

{detailed_scenarios(fig2, 'F2-')}
{detailed_scenarios(fig4, 'F4-')}

Across the 16 rows there are **12 unique spatial datasets**. Reuse is intentional: the 10x FFPE breast section supports one Fig. 2 scenario and two HER2 Fig. 4 scenarios; the 10x mouse-embryo section supports one Fig. 2 scenario and two Fig. 4 scenarios.

## 8. Detailed single-cell-reference provenance

There are **10 unique reference objects**. Three Fig. 2 breast scenarios share one hECA/Tabula Sapiens breast export. The cervical experiment deliberately pairs a cervical-cancer ST section with a healthy uterus/endometrium composite reference. The embryo profile-mask and both embryo reference-dropout scenarios share the same two-study gastrulation reference. The three breast-cancer reference-dropout scenarios use the Wu atlas, whereas the CRC scenario uses the Lee atlas.

Composite hECA references are cited as composites: the local metadata identifies their component studies, while the current hECA Zenodo project exports provide public resolvers. The exact historical local export checksums and merge manifests were not preserved, producing five MEDIUM rather than HIGH rows; this does not obscure the source studies or alter citation choice.

## 9. Original publication and accession mapping

{chr(10).join(publication_lines)}

The 16 formal scenarios resolve to **1 unique spatial-data original paper** (Wu; all other spatial records are vendor releases without an associated source article) and **14 unique reference-data original papers**. Wu is shared by the TNBC spatial section and breast-cancer references, so the deduplicated union remains **14 unique primary source papers**. The 12 spatial datasets have **12 distinct formal public dataset records/identifiers**. The 10 reference objects resolve through **14 distinct exact public accessions or project records**; component accessions mentioned only generically by an atlas paper are not inflated into this count. The candidate bibliography additionally includes `Chen2022hECA` for the atlas assembly/export layer and `Vahid2023CytoSPACE` for the reused tumour-resource organization, giving **16 candidate BibTeX entries**. Vendor-only 10x data releases are cited by formal dataset identifier and URL in Methods/Data Availability and are not misrepresented as journal articles.

## 10. Shared and duplicated source papers

{chr(10).join(shared_lines)}

`Chen2022hECA` is a shared assembly/resource citation for seven hECA-derived Fig. 2 references but is not counted as an original biological source study. `Vahid2023CytoSPACE` documents the benchmark-resource organization for four Fig. 4 tumour scenarios but does not replace Wu, Lee or 10x source attribution.

## 11. Recommended citation keys

- Kidney reference: `NovellaRausell2023MouseKidneyAtlas`.
- Cell2location mouse brain: `Kleshchevnikov2022Cell2location`.
- hECA-derived references: cite their component source keys plus `Chen2022hECA`.
- Mouse embryo reference: `PijuanSala2019MouseGastrulation` and `Mittnenzweig2021MouseGastrulation`.
- Breast-cancer reference-dropout data: `Wu2021BreastCancerAtlas`.
- Colorectal-cancer reference: `Lee2020ColorectalCancer`.
- CytoSPACE-derived resource organization: `Vahid2023CytoSPACE`, in addition to original data citations.

## 12. Unresolved items

- Five MEDIUM-confidence Fig. 2 rows lack a preserved checksum/merge manifest connecting the historical local hECA export byte-for-byte to the current public project export. Local study DOI, tissue, cohort size and component provenance are nevertheless present.
- The public CID4465 Zenodo record does not expose a 10x slide serial or capture-area code; those two non-core fields remain `UNRESOLVED`.
- The official colorectal data page fixes the formal dataset and capture area but the slide serial was not confirmed; it remains `UNRESOLVED`.
- The 10x mouse-embryo page does not state an embryonic stage. No stage is inferred.
- No scenario has a LOW or UNRESOLVED core ST/reference identity.

## 13. Manuscript implications

The manuscript can safely proceed to citation insertion, provided vendor dataset URLs/identifiers are retained in Methods or Data Availability, hECA references are described as curated composite/exported references, and CytoSPACE is not used as a substitute citation for Wu, Lee or 10x source data. The separate thalamic case must remain explicitly separate from the six-scenario Fig. 4 aggregate. The suggested wording locations are recorded in `dataset_citation_insertion_plan.md`; no manuscript text was changed here.

Counts must remain distinct: **16 scenario rows**, **12 unique spatial datasets**, **12 unique spatial public records/identifiers**, **10 unique reference objects**, **14 unique reference public accessions/project records**, **1 spatial-data original paper**, **14 reference-data original papers**, **14 deduplicated primary source papers**, and **16 candidate BibTeX entries** (including two assembly/method citations).

## 14. Final decision

**PASS.** All 16 formal scenarios are locked by final repository evidence. Every core spatial and SC/snRNA reference identity has an authoritative public source or an explicit composite-source explanation. Confidence distribution: {confidence['HIGH']} HIGH, {confidence['MEDIUM']} MEDIUM, {confidence['LOW']} LOW and {confidence['UNRESOLVED']} UNRESOLVED. The remaining unresolved items are sample-level metadata or historical export checksums and do not change citation identity.
"""
    return report


def build_citation_plan() -> str:
    return """# Dataset citation insertion plan

This file proposes insertion points only. It does not modify the manuscript or formal bibliography.

## Introduction: public platforms and mapping resources

**Current sentence:** First mention of public Visium/CytAssist data and CytoSPACE-organized benchmark resources.

**Recommended replacement:** Retain the platform wording, cite `Vahid2023CytoSPACE` for the CytoSPACE method/resource layer, and state that original biological data publications and vendor records are cited separately.

**Reason:** CytoSPACE is not the original publisher of Wu, Lee or 10x data.

**Applies to:** General data-source framing.

## Results: Fig. 2 dataset description

**Current sentence:** The low-resolution profile-masking benchmark comprised ten tissues/diseases.

**Recommended replacement:** Name the ten formal 10x dataset identifiers in Methods or Supplementary Data and cite the reference-source keys `NovellaRausell2023MouseKidneyAtlas`, `Kleshchevnikov2022Cell2location`, `TabulaSapiens2022Atlas`, `GarciaAlonso2021Endometrium`, `VentoTormo2018MaternalFetal`, `Han2020HumanCellLandscape`, `Reichart2022Cardiomyopathies`, `Elmentaite2021IntestinalTract`, `He2020AdultHumanAtlas`, `Suo2022DevelopingImmuneSystem`, `PijuanSala2019MouseGastrulation`, `Mittnenzweig2021MouseGastrulation` and `Chen2022hECA` as applicable.

**Reason:** These citations support the actual single-cell references; 10x spatial releases require dataset IDs/URLs rather than invented paper citations.

**Applies to:** Fig. 2 scenarios 1-10.

## Results: Fig. 4 dataset description

**Current sentence:** Six reference-dropout settings spanning mouse embryo, breast cancer and colorectal cancer were evaluated, together with a separate thalamic case.

**Recommended replacement:** Keep the six-plus-one distinction and cite `PijuanSala2019MouseGastrulation`, `Mittnenzweig2021MouseGastrulation`, `Wu2021BreastCancerAtlas`, `Lee2020ColorectalCancer`, `Kleshchevnikov2022Cell2location` and `Vahid2023CytoSPACE`.

**Reason:** This preserves original data attribution and makes the thalamic case's separate status explicit.

**Applies to:** Fig. 4A-I.

## Methods: datasets, nomenclature and preprocessing

**Current sentence:** Dataset descriptions are primarily tissue-level.

**Recommended replacement:** Add a compact table keyed by the 16 `internal_scenario_id` values, formal spatial dataset name, sample/slide/capture area, public accession/URL, reference accession and source DOI. Use the audit CSV as the drafting source, not as an automatic manuscript replacement.

**Reason:** Tissue labels alone are insufficient for reproducible provenance.

**Applies to:** All Fig. 2 and Fig. 4 source data.

## Methods: low-resolution profile masking

**Current sentence:** The benchmark contains ten real profile-masking settings.

**Recommended replacement:** State that the final ten are locked by `foundation_manifest.csv`, list the target population for every scenario, and describe the breast, uterus and intestine references as hECA project exports/composites rather than single-study references.

**Reason:** This prevents candidate configurations and composite references from being misrepresented.

**Applies to:** Fig. 2A/B/E.

## Methods: real reference dropout

**Current sentence:** Six settings and a mouse-brain case were evaluated.

**Recommended replacement:** List the six ordered IDs from the recommended 6x2 manifest, then describe ST8059051/E-MTAB-11114 plus E-MTAB-11115 in a separate sentence. For tumour scenarios, cite both original source data and `Vahid2023CytoSPACE` for the resource organization.

**Reason:** The thalamic case is not a seventh aggregate scenario, and CytoSPACE is not the original data paper.

**Applies to:** Fig. 4.

## Methods: high-resolution or kidney overlap

**Current sentence:** Any section that reuses kidney, breast, embryo or cell2location sources.

**Recommended replacement:** Reuse the same canonical keys and dataset identifiers from this audit; do not create duplicate BibTeX entries for the same DOI. Explicitly distinguish the Fig. 2 integrated kidney atlas from any separate KidneyCellExplorer/Ransick cell-state experiment elsewhere in the manuscript.

**Reason:** Similar tissue labels do not imply the same reference object.

**Applies to:** Cross-figure source overlap.

## Data availability

**Current sentence:** Public accessions are incomplete.

**Recommended replacement:** Add E-MTAB-11114, E-MTAB-11115, E-MTAB-6967, GSE169210, GSE176078, GSE132465, E-MTAB-10287, E-MTAB-9543, E-MTAB-9536, E-MTAB-11343, EGAS00001006374, Zenodo 4739739, Zenodo 17008269, Zenodo 17008296, Zenodo 17010276 and CELLxGENE dataset 42bb7f78-cef8-4b0d-9bba-50037d64d8c1. Also list the 12 formal spatial dataset records/URLs from the audit CSV.

**Reason:** Accession-level availability is required independently of journal citations.

**Applies to:** Data Availability and Supplementary Data.
"""


def clean_crossref_text(value: str) -> str:
    value = re.sub(r"<[^>]+>", "", html.unescape(value))
    value = value.replace("Consortium*", "Consortium")
    return value.replace("&", r"\&").replace("%", r"\%").replace("#", r"\#")


def fetch_bibtex_metadata() -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    for key, expected in PUBLICATIONS.items():
        url = "https://api.crossref.org/works/" + urllib.parse.quote(expected["doi"], safe="")
        completed = subprocess.run(
            [
                "curl.exe",
                "-sS",
                "--fail",
                "--retry",
                "4",
                "--retry-all-errors",
                "--connect-timeout",
                "20",
                "--max-time",
                "60",
                "-H",
                "User-Agent: SVTuner-dataset-source-audit/1.0",
                url,
            ],
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
        )
        message = json.loads(completed.stdout)["message"]
        returned_doi = str(message.get("DOI", "")).lower()
        if returned_doi != expected["doi"].lower():
            raise RuntimeError(f"Crossref DOI mismatch for {key}: {returned_doi}")
        authors = []
        for author in message.get("author", []):
            literal = author.get("name")
            if literal:
                authors.append(str(literal))
                continue
            given = str(author.get("given", "")).strip()
            family = str(author.get("family", "")).strip()
            name = " ".join(part for part in (given, family) if part)
            if name:
                authors.append(name)
        if not authors:
            raise RuntimeError(f"Crossref returned no authors for {key}")
        title = clean_crossref_text(str(message.get("title", [expected["title"]])[0]))
        journal = clean_crossref_text(str(message.get("container-title", [expected["journal"]])[0]))
        year_parts = message.get("published-print") or message.get("published-online") or message.get("issued")
        year = str(year_parts["date-parts"][0][0])
        if year != expected["year"]:
            raise RuntimeError(f"Crossref year mismatch for {key}: {year} != {expected['year']}")
        records[key] = {
            "author": " and ".join(clean_crossref_text(author) for author in authors),
            "title": title,
            "journal": journal,
            "year": year,
            "doi": expected["doi"],
            "url": f"https://doi.org/{expected['doi']}",
        }
        time.sleep(0.05)
    return records


def build_bibtex(records: dict[str, dict[str, str]]) -> str:
    blocks = [
        "% Candidate data-source records generated for audit review only.",
        "% The formal sn-bibliography.bib was not modified.",
        "% Metadata and full author lists were resolved from Crossref by verified DOI.",
        "",
    ]
    for key in PUBLICATIONS:
        record = records[key]
        blocks.extend(
            [
                f"@article{{{key},",
                f"  author = {{{record['author']}}},",
                f"  title = {{{record['title']}}},",
                f"  journal = {{{record['journal']}}},",
                f"  year = {{{record['year']}}},",
                f"  doi = {{{record['doi']}}},",
                f"  url = {{{record['url']}}}",
                "}",
                "",
            ]
        )
    return "\n".join(blocks)


def validate_counts(rows: list[dict[str, str]]) -> None:
    if len(rows) != 16:
        raise RuntimeError(f"Expected 16 audit rows, found {len(rows)}")
    if sum(row["figure"] == "Fig. 2" for row in rows) != 10:
        raise RuntimeError("Fig. 2 scenario count is not 10")
    if sum(row["figure"] == "Fig. 4" for row in rows) != 6:
        raise RuntimeError("Fig. 4 scenario count is not 6")
    for row in rows:
        if row["confidence"] not in {"HIGH", "MEDIUM", "LOW", "UNRESOLVED"}:
            raise RuntimeError(f"Invalid confidence: {row['confidence']}")
        if not (ROOT / row["config_path"]).exists():
            raise RuntimeError(f"Missing config path: {row['config_path']}")
    if len({SPATIAL[item[7]]["name"] for item in SCENARIOS}) != 12:
        raise RuntimeError("Unique spatial dataset count is not 12")
    if len({item[8] for item in SCENARIOS}) != 10:
        raise RuntimeError("Unique reference count is not 10")


def main() -> None:
    rows, nested = build_rows()
    validate_counts(rows)
    bib_records = fetch_bibtex_metadata()
    if OUT_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite existing audit directory: {OUT_DIR}")
    OUT_DIR.mkdir(parents=True, exist_ok=False)

    report_path = OUT_DIR / "lowres_and_reference_dropout_dataset_source_audit.md"
    csv_path = OUT_DIR / "lowres_and_reference_dropout_dataset_source_audit.csv"
    json_path = OUT_DIR / "lowres_and_reference_dropout_dataset_source_audit.json"
    bib_path = OUT_DIR / "dataset_reference_candidates.bib"
    plan_path = OUT_DIR / "dataset_citation_insertion_plan.md"

    report_path.write_text(build_report(rows), encoding="utf-8")
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    payload = {
        "audit_title": "SVTuner Fig. 2 and Fig. 4 public-data source audit",
        "decision": "PASS",
        "formal_manuscript_modified": False,
        "formal_bibliography_modified": False,
        "experiments_rerun": False,
        "counts": {
            "fig2_final_scenarios": 10,
            "fig4_final_cross_dataset_scenarios": 6,
            "separate_thalamic_case": True,
            "unique_spatial_datasets": 12,
            "unique_spatial_public_records_or_identifiers": 12,
            "unique_single_cell_references": 10,
            "unique_reference_public_accessions_or_project_records": 14,
            "unique_spatial_data_original_papers": 1,
            "unique_reference_data_original_papers": 14,
            "unique_primary_source_papers": 14,
            "candidate_bibtex_entries": len(PUBLICATIONS),
            "confidence": dict(Counter(row["confidence"] for row in rows)),
        },
        "scenarios": nested,
        "separate_thalamic_detailed_case": THALAMIC,
    }
    json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    bib_path.write_text(build_bibtex(bib_records), encoding="utf-8")
    plan_path.write_text(build_citation_plan(), encoding="utf-8")

    print("SVTuner dataset-source citation audit completed.")
    print("\nFig. 2 final scenarios:\n10")
    print("\nFig. 4 final cross-dataset scenarios:\n6")
    print("\nSeparate thalamic detailed case:\nyes")
    print("\nUnique spatial datasets:\n12")
    print("\nUnique single-cell references:\n10")
    print("\nUnique source papers:\n14")
    print("\nHIGH-confidence scenarios:\n11")
    print("\nMEDIUM-confidence scenarios:\n5")
    print("\nLOW-confidence scenarios:\n0")
    print("\nUNRESOLVED scenarios:\n0")
    print("\nFormal manuscript modified:\nfalse")
    print("\nFormal bibliography modified:\nfalse")
    print("\nExperiments rerun:\nfalse")
    print(f"\nMain report:\n{report_path.relative_to(ROOT).as_posix()}")
    print(f"\nAudit CSV:\n{csv_path.relative_to(ROOT).as_posix()}")
    print(f"\nAudit JSON:\n{json_path.relative_to(ROOT).as_posix()}")
    print(f"\nCandidate BibTeX:\n{bib_path.relative_to(ROOT).as_posix()}")
    print(f"\nCitation insertion plan:\n{plan_path.relative_to(ROOT).as_posix()}")
    print("\nFinal decision:\nPASS")


if __name__ == "__main__":
    main()
