# Secondary mito delivery
This project is ment to check the capacity of a cell to do secondary mitocondria delivery, and to find the best cell that can do that in evry tissue.

The logic is this:

We check if the cell have the ability to express gene associated with mitocondria donation, and safety of this cell and donation, and have a score.

# Reasoning

In the first model, the gene are chose manualy:

Ability | Why Important for Mito Transfer | Associated Genes
---|---|---
**Ability to make Vesicle** | Package mitochondria into membrane-bound carriers for safe transport between cells without direct contact | TSG101, ALIX, VPS4A, VPS4B, SDCBP, CD63, CD9, CD81, RAB27A, RAB27B, PDCD6IP
**Ability to make Nanotube (TNT)** | Form direct cytoplasmic bridges for high-efficiency, contact-dependent mitochondrial transfer via intercellular tunnels | RHOA, CDC42, ACTA2, INF2, FMNL1, FMNL2, ARPC2, ARPC3, MYH9
**Rab Traffic** | GTPases that dock vesicles at the plasma membrane and ensure directional secretion of exosomes toward recipient cells | RAB27A, RAB27B, RAB11A, RAB11B, RAB35, RAB5A, RAB7A, RAB8A
**Adhesion** | Enable cell-cell contact and homing to damaged tissue; prerequisite for both vesicle transfer and nanotube formation | CXCR4, CCR7, ITGA4, ITGB1, ICAM1, SELP, SELE, L1CAM
**Metabolism** | High metabolic activity provides ATP and cofactors to power mitochondrial packaging, vesicle biogenesis, and transfer machinery | PFKFB3, PFKFB4, SLC2A1, LDHA, CPT1A, HIF1A, EPAS1
**Biogenesis** | Continuously generate new, healthy mitochondria via PGC-1α-driven transcription; ensures donor pool remains robust | TFAM, PPARGC1A, PPARGC1B, NRF1, NRF2, ESRRA
**Fusion** | Merge mitochondrial membranes to maintain cristae integrity and respiratory capacity; fused mito are more stable donors | OPA1, MFN1, MFN2, MFUS, MAOB
**Fission** | Fragment mitochondria into transfer-competent units that fit into vesicles or nanotube lumens | DNM1L (DRP1), FIS1, IMMT, MIEF1, MIEF2
**OXPHOS** | Ensure donor mitochondria produce ATP; non-functional OXPHOS = useless gift to recipient cell | ATP5F1A, NDUFA1, COX7A1, CYB, UQCRC1, SDHA
**Remove damage mitochondria** | Remove damaged mitochondria before donation via mitophagy; only healthy mito reach recipient | PINK1, PRKN, BNIP3, NIX, OPTN, AMBRA1
**Calcium** | Regulate intracellular Ca²⁺; Ca²⁺ imbalance triggers mitochondrial damage or apoptosis during transfer | MCU, MICU1, MICU2, EFHA1, SLC8B1, CALM
**Transport** | Move mitochondrial cargo (lipids, metabolites) into vesicles; substrate-specific transporters for different mito membrane compartments | SLC25A1, SLC25A3, SLC25A4, SLC25A5, SLC25A10, SLC25A24, CPT1A
**Structure** | Maintain cristae organization and protein folding capacity; collapsed cristae = loss of electron transport chain | OPA1, IMMT, MIC60, CHCHD3, CHCHD6, MTERF3
**Safety** | Inflammatory, Antigen presentation, pro-apoptic, danger signal | TNF, IL1B, IL6, IFNG, CD69, HLA-DR, BAX,, HMGB1

After **thinking more** instead of relying on assumptions, we use MitoCarta and PNAS data to identify which genes are actually present and validated.

| Feature | Associated Genes |
|---|---|
| Exosome Machinery | TSG101, PDCD6IP, VPS4A, VPS4B, CD63, CD9, CD81, SDCBP |
| Donor Fitness (Mito Health) | TFAM, OPA1, MFN1, MFN2, DNM1L, PINK1, PRKN, PPARGC1A, ATP5F1A |
| Rab Trafficking | RAB27A, RAB27B, RAB11A, RAB11B, RAB35, RAB5A, RAB7A |
| MDV Pathway | VPS35, SNX9, ATG9A, MUL1, RAB9A, OPTN |
| TNT Capacity | RHOA, CDC42, ACTA2, INF2, FMNL1, FMNL2, ARPC2, ARPC3 |
| Adhesion & Migration | CXCR4, CCR7, ITGA4, ITGB1, ICAM1 |
| Metabolic Vigor | PFKFB3, SLC2A1, LDHA, CPT1A, HIF1A |
| Microvesicle Formation | ARF6, ROCK1, ROCK2, MYH9, ARRDC1, FLOT1, FLOT2, PLD2 |
| Immune Shedding (Safety+) | CD63, LAMP1, TSG101 |
| Activation Markers (Safety−) | IL2RA, CD69, TNF, IL1B, IFNG, HLA-DRA |
| Viability / Survival | BCL2, MCL1, BIRC2, XIAP, BAX, BAK1 |

# Scoring Methodology
 
The scoring process follows these steps:
 
1. Identify the immune cell type
2. Check if the associated genes are found and expressed in the dataset
3. Group cells by cell type and tissue
4. Calculate normalized feature scores for each group

 
The Über Score is lovely called "Über" because the cell essentially acts as an Uber driver for mitochondria :)
 
```
Über Score = 0.6 × [Σ(w_f × score_f) + 0.10×shedding − 0.25×activation]
           + 0.4 × viability
```
 
Where:
- `w_f` = weight of feature f (e.g., Exosome = 0.20, Mito Health = 0.18)
- `score_f` = min-max normalized score (0→1) across all (cell_type, tissue) groups
- `shedding` = bonus for controlled immune shedding
- `activation` = penalty for pro-inflammatory markers
- `viability` = anti-apoptotic capacity (BCL2, MCL1, XIAP, BIRC2) minus pro-apoptotic markers (BAX, BAK1)

# Data Source
 
The project uses the Human Immune Health Atlas dataset, but any single-cell RNA-seq dataset with immune cells and gene expression data can be used.
 
# Results and Output

1. **Cell Type × Tissue Heatmap**: Rows are cell types, columns are tissues. Color intensity represents score strength (red = high, white = absent).
(!image)()
2. **Top-Ranked Cells per Tissue**: Displays the best candidate cell type for each tissue with its score breakdown.
(!image)()
3. **Feature × Cell Type Heatmap**: Shows expression of individual features (rows) across cell types (columns), stratified by tissue.
(!image)()


Example Result
 
**Best cell for anterior segment of eyeball (blood):**
 
Neutrophil
- Über Score = 0.6577
- Strong in:
  - TNT Capacity: 0.987 (excellent nanotube formation)
  - Metabolism: 1.078 (maximum energy)
  - Exosome: 0.571 (moderate exosome production)
- Viability: 0.891 (healthy cell)


# TNT vs. Vesicle Delivery Route

For each top-scoring cell type, the optimal delivery mechanism is determined based on feature scores:

example

TOP 5 TISSUES FOR VESICLE TRANSFER:
  submandibular gland                  intermediate monocyte         vesicle=0.614
  buccal mucosa                        T cell                        vesicle=0.538
  crista ampullaris                    mast cell                     vesicle=0.509
  cardiac atrium                       monocyte                      vesicle=0.480
  inguinal lymph node                  hematopoietic precursor cell  vesicle=0.465

TOP 5 TISSUES FOR TNT TRANSFER:
  uterus                               myeloid dendritic cell        TNT=0.852
  endocrine pancreas                   monocyte                      TNT=0.822
  brown adipose tissue                 platelet                      TNT=0.811
  exocrine pancreas                    intermediate monocyte         TNT=0.741
  sublingual gland                     leukocyte                     TNT=0.728

# Running the Pipeline
 
```bash
python score_vesicle_nanotube.py       # Compute Über Scores
python visualize_sec_transfer.py       # Generate heatmaps and figures
```
 
# Disclaimer
 
Grammar and code review performed with assistance from Claude (Anthropic AI).
