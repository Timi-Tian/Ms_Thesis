# Amyloid-β Seeding in Early Disease Stages: Effects of early-life infection and sex differences on glial cells in Alzheimer’s Disease 

This repository contains analysis code for single-cell RNA-seq from WT and 5xFAD mice under vehicle or LPS treatment at P9, focusing on glial cells, microglia, astrocytes, oligodendrocytes from mice brain and ligand–receptor (LR) communication (LIANA). 

## Repository structure
Timi_Ms_Thesis/
├─ ACs.ipynb # Astrocyte analysis
├─ MG.ipynb # Microglia analysis
├─ MOLs.ipynb # Mature Oligodendrocyte analysis
├─ adata_batch_qc.ipynb # Integrate all cells; then subset per cell type for focused analyses
├─ rawdata_qc.ipynb # Load batches; HashSolo to remove doublets and batch correction; concatenate batches
├─ LIANA.ipynb # filter and post-process LR interactions
├─ Conditional.R # Chord plot utilities (circlize/ComplexHeatmap)
├─ environment.yml # Conda environment (reproducible dependencies)
└─ results/ # Figures & tables (output)

## Data
- **Input:** processed and subsetted adatas (e.g., `.h5ad`).
- **Outputs:** figures/csvs are written to `results/`.

## Environment (dependencies)
Create the conda environment:environment.yml

## QC & integration
Run rawdata_qc.ipynb first : load in each batch file, containing features, matrix and barcodes files. Then remove doublets with HashSolo and do quality control within batches, lastly concatenate all batches into an adata.
Run adata_batch_qc.ipynb, loading a combined AnnData, then subset each cell type for further focused analyses.
Output: ACs/MG/MOL.h5ad files used as inputs for downstream steps.

## Cell-type analyses
Run ACs.ipynb, MG.ipynb, and MOLs.ipynb to perform clustering, annotation, and differential/pathway/compositional analyses, gaining ligands and receptors based on differential conditions.
Investigating the effect of AD and LPS on each cell type by comparing the differential geens bwtween WT_LPs vs.WT_Veh, 5XFAD_Veh vs.WT_Veh and 5XFAD_LPS vs.5XFAD_Veh.
Output: per–cell-type figures/tables and optional trimmed objects.

## Ligand–receptor (LIANA)
Run LIANA.ipynb to predict possible interactions between each cell type and under specific conditions, generating LR result tables for your contrasts/conditions, WT_LPs vs.WT_Veh, 5XFAD_Veh vs.WT_Veh and 5XFAD_LPS vs.5XFAD_Veh.
Extract interaction interested based on interested ligands/
Output: filtered LR tables of a pair of celltypes used for the chord plotting.

## Chord diagrams (R)
In R, source("Conditional.R"), set input/output paths and the list of conditions, then run the plotting functions/loop to produce Circos chord PDFs
Output: _Chord.pdf of each pair of celltypes
