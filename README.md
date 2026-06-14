# DNA Methylation Deconvolution Reference Map

## Overview
For the deconvolution task, a reference map is required. In this project, we construct a DNA methylation deconvolution reference map using DNA methylation density as the main feature.

---

## DNA Methylation Density Definition
DNA methylation density is defined as:

$$
D_{block} = \frac{1}{N}\sum_{i=1}^{N} \beta_i
$$

---

## Data Sources

We use two public datasets to construct the reference map:

1. **Main methylation atlas**  
   *A DNA Methylation Atlas of Normal Human Cell Types*  
   - Covers almost all major human cell types.

2. **Megakaryocyte dataset**  
   *DNA Methylation Dynamics of Human Hematopoietic Stem Cell Differentiation*  
   - Used to supplement Megakaryocyte methylation profiles.

---

## Reference Map Construction Pipeline

1. Convert public methylation datasets into DNA methylation density format.
2. Combine datasets using:  
   `combined.py`  

3. Process the combined matrix with the following steps:

   - Drop missing values (NA)
   - Compute mean methylation density grouped by cell/tissue type
   - Filter blocks with CpG count < 5 using:  
     `hg38_500bp_blocks_with_CpG_dyad.bed`
   - Filter blocks with coefficient of variation (CV) < 0.35
   - Compute:
     - Hyper score: `this_value - other_max`
     - Hypo score: `other_min - this_value`
   - Keep only positive hyper/hypo scores
   - Rank blocks by absolute score values
   - Select top 50 blocks per cell type/tissue
   - Remove duplicated blocks

---

## Final Outputs

- Reference map:  
  `reference_map_df_CpG5_CV0.35_Top50_NoDis.parquet`
form:
```
reference_map_df.iloc[0:5, 0:5]
| Block | Adipocytes | Aorta-Endothel | Aorta-Smooth-Muscle | Bladder-Epithelial | Bladder-Smooth-Muscle |
|------|------------|----------------|---------------------|--------------------|----------------------|
| chr18_79358000_79358500 | 67.405333 | 20.8330 | 27.083 | 29.0000 | 8.333 |
| chr5_181053500_181054000 | 21.302000 | 5.2405 | 3.851 | 5.5152 | 4.864 |
| chr8_127391500_127392000 | 16.136667 | 46.1000 | 31.464 | 46.8748 | 24.076 |
| chr22_38317000_38317500 | 11.347667 | 2.3370 | 4.028 | 3.5728 | 2.410 |
| chr1_153774500_153775000 | 13.373667 | 26.5775 | 24.642 | 41.6622 | 22.517 |

```


- Annotation file:  
  `final_marker_annotation_CpG5_CV0.35_Top50_NoDis.parquet`  
  (Contains hyper/hypo methylation labels and CpG counts of each blocks)
form:
```
final_marker_annotation.iloc[0:5, 0:5]
| block | tissue | marker_type | cpg_count |
|------|--------|-------------|----------|
| chr10_100342000_100342500 | Cortex-Neuron | hyper | 5 |
| chr10_100659500_100660000 | Kidney-Glomerular-Epithelial | hyper | 32 |
| chr10_100881000_100881500 | Small-int-Epithelial | hypo | 5 |
| chr10_101113000_101113500 | Breast-Basal-Epithelial | hyper | 7 |
| chr10_101123000_101123500 | Liver-Hepatocytes | hyper | 28 |

```

---

## Important Notes

Because raw matrices are large, only the final reference files are provided.

---

## Genome & CpG Reference Consistency

To ensure compatibility when using this reference map, you MUST use the same genome and CpG definitions.

### Genome
- hg38 reference genome:  
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

### CpG Definitions
- CpG dyad file:  
  `hg38_CpG_dyad_clean.sorted.bed`

- Generated using `find_CpG.slurm`  
- "clean" means CpGs located in standard autosomal chromosomes  
- File is sorted

### Genomic Blocks
- Block file:  
  `hg38_500bp_blocks_with_CpG_dyad.bed`

- Generated from:
  - `hg38_500bp_four_col_blocks.bed` (via `block_generate.sh`)
  - processed using `find_CpG.slurm`

Each block records:
- genomic coordinates
- CpG counts

---

## Methylation Level Calculation

### Strand-aware case
If Watson and Crick strand depths are available:

$$
\beta = \frac{mC_{Watson} + mC_{Crick}}{depth_{Watson} + depth_{Crick}}
$$

### Strand-agnostic case
If strand-specific depth is not available:

$$
\beta = \frac{\beta_{Watson} + \beta_{Crick}}{2}
$$

---

## Methylation Density Input Format

Your input file should follow this format:

```
chrom  start   end     methylation_density
```

Example:

```
chr1    12500   13000   7.143
chr1    13000   13500   51.000
chr1    13500   14000   30.000
```

---


How i process Main methylation atlas and Megakaryocyte dataset

## ENCODE BigWig Processing Pipeline

Scripts used:

- BigWig → BED conversion:  
  `ToBed.slurm`

- Methylation density per tissue/cell type:  
  `MeDen_per_tissue_cell.slurm`

Big wig input
```
bigWigToBedGraph GSM5652176_Adipocytes-Z000000T7.hg38.bigwig  stdout | head -10
chr1	10468	10470	0.667
chr1	10470	10472	0.333
chr1	10483	10485	0.667
chr1	10488	10490	1
chr1	10492	10494	0.5
chr1	10496	10498	0.5
```



Example outputs:

### Sorted BED
```
chr1    10468   10470   0.667
chr1    10470   10472   1
```

### Methylation density BED
```
chr1    10000   10500   86.117
chr1    10500   11000   59.306
```

---

## Single-cell Megakaryocyte Data Processing

Pipeline:

- Source directory:  
  `mk.slurm`

- Merge all samples from one individual
- Convert into CpG dyad-based methylation density

Example raw data:
```
chr15   17005692   1   1
```

Final output:
```
chr1    10000   10500   50.000
chr1    10500   11000   8.537
```

---

## Tissue Grouping Strategy (Optional)

To reduce granularity, similar tissues/cell types can be grouped as follows:

```python
tissue_mapping = {
    'Macrophages/Monocytes': [
        'Lung-Alveolar-Macrophages', 'Lung-Interstitial-Macrophages', 'Liver-Macrophages', 
        'Colon-Macrophages', 'Blood-Monocytes'
    ],
    'T cells': [
        'Blood-T-Naive-CD8', 'Blood-T-Naive-CD4', 'Blood-T-EffMem-CD8', 'Blood-T-EffMem-CD4',
        'Blood-T-Eff-CD8', 'Blood-T-CD8', 'Blood-T-CenMem-CD4', 'Blood-T-CD4', 'Blood-T-CD3'
    ],
    'B cells': ['Blood-B-Mem', 'Blood-B'],
    'NK cells': ['Blood-NK'],
    'Granulocytes': ['Blood-Granulocytes'],
    'Megakaryocytes': ['Megakaryocytes'],
    'Bone': ['Bone-Osteoblasts'],
    'Erythrocyte progenitor': ['Bone_marrow-Erythrocyte_progenitors'],
    'Liver': ['Liver-Hepatocytes'],
    'Brain': ['Neuron', 'Cerebellum-Neuron', 'Cortex-Neuron', 'Oligodendrocytes'],
    'Lung': [
        'Lung-Pleura', 'Lung-Bronchus-Epithelial',
        'Lung-Alveolar-Epithelial'],
    'Kidney': [
        'Kidney-Tubular-Epithelial', 'Kidney-Glomerular-Podocytes', 'Kidney-Glomerular-Epithelial'
    ],
    'Thyroid': ['Thyroid-Epithelial'],
    'Larynx': ['Larynx-Epithelial'],
    'Pharynx': ['Pharynx-Epithelial'],
    'Tongue': ['Tongue-Epithelial', 'Tongue_base-Epithelial'],
    'Tonsil': ['Tonsil-Palatine-Epithelial', 'Tonsil-Pharyngeal-Epithelial'],
    'Gastric': [
        'Gastric-antrum-Endocrine', 'Gastric-antrum-Epithelial', 'Gastric-fundus-Epithelial',
        'Gastric-body-Epithelial'
    ],
    'Small intestine': ['Small-int-Epithelial', 'Small-int-Endocrine'],
    'Colon': [
        'Colon-Right-Endocrine', 'Colon-Right-Epithelial', 'Colon-Fibroblasts',
        'Colon-Left-Endocrine', 'Colon-Left-Epithelial'
    ],
    'Esophagus': ['Esophagus-Epithelial'],
    'Gallbladder': ['Gallbladder-Epithelial'],
    'Pancreas': [
        'Pancreas-Acinar', 'Pancreas-Alpha', 'Pancreas-Delta', 'Pancreas-Beta', 'Pancreas-Duct'
    ],
    'Adipocytes': ['Adipocytes'],
    'Bladder': ['Bladder-Epithelial'],
    'Skin': ['Epidermal-Keratinocytes', 'Dermal-Fibroblasts'],
    'Endometrium': ['Endometrium-Epithelial'],
    'Fallopian': ['Fallopian-Epithelial'],
    'Ovary': ['Ovary-Epithelial'],
    'Breast': ['Breast-Basal-Epithelial', 'Breast-Luminal-Epithelial'],
    'Prostate': ['Prostate-Epithelial'],
    'Heart': ['Heart-Cardiomyocyte', 'Heart-Fibroblasts'],
    'Muscle': ['Skeletal-Muscle'],
    'Endothelial cells': ['Saphenous-Vein-Endothel','Aorta-Endothel', 'Pancreas-Endothel', 'Pancreas-Islet-Endothel', 'Lung-Alveolar-Endothel', 'Kidney-Glomerular-Endothel','Kidney-Tubular-Endothel', 'Liver-Endothelium'],
    'Smooth muscle cells': ['Aorta-Smooth-Muscle', 'Coronary-Artery-Smooth-Muscle', 'Prostate-Smooth-Muscle', 'Bladder-Smooth-Muscle', 'Lung-Bronchus-Smooth-Muscle']
}

```

---

## Notes

- Classification granularity can be adjusted depending on downstream analysis needs.
- The provided mapping is not fixed and can be modified.
