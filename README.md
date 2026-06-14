# DNA Methylation Deconvolution Reference Map

## Overview
For the deconvolution task, a reference map is required. In this project, we construct a DNA methylation deconvolution reference map using DNA methylation density as the main feature.

---

## How to Use

1. Generate `hg38_CpG_dyad_clean.sorted.bed` and `hg38_500bp_blocks_with_CpG_dyad.bed`.
2. Download reference map and reference map annotation.
3. Reorganize your methylation data into this form:

```
chr1    10468   10470   0.667
chr1    10470   10472   1
```

4. Use `MeDen_per_tissue_cell.sh` to convert BED file to methylation density file.

Example:
```
chr1    12500   13000   7.143
chr1    13000   13500   51.000
chr1    13500   14000   30.000
```

5. Use `run_nnls.ipynb` to perform deconvolution on your sample.



---

## DNA Methylation Density Definition
DNA methylation density is defined as:

$$
D_{block} = \frac{1}{N}\sum_{i=1}^{N} \beta_i
$$

---


## Methylation Level (Beta) Calculation

Note: We do not consider strand-specific DNA methylation levels. Instead, we treat:

- 5' C–G 3'
- 3' G–C 5'

as the same CpG site.

Below are three cases for calculating beta values.

---

### Strand-aware Case 1

If Watson and Crick strand depths are available:

Example:
```
#CHROM  START   END     MOD_LEVEL   MOD UNMOD REF ALT SPECIFIC_CONTEXT CONTEXT SNV TOTAL_DEPTH
chr1    48150   48151   1.0         5   0     C   T   CGG              CpG     No  29
chr1    48151   48152   0.895       17  2     G   A   CGC              CpG     No  29
chr1    48172   48173   1.0         4   0     C   T   CGT              CpG     No  33
chr1    48173   48174   0.867       13  2     G   A   CGA              CpG     No  33
```

$$
\beta = \frac{mC_{Watson} + mC_{Crick}}{depth_{Watson} + depth_{Crick}}
$$

Example:
```
beta_48150-48152 = (5 + 17) / (5 + 0 + 17 + 2)
```

---

### Strand-aware Case 2

If strand-specific depth is not available:

$$
\beta = \frac{\beta_{Watson} + \beta_{Crick}}{2}
$$

---

### Non-strand-aware Case

```
chr1    10468   10470   0.667
chr1    10470   10472   1
chr1    10483   10485   0.667
chr1    10488   10490   1
```

In this case, strand information is not distinguished, so no merging is required.

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
`generate_test_ref_mat.ipynb`

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
  
Example (`reference_map_df.iloc[0:5, 0:5]`):

| Block | Adipocytes | Aorta-Endothel | Aorta-Smooth-Muscle | Bladder-Epithelial | Bladder-Smooth-Muscle |
|------|------------|----------------|---------------------|--------------------|----------------------|
| chr18_79358000_79358500 | 67.405333 | 20.8330 | 27.083 | 29.0000 | 8.333 |
| chr5_181053500_181054000 | 21.302000 | 5.2405 | 3.851 | 5.5152 | 4.864 |
| chr8_127391500_127392000 | 16.136667 | 46.1000 | 31.464 | 46.8748 | 24.076 |
| chr22_38317000_38317500 | 11.347667 | 2.3370 | 4.028 | 3.5728 | 2.410 |
| chr1_153774500_153775000 | 13.373667 | 26.5775 | 24.642 | 41.6622 | 22.517 |


- Annotation file:  
  `final_marker_annotation_CpG5_CV0.35_Top50_NoDis.parquet`  
  (Contains hyper/hypo methylation labels and CpG counts of each blocks)

Example (`final_marker_annotation.iloc[0:5, 0:4]`):

| block | tissue | marker_type | cpg_count |
|------|--------|-------------|----------|
| chr10_100342000_100342500 | Cortex-Neuron | hyper | 5 |
| chr10_100659500_100660000 | Kidney-Glomerular-Epithelial | hyper | 32 |
| chr10_100881000_100881500 | Small-int-Epithelial | hypo | 5 |
| chr10_101113000_101113500 | Breast-Basal-Epithelial | hyper | 7 |
| chr10_101123000_101123500 | Liver-Hepatocytes | hyper | 28 |

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

```
head -10 hg38_CpG_dyad_clean.sorted.bed
chr1    10468   10470   CpG_1
chr1    10470   10472   CpG_2
chr1    10483   10485   CpG_3
chr1    10488   10490   CpG_4
chr1    10492   10494   CpG_5
chr1    10496   10498   CpG_6
chr1    10524   10526   CpG_7
chr1    10541   10543   CpG_8
chr1    10562   10564   CpG_9
chr1    10570   10572   CpG_10
```


### Genomic Blocks
- Block file:  
  `hg38_500bp_blocks_with_CpG_dyad.bed`

- Generated from:
  - `hg38_500bp_four_col_blocks.bed` (via `block_generate.sh`)
  - processed using `find_CpG.slurm`

```
head -10 hg38_500bp_blocks_with_CpG_dyad.bed
chr1    10000   10500   block_21        6
chr1    10500   11000   block_22        82
chr1    11000   11500   block_23        36
chr1    11500   12000   block_24        9
chr1    12000   12500   block_25        8
chr1    12500   13000   block_26        14
chr1    13000   13500   block_27        5
chr1    13500   14000   block_28        10
chr1    14000   14500   block_29        7
chr1    14500   15000   block_30        16
```

Each block records:
- genomic coordinates
- CpG counts

---



## Testing Results

In the last part of the notebook, I used ENCODE pancreas WGBS data and breast cancer TAPS data to test my reference map.

The results are shown below. For the ENCODE pancreas deconvolution, our results are highly consistent with those reported in the paper "A DNA Methylation Atlas of Normal Human Cell Types". In Figure 6i of that paper, the authors applied their own method to deconvolve the same ENCODE pancreas WGBS dataset (you can verify this).

For the breast cancer deconvolution, I did not compare our results with any published papers from other researchers.


```
============================================================
ENCODE_Pancreas1
共同 block 数: 4107
residual: 437.7591038289946
                                 tissue  coefficient  proportion
0                       Pancreas-Acinar     0.653705    0.670674
1                         Pancreas-Duct     0.206556    0.211918
2                     Pancreas-Endothel     0.032923    0.033778
3                        Megakaryocytes     0.030332    0.031119
4                   Aorta-Smooth-Muscle     0.024442    0.025077
5         Coronary-Artery-Smooth-Muscle     0.015831    0.016242
6          Kidney-Glomerular-Epithelial     0.003183    0.003266
7                         Pancreas-Beta     0.002719    0.002789
8                     Colon-Macrophages     0.002588    0.002655
9                     Cerebellum-Neuron     0.001403    0.001439
10              Pancreas-Islet-Endothel     0.000549    0.000563
11               Prostate-Smooth-Muscle     0.000396    0.000407
12               Endometrium-Epithelial     0.000072    0.000073
13                              Blood-B     0.000000    0.000000
14                          Blood-B-Mem     0.000000    0.000000
15                   Bladder-Epithelial     0.000000    0.000000
16                Bladder-Smooth-Muscle     0.000000    0.000000
17                           Adipocytes     0.000000    0.000000
18                       Aorta-Endothel     0.000000    0.000000
19                    Blood-T-Naive-CD8     0.000000    0.000000
20                    Blood-T-Naive-CD4     0.000000    0.000000
21                   Blood-T-EffMem-CD8     0.000000    0.000000
22                   Blood-T-EffMem-CD4     0.000000    0.000000
23                      Blood-T-Eff-CD8     0.000000    0.000000
24                   Blood-T-CenMem-CD4     0.000000    0.000000
25                          Blood-T-CD8     0.000000    0.000000
26                          Blood-T-CD4     0.000000    0.000000
27                          Blood-T-CD3     0.000000    0.000000
28                             Blood-NK     0.000000    0.000000
29                      Blood-Monocytes     0.000000    0.000000
30                   Blood-Granulocytes     0.000000    0.000000
31                Colon-Left-Epithelial     0.000000    0.000000
32            Breast-Luminal-Epithelial     0.000000    0.000000
33              Breast-Basal-Epithelial     0.000000    0.000000
34                     Bone-Osteoblasts     0.000000    0.000000
35              Epidermal-Keratinocytes     0.000000    0.000000
36                 Esophagus-Epithelial     0.000000    0.000000
37               Gallbladder-Epithelial     0.000000    0.000000
38                 Fallopian-Epithelial     0.000000    0.000000
39            Gastric-antrum-Epithelial     0.000000    0.000000
40              Gastric-body-Epithelial     0.000000    0.000000
41            Gastric-fundus-Epithelial     0.000000    0.000000
42             Gastric-antrum-Endocrine     0.000000    0.000000
43                  Heart-Cardiomyocyte     0.000000    0.000000
44                 Colon-Left-Endocrine     0.000000    0.000000
45                    Colon-Fibroblasts     0.000000    0.000000
46                Colon-Right-Endocrine     0.000000    0.000000
47               Colon-Right-Epithelial     0.000000    0.000000
48                        Cortex-Neuron     0.000000    0.000000
49                   Dermal-Fibroblasts     0.000000    0.000000
50  Bone_marrow-Erythrocyte_progenitors     0.000000    0.000000
51                    Liver-Endothelium     0.000000    0.000000
52                    Larynx-Epithelial     0.000000    0.000000
53            Kidney-Tubular-Epithelial     0.000000    0.000000
54              Kidney-Tubular-Endothel     0.000000    0.000000
55          Kidney-Glomerular-Podocytes     0.000000    0.000000
56           Kidney-Glomerular-Endothel     0.000000    0.000000
57                    Heart-Fibroblasts     0.000000    0.000000
58                    Liver-Hepatocytes     0.000000    0.000000
59        Lung-Interstitial-Macrophages     0.000000    0.000000
60          Lung-Bronchus-Smooth-Muscle     0.000000    0.000000
61             Lung-Bronchus-Epithelial     0.000000    0.000000
62                          Lung-Pleura     0.000000    0.000000
63            Lung-Alveolar-Macrophages     0.000000    0.000000
64                    Liver-Macrophages     0.000000    0.000000
65               Lung-Alveolar-Endothel     0.000000    0.000000
66             Lung-Alveolar-Epithelial     0.000000    0.000000
67                       Pancreas-Alpha     0.000000    0.000000
68                     Ovary-Epithelial     0.000000    0.000000
69                               Neuron     0.000000    0.000000
70                     Oligodendrocytes     0.000000    0.000000
71                       Pancreas-Delta     0.000000    0.000000
72                   Pharynx-Epithelial     0.000000    0.000000
73                  Prostate-Epithelial     0.000000    0.000000
74              Saphenous-Vein-Endothel     0.000000    0.000000
75                      Skeletal-Muscle     0.000000    0.000000
76                  Small-int-Endocrine     0.000000    0.000000
77                 Small-int-Epithelial     0.000000    0.000000
78                   Thyroid-Epithelial     0.000000    0.000000
79                    Tongue-Epithelial     0.000000    0.000000
80               Tongue_base-Epithelial     0.000000    0.000000
81           Tonsil-Palatine-Epithelial     0.000000    0.000000
82         Tonsil-Pharyngeal-Epithelial     0.000000    0.000000
============================================================
ENCODE_Pancreas2
共同 block 数: 4106
residual: 468.79092580478743
                                 tissue  coefficient  proportion
0                       Pancreas-Acinar     0.583917    0.584457
1                         Pancreas-Duct     0.261615    0.261858
2                     Pancreas-Endothel     0.065720    0.065780
3                         Pancreas-Beta     0.041929    0.041968
4                        Megakaryocytes     0.013395    0.013408
5                                Neuron     0.009950    0.009960
6                   Heart-Cardiomyocyte     0.009813    0.009822
7          Kidney-Glomerular-Epithelial     0.003979    0.003983
8                Gallbladder-Epithelial     0.003851    0.003854
9              Gastric-antrum-Endocrine     0.003486    0.003489
10                       Pancreas-Delta     0.001420    0.001421
11                   Bladder-Epithelial     0.000000    0.000000
12                Bladder-Smooth-Muscle     0.000000    0.000000
13                              Blood-B     0.000000    0.000000
14                          Blood-B-Mem     0.000000    0.000000
15                   Blood-Granulocytes     0.000000    0.000000
16                           Adipocytes     0.000000    0.000000
17                       Aorta-Endothel     0.000000    0.000000
18                  Aorta-Smooth-Muscle     0.000000    0.000000
19                    Blood-T-Naive-CD8     0.000000    0.000000
20                    Blood-T-Naive-CD4     0.000000    0.000000
21                   Blood-T-EffMem-CD8     0.000000    0.000000
22                   Blood-T-EffMem-CD4     0.000000    0.000000
23                      Blood-T-Eff-CD8     0.000000    0.000000
24                   Blood-T-CenMem-CD4     0.000000    0.000000
25                          Blood-T-CD8     0.000000    0.000000
26                          Blood-T-CD4     0.000000    0.000000
27                          Blood-T-CD3     0.000000    0.000000
28                             Blood-NK     0.000000    0.000000
29                      Blood-Monocytes     0.000000    0.000000
30                    Colon-Macrophages     0.000000    0.000000
31                     Bone-Osteoblasts     0.000000    0.000000
32  Bone_marrow-Erythrocyte_progenitors     0.000000    0.000000
33              Breast-Basal-Epithelial     0.000000    0.000000
34            Breast-Luminal-Epithelial     0.000000    0.000000
35              Epidermal-Keratinocytes     0.000000    0.000000
36               Endometrium-Epithelial     0.000000    0.000000
37                   Dermal-Fibroblasts     0.000000    0.000000
38                        Cortex-Neuron     0.000000    0.000000
39                 Esophagus-Epithelial     0.000000    0.000000
40            Gastric-antrum-Epithelial     0.000000    0.000000
41              Gastric-body-Epithelial     0.000000    0.000000
42                 Fallopian-Epithelial     0.000000    0.000000
43            Gastric-fundus-Epithelial     0.000000    0.000000
44                Colon-Right-Endocrine     0.000000    0.000000
45               Colon-Right-Epithelial     0.000000    0.000000
46        Coronary-Artery-Smooth-Muscle     0.000000    0.000000
47                    Cerebellum-Neuron     0.000000    0.000000
48                    Colon-Fibroblasts     0.000000    0.000000
49                 Colon-Left-Endocrine     0.000000    0.000000
50                Colon-Left-Epithelial     0.000000    0.000000
51                    Liver-Endothelium     0.000000    0.000000
52                    Larynx-Epithelial     0.000000    0.000000
53            Kidney-Tubular-Epithelial     0.000000    0.000000
54              Kidney-Tubular-Endothel     0.000000    0.000000
55          Kidney-Glomerular-Podocytes     0.000000    0.000000
56           Kidney-Glomerular-Endothel     0.000000    0.000000
57                    Heart-Fibroblasts     0.000000    0.000000
58                    Liver-Hepatocytes     0.000000    0.000000
59        Lung-Interstitial-Macrophages     0.000000    0.000000
60          Lung-Bronchus-Smooth-Muscle     0.000000    0.000000
61             Lung-Bronchus-Epithelial     0.000000    0.000000
62            Lung-Alveolar-Macrophages     0.000000    0.000000
63                     Oligodendrocytes     0.000000    0.000000
64                    Liver-Macrophages     0.000000    0.000000
65               Lung-Alveolar-Endothel     0.000000    0.000000
66             Lung-Alveolar-Epithelial     0.000000    0.000000
67                       Pancreas-Alpha     0.000000    0.000000
68                     Ovary-Epithelial     0.000000    0.000000
69                          Lung-Pleura     0.000000    0.000000
70              Pancreas-Islet-Endothel     0.000000    0.000000
71                   Pharynx-Epithelial     0.000000    0.000000
72                  Prostate-Epithelial     0.000000    0.000000
73               Prostate-Smooth-Muscle     0.000000    0.000000
74              Saphenous-Vein-Endothel     0.000000    0.000000
75                      Skeletal-Muscle     0.000000    0.000000
76                  Small-int-Endocrine     0.000000    0.000000
77                 Small-int-Epithelial     0.000000    0.000000
78                   Thyroid-Epithelial     0.000000    0.000000
79                    Tongue-Epithelial     0.000000    0.000000
80               Tongue_base-Epithelial     0.000000    0.000000
81           Tonsil-Palatine-Epithelial     0.000000    0.000000
82         Tonsil-Pharyngeal-Epithelial     0.000000    0.000000
============================================================
breast_tumor
共同 block 数: 4005
residual: 841.9301245656545
                                 tissue  coefficient  proportion
0             Breast-Luminal-Epithelial     0.534097    0.480259
1                           Blood-B-Mem     0.108799    0.097832
2                    Bladder-Epithelial     0.087811    0.078959
3                   Small-int-Endocrine     0.056872    0.051139
4                     Colon-Fibroblasts     0.040402    0.036330
5                        Pancreas-Alpha     0.039137    0.035192
6                       Blood-T-Eff-CD8     0.035234    0.031682
7                           Lung-Pleura     0.033253    0.029901
8           Lung-Bronchus-Smooth-Muscle     0.031549    0.028368
9                    Blood-T-EffMem-CD8     0.029177    0.026236
10                          Blood-T-CD8     0.025870    0.023262
11               Endometrium-Epithelial     0.023275    0.020929
12                Colon-Left-Epithelial     0.018298    0.016454
13        Lung-Interstitial-Macrophages     0.010543    0.009480
14                    Liver-Endothelium     0.010271    0.009236
15                        Pancreas-Beta     0.008753    0.007871
16                      Pancreas-Acinar     0.008323    0.007484
17                   Dermal-Fibroblasts     0.006321    0.005684
18             Gastric-antrum-Endocrine     0.003578    0.003218
19                       Megakaryocytes     0.000537    0.000483
20                    Blood-T-Naive-CD4     0.000000    0.000000
21                    Blood-T-Naive-CD8     0.000000    0.000000
22                   Blood-T-CenMem-CD4     0.000000    0.000000
23                          Blood-T-CD4     0.000000    0.000000
24                      Blood-Monocytes     0.000000    0.000000
25                   Blood-T-EffMem-CD4     0.000000    0.000000
26                          Blood-T-CD3     0.000000    0.000000
27                              Blood-B     0.000000    0.000000
28                Bladder-Smooth-Muscle     0.000000    0.000000
29                   Blood-Granulocytes     0.000000    0.000000
30                  Aorta-Smooth-Muscle     0.000000    0.000000
31                       Aorta-Endothel     0.000000    0.000000
32                           Adipocytes     0.000000    0.000000
33                             Blood-NK     0.000000    0.000000
34  Bone_marrow-Erythrocyte_progenitors     0.000000    0.000000
35              Epidermal-Keratinocytes     0.000000    0.000000
36                        Cortex-Neuron     0.000000    0.000000
37                 Fallopian-Epithelial     0.000000    0.000000
38                 Esophagus-Epithelial     0.000000    0.000000
39            Gastric-antrum-Epithelial     0.000000    0.000000
40              Gastric-body-Epithelial     0.000000    0.000000
41            Gastric-fundus-Epithelial     0.000000    0.000000
42               Gallbladder-Epithelial     0.000000    0.000000
43                 Colon-Left-Endocrine     0.000000    0.000000
44                    Cerebellum-Neuron     0.000000    0.000000
45                    Colon-Macrophages     0.000000    0.000000
46                Colon-Right-Endocrine     0.000000    0.000000
47        Coronary-Artery-Smooth-Muscle     0.000000    0.000000
48               Colon-Right-Epithelial     0.000000    0.000000
49              Breast-Basal-Epithelial     0.000000    0.000000
50                     Bone-Osteoblasts     0.000000    0.000000
51                    Larynx-Epithelial     0.000000    0.000000
52            Kidney-Tubular-Epithelial     0.000000    0.000000
53              Kidney-Tubular-Endothel     0.000000    0.000000
54          Kidney-Glomerular-Podocytes     0.000000    0.000000
55           Kidney-Glomerular-Endothel     0.000000    0.000000
56         Kidney-Glomerular-Epithelial     0.000000    0.000000
57                    Heart-Fibroblasts     0.000000    0.000000
58                  Heart-Cardiomyocyte     0.000000    0.000000
59            Lung-Alveolar-Macrophages     0.000000    0.000000
60             Lung-Bronchus-Epithelial     0.000000    0.000000
61               Lung-Alveolar-Endothel     0.000000    0.000000
62             Lung-Alveolar-Epithelial     0.000000    0.000000
63                    Liver-Macrophages     0.000000    0.000000
64                               Neuron     0.000000    0.000000
65                     Oligodendrocytes     0.000000    0.000000
66                    Liver-Hepatocytes     0.000000    0.000000
67                     Ovary-Epithelial     0.000000    0.000000
68                       Pancreas-Delta     0.000000    0.000000
69                    Pancreas-Endothel     0.000000    0.000000
70                        Pancreas-Duct     0.000000    0.000000
71                   Pharynx-Epithelial     0.000000    0.000000
72                  Prostate-Epithelial     0.000000    0.000000
73               Prostate-Smooth-Muscle     0.000000    0.000000
74              Pancreas-Islet-Endothel     0.000000    0.000000
75              Saphenous-Vein-Endothel     0.000000    0.000000
76                      Skeletal-Muscle     0.000000    0.000000
77                 Small-int-Epithelial     0.000000    0.000000
78                   Thyroid-Epithelial     0.000000    0.000000
79                    Tongue-Epithelial     0.000000    0.000000
80               Tongue_base-Epithelial     0.000000    0.000000
81           Tonsil-Palatine-Epithelial     0.000000    0.000000
82         Tonsil-Pharyngeal-Epithelial     0.000000    0.000000

```



## How I process Main methylation atlas and Megakaryocyte dataset

### ENCODE BigWig Processing Pipeline

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

#### Sorted BED
```
chr1    10468   10470   0.667
chr1    10470   10472   1
```

#### Methylation density BED
```
chr1    10000   10500   86.117
chr1    10500   11000   59.306
```

---

### Single-cell Megakaryocyte Data Processing

Pipeline:

- Source directory:  
  `mk.slurm`

- Merge all samples from one individual
- Convert into CpG dyad-based methylation density

#### Example raw data:
```
chr15   17005692   1   1
```

#### Final output:
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


## if you use my ref map, please cite me:
https://github.com/pandas-zwx/DNA_methylation_deconvolution
