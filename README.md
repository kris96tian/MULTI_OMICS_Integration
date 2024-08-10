# 

# MULTI-OMICS INTEGRATION APPROACH TO UNDERSTAND GENE-METABOLITE INTERACTIONS IN HFPEF

## Introduction

Heart failure with preserved ejection fraction (HFpEF) is increasingly  
recognized as a significant public health issue, poised to become the  
predominant form of heart failure globally.  
The disease's underlying molecular mechanisms remain elusive, and current  
therapeutic options are limited.  
A systems biology approach, leveraging high-throughput sequencing  
technologies and integrating multi-omics datasets, offers a promising avenue  
for elucidating these mechanisms.  
By examining the interplay between genes, metabolites, and signaling  
pathways, this approach aims to uncover novel biomarkers and therapeutic  
targets for HFpEF.

## Methodology

The project employs an integrative analysis of transcriptomic, metabolomic, and phosphoproteomic datasets. The transcriptomic data is preprocessed as a matrix of log2 fold changes, while metabolomic data is prepared for COSMOS input analysis. Phosphoproteomic data is presented as a matrix of negative log10 adjusted pvalues.

## 

## 

# Stattifically signigicant

### cosmosR

Analysis COSMOS was used for both signaling-to-metabolism and metabolism-tosignaling analyses. This involved preprocessing steps lie network compression and filtering to manage missing data points and partial matches between reactions. 1\. Signaling to Metabolism:(test\_for) The signaling-to-metabolism analysis aimed to connect signaling events to metabolic outcomes. The preprocessing step involved removing unexpressed nodes and pruning the network to retain only relevant interactions. 2.Metabolism to Signaling: (test\_back) The metabolism-to-signaling analysis hasthe goal to study how metabolic changes influence signaling pathways. Similar preprocessing steps were applied, with a focus on retaining nodes and interactions relevant to the metabolic data. But here I faced most technical challenges, due to the program terminatin.

Pathways such as oxidative **phosphorylation**, **glycolysis**, and **fatty acid metabolism** are interconnected and interact with other pathways involved in inflammation, **hypoxia**, and **apoptosis**.   
This could suggest that metabolic dysregulation, including impaired energy productionand altered substrate utilization, may contribute to the development and progression of HFpEF. The network highlights the role of growth factor signaling pathways in HFpEF. Pathways such as PI3K-Akt-mTOR signaling, E2F targets, and Notch signaling are interconnected and interact with other pathways involved in cell cycle regulation, apoptosis, and inflammation. This could suggest that dysregulation of growth factor signaling may contribute to the pathogenesis of HFpEF by promoting maladaptive remodeling, hypertrophy, and fibrosis. Overall, the network reveals the interplay between various biological processes, including inflammation,immune response,metabolism, and growth factor signaling, and highlights the central role of the p53 pathway in HFpEF pathogenesis

!

* the pronounced connection between **MYOD1 and GLI3** hints at a coordinated regulation of genes implicated in m**uscle development and cellular differentiation**  
* The multiple connections of **TP53 with other factors like JUN and SMAD4** underscore its involvement in a wide array of cellular processes like **DNA repair, apoptosis, and cell cycle regulation.**   
* The connection between **NR3C1, a glucocorticoid receptor, and SMAD4, a key component of the TGF-beta signaling pathway.**

* **MYOD1:** *Regulator of muscle differentiation*  
  MYOD1's interactions MEF2C and MYOG suggest a coordinated regulation of genes involved in myogenesis.   
* **TP53:** *‘Guardian of the genome’*  
  Its interactions with DNA damage response factors, cell cycle regulators, and apoptotic proteins reveal its multifaceted role in maintaining genomic stability and preventing tumorigenesis.  
* **NR3C1:** *Glucocorticoid receptor;*   
  Interactions with **SMAD4 and NFKB1**, indicate its involvement in ***immune response*****, *inflammation,* and *stress response*.** 

* **UPREGULATED PATHWAYS:**  
  * **mTORC1 signaling:** This pathway is involved in cell growth, proliferation, and metabolism.  
  * **Adipogenesis:** Markedly upregulated, indicating increased fat cell formation or lipid accumulation  
  * **MYC targets V1:** Known to regulate cell growth, proliferation, and apoptosis.  
  * **Unfolded protein response:** Indicates cellular stress response to misfolded proteins, that could be triggered by various factors, including inflammation, oxidative stress, or metabolic imbalances.  
      
* **Moderately Upregulated Pathways:**  
  * **Androgen Response, Oxidative Phosphorylation, Protein Secretion, UV Response DN:** These pathways show moderate activation.  
    

**Linear Discriminant Analysis (LDA) Effect Size (LEfSe) Analysis of Interacting Molecules**  
The LEfSe analysis reveals the most differentially active molecular interactions between two groups (e.g., disease vs. control, treated vs. untreated).

**UPREGULATED INTERACTIONS:**

* **PLA2G2A & EGFR:**   
  PLA2G2A **(Phospholipase A2 Group IIA)** and EGFR **(Epidermal Growth Factor Receptor)** suggests **increased release of arachidonic acid**, a precursor to inflammatory molecules, and potential dysregulation of cell growth and survival pathways.  
    
* **GNA12 & PTPRU:**  
  GNA12 **(G Protein Subunit Alpha 12\)** and PTPRU **(Protein Tyrosine Phosphatase Receptor Type U).** This dysregulation may contribute to processes such as wound healing, immune response, and tumor metastasis.  
    
* **ADM & RAMP1:**  
  ADM **(Adrenomedullin)** and RAMP1 **(Receptor Activity-Modifying Protein 1\)** implies augmented adrenomedullin signaling, potentially leading to beneficial **vasodilatory effects in conditions like heart failure or hypertension.**  
    
* **HLA-DPB1 & CD4:**   
  HLA-DPB1 **(Human Leukocyte Antigen \- DPB1)** and **CD4** (Cluster of Differentiation 4\) indicate a potentially **enhanced immune response.** 


**DOWNREGULATED INTERACTIONS:**

* **EFNB2 & GRM1:**  
  EFNB2 **(Ephrin B2)** and GRM1 **(Metabotropic Glutamate Receptor 1\)** may result in **impaired cell-cell communication and synaptic plasticity.**   
* **CCL21 & ACKR4:**   
  CCL21 **(C-C Motif Chemokine Ligand 21\)** and ACKR4 **(Atypical Chemokine Receptor 4\)** suggest **dysregulated immune cell migration and accumulation**.   
* **YBX1 & NOTCH1:**   
  YBX1 **(Y-Box Binding Protein 1\)** and NOTCH1 **(Notch Receptor 1\)** might affect cell fate decisions and differentiation.  
* **RGMA & BMPR1B:**   
  RGMA **(Repulsive Guidance Molecule A)** and BMPR1B **(Bone Morphogenetic Protein Receptor Type 1B)** may impair neuronal development and axonal guidance.

**UPREGULATED PATHWAYS:**

* **mTORC1 signaling:**   
  This pathway is involved in cell growth, proliferation, and metabolism.  
  * **Adipogenesis:**  
    Markedly upregulated, indicating increased fat cell formation or lipid accumulation.  
  * **MYC targets V1:**  
    MYC is known to regulate cell growth, proliferation, and apoptosis.  
  * **Unfolded protein response:**  
    Cellular stress response to misfolded proteins. This could be triggered by various factors, including inflammation, oxidative stress, or metabolic imbalances.  
    

**MODERATELY UPREGULATED PATHWAYS:**

* **Androgen Response, Oxidative Phosphorylation, Protein Secretion, UV Response DN:** 


* **Central Role:** p53 pathway  
* **THERAPEUTIC TARGETS:** INFLAMMATION, METABOLISM  
    
    
    
    
    
  

The network highlights the involvement of inflammatory and immune responses in HFpEF. 

Pathways such as **TNF-alpha signaling via NF-kB, interferon alpha/gamma response**, and inflammatory response are **interconnected** and interact with other **pathways involved in oxidative phosphorylation, hypoxia, and apoptosis.** 

* ***This could suggest that chronic inflammation and immune activation may contribute to the pathogenesis of HFpEF by triggering oxidative stress, impairing mitochondrial function, and promoting cell death.***

Pathways such as **oxidative phosphorylation, glycolysis, and fatty acid metabolism** are **interconnected** and interact with other pathways involved in **inflammation, hypoxia, and apoptosis.** 

* ***This could suggest that metabolic dysregulation, including impaired energy production and altered substrate utilization, may contribute to the development and progression of HFpEF.***


The network **highlights** the role of **growth factor signaling pathways** in HFpEF. Pathways such as **PI3K-Akt-mTOR signaling,** **E2F targets,** and **Notch signaling** are interconnected and interact with other pathways involved in **cell cycle regulation, apoptosis, and inflammation.** 

* ***This could suggest that dysregulation of growth factor signaling may contribute to the pathogenesis of HFpEF by promoting maladaptive remodeling, hypertrophy, and fibrosis.*** 

Overall, the network reveals the interplay between various biological processes, including *inflammation*, *immune response*, **metabolism**, and growth factor signaling, and highlights the central role of the p53 pathway in HFpEF pathogenesis.
