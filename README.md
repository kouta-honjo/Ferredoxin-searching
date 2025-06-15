# Protocol for *Moorella thermoacetica* Gene Deletion Strain Construction Support

## 1. Overview

This protocol is a semi-automated bioinformatics tool for constructing gene deletion strains of specific proteins in the acetogen *Moorella thermoacetica* ATCC 39073. By integrating bioinformatics for understanding the role of ferredoxins in various biological processes such as electron transfer, metabolic pathways, and cellular structures, with molecular biology for plasmid construction for deletion cloning, we aim to contribute to rapid Proof-of-Concept (PoC) implementation at the lab scale.

**Target Users:** Researchers and students with basic knowledge of genetic manipulation and genome analysis
**Purpose:** To support rapid PoC implementation in the laboratory

---

## 2. Pre-confirmation of Applicability

### Recommended Targets
- Microbial species belonging to the acetogenic category (refer to Rosenbaum & V. Müller (2023))
- Primary objective is the analysis of FeS cluster-containing proteins

### Application Limitations
- **Step i)** The first stage of selection is assumed to be FeS cluster-based.
- **Step ii)** Only supports plasmid construction for homologous recombination.

**→ If the above conditions are not met, this protocol may not be suitable for use.**

---

## 3. System Configuration

### Step i) FeS Cluster Analysis System
- **Function:** Clustering by FeS cluster units, search and visualization based on function.
- **Data Source:** F.P. Rosenbaum & V. Müller (2023)
- **Constraints:** Depends on metabolic pathway data of acetogens.

### Step ii) Plasmid Design Automation System
- **Function:** Automation of plasmid design for homologous recombination.
- **Constraints:** Only supports homologous recombination (other gene introduction methods are not supported).

---

## 4. Implementation Requirements

### Technical Environment
- **Recommended IDE:** VSCode
- **Python:** 3.8.8 or higher

### Prerequisite Knowledge/Skills
- Basic genetic manipulation techniques
- Experience with plasmid construction
- Python execution environment setup

### Preparation Checklist
- [ ] Genome sequence data of the target strain
- [ ] Sequence information of the target protein
- [ ] Python execution environment setup
- [ ] Installation of necessary libraries
---

## 5. Example Usage Scenarios

1.  **Construction of plasmids for metabolic modification of acetogens**
    - Clustering by FeS motif search
    - Comparative analysis with similar proteins

2.  **Rapid creation of gene deletion strains**
    - Selection of target gene
    - Generation of construct for homologous recombination

---

## 6. Methods

### i) FeS Cluster Analysis System

This system enables **clustering by FeS cluster units** and **search and visualization based on their function**.

#### 1. First Stage: Primary Screening by FeS Motif Search

**Function**: Extracts genes encoding iron-sulfur clusters and outputs them in TSV format.

**Required**: Whole genome information of the microorganism to be analyzed (please download in GenBank format).

**How to Use**:

1.  Place `protein_info_extraction.py` and the whole genome file (e.g., `Moorella_thermoacetica_ATCC_39073_complete_genome.gb`) in the same folder.

    ![](/image_fig/1_1.png)

2.  Set the file paths within the script. Change the input and output file names as appropriate.

    ![](/image_fig/1_2.png)
3.  Execute the python script in bash.
    ```bash
    python protein_info_extraction.py
    ```

4.  After execution, confirm that the primary screened TSV file (e.g., `Moorella_protein_feature.tsv`) has been outputted. Locus_tag, product_protein_id, and translation are outputted, but you can change the code as needed.

    ![](/image_fig/1_content.png)

#### 2. Second Stage: Sorting by FeS Cluster Sequence Motif

**Function**: Sorts only genes encoding FeS cluster sequence motifs from the information extracted in the first stage.

**Required**: Output file from the first stage (e.g., `Moorella_protein_feature.tsv`).

**How to Use**:

1.  Place `FeS_screening.py` in the folder set as your current directory.

2.  Set the file paths within the script. Change the input and output file names as appropriate.
    ![](/image_fig/2_2.png)

3.  Execute the python script in bash.
    ```bash
    python Fes_screening.py
    ```

4.  After execution, confirm that a TSV file (e.g., `Moorella_FeS_annotated.tsv`) containing only sequences encoding FeS clusters has been outputted.

    ![](/image_fig/2_result.png)

#### 3. Third Stage: Functional Screening

**Function**: From the information extracted in the second stage, sorts proteins based on function, such as metabolism, cellular organization, and coenzymes. The sorted protein information is outputted as a TSV file, and simultaneously, each sequence including 10 surrounding genes is outputted in gb format for visualization.
**Required**: Output files from the second stage (e.g., `Moorella_protein_feature.tsv`, `Moorella_FeS_annotated.tsv`).

**How to Use**:

1.  Place `Function_based_screening.py` in the folder set as your current directory.

2.  Set the file paths within the script. Change the input and output file names as appropriate.
    ![](/image_fig/3_2.png)

3.  Execute the python script in bash.
    ```bash
    Function_based_screening.py
    ```

4.  After execution, a list of searchable functions will be displayed. Enter the desired function in the comment box below.

    ![](/image_fig/3_3.png)

5.  A file "Moorella_Fd_related_to_Ferredoxin.tsv" containing the functional information of the sorted genes will be outputted.
    ![](/image_fig/3_result1.png)

5.  A file "Moorella_Fd_related_to_Ferredoxin.tsv" containing the functional information of the sorted genes will be outputted.
    ![](/image_fig/3_result2.png)
```