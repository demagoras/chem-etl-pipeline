# Cheminformatics Data Validation & Screening Pipeline

A cheminformatics pipeline built using **RDKit** and **Pandas**. It automates raw input, canonization, sanitization, deduplication, and physicochemical filtering of SMILES strings using **SMARTS**, isolating drug-like chemical structures based on **Lipinski’s Rule of Five** while flagging aromatic features and sub-structures. A clean, analysis-ready dataset is produced alongside a complete exclusion log documenting all rejected compounds.

---

## Pipeline Architecture

The script processes data sequentially through a four-tier filtering architecture to ensure data integrity before analytical visualization:

```
[ Raw SMILES Input ]
         │
         ▼
 1. SMILES Sanitization      ──► [ Fails Parse ] ──► (Exclusion Log)
         │
         ▼
 2. Canonical Deduplication  ──► [ Duplicate ]   ──► (Exclusion Log)
         │
         ▼
 3. Lipinski Ro5 Filtering   ──► [ >1 Violation] ──► (Exclusion Log)
         │
         ▼
 4. Substructure Profiling   ──► [ Valid Compounds ] ──► [ PNG Image Generation ]
```
## Features

### 1. SMILES Validation

Invalid molecules are automatically detected and moved to an exclusion log.

Examples:
* Invalid syntax
* Unclosed rings
* Impossible valence states
* Non-SMILES input

### 2. Structure Standardization & Deduplication

Molecules are converted to canonical SMILES representations before comparison. Any duplicates are detected and only one copy is retained, the other moved to the exclusion log.

Example:
```text
CCO
OCC
```

### 3. Functional Group Detection

The pipeline performs substructure searching using **SMARTS**.

Current implementation:

* Phenol detection

```python
c1ccccc1[OH]
```

Additional SMARTS patterns can be easily added.

### 4. Nitrogen Analysis

Each molecule is scanned for total nitrogen count and their formal charges.

### 5. Physicochemical Property Calculation

Using RDKit descriptors:

| Property         | Description                         |
| ---------------- | ----------------------------------- |
| MW               | Molecular Weight                    |
| LogP             | Octanol-water partition coefficient |
| H-Bond Donors    | Lipinski donor count                |
| H-Bond Acceptors | Lipinski acceptor count             |

### 6. Drug-Likeness Filtering

Compounds are evaluated using **Lipinski's Rule of Five**, so molecules with more than one violation are excluded and recorded separately.

### 7. Molecular Visualization

Accepted molecules are rendered as a grid image.

Features:
* Automatic 2D coordinate generation
* Aromatic ring highlighting
* Per-molecule annotations showing:
  * Molecular weight
  * Phenol detection
  * Nitrogen count

Example output:

```text
MW: 180.16 | Phenol: NO | N: 0
```

---

## Example Output

**Note**: only three molecules are shown in each table.

### Clean Dataset

| SMILES                | MolWt  | LogP | H-Donors | H-Acceptors | Lipinski_Violations | Contains_Phenol | Nitrogen_Count |
| --------------------- | ------ | ---- | -------- | ----------- | ------------------- | --------------- | -------------- |
| CC(=O)Oc1ccccc1C(=O)O | 180.16 | 1.31 | 1        | 3           | 0                   | NO              | 0              |
| Oc1ccccc1             | 94.11  | 1.39 | 1        | 1           | 0                   | YES             | 0              |
| CC(=O)Nc1ccccc1       | 135.17 | 1.64 | 1        | 1           | 0                   | NO              | 1              |

### Exclusion Log

| Identifier          | Pipeline Stage | Reason                       |
| ------------------- | -------------- | ---------------------------- |
| NOT_A_SMILES_STRING | Sanitization   | Invalid SMILES Syntax        |
| OCC                 | Deduplication  | Duplicate Structure          |
| Cyclosporine        | Filtering      | Failed Lipinski Rule of Five |

### Image ``

<img src="https://github.com/demagoras/chem-etl-pipeline/blob/main/example_grid.png" width="650" height="650">

---

## Technologies Used

* Python (in Jupyter)
* RDKit
* Pandas

To run install required dependencies:

```bash
pip install pandas rdkit
```

---

## Future Improvements

Potential extensions:

* Implementation of PAINS
* Removal of salt forms (important for pharma)
* Molecular fingerprints (Morgan/ECFP)
* Expand SMARTS patterns beyond phenols
* SDF/CSV export
* Machine learning (QSAR)
* Interactive dashboard with user input
* Performance optimisation
* Similarity search and clustering

---

## Why This Project Matters

Raw datasets sourced from external vendors or public databases frequently contain invalid SMILES strings, duplicate structures, complex non-drug-like molecules. Feeding this raw, uncurated data into a QSAR model or machine learning pipeline leads to a classic GIGO problem.

This project simulates the crucial pre-processing and curation stage of data, demonstrating how this raw and messy data can be automatically and very quickly transformed into structured data.

In any production pipeline, this automated standardization is essential for ensuring data integrity and optimising downstream resources to ensure model readiness.
