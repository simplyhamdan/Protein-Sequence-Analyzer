# Protein Sequence Analyzer

A Python-based tool for analysing protein sequences and calculating basic physicochemical properties.

## Features

* Protein sequence validation
* FASTA and plain-text sequence input
* Multiline FASTA support
* Amino acid counting and composition analysis
* Protein sequence length calculation
* Molecular weight estimation
* Isoelectric point (pI) estimation
* Average hydropathy index calculation
* Hydrophobic/hydrophilic interpretation
* Amino acid group composition
* Acidic/basic charge tendency

## How to Run

Make sure Python 3 is installed.

Run the program:

```bash
python main.py
```

Enter a protein sequence in plain-text or FASTA format.

Press **Enter on an empty line** when you have finished entering the sequence.

## Example Input

```text
>my_protein
ACDEFGHIKL
MNPQRSTVWY
```

## Example Output

```text
Valid protein sequence!
Length of the sequence: 20
Molecular Weight: 2738.01 Da
Isoelectric Point: 7.11
Hydropathy Index: -0.49
Hydropathy Interpretation: Hydrophilic
Charge Tendency: Basic

Amino Acid Counts:
A: 1 (5.00%)
C: 1 (5.00%)
D: 1 (5.00%)
...
Y: 1 (5.00%)

Amino Acid Groups:
Nonpolar: 11 (55.00%)
Polar: 4 (20.00%)
Acidic: 2 (10.00%)
Basic: 3 (15.00%)
```

## Technologies

* Python 3
* Functions
* Dictionaries
* Lists and sets
* Loops
* Conditional statements
* Basic numerical calculations
* FASTA sequence handling

## What I Learned

This project provided practice in Python programming while applying concepts from bioinformatics and molecular biology.

Key concepts included:

* Working with biological sequence data
* Sequence validation and cleaning
* Amino acid composition analysis
* Protein physicochemical properties
* pKa and charge calculations
* Binary search for pI estimation
* Hydropathy analysis
* FASTA format processing
* Structuring a Python project using functions

## Notes

The physicochemical calculations use simplified models and approximate values for educational purposes. They are not intended to replace specialized bioinformatics tools.

## Project Series

Part of a series of progressively more advanced bioinformatics projects.

**Project 1:** DNA Sequencing Toolkit
**Project 2:** Protein Sequence Analyzer
