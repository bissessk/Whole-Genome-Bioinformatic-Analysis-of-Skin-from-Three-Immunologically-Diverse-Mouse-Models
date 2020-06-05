# Whole-Genome-Bioinformatic-Analysis-of-Skin-from-Three-Immunologically-Diverse-Mouse-Models
Performed a differential gene expression analysis with RNA-seq on three immunologically diverse mouse models. Differentially expressed genes were identified and a GO term enrichment analysis was done to determine their roles.

## Experimental Context
In the experiment, Whole Genome Analysis of Skin from Three Immunologically Diverse Mouse Models by
Luckett-Chastain and Gallucci, expression profiling by high throughput sequencing was done on three mouse
models that differ immunologically. The first mouse model, C57BL/6 was used to represent the properties of Th1.
The second mouse model, BALB/c, was used to demonstrate the properties of Th2. The third mouse model was in
the same family as the C57BL/6 mouse model but had its IL-6 gene knocked out. The goal of profiling each mouse
model’s differential gene expression through RNA sequencing was to highlight possible influences of the Th1/Th2
balance and provide further insight into the role of IL-6 in overall basal immunity.

## Analysis
The analysis involved processing the unique reads of the data to construct a DESeq dataset using DESeq2. After filtering out low read counts and normalizing the data, three comparison objects were created using DESeq2’s results function which were then
processed to determine lists of differentially expressed genes. After performing a GO Term Enrichment Analysis to
determine the roles of each differentially expressed gene, conclusions about the three mouse models were made.

## Conclusions
There were many findings for each mouse model. The C57BL/6 mouse model showed Th1 production to be involved
with the expression of genes that are involved with the regulation of leukocyte migration, regulation of neutrophil
chemotaxis, and the regulation of T Cell Migration. The IL-6 KO mouse model showed that when the IL-6 gene was
knocked out there was a higher expression of genes involved with the epidermis and the cell cycle. It also showed
that when knocked out, muscle-related genes are expressed more as well with transport-related genes and
metabolism related genes. The BALB/c mouse model showed Th2 production to be involved with the expression of
genes involved with Antigen processing, T cell regulation, apoptosis regulation, and peptide crosslinking.

**Dataset Used:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119202
