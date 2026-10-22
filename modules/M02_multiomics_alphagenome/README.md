# M02: Multiomics with AlphaGenome

**Time:** 09:45-10:30 | **Type:** Hands-on

From DNA sequence to regulatory activity, gene expression, chromatin architecture and non-coding variant prioritization with AlphaGenome.
Case study: oncogenic non-coding variants at the *TAL1* locus in T-cell acute lymphoblastic leukemia (T-ALL).

| Notebook | Content |
|----------|---------|
| `notebooks/00_setup.ipynb` | Install, API key, smoke test |
| `notebooks/01_explore_locus_TAL1.ipynb` | Exercise 1: multi-modal exploration of the TAL1 locus |
| `notebooks/02_score_variant_TAL1.ipynb` | Exercise 2: scoring and validating a TAL1 enhancer variant |

_Content in preparation (see milestone "M02 AlphaGenome")._

## Learning goals

By the end of this module, participants will be able to:

- Explore the TAL1 locus using multiple types of AlphaGenome output.
- Relate regulatory activity and chromatin information to gene expression.
- Score and evaluate a non-coding variant in the context of a candidate enhancer.

## Session plan

| Activity | Time |
|---|---:|
| Theory and AlphaGenome overview | 15 min |
| Exercise 1: Explore the TAL1 locus | 12 min |
| Exercise 2: Score and validate a variant | 12 min |
| Discussion | 6 min |
| **Total** | **45 min** |

## Case study: TAL1 MuTE insertion

### Biological context

TAL1 is a transcription factor that is frequently dysregulated in T-cell acute lymphoblastic leukemia (T-ALL). This case study examines a non-coding insertion near TAL1 that was reported in Jurkat cells and is associated with the formation of an oncogenic enhancer.

### Reference assembly and gene annotation

- **Reference assembly:** GRCh38/hg38.
- **Locus window used for visualization:** `chr1:47209255-47242023` (AlphaGenome example workflow; TAL1 is on the negative strand).
- **Gene annotation:** GENCODE v46, using the protein-coding MANE Select transcript.
- **TAL1 identifiers:** `ENSG00000162367`; MANE Select transcript `ENST00000691006.1`, corresponding to RefSeq `NM_001290403.2`.
- **Gene coordinates cross-check:** `chr1:47216290-47232335`, negative strand. This gene span is distinct from the larger locus window used in the example workflow.

The GENCODE v46 annotation and MANE Select transcript are used to display the gene model in the analysis. The locus window provides context for visualizing the variant and nearby regulatory regions.

### Main variant

The primary variant is the Jurkat MuTE insertion reported by Mansour et al. (2014):

- **Variant:** `chr1:47239296 C>CCGTTTCCTAACC` (GRCh38/hg38)
- **Inserted sequence:** `CGTTTCCTAACC` (12 bp)
- **Source:** Mansour et al. (2014), PMID: [25394790](https://pubmed.ncbi.nlm.nih.gov/25394790/)

The initial `C` in the alternate allele is the reference base used to anchor the insertion representation.

### Matched synthetic control

For a simple matched comparison, use a synthetic 12-bp insertion at the same position:

- **Control variant:** `chr1:47239296 C>CGATCCAGTACGA`
- **Inserted control sequence:** `GATCCAGTACGA` (12 bp)

This is a length- and position-matched synthetic control for the exercise. It is not an experimentally observed allele and should not be interpreted as a validated biologically neutral variant. A single control illustrates the comparison; multiple matched random insertions would provide a better background distribution.

### Biological contexts

The exercise compares a hematopoietic context with an unrelated tissue context, with T cells as an optional lineage context:

| Context | Ontology term | Description | Use |
| --- | --- | --- | --- |
| Hematopoietic progenitor | `CL:0001059` | Common myeloid progenitor, CD34-positive | Matches the CD34-positive hematopoietic context used in the AlphaGenome TAL1 example |
| T cell | `CL:0000084` | T-cell | Optional comparison with the lineage of origin of T-ALL |
| Liver | `UBERON:0002107` | Liver | Unrelated tissue comparison |

All three terms were checked against the AlphaGenome human output metadata (`model.output_metadata`). Each has RNA-seq, DNase and histone ChIP-seq tracks (including H3K27ac and H3K4me3). None has contact map tracks.

### Recommended background readings

- Mansour MR, Abraham BJ, Anders L, et al. **An oncogenic super-enhancer formed through somatic mutation of a noncoding intergenic element.** *Science*. 2014;346(6215):1373–1377. (https://doi.org/10.1126/science.1259037)
- Liu Y, Li C, Shen S, et al. **Discovery of regulatory noncoding variants in individual cancer genomes by using cis-X.** *Nature Genetics*. 2020;52(8):811–818. (https://doi.org/10.1038/s41588-020-0659-5)
- Smith C, Goyal A, Weichenhan D, et al. **TAL1 activation in T-cell acute lymphoblastic leukemia: a novel oncogenic 3′ neo-enhancer.** *Haematologica*. 2023;108(5):1259–1271. (https://doi.org/10.3324/haematol.2022.281583)

### Analysis reference

The variant representation, locus window, GENCODE v46 annotation workflow, and CD34-positive progenitor context are based on the [AlphaGenome TAL1 example workflow](https://github.com/google-deepmind/alphagenome/blob/main/colabs/v1/example_analysis_workflow.ipynb).
