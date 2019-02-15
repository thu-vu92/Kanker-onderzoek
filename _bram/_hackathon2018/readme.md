# MoMa 

Multi-Omic Method aggregator 


# Background..


Existing methods and frameworks


* effectivity immuno for mela/lunga/..
* difference between metastasis/primair for ..
* effectivitiy chemo for mela/lung/..
* most important phenotypic data..
* multi omic analysis
 ** pathways
 ** drivers per layer

 Output:
 * biomarkers
 * diagnostic profiles
 * pathways
 * methods to find pathways
 * library to perform multi-omic analysis/prediction


# Data sets

We have two collections of data sets, one for melanoma tumors and one for lung cancer 
tumors. We will focus our efforts on the melanoma tumors. Although we will 
produce melanoma-specific results, the underlying methods should be equally applicable to the 
lung cancer data set as the data schema's are identical.

For both the melanoma and the lung cancer data set we have:
* DNA, mutations: genomics
* DNA, Copy number variations: genomics
* DNA, Methylation: epigenomics
* RNA, genomic expression: epigenomics
* RNA, miRNA: transcriptomics
* Proteins: proteomics

## DNA, [Mutation](https://ghr.nlm.nih.gov/primer/mutationsanddisorders/possiblemutations)

Literally, per genome and chromosome the change in the pair compared 
to a normal reference. Remember we have (Adenine,Thymine) and (Guanine,Cytosine) as the base pairs.

The types of mutations include (taken [from here]((https://ghr.nlm.nih.gov/primer/mutationsanddisorders/possiblemutations)):

**Missense mutation**: This type of mutation is a change in one DNA base pair that results in the substitution of one amino acid for another in the protein made by a gene. 

**Nonsense mutation**: is also a change in one DNA base pair. Instead of substituting one amino acid for another, however, the altered DNA sequence prematurely signals the cell to stop building a protein. This type of mutation results in a shortened protein that may function improperly or not at all.

**Insertion**: An insertion changes the number of DNA bases in a gene by adding a piece of DNA. As a result, the protein made by the gene may not function properly.

**Deletion**: A deletion changes the number of DNA bases by removing a piece of DNA. Small deletions may remove one or a few base pairs within a gene, while larger deletions can remove an entire gene or several neighboring genes. The deleted DNA may alter the function of the resulting protein(s).

**Duplication**: A duplication consists of a piece of DNA that is abnormally copied one or more times. This type of mutation may alter the function of the resulting protein.

**Frameshift mutation**: This type of mutation occurs when the addition or loss of DNA bases changes a gene's reading frame. A reading frame consists of groups of 3 bases that each code for one amino acid. A frameshift mutation shifts the grouping of these bases and changes the code for amino acids. The resulting protein is usually nonfunctional. Insertions, deletions, and duplications can all be frameshift mutations.

**Repeat expansion**: Nucleotide repeats are short DNA sequences that are repeated a number of times in a row. For example, a trinucleotide repeat is made up of 3-base-pair sequences, and a tetranucleotide repeat is made up of 4-base-pair sequences. A repeat expansion is a mutation that increases the number of times that the short DNA sequence is repeated. This type of mutation can cause the resulting protein to function improperly.

### DATA FIELDS, shape (422553, 11)
```ID      |  Location        | Change     |  Gene   | Mutation type|  Var.Allele.Frequency  | Amino acid```

```SampleID,| Chr, Start, Stop|  Ref, Alt  | Gene    |    Effect    |  DNA_VAF, RNA_VAF      | Amino_Acid_Change```

```string   |string, int, int | char, char | string  |    string    |  float, float          |  string```

NOTE: this gives us direct insight in how genetic mutations lead to changes in amino-acids.

## Copy Number Variations

A copy number variation (CNV) is when the number of copies of a particular gene varies from one individual to the next.

### DATA FIELDS, shape (24802, 372)
```Gene      | Chr, Start, Stop | Strand     |   SampleID 1..SampleID N```

```string    |string, int, int  | int        |  int..int```


## Methylation, gene expression regulation

Degree of [methylation](https://en.wikipedia.org/wiki/DNA_methylation)
indicates addition of Methyl groups to the DNA. Increased methylation is associated with less transcription of the DNA:
Methylated means the gene is switched OFF, Unmethylated means the gene is switched ON.

Alterations of DNA methylation have been recognized as an important component of cancer development.


### DATA FIELDS, shape (485577, 483) 
```probeID   | Chr, Start, Stop | Strand  | Gene   |  Relation_CpG_island | SampleID 1..SampleID N```

```string    |string, int, int  | int     | string |   string             | float..float```


## RNA, gene expression

Again four building blocks; Adenosine (A), Uracil (U), Guanine (G), Cytosine (C).

(DNA) --> (RNA)

A --> U 

T --> A

C --> G

G --> C

Gene expression profiles, continuous values resulting from the normalisation of counts.

### DATA FIELDS, shape (60531, 477)
```Gene      | Chr, Start, Stop | Strand  | SampleID 1..SampleID N```

```string    |string, int, int  | int     |  float..float```


## miRNA, transcriptomics

The connection between the RNA production and protein creation. I.e. perhaps miRNA expression values can be associated with specific proteins.

### DATA FIELDS, shape (2220, 458)
```MIMATID  | Name   | Chr, Start, Stop | Strand  | SampleID 1..SampleID N```

```string   | string |string, int, int  | int     |  float..float```


## Proteomes

Proteine expression profiles, ditto, continuous values resulting from the normalisation of counts


### DATA FIELDS, shape (282, 355)
```ProteinID  | SampleID 1..SampleID N```

```string     | float..float```

### QUIZ, identify our data sets in the following image!


![image.png](_images/overview.png)


## GOAL

Some degree of multi-omic analysis and identification of pathways.

![image.png](_images/multi_omic.png)

# Targets

First we need get a picture of what a "signature" actually means in this context. We basically have hierarchically dependent data with "pathways" going through those layers, those pathways are connected by mutations on the one end (DNA) and proteins on the other end. How to find those pathways is the main question, because once we can do that, we only have to identify either which pathways are typical for people that do or do not respond well to immunotherapy, or what part of the pathway is typically different for those patients.

So what is a pathway? A pathway is a chain of molecular changes that leads to, in our case, the production of certain proteines, 
or (since we don't have many proteomic measurements) certain RNA codes. In the simplest form it is a sequence, but it 
is more likely similar to a bi-directed graph, in it's simplest form; DNA mutation <--> RNA <--> mRNA <--> proteines.
The collection of molecular regulators that govern the genomic expression levels of mRNA and proteins is called the 
[gene regulatory network (GRN)](https://en.wikipedia.org/wiki/Gene_regulatory_network). So, instead of chain, it is better
to say a network of molecular changes.

![GRN](_images/Gene_Regulatory_Network_2.jpg)

![GRN](_images/Gene_Regulatory_Network.jpg)

![GRN](_images/DG_Network_in_Hybrid_Rice.png)

Specifically, given that we find a pathway, the genomes that, given a mutation, will lead to a proto-oncogenic or an inhibiting effect on tumor development
are denoted as proto-oncogenics and inhibitors, the former promotes tumor growth, the latter slows it down.

One such pathway is the Mitogen Activated Protein Kinase (MAPK) pathway. This pathway connects certain mutations of the
[BRAF](https://en.wikipedia.org/wiki/BRAF_(gene)) oncogene (i.e. DNA) to the generation of certain proteins that lead to the promotion of cell growth. This 
[MAPK pathway](https://en.wikipedia.org/wiki/MAPK/ERK_pathway) looks as follows:

![MAPK](_images/MAPKpathway_diagram.svg.png)

Just to show the complexity, the [PI3K/AKT/mTOR pathway](https://en.wikipedia.org/wiki/PI3K/AKT/mTOR_pathway) looks as follows:
![PI3k](_images/m836px-MTOR-pathway-v1.7.svg.png)

![mind blown!](_images/mind_blown.jpg)

Don't worry, this is not expected from us..although I can produce this stuff in paint, hands down (I can actually use my 
nose to draw..). Anyways, back to reality: we only have a few thousand protein measurements, and we do not have any
time series data so it is practically impossible to extract any feedback effects from downstream changes. 
Ooofff, that leaves us with a top down approach, from instruction/mutation to RNA and in some cases, proteins.
Finding any feedback effect is secondary, and perhaps for continued work after the hackathon.

Now, given such a pathway we can frustrate the signal anywhere on the chain, as long as it prevents
the cell growth stimulation.

We should keep in mind that the inhibitors/proto-oncogenes are likely specific to
the type of melanoma's, we distinguish at least the following by their genomic mutations:
* (proto-oncogenes) BRAF wild-type
* Triple Wild-type
* NF1 
* KIT
* MITF
* RAS
* (inhibitor) PD-1/PD-L1

The current inhibitor, the one they likely use in the immunotherapy is PD-L1. It is called a checkpoint inhibitor, retrieving this 
from our models would be a good validator. In fact, PD-L1 itself is a proto-oncogene, but apparently
it can be inhibited. So, instead of searching for inhibitors specifically, we should be looking for
proto-oncogenes.


The main questions:
* Why do some patients respond to immunotherapy and others not? --> can we predict who will not respond?
* What are the pathways related to melanoma (firstly), to immunotherapy response (secondly)


"Official" supporting questions:
* Can you integrate all the data to make more accurate predictions for each patient than you would by only looking at one data source?
* As melanoma is a set of diverse diseases, can you stratify the patients based on all the data in to subgroups?
* Can you show and visualize the correlations and concepts between the different datasets?
* Can you select a list of most informational variables that drive the predictions?
* Can you select a list of most informational variables distinctive for each patient subgroup?

* What is the difference between patients that do and do not respond
* Can you identify a signature based on an integrative approach that can predict response to immunotherapy?
* Can you identify a signature that correlates with the prognosis of immunotherapy?

Basic hypotheses that would be nice to confirm (i.e. nice to haves, feel free to ignore, but don't ;))
* T(tumor), increased Bresow-thickness correlates with more malignancy (Tis, T1a/b, T2a/b, T3a/b, T4a/b), i.e. decreasing survival rate
* N(nodal stage), Local spread correlates with more malignancy (N0, N1a/b, N2a/b/c, N3)
* M(metastasis location), distant metastasis (beyond regional lymph nodes)  corresponds with mmore malignancy (M0, M1)
* BRAF proto-oncogenic mutations should occcur in about 50% of all cutaneous melanomas.
* Gene, if inactivated by mutation then the deactivation by the methylation will be less prominent
* CNV says something about visibility for the immune system
* Spink5 not in mestastised tumor
* We should be able to identify 4 subtypes of cutaneous melanomas: BRAS/RAS(N/H/K)/NF1/Triple-WT
* order 3 clusters in the mRNA profiles of the most variant genes (keratin, immune, MITF-low)
* inhibitor: PD-L1/PD-1, our method should be able to retrieve this specific mutation as an inhibitor
* inhibitor: MEK, our method should be able to retrieve this specific mutation as an inhibitor for BRAF wild-type/NF1 mutant melanoma's
* inhibitor: PTEN/TP53/APC, our method should be able to retrieve this specific mutation as an inhibitor
* inhibitor: CTLA4, APRIL
* proto-oncogenic: BRAF, our method should be able to retrieve this specific mutation as a proto-oncogene
* LCK protein expression: correlates positively with patient survival
* genetic markers for melatonine may be proxy for higher risk of melanoma
* proto-oncogenes with higher copy number (in absolute sense) are more visible for the immune system
* Are the strongest proto-oncogenes mutually exclusive: BRAF, NRAS, MAP2K1, KIT, CTNNB1, GNA11, GNAQ?



Basic information:
* Moles spatially near eachother, in combination with discolouring is indicative for a higher likelihood
of metastasis 
* (Breslow depth) under skin correlates with more malignancy 
* RNA, each 3 letter combo is associated with a specific amino-acid --> an amino-acid is associated with multiple 
combinations --> i.e. not every change in RNA coding leads to a delta amino acid.
* mutations may lead to a change in protein function, and a change in genetic expression of other genes
* The multi-omic information we have in our hands says nothing about the environmental factors influencing it..but the clinical 
information may give us a hint; age is likely important, also, someone who has stage 4 cancer is likely to have a significantly 
different level of immune function.

Questions to people from Erasmus:
* can we get mappings from genomes to genome groups ? I added a folder to the google drive: mapping_data, which has
a mapping information from RNA probesets to genomes.
* how can we couple miRNA to proteomics proteine/miRNA-wise (so not sample wise..)?

# Things biologists like

Hierarchical cluster diagrams, linked graphs, flow diagrams and simple tables/heatmaps with 
the most important genomes:

![Heat maps!](_images/ex_plot.png)

![Heat maps!](_images/ex_plot4.png)

![Importances](_images/ex_plot2.png)

![Importances](_images/ex_plot3.png)

![Circle](_images/ex_plot5.png)

![Circle](_images/ex_plot6.png)


Biologists like to understand the results, not very strange since they will base their laboratory 
work on it, medications will be derived from it and it will be applied to real patients. 
This is very important to keep in mind since it (sort of) **excludes a neural-network-only approach**.

# Suggested approaches 

Please add ideas with your name in the section header (change the readme.md file in the `_doc` folder, then use `grip readme.md`, to install just do `pip3 install grip`)

Overall, I see three paths:
1. unsupervised learning and general exploratory data analysis to identify promising target variables, plus worth while hypotheses
2. feature engineering --> transposition of tables --> dimension reduction --> normalisation --> classification --> viz
3. graph generation and identification of common paths and graph clusters per classification  --> viz

## Pre-processing

Determine feature importances bottom-up: from the proteins->RNA->DNA
Determine feature importances top-down: from the DNA->RNA->proteins

Omics --> sub-omics
*	split the omics 
*	apply tf-idf of sorts to all count-based features
*	normalise if appropriate 
*	reduce if necessary (preferably keeping the features intact)
* 	find self-similar features and combine

## Clusters per layer

For the non-graphs, use some density-based clustering algorithm like HDBSCAN and 
lower dimension embedding like t-SNE

For the graphs, assuming we can construct them we can try to find 
* communities
* exemplars 
* cliques

Suggested algorithms/tools are :
* Sparse Affinity Propagation, for exemplars and communities
* Markov Clustering for cliques
* t-SNE (in sklearn but not hierarchical): [multicore](https://github.com/DmitryUlyanov/Multicore-TSNE), [multicore2](https://github.com/danielfrg/tsne)
* instead of t-SNE use UMAP
* HDBSCAN (you can find that [here](http://hdbscan.readthedocs.io/en/latest/how_hdbscan_works.html))

Interesting candidates for clustering are:
* similarity between miRNA and proteine
* similarity between RNA and CNV
* self-similarities

## Dimension reduction

I would suggest the golden oldies, because they work :D
* PCA
* LDA
* FDR with ANOVA

If anyone can whip up an autoencoder that we would be cool but likely the above methods will do fine.
For the technical jury we should at least mention that we looked at using autoencoders though.

We can also use the embedding technique UMAP that, unlike t-SNE, preserves the actual distances in higher dimensions.


## Classifications

This should be easy to do. First we should define the targets that are relevant to our end goal,
which is to recognize pathways, and inhibitors on those pathways. I.e. we need to be able to 
predict the **level of malignancy**, the **survival rate** and the **response** to immunotherapy.

We choose (for the sake of time):
*   therapy response: yes/no
*   tumor type: primary/metastasis
*   tumor stage: ```T0, T1, T2, N0, N1, N2```

This generates weights/importances per feature and gives the predictive power of each layer.

A novel thing to do is to apply a tool like Quiver to visualise which genomes are important 
per classification. For this we would need to transform the 1-dimensional tensors 
into 2-dimensional tensors. This would only work per patient, but as transparency 
is one of the key-ingredients of personalised medicine, and in my view of ML applied to
health care in general, it would be a nice touch. 

Suggested algorithms:
* lightGBM, boosting type = GBDT, or DART, or GOSS in case of scaling problems
* CNN in Keras, we already have something laying around from last year

In general I see two approaches here:
* classification per layer plus stacking of the classifiers, by whatever means.
* classification of the merged layers [CNV, mutation, methylation and RNA expression]+[miRNA, proteins], where the layers are filtered a priori by whatever means.

In either case, greedily collect the most important features and visualise the overlap .

### Combining the layers

Below we describe a concatenated omic integration
1. **model based**: proba's per layer as input for a new classifier 
2. **model based nt**: weighted average of the proba's using prediction uncertainty and CV accuracy, or majority vote
3. **reduced concatenated**: merge different datasets based on most important features per sub-omic 
4. **model-based inter-omic**: apply heuristics-based models to connect the different sub-omics
5. **similarity-based path constructor**:  merge different sub-omics based on inter-omics similarities

Point 3. needs elaboration:
*   we need to reduce the omics sets **first**
    *   create sub-omic sets; split based on meta-descriptors: such as ```strand, Mutation Effect```
    *   reduce number of features; non-parametric statistical test over the classifications
*   Sub-omic sets can be combined first within omic's, then between omics
*   Between omics by:
    *   simply concatening columns on SampleID --> we lose the relationship between datasets?
    *   combining by model-based coupling between datasets:```f(x_mut, x_methyl, x_cvn,..)```

**These classifications can provide use with input for finding the pathways!** 
Obviously by looking at the feature importances. 

In particular this will say per layer which features are important, which for the miRNA and the protein 
data may be key for identifying the final pieces of the pathway puzzle.

Point 4. Heuristics are for instance: RNA expression as weight for mutations, methylation as weight for RNA

Point 5. Use known proto-oncogenes as starting point initially; BRAF, NRAS, NF1, Triple Wild-Type
Point 5. Search for mutually exclusive/inclusive gene's/mutations --> intra-omic self-similarity



## Correlations between different layers, 

Extract layer pairs: DNA-RNA, RNA-mRNa, mRNa-proteomics, using 
* raw features
* raw, filtered features (using variance over the different classes)
* reduced features per layer (PCA/LDA whatever)

-> generate new cross-layers that combine all possible layer pairs into single layer,
train classifier and characterise the multilayer pairs.

I.e. Correlated features link the layers pairwise, after which the layers can be connected into 
a single layer. 

## Bayesian Networks

Given potential pathways we can infer Bayesian Networks as approximations for the GRN and visualize them
with some graph viz. tool. [Watch](https://www.youtube.com/watch?v=TuGDMj43ehw),
 [Read the wiki :D](https://en.wikipedia.org/wiki/Bayesian_network)


![Connected layers](_images/Network_layers.png)

## Markov Networks

Train a Conditional Random Fields to describe the inter-omic network.


## Graphs, from the ground up 

Per patient we have a graph connecting the CNV, mutation, methylation and 
RNA expression data using (Gene, Chr, start, stop). 
When looking at the gene-connectivity (i.e. counting occurrence of chr, start, stop, amino-acid change, 
type of mutation), this graph will mostly be similar per patient in terms of the adjacency matrix 
but dissimilar in terms of the similarity matrix. This opens up some possibilities: 
We can 
* determine clusters per patient graph: exemplars, communities, cliques. Then determine cluster overlap
per target label.
* create multi-layer graph per target label, count edges (or sum edge weights), normalise edge sums. 
Flatten multi-layer graph and interpret normalised sums as edge weights. Determine characteristics clusters
per target label. I.e. (N, N, m) --> (N*, N*, 1)

The resulting clusters, and their characteristics can be used to feed a predictor. This has the benefit of 
* transparency: it is clear why a target value is predicted
* compatibility: compared to simply merging the data into one matrix we have more guarantee to obtain biologically
sound estimations
* it looks fucking nice 
Suggested tools/algo's are:
* Markov Clustering, 
* Neo4j-Bloom, networkx

# Sources

Hardware: Google collab and custom Google compute engines

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4731297/
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5496318/#MOESM8

https://www.nature.com/articles/nature13385

https://www.ncbi.nlm.nih.gov/pubmed/26091043
https://www.ncbi.nlm.nih.gov/pubmed/22960745