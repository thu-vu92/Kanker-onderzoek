# Continuation..

Main research questions:
* (how) can we identify patient who will develop metastasis?
* (how) can we predict the efficacy of immunotherapy for melanoma and lung cancer
* (how) can we predict the efficacy of chemotherapy for melanoma and lung cancer



## Background cont.'d

Tumors try to fool the immune system but adopting ligands that mimick
* co-inhibitory receptors
* co-stimulatory receptors

..and create 
* cytokines/chemokines that make the immune system think that it is already 'treated' by T-cells

Cancer cells evolve: they change because they are attacked by the chemo/immuno

Going from DNA mutation to methylation to (mi)RNA expression to proteomics is less and less susceptible to noise from externalities.


# Paper 1

* Heatmap with {samples (grouped by xx), methylation (agg'ed)}
* Visualisation of phenotypical groups using dimensionality reduction
* Methylation versus RNA expression
* cohort-bias removal (ask for raw data from Rogier), use [Combat](https://www.bu.edu/jlab/wp-assets/ComBat/Usage.html), and the L/S method
* check for gezonde samples (zie ook 1000genome project)

Otherwise: 
* document/knowledge sharing

# Paper 2

* Simple multi-omic predictor: omic-wise normalisation and dimensionality reduction