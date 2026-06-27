# Background
The code and data associated with this repository may be used to recreate the results from two related research products:
1.  A research presentation at the American Association of Biological Anthropology:
> Wolfe, C. A., E. Y. Chu, L. K. Corron, M. H. Price, and K. E. Stull. 2022. "Advances in subadult age estimation: using information theory to explore the relationship between growth indicators and age." AMERICAN JOURNAL OF BIOLOGICAL ANTHROPOLOGY, vol. 177, pp. 199-199.

2.  A publication submitted for review in July 2026:
> Stull, K. E., L. K. Corron, E. Y. Chu, M. H. Price, and C. A. Wolfe. 2026. "Using information theory to evaluate the predictive capability of subadult age indicators throughout ontogeny." In Review. 

All necessary data files can be found in the data folder. If one were to use the data in this work, please refer to and cite the initial source and publication:

>Stull, K. E., and L. K. Corron. 2021. “SVAD_US (1.0.0) [Data Set].” Zenodo, ahead of print. https://doi.org/10.5281/zenodo.5193208.

>Stull, K. E., and L. K. Corron. 2022. “The Subadult Virtual Anthropology Database (SVAD): An Accessible Repository of Contemporary Subadult Reference Data.” Forensic Science 2: 20–36. https://doi.org/10.3390/forensicsci2010003%2520doi.

The results of this work are a derivation of research first described in the citation below. All code associated with said work can be found at the [yada](https://github.com/MichaelHoltonPrice/yada/tree/dev) repository:

> Stull, K. E., E. Y. Chu, L. K. Corron, and M. H. Price. 2022. “Subadult Age Estimation Using the Mixed Cumulative Probit and a Contemporary United States Population.” Forensic Science 2: 741–79. https://doi.org/10.3390/forensicsci2040055.

The code associated with th

## Summary of Files

1. The `data` folder contains all requisite data and associated files.
2. The `AABA2022` folder contains initial code used in the conference presentation. Note, this code is nearly identical to that described below with changes only related to versioning in the underlying `yada` package.
3. `S1...` installs all packages necessary for the analysis.
4. `S2...` runs all analyses including prior sampling, likelihood modeling, and quantification of information gain.
5. `S3...` imports the results from `S2...` and visualizes information gain across age. Note, to fully recreate the figures, a user will need to combine individual curves from `S3... `.
6. `S4...` creates Figure 3 in the text, visually demonstrating the information gain from prior to posterior. This code was created by co-author Michael Holton Price.

To cite this software repository please include:

>Wolfe, C. A., Stull, K. E., Price, M. H., Chu, E. Y., & Corron, L. K. (2026). StullEtAl_2026_InfoGain (Version 2.0.0) [Computer software]. https://doi.org/10.5281/zenodo.20975429
