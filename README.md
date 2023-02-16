# Valega-Mackenzie, Wencel et al. (2023) (CanVaxDOpt_Model)

The Jupyter Notebook, CanVaxDOpt, contains the code to reproduce the graphs of activated T cells (Fig 3 and 4), tumor diameter (Fig 5), optimal vaccine doses (Fig 6 and 7) and the vaccine dose-response curves (Fig 8 and 9) from the main paper. The easiest way to reproduce the results is to run each cell from the top down to avoid unexpected errors. Running the cells to find the optimal solutions can take approximately 20-30 minutes. All other cells should run relatively quickly. 

The excel file, Patients_HLA_BA, has the binding affinities of MHC-I/II to specific alleles associated to the set of immunogenic peptides of each patient. These binding affinities were used to compute the dissociation constant, $K_{D,j}$. All the necessary information in the spreadsheet has been merged to the CanVaxDOpt Notebook.
