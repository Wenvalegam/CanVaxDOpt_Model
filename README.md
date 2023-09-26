# Dose optimization of an adjuvanted peptide-based personalized neoantigen melanoma vaccine (2023)
Authors: Wencel Valega-Mackenzie, Marisabel Rodriguez Messan, Osman N. Yogurtcu, Ujwani Nukala, Zuben E. Sauna, and Hong Yang 
# Dose Optimization of an Adjuvanted Peptide-Based Personalized Neoantigen Melanoma Vaccine (2023)

**Authors**: Name1, Name2, ...

## Overview

This repository provides the code to reproduce the numerical results presented in our paper [insert citation here]. It includes a Jupyter Notebook named CanVaxDOpt, which guides users through the process of simulating optimal vaccine doses, immune responses, and tumor responses for each patient.

### Key Features:

- Numerical solution of the optimization problem using the Forward Backward Sweep method.
- Selection of the optimal peptide dose with the most clinical benefit.
- Reproduction of figures from our paper: Optimal vaccine doses (Fig 5), Activated T cells (Fig 6), Tumor diameter (Fig 7).
- Comparison of clinical benefits using the predicted optimal vaccine dose against clinical trial doses or other tested doses.
- Tables presenting optimal vaccine dose concentrations of peptides and adjuvants per vaccination day (Supplemental Information).

### Usage:

To reproduce the results, execute each cell in the Jupyter Notebook from top to bottom. Note that running cells to find optimal solutions may take approximately 20-30 minutes per subset $\mathcal{V_i}$ for each patient. Other cells should run quickly.

### Dependencies:

Ensure that you have all the necessary dependencies installed. [Provide installation instructions if required.]

### Excel File:

The Excel file [\Patients_HLA_BA.xlsx](Patients_HLA_BA) contains binding affinities of MHC-I/II to specific alleles associated with the immunogenic peptides of each patient. These affinities were used to compute the dissociation constant, $K_{D,j}$. The relevant information from this spreadsheet has been integrated into the CanVaxDOpt Notebook.

### Citation:

If you use this code in your research, please cite our paper [citation here].


# Dose optimization of an adjuvanted peptide-based personalized neoantigen melanoma vaccine (2023)
## Overview 
This repository contains the codes to reproduce the numerical results in our paper [citation goes here]. It provides a Jupyter Notebook with instructions to simulate the optimal vaccine doses along with their corresponding immune and tumor responses for each patient. 

The Jupyter Notebook, CanVaxDOpt, starts with the code to run the Forward Backward Sweep method to numerically solve the optimization problem using the optimiality system, that is, the state equations, adjoint equations and optimal control characterization. We present the code to solve the optimization problem in each subset $\mathcal{V_i} \subseteq \mathcal{V}$, and how to select the peptide dose with most clinical benefit than any other tested dose. The notebook also contains the code to reproduce the graphs of optimal vaccine doses (Fig 5), activated T cells (Fig 6), tumor diameter (Fig 7) from our paper. Moreover, it provides the tables to compare the clinical benefits of using the predicted optimal vaccine dose when compared to the clinical trial dose (or any other tested dose). The tables at the end correspond to the optimal vaccine dose concentration of peptides and adjuvant per vaccination day as shown in the Suplemental Information. 

The easiest way to reproduce the results is to run each cell from the top down to avoid unexpected errors. Running the cells to find optimal solutions can take approximately 20-30 minutes per subset $\mathcal{V_i}$ for each patient. All other cells should run relatively quickly. 

The excel file, Patients_HLA_BA, has the binding affinities of MHC-I/II to specific alleles associated with the set of immunogenic peptides of each patient. These binding affinities were used to compute the dissociation constant, $K_{D,j}$. All the necessary information in the spreadsheet has been merged to the CanVaxDOpt Notebook.
