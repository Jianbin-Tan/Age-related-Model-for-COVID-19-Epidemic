# tf-fsvd
TensorFlow Implementation of Functional Singular Value Decomposition for paper
[**Fast Graph Learning with Unique Optimal Solutions**](https://arxiv.org/abs/2102.08530)

## Cite
If you find our code useful, you may cite us as:

    @inproceedings{haija2021fsvd,
      title={Fast Graph Learning with Unique Optimal Solutions},
      author={Sami Abu-El-Haija AND Valentino Crespi AND Greg Ver Steeg AND Aram Galstyan},
      year={2021},
      booktitle={arxiv:2102.08530},
    }
    
---
## 1. Data
### 1.1 Abstract

- Information was collected from all confirmed cases during this wave, totaling 1,342 individuals. The data includes individual-level symptom onset dates, confirmation dates, and ages for both symptomatic and asymptomatic patients. This dataset was instrumental in distinguishing between asymptomatic and pre-symptomatic cases, noting that all were infected by early strains of the SARS-CoV-2 virus.

- Additional data comprising contact matrices, susceptibility parameters, and demographic details of different age groups in Zhejiang province were sourced from existing literature.

### 1.2 Availability
The datasets necessary to reproduce our findings are available.

### 1.3 Data Dictionary
The dataset "ALL_dat.rda" within the "Data" folder contains contact matrices, susceptibility parameters, demographic details of different age groups, initial cases for all compartments, daily cases for seven age groups, and observed periods of symptomatic transmission across all age groups.

---
## 2. Code
### 2.1 Abstract
The age-related model, along with two fundamental models referenced in the paper, were utilized for simulations and data analysis.

### 2.2 Reproducibility
- The simulations detailed in Section 3 of the paper (including Table 1 and Figure 3) were generated using "Simulation.R" with settings from "Sim_setting.R".
- Model estimations and comparisons in Section 4 were performed using "Running_case_study.R".
- Sensitivity analyses (incorporating Figures 2-5 from the Supporting Information) were conducted using "Sensitivity_analysis.R".
- Additional figures and tables found in the main text and Supporting Information were produced using "Plot_case_study.R".
