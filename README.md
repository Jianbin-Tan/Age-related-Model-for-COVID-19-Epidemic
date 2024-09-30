# Age-Related Model for Estimating the Symptomatic and Asymptomatic Transmissibility of COVID-19 Patients

This README accompanies the paper titled "Age-Related Model for Estimating the Symptomatic and Asymptomatic Transmissibility of COVID-19 Patients," authored by Jianbin Tan, Ye Shen, Yang Ge, Leonardo Martinez, and Hui Huang. The paper is published in [Biometrics](https://academic.oup.com/biometrics/article/79/3/2525/7513834?login=false).

## Citation
Please cite the following if you utilize this code, data, or methodology in your research:

```bibtex
@article{tan2023age,
  title={Age-related model for estimating the symptomatic and asymptomatic transmissibility of COVID-19 patients},
  author={Tan, Jianbin and Shen, Ye and Ge, Yang and Martinez, Leonardo and Huang, Hui},
  journal={Biometrics},
  volume={79},
  number={3},
  pages={2525--2536},
  year={2023},
  publisher={Oxford University Press}
}

---
## Data
### 1.1 Abstract

- Information was collected from all confirmed cases during this wave, totaling 1,342 individuals. The data includes individual-level symptom onset dates, confirmation dates, and ages for both symptomatic and asymptomatic patients. This dataset was instrumental in distinguishing between asymptomatic and pre-symptomatic cases, noting that all were infected by early strains of the SARS-CoV-2 virus.

- Additional data comprising contact matrices, susceptibility parameters, and demographic details of different age groups in Zhejiang province were sourced from existing literature.

### 1.2 Availability
The datasets necessary to reproduce our findings are available.

### 1.3 Data Dictionary
The dataset "ALL_dat.rda" within the "Data" folder contains contact matrices, susceptibility parameters, demographic details of different age groups, initial cases for all compartments, daily cases for seven age groups, and observed periods of symptomatic transmission across all age groups.

---
## Code
### 2.1 Abstract
The age-related model, along with two fundamental models referenced in the paper, were utilized for simulations and data analysis.

### 2.2 Reproducibility
- The simulations detailed in Section 3 of the paper (including Table 1 and Figure 3) were generated using "Simulation.R" with settings from "Sim_setting.R".
- Model estimations and comparisons in Section 4 were performed using "Running_case_study.R".
- Sensitivity analyses (incorporating Figures 2-5 from the Supporting Information) were conducted using "Sensitivity_analysis.R".
- Additional figures and tables found in the main text and Supporting Information were produced using "Plot_case_study.R".
