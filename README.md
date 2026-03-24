# Surveillance-as-Governance.Paper-Repo

## Repository Structure

This repository contains the data and code used in the analysis for the paper *“Surveillance as Governance: Policing Effectiveness and Political Control in the Moscow AI Experiment”* by Dmitrii Serebrennikov,	Eleonora Minaeva and	Sergey Ross.

The project is organized into two main directories: `code/` and `data/`:
- `code/` — scripts for data processing and statistical modeling  
- `data/` — datasets used in the analysis  

---

### `code/`

This folder contains R scripts implementing the empirical analysis presented in the paper. The scripts are organized by analytical components:

- **Facial Recognition, CCTV, and Crime Clearance. Models**  
  Scripts for estimating panel regression models of crime clearance rates, including:
  - balanced panel models with multiple imputation  
  - unbalanced models for specific crime categories (e.g., theft, public-space crimes)  

- **Protest Activity and the Expansion of Surveillance Infrastructure. Models**  
  Scripts for estimating spatial and temporal models of CCTV expansion, including:
  - negative binomial regression models  with spatial lag and other specifications for protest activity  

---

### `data/`

This folder contains processed datasets used in the analysis. All files are stored in `.RDS` format for efficient loading in R.

- `cctv_crimes_data_balanced.RDS`  
  Balanced district-level panel dataset used for the main crime clearance models.

- `cctv_crimes_data_imputed.RDS`  
  Dataset with multiple imputation via MICE approach,, including auxiliary variables used in the imputation procedure.

- `cctv_crimes_data_unbalanced.RDS`  
  Unbalanced panel dataset with observed data only, used for robustness checks and models on specific crime categories.

- `cctv_protests_data_long.RDS`  
  Long-format dataset at the hexagon–halfyear level, combining:
  - number of newly installed CCTV  
  - protest activity (data from OVD-Info)  
  - spatial and temporal covariates  

---

### Suggested Workflow

1. Run scripts in `code/` depending on the analysis:
   - clearance models  
   - protest–CCTV models
2. The datasets automatically will be loaded from the `data/` directory  
3. Reproduce tables and results reported in the paper  
