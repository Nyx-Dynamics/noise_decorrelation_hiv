# Data Directory Reference

## Quick Navigation

```
data/
├── extracted/      → Group-level MRS statistics (n=3 model)
├── individual/     → Patient-level data (validation)
├── raw/           → Original source files (READ-ONLY)
├── master/        → Consolidated databases
├── processed/     → Cleaned/merged datasets (for future use)
├── analysis_outputs/ → Cross-analysis results (for future use)
├── documentation/ → Extraction methodology notes
├── figures/       → Data visualizations
└── papers/        → PDF copies of source publications
```

---

## Directory Purposes

### `extracted/` - Group-Level Statistics
Clean CSV files containing **means, SDs, and n** for each study group.
Used as input for Bayesian v3.6 model (n=3: control, chronic, acute).

**Files:**
- `SAILASUTA_2012_ACUTE_DATA.csv` - RV254 acute cohort
- `VALCOUR_2015_REGIONAL_SUMMARY.csv` - Regional MRS data
- `YOUNG_2014_CROSS_SECTIONAL_DATA.csv` - Baseline measurements
- `YOUNG_2014_LONGITUDINAL.csv` - Follow-up data
- `CHANG_2002_EXTRACTED.csv` - Chronic HIV cohort
- `CRITICAL_STUDIES_COMPLETE_DATA.csv` - Multi-study summary

### `individual/` - Patient-Level Data
Individual patient measurements for validation and statistical testing.

**Files:**
- `VALCOUR_2015_INDIVIDUAL_PATIENTS.csv` - n=62 acute patients (KEY FILE)
- `VALCOUR_2015_DATA_FOR_MASTER.csv` - Formatted for master database
- `VALCOUR_2015_ACUTE_DATA.csv` - Acute phase subset

**Usage:** Individual-level t-tests, correlation analyses, subgroup validation

### `raw/` - Original Source Materials
**READ-ONLY.** Never modify these files - they are the source of truth.

Contains:
- Excel supplementary tables from publications
- Word documents with data tables
- Original data exports

**Citation:** Reference these files in manuscript methods for data provenance.

### `master/` - Consolidated Databases
Merged datasets combining multiple studies.

**Files:**
- `MASTER_HIV_MRS_DATABASE_v2.csv` - Comprehensive multi-study database

### `documentation/` - Extraction Notes
Markdown files documenting how data was extracted from papers.

Contains:
- Extraction methodology
- Data transformations applied
- Quality control notes
- Study inclusion/exclusion criteria

### `figures/` - Data Visualizations
Plots and visualizations of raw data (not model results).

### `papers/` - Source Publications
PDF copies of papers from which data was extracted.

---

## Data Flow

```
raw/ (Excel/Word)
  ↓ [Manual extraction]
extracted/ (clean CSVs)
  ↓ [Bayesian inference]
../results/bayesian_v3_6/ (model outputs)

individual/ (patient data)
  ↓ [Statistical tests]
Validation of model predictions
```

---

## Important Notes

1. **All files in `extracted/` are derived from `raw/`** - full audit trail
2. **`individual/` data enables validation** - not used for model training
3. **Never modify `raw/` or `extracted/`** - these are data of record
4. **Use `processed/` for any data transformations** - keep originals pristine
5. **`master/` is for reference** - not used directly in primary analyses

---

## Data Standards

### File Naming Convention
`AUTHOR_YEAR_DESCRIPTION.csv`

Examples:
- `SAILASUTA_2012_ACUTE_DATA.csv`
- `VALCOUR_2015_REGIONAL_SUMMARY.csv`

### CSV Format
- UTF-8 encoding
- Comma-separated
- Headers in first row
- Missing data: `NA` or empty cell
- Standard column names where possible:
  - `study`, `subject_id`, `group`, `region`
  - `NAA`, `Cho`, `Cr`, `Glu` (metabolite names)
  - `NAA_mean`, `NAA_sd`, `n` (for group data)

---

## For Manuscript Preparation

**Methods - Data Sources:**
```
Data were extracted from published studies [citations].
Original source materials are available in data/raw/.
Extracted data available in data/extracted/ and data/individual/.
```

**Methods - Data Processing:**
```
See data/documentation/ for complete extraction methodology.
All data transformations documented and reproducible.
```

**Data Availability Statement:**
```
All data used in this study are available in the data/
directory of the project repository. Original source
materials (data/raw/) are from published studies [citations].
```

---

*See PROJECT_STRUCTURE.md for complete repository organization*
