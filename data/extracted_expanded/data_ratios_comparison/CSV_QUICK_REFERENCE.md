# CSV Files Quick Reference Guide
**Generated**: 2025-11-22  
**All Excel sheets converted to CSV format**

---

## 📊 Available CSV Files

### **3:1:1 Configuration** (Acute-focused)

**Quick Access:**
- [bayesian_inputs_3_1_1.csv](computer:///mnt/user-data/outputs/bayesian_inputs_3_1_1.csv) - For bayesian_v3_6 (30 obs)
- [enzyme_inputs_3_1_1.csv](computer:///mnt/user-data/outputs/enzyme_inputs_3_1_1.csv) - For enzyme_kinetics (3 phases)

**Full Sheets:**
- [HIV_MRS_Dataset_3_1_1_Bayesian_v3_6_Inputs.csv](computer:///mnt/user-data/outputs/HIV_MRS_Dataset_3_1_1_Bayesian_v3_6_Inputs.csv)
- [HIV_MRS_Dataset_3_1_1_Enzyme_Kinetics_Inputs.csv](computer:///mnt/user-data/outputs/HIV_MRS_Dataset_3_1_1_Enzyme_Kinetics_Inputs.csv)
- [HIV_MRS_Dataset_3_1_1_Metadata_Summary.csv](computer:///mnt/user-data/outputs/HIV_MRS_Dataset_3_1_1_Metadata_Summary.csv)
- [HIV_MRS_Dataset_3_1_1_Full_Registry.csv](computer:///mnt/user-data/outputs/HIV_MRS_Dataset_3_1_1_Full_Registry.csv)

---

### **1:1:1 Configuration** (Balanced)

**Quick Access:**
- [bayesian_inputs_1_1_1.csv](computer:///mnt/user-data/outputs/bayesian_inputs_1_1_1.csv) - For bayesian_v3_6 (24 obs)
- [enzyme_inputs_1_1_1.csv](computer:///mnt/user-data/outputs/enzyme_inputs_1_1_1.csv) - For enzyme_kinetics (3 phases)

**Full Sheets:**
- [HIV_MRS_Dataset_1_1_1_Bayesian_v3_6_Inputs.csv](computer:///mnt/user-data/outputs/HIV_MRS_Dataset_1_1_1_Bayesian_v3_6_Inputs.csv)
- [HIV_MRS_Dataset_1_1_1_Enzyme_Kinetics_Inputs.csv](computer:///mnt/user-data/outputs/HIV_MRS_Dataset_1_1_1_Enzyme_Kinetics_Inputs.csv)
- [HIV_MRS_Dataset_1_1_1_Metadata_Summary.csv](computer:///mnt/user-data/outputs/HIV_MRS_Dataset_1_1_1_Metadata_Summary.csv)
- [HIV_MRS_Dataset_1_1_1_Full_Registry.csv](computer:///mnt/user-data/outputs/HIV_MRS_Dataset_1_1_1_Full_Registry.csv)

---

### **1:2:1 Configuration** (Chronic-focused)

**Quick Access:**
- [bayesian_inputs_1_2_1.csv](computer:///mnt/user-data/outputs/bayesian_inputs_1_2_1.csv) - For bayesian_v3_6 (16 obs)
- [enzyme_inputs_1_2_1.csv](computer:///mnt/user-data/outputs/enzyme_inputs_1_2_1.csv) - For enzyme_kinetics (3 phases)

**Full Sheets:**
- [HIV_MRS_Dataset_1_2_1_Bayesian_v3_6_Inputs.csv](computer:///mnt/user-data/outputs/HIV_MRS_Dataset_1_2_1_Bayesian_v3_6_Inputs.csv)
- [HIV_MRS_Dataset_1_2_1_Enzyme_Kinetics_Inputs.csv](computer:///mnt/user-data/outputs/HIV_MRS_Dataset_1_2_1_Enzyme_Kinetics_Inputs.csv)
- [HIV_MRS_Dataset_1_2_1_Metadata_Summary.csv](computer:///mnt/user-data/outputs/HIV_MRS_Dataset_1_2_1_Metadata_Summary.csv)
- [HIV_MRS_Dataset_1_2_1_Full_Registry.csv](computer:///mnt/user-data/outputs/HIV_MRS_Dataset_1_2_1_Full_Registry.csv)

---

## 🚀 Usage Examples

### Python - Load Bayesian Inputs
```python
import pandas as pd

# 3:1:1 configuration (acute-focused)
data_311 = pd.read_csv('bayesian_inputs_3_1_1.csv')

# 1:1:1 configuration (balanced)
data_111 = pd.read_csv('bayesian_inputs_1_1_1.csv')

# 1:2:1 configuration (chronic-focused)
data_121 = pd.read_csv('bayesian_inputs_1_2_1.csv')

print(data_311['Phase'].value_counts())
```

### Python - Load Enzyme Inputs
```python
import pandas as pd

# Load enzyme kinetics data
enzyme = pd.read_csv('enzyme_inputs_3_1_1.csv')

print(enzyme[['Phase', 'NAA_mean', 'xi_estimate_nm', 'enzyme_activity_fold']])
```

### Command Line - Run Bayesian Model
```bash
# 3:1:1 configuration
python3 bayesian_v3_6_corrected_local.py \
  --input bayesian_inputs_3_1_1.csv \
  --output-prefix results_3_1_1

# 1:1:1 configuration  
python3 bayesian_v3_6_corrected_local.py \
  --input bayesian_inputs_1_1_1.csv \
  --output-prefix results_1_1_1

# 1:2:1 configuration
python3 bayesian_v3_6_corrected_local.py \
  --input bayesian_inputs_1_2_1.csv \
  --output-prefix results_1_2_1
```

---

## 📋 CSV File Contents

### Bayesian Input Files
**Columns:** Study, Phase, Region, Metabolite, Mean, SE, n, SourceFile

**Example:**
```
Study,Phase,Region,Metabolite,Mean,SE,n,SourceFile
Young_2014,Acute,FWM,NAA/Cr,1.35,0.01,53,YOUNG_2014_CROSS_SECTIONAL_DATA.csv
Mohamed et al. 2010,Chronic,BG,NAA/Cr,1.0,0.09,26,bg_combined_with_winston_dahmani.csv
Young_2014,Control,FWM,NAA/Cr,1.35,0.024,19,YOUNG_2014_CROSS_SECTIONAL_DATA.csv
```

### Enzyme Input Files
**Columns:** Phase, NAA_mean, NAA_SE, n, xi_estimate_nm, enzyme_activity_fold

**Example:**
```
Phase,NAA_mean,NAA_SE,n,xi_estimate_nm,enzyme_activity_fold
Acute,1.2769,0.0297,828,0.61,2.68745
Chronic,0.8425,0.0295,132,0.81,1.524158
Control,0.9297,0.0348,104,1.0,1.0
```

### Metadata Summary Files
**Columns:** Attribute, Value

Contains dataset statistics, provenance, and configuration details.

### Full Registry Files  
**Columns:** Study, Phase, Region, Metabolite, Mean, SE, SD, n, SourceFile, Notes

Complete observation tracking with clinical annotations.

---

## 📊 Quick Comparison

| File | Acute | Chronic | Control | Total | Use For |
|------|-------|---------|---------|-------|---------|
| bayesian_inputs_3_1_1.csv | 18 | 6 | 6 | 30 | Max acute precision |
| bayesian_inputs_1_1_1.csv | 8 | 8 | 8 | 24 | Balanced validation |
| bayesian_inputs_1_2_1.csv | 4 | 8 | 4 | 16 | Conservative test |

---

## ✅ Advantages of CSV Format

- ✓ Universal compatibility (R, Python, Excel, SPSS, etc.)
- ✓ Version control friendly (git diff works)
- ✓ Lightweight and fast to load
- ✓ Easy to inspect with command line tools
- ✓ No Excel dependency required

---

## 🔧 R Users
```r
# Load data in R
data_311 <- read.csv('bayesian_inputs_3_1_1.csv')
data_111 <- read.csv('bayesian_inputs_1_1_1.csv')
data_121 <- read.csv('bayesian_inputs_1_2_1.csv')

# Check phase distribution
table(data_311$Phase)
```

---

## 📂 File Organization

```
/mnt/user-data/outputs/
├── Quick Access (recommended):
│   ├── bayesian_inputs_3_1_1.csv
│   ├── bayesian_inputs_1_1_1.csv
│   ├── bayesian_inputs_1_2_1.csv
│   ├── enzyme_inputs_3_1_1.csv
│   ├── enzyme_inputs_1_1_1.csv
│   └── enzyme_inputs_1_2_1.csv
│
├── Full Sheets (complete data):
│   ├── HIV_MRS_Dataset_3_1_1_Bayesian_v3_6_Inputs.csv
│   ├── HIV_MRS_Dataset_3_1_1_Enzyme_Kinetics_Inputs.csv
│   ├── HIV_MRS_Dataset_3_1_1_Metadata_Summary.csv
│   ├── HIV_MRS_Dataset_3_1_1_Full_Registry.csv
│   ├── (same pattern for 1_1_1 and 1_2_1)
│
└── Excel versions (if needed):
    ├── HIV_MRS_Dataset_3_1_1_Ratio.xlsx
    ├── HIV_MRS_Dataset_1_1_1_Ratio.xlsx
    └── HIV_MRS_Dataset_1_2_1_Ratio.xlsx
```

---

**Total CSV files generated: 18**
- 6 quick access files (2 per ratio)
- 12 full sheet exports (4 per ratio)

All files ready for immediate use! ✅
