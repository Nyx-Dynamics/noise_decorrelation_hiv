# HIV MRS Data Inventory

## Study-by-Study Breakdown

### 1. Sailasuta 2012 (SEARCH 010/RV254)
- **Type**: Aggregated (means ± SD)
- **Groups**: Control (n=10), Acute HIV (n=31, 14d post-infection), Chronic (n=26)
- **Regions**: BG, OGM
- **Metabolites**: NAA, Cho, mI, Glx
- **Format**: Ratio to Creatine
- **Key finding**: NAA preserved in acute phase despite peak viremia

### 2. Valcour 2015 (RV254)
- **Type**: INDIVIDUAL patient data (n=62)
- **Time point**: Acute HIV (Fiebig I-V)
- **Regions**: FGM, FWM, BG, PGM
- **Metabolites**: NAA, Cho (absolute concentrations)
- **Also includes**: Age, CD4, log viral load, NPZ4 cognitive scores
- **Key finding**: +7.7% NAA elevation (p=0.0317) in BG

### 3. Young 2014
- **Type**: Aggregated
- **Groups**: Control (n=19), PHI ~6mo (n=9), Chronic ~10y (n=18)
- **Region**: BG
- **Format**: Ratio to Creatine
- **Key finding**: Intermediate NAA in PHI group

### 4. Chang 2002
- **Type**: Aggregated
- **Groups**: Control (n=25), Early HIV ~2y (n=45)
- **Regions**: BG, FWM, FGM
- **Metabolites**: NAA, Cho, mI
- **Format**: Absolute concentrations (mM)

### 5. Sailasuta 2016
- **Type**: Aggregated longitudinal
- **Groups**: ART-naive baseline, 6mo, 12mo
- **Regions**: BG, FWM
- **Intervention**: Treatment response to cART

### 6. Dahmani 2021
- **Type**: Meta-analysis
- **Acute**: 146 patients, 9 studies
- **Chronic**: 943 patients, 60 studies
- **Result**: No significant NAA change in acute (p=NS), significant decrease in chronic (p<0.001)

## Data Quality Assessment

### High Quality (usable for modeling)
✓ Sailasuta 2012 - Full statistics, multiple metabolites
✓ Valcour 2015 - Individual data (gold standard)
✓ Young 2014 - Clear group definitions
✓ Chang 2002 - Absolute concentrations

### Medium Quality (usable with caveats)
~ Sailasuta 2016 - Missing SD/SE, longitudinal confounds

### Meta-analytic Only
✓ Dahmani 2021 - Validates overall pattern, not for mechanistic modeling

## Current Model Inputs

### v3.6 Model (VALIDATED)
**Input**: `NAA_DATA_FOR_MODEL.csv`
- Control: 8.76 ± 0.79 mM (Chang 2002, BG)
- Acute: 10.0 (Sailasuta 2016, BG, ART-naive)
- Chronic: 7.96 ± 0.91 mM (Chang 2002, BG)

**Result**: P(ξ_acute < ξ_chronic) > 0.999

### Future v4 Model Options
**Option A**: Same aggregated data with enzyme kinetics
**Option B**: Valcour individual data with hierarchical structure
**Option C**: Combined approach (requires careful likelihood specification)

## Missing Data / Opportunities

### What We Need
1. More acute phase individual data (not just Valcour)
2. Longitudinal follow-up: acute → chronic same patients
3. Viral load + metabolite + cognitive scores
4. Regional specificity in acute phase

### Who Might Have It
- Dr. Serena Spudich (Yale) - RV254 extended analyses
- Dr. Sailasuta (?)
- Dr. Valcour - may have additional unpublished data

## Data Access Notes

- Most published data is aggregated (authors report means, not raw data)
- Individual data is rare because:
  * Privacy concerns
  * Not typically required for publication
  * Authors may not share
- Valcour 2015 individual data = GOLD MINE
- Consider reaching out for expanded datasets after first publication
