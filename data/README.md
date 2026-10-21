# Fluorescent protein maturation data

Single-cell fluorescence time series of cells expressing fluorescent
proteins (FPs) from an EL222-inducible (`EL222_AQTrip`) construct, measured
by time-lapse microscopy in a microfluidic device. These data are analysed in
`experiments/inference_experimental.m` to infer FP maturation rate
distributions, see [1] in the main README for the experimental details.

## Files

Each `EL222_AQTrip*_<FP>_<imaging settings>_<replicate>.xlsx` file contains
one row per cell and time point. The columns used by the analysis are:

| Column | Meaning |
| --- | --- |
| `TimeID` | Time point index |
| `Time` | Time (min); frames are 5 min apart |
| `Cell_pos` | Cell identifier (imaging position and cell number) |
| `CFP_total`, `GFP_total`, `YFP_total`, `RFP_total` | Total cell fluorescence in the channel matching the FP |
| `volume` | Cell volume; concentrations are computed as total fluorescence / volume |

`FP_dilution_rates.xlsx` lists the dilution rates `kdil` (min^-1) used for
each FP; these values are hard-coded in `inference_experimental.m`.

Files used in the paper (with the FP channel used):

| FP | File | Channel |
| --- | --- | --- |
| mTurquoise2 | `EL222_AQTrip-mTurquoise2_microfluidic_CFP_15x300ms_rep1.xlsx` | CFP |
| mTFP1 | `EL222_AQTrip_mTFP1_microfluidic_GFP(new)_12x300ms_rep1.xlsx` | GFP |
| mScarlet-I | `EL222_AQTrip-mScarletI_microfluidic_CFP_1x1s_rep1.xlsx` | RFP |
| mCherry | `EL222_AQTrip-mCherry_microfluidic_CFP_1x1s_rep2.xlsx` | RFP |
| tdTomato | `EL222_AQTrip_tdTomato_microfluidic_CFP_1x1s_rep1.xlsx` | RFP |
| mKate2 | `EL222_AQTrip-mKate2_microfluidic_CFPnew_1x1s_rep2.xlsx` | RFP |
| CFP | `EL222_AQTrip-CFP_microfluidic_CFP_15x300ms_rep1.xlsx` | CFP |
| sfGFP | `EL222_AQTrip-sfGFP_microfluidic_CFP_1x1s_rep1.xlsx` | GFP |
| pHtdGFP | `EL222_AQTrip-pHtdGFP_microfluidic_GFP(new)_12x300ms_rep1.xlsx` | GFP |
| mVenus | `EL222_AQTrip-mVenus_microfluidic_CFP_1x1s_rep2.xlsx` | YFP |
| mCitrine | `EL222_AQTrip-mCitrine_microfluidic_CFP_1x1s.xlsx` | YFP |
| mNeonGreen | `EL222_AQTrip-mNeonGreen_microfluidic_CFP_1x1s_rep1.xlsx` | YFP |

The files `EL222_AQTrip-sfGFP_microfluidic_CFP_1x1s_rep2.xlsx` and
`EL222_AQTrip_sfGFP_YFPimaging_rep1.xlsx` are additional replicates that
are not used in the analysis.

Preprocessing (selection of cells and time window, interpolation to a
regular grid, background subtraction) is done in
`inference_experimental.m`.
