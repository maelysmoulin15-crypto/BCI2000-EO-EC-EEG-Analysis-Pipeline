BCI2000-EO-EC-EEG-Analysis-Pipeline

Reproducible MATLAB pipeline for the PhysioNet BCI2000 EEG dataset
(Eyes Open R01 vs Eyes Closed R02)

Overview : 
This repository provides a reproducible MATLAB pipeline for EEG preprocessing and nonlinear analysis, combining spectral (Welch) and nonlinear (RQA, SampEn) metrics on the PhysioNet BCI2000 dataset.

The pipeline computes:
- Time delay (Average Mutual Information minimum & ACF at 1/e) on a 20 s common window
- RQA with data-driven radius selection
- RQA with fixed Recurrence Rate (RR = 0.10 via iterative radius adjustment)
- Power Spectral Density (Welch): absolute and relative band powers (α, β, γ) and mean frequency
- Sample Entropy (SampEn) for each electrode and condition

FieldTrip is used for EDF I/O and filtering.
All outputs are written as plain TXT/TSV tables and PNG/PDF/FIG figures.

Main Features : 
Common 20 s window at the beginning of each recording (EO/EC)
→ High-pass = 1 Hz, Low-pass = 70 Hz, optional 50 Hz notch if Nyquist > 50 Hz
RQA parameters:
→ Data-driven radius (radsel) or RR-fixed (RR = 0.10)
→ Dimensions = 1, 3, 5; Lmin = 2, 3, 4; τ = 4; Theiler = 1
PSD (Welch): absolute/relative power + mean frequency
→ Alpha (8–13 Hz), Beta (14–30 Hz), Gamma (30–80 Hz)
→ 64-channel per-subject plots + group mean plots (EO vs EC)
SampEn: default parameters m = 2, r = 0.2·σ

Automatic label cleaning and optional exclusion of non-EEG channels (EOG, ECG, REF, M1, M2, …)

Dependencies (not bundled)
- MATLAB R2024b
- Signal Processing Toolbox (for pwelch, findpeaks)
- FieldTrip (for EDF I/O and filtering): https://www.fieldtriptoolbox.org/

Custom Helper Functions (not redistributed)
1. phasespace1.m
Phase-space reconstruction via the Method of Delays (Takens embedding).
Reference : Leontitsis A. (2001) – phasespace1.m, University of Ioannina. Based on : Takens F. (1981). Detecting strange attractors in turbulence, Lecture Notes in Mathematics 898, Springer.

2. radsel_mai2024.m
Data-driven radius selection for nonlinear measures (e.g. correlation dimension, KS entropy) using Kernel Density Estimation (KDE).
Reference : Medrano J., Kheddar A., Lesne A., Ramdani S. (2021). Radius selection using kernel density estimation for the computation of nonlinear measures, Chaos 31(8): 083131. DOI: https://doi.org/10.1063/5.0055797/

3. sampen.m
Sample Entropy (SampEn) with options for standardization, C-acceleration and SE estimates.
Reference : Lake D.E., Richman J.S., Moorman J.R., Goldberger A.L. (2002). Sample entropy analysis of neonatal heart rate variability, AJP – Regulatory, Integrative and Comparative Physiology 283(3): R789–R797. DOI: https://doi.org/10.1152/ajpregu.00069.2002

4. ami.m
Average Mutual Information for univariate/bivariate time series.
Reference : Durga Lal Shrestha (2025). ami and correlation (https://fr.mathworks.com/matlabcentral/fileexchange/7936-ami-and-correlation), MATLAB Central File Exchange. Extrait(e) le octobre 31, 2025. 

5. RQAmac.m → rp_mac64
MATLAB wrapper for the CRP/rp binary (computes RR, DET, Lmean, Lmax, ENTR, LAM, Vmax).
Reference : Marwan N. (2007). Recurrence plots for the analysis of complex systems, Physics Reports 438: 237–329.
Software : Commandline Recurrence Plots (rp) v1.13z2 (2006-03-09), CRP Toolbox – PIK Potsdam.

Key Parameters (defaults)
- Windowing
     COMMON_MAX_SECONDS = 20
     MIN_COMMON_SECONDS = 5
- Filtering
     HP_FREQ = 1
     LP_FREQ = 70
     NOTCH_BAND = [49.9 50.1]
- RQA
     DIMENSIONS = [1 3 5]
     LMINS = [2 3 4]
     TAU_FIXED = 4
     THEILER_WIN = 1

- For RR-fixed: TARGET_RR = 0.10, RR_TOL = 1e-3, RR_MAX_ITER = 1500
- PSD (Welch)
      WELCH_WIN_S = 1.0
      WELCH_OVERLAP = 0.5
      Bands: Alpha, Beta, Gamma
- SampEn
      m = 2
      r_frac = 0.20
      dist = 1
      norm = 0
      vec = 1

Outputs
- All outputs are plain text (TSV/TXT) with ASCII headers.
- Time Delay : Figures_Time_delay.pdf (demo + global histograms) & TimeDelay_Medians_perElectrode.tsv
- RQA (radius-optimal) : Processing_RQA_Radjust_m{m}_td4_lmin{L}.txt
- RQA (RR-fixed) : Processing_RQA_RRfixed_m{m}_td4_lmin{L}.txt
- PSD (Welch) : Processing_PSD.txt & PSD/PSD_S%03d_64ch.fig or .png (per subject) & PSD/PSD_mean_64ch.fig or .png (group mean)
- Sample Entropy : Processing_SampEn.txt

If you use this pipeline, please cite:

Code / Methods
Moulin M. et al. (2025). Recurrence Quantification Analysis of EEG signals: Parameters Selection and Correlations with Spectral Features. (venue TBD)

Dataset
Goldberger A.L. et al. (2000). PhysioNet/PhysioBank BCI2000 EEG Dataset.
https://physionet.org/
