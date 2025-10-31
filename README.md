BCI2000-EO-EC-EEG-Analysis-Pipeline

Reproducible MATLAB pipeline for the PhysioNet BCI2000 EEG dataset (Eyes Open R01 vs Eyes Closed R02).

It computes:
- Time delay (AMI first minimum & ACF at 1/e) on a common 20 s window,
- RQA with data-driven radius,
- RQA with fixed Recurrence Rate (RR = 0.10 via robust bisection on radius),
- PSD (Welch) band summaries + subject & group figures,
- Sample Entropy (SampEn).

FieldTrip is used for EDF I/O and filtering; custom helper functions are referenced (not redistributed).
Outputs are plain TXT/TSV tables and PNG/PDF/FIG figures.

Features
- Common 20 s window at recording start (EO/EC) with identical filtering: HP 1 Hz, LP 70 Hz, optional 50 Hz notch (if Nyquist > 50).
- RQA: radius-optimal and RR-fixed (RR = 0.10), m ∈ {1, 3, 5}, Lmin ∈ {2, 3, 4}, τ = 4, Theiler = 1.
- PSD (Welch): absolute/relative power and mean frequency for Alpha (8–13 Hz), Beta (14–30 Hz), Gamma (30–80 Hz);
per-subject 64-channel plots + group mean plots (EO vs EC).
- SampEn: defaults m = 2, r = 0.2 · σ (adapt to your sampen function signature if needed).
- Label cleaning + optional exclusion of non-EEG channels for RR-fixed analyses (EOG|ECG|EMG|REF|M1|M2|…).

Dependencies (not bundled)
- MATLAB R20XXa/b (⚑ specify your version)
- Toolboxes : Signal Processing Toolbox (for pwelch, findpeaks)
- FieldTrip (for EDF I/O and filtering): https://www.fieldtriptoolbox.org/

Custom helpers (to provide or reference):
- phasespace1.m — phase-space embedding
- radsel_mai2024.m — radius heuristic
- RQAmac.m — RQA metrics (RR, DET, Lmean, Lmax, ENTR, LAM, Vmax)
- ami.m — Average Mutual Information
- sampen.m — Sample Entropy

Default Folder Layout; By default, the master script expects:
- BASE_DIR = /.../Dataset_109
- DATA_DIR = BASE_DIR/Data
- RESULT_DIR = BASE_DIR/Results

Each subject folder must contain both runs:

Data/
 ├─ S001/S001R01.edf   # Eyes Open
 └─ S001/S001R02.edf   # Eyes Closed

Key Parameters (default values)
- All parameters can be edited at the top of BCI2000_EEG_Pipeline_RunAll.m.
- Windowing : COMMON_MAX_SECONDS = 20 (use Inf for full common duration); MIN_COMMON_SECONDS = 5
- Filtering : HP_FREQ = 1, LP_FREQ = 70, NOTCH_BAND = [49.9 50.1] --> Filters applied only if within Nyquist range.
- RQA : DIMENSIONS = [1 3 5], LMINS = [2 3 4], TAU_FIXED = 4, THEILER_WIN = 1; RR-fixed: TARGET_RR = 0.10, RR_TOL = 1e-3, RR_MAX_ITER = 1500
- PSD (Welch): WELCH_WIN_S = 1.0, WELCH_OVERLAP = 0.5 --> Frequency bands : {Alpha, Beta, Gamma}
- SampEn : m = 2, r_frac = 0.20, dist = 1, norm = 0, vec = 1

Outputs : All result tables are tab-separated (TSV/TXT) with plain ASCII headers.
- Time Delay
      Figures_Time_delay.pdf — demo pages + global histograms
      TimeDelay_Medians_perElectrode.tsv
- RQA (radius-optimal)
      Processing_RQA_Radjust_m{m}_td4_lmin{L}.txt
- RQA (RR-fixed)
      Processing_RQA_RRfixed_m{m}_td4_lmin{L}.txt
- PSD (Welch)
      Processing_PSD.txt
      PSD/PSD_S%03d_64ch.fig|.png — per-subject
      PSD/PSD_mean_64ch.fig|.png — group mean
- Sample Entropy
      Processing_SampEn.txt

If you use this pipeline, please cite:

Code / Methods
Moulin, M. et al. “Recurrence Quantification Analysis of EEG signals: Parameters Selection and Correlations with Spectral Features.” (2025, venue TBD)

Dataset
Goldberger A. L. et al., PhysioNet/PhysioBank BCI2000 EEG dataset, https://physionet.org/
