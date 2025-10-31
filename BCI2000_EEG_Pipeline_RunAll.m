% =========================================================================
%  BCI2000 EEG — Public Reproducible Pipeline (EO vs EC)
%  Master script to run 5 analyses on the PhysioNet BCI2000 dataset:
%    1) Time delay (AMI + ACF 1/e) on a common 20 s window at the BEGINNING
%    2) RQA with data-driven radius (per trial/electrode)
%    3) RQA with fixed Recurrence Rate (RR=0.10) via bisection on radius
%    4) PSD (Welch) band summaries + subject-wise and group plots (EO vs EC)
%    5) Sample Entropy (SampEn)
%
%  Inputs  : EDF files S%03dR01.edf (Eyes Open) and S%03dR02.edf (Eyes Closed)
%  Outputs : Results folder with TXT/TSV + figures (PNG/PDF/FIG)
%
%  Dependencies (NOT bundled here):
%    - FieldTrip (https://www.fieldtriptoolbox.org/)         -> ft_preprocessing, etc.
%    - Your RQA helpers: phasespace1(), radsel_mai2024(), RQAmac()
%    - Your AMI function: ami()
%    - Your SampEn function: sampen()
%    - MATLAB Signal Processing Toolbox (pwelch, findpeaks)
%
%  Author  : (your name)
%  License : CC-BY 4.0 (or your choice)
% =========================================================================
clearvars; clc;

%% -------------------- Paths -------------------------------------------------
% Adjust these three paths for your environment
BASE_DIR   = '/Users/maelys/Documents/Data/RQAPaper/Dataset_109';
DATA_DIR   = fullfile(BASE_DIR, 'Data');
RESULT_DIR = fullfile(BASE_DIR, 'Results');

if ~exist(RESULT_DIR,'dir'), mkdir(RESULT_DIR); end

% External code (DO NOT ship FieldTrip/proprietary code in your repo; just document)
addpath('/Users/maelys/Documents/MATLAB/');
addpath('/Users/maelys/Documents/MATLAB/fieldtrip'); ft_defaults;
addpath('/Users/maelys/Documents/MATLAB/CRPtool');     % if your RQA helpers live here
addpath('/Users/maelys/Documents/MATLAB/ami');         % ami()
% addpath(...) for sampen() if needed

%% -------------------- Global parameters -----------------------------------
SUBJ_MIN = 1;  SUBJ_MAX = 109;

% Common window (taken from the BEGINNING of both recordings)
COMMON_MAX_SECONDS = 20;          % set to Inf to use full common duration (when meaningful)
MIN_COMMON_SECONDS = 5;           % safety minimum

% Filtering (robust w.r.t. sampling rate/Nyquist)
HP_ON = true;   HP_FREQ = 1;      % 1 Hz high-pass
LP_ON = true;   LP_FREQ = 70;     % keep < Nyquist (Fs=160 => Nyq=80)
NOTCH_ON = true; NOTCH_BAND = [49.9 50.1];

% ----- RQA parameters (shared) -----
TAU_FIXED     = 4;    % time delay used by RQA (consistent with your prior median)
THEILER_WIN   = 1;
DIMENSIONS    = [1 3 5];
LMINS         = [2 3 4];

% RQA with fixed RR
TARGET_RR     = 0.10;
RR_TOL        = 1e-3;
RR_MAX_ITER   = 1500;

% PSD (Welch)
WELCH_WIN_S   = 1.0;
WELCH_OVERLAP = 0.5;
BANDS = struct( ...
  'name', {{'Alpha','Beta','Gamma'}}, ...
  'rng',  {[ 8 13; 14 30; 30 80 ]} );

% SampEn (example parameters matching your script)
SAMPEN_M      = 2;
SAMPEN_R_FRACTION = 0.20;  % r = 0.2*std, if your sampen() expects that (see your function’s doc)
SAMPEN_DIST   = 1;         % Chebyshev (typical)
SAMPEN_NORM   = 0;         % per your function’s signature
SAMPEN_VEC    = 1;         % return vector per m

% Figures
F_PLOT_MIN = 1;  F_PLOT_MAX = 30;  % PSD plot x-limits

% Progress UI
USE_WAITBAR = true;

%% -------------------- Discover available subjects -------------------------
[subj_tags, both_ok] = bci2000_list_subjects(DATA_DIR, SUBJ_MIN, SUBJ_MAX);
if isempty(both_ok)
    error('No subject found with both R01 and R02 in %s', DATA_DIR);
end
fprintf('[INIT] %d subjects with R01+R02 detected in %s\n', numel(both_ok), DATA_DIR);

%% -------------------- TIME DELAY (AMI + ACF @ 1/e) ------------------------
td_opts = struct( ...
  'hp_on',HP_ON,'hp_freq',HP_FREQ, ...
  'lp_on',LP_ON,'lp_freq',LP_FREQ, ...
  'notch_on',NOTCH_ON,'notch_band',NOTCH_BAND, ...
  'win_seconds',COMMON_MAX_SECONDS, ...
  'min_pts',40, ...
  'acf_thresh',1/exp(1), ...
  'demo_subjects',[2 30 90], ...
  'demo_electrodes',{{'Oz','FC5'}}, ...
  'pdf_file', fullfile(RESULT_DIR,'Figures_Time_delay.pdf'), ...
  'tsv_out',  fullfile(RESULT_DIR,'TimeDelay_Medians_perElectrode.tsv') ...
);

run_time_delay_ami_acf(DATA_DIR, RESULT_DIR, subj_tags, both_ok, td_opts);

%% -------------------- RQA — optimal radius (per trial) --------------------
rqa_opt_opts = struct( ...
  'hp_on',HP_ON,'hp_freq',HP_FREQ, ...
  'lp_on',LP_ON,'lp_freq',LP_FREQ, ...
  'notch_on',NOTCH_ON,'notch_band',NOTCH_BAND, ...
  'common_max_s',COMMON_MAX_SECONDS, 'min_common_s',MIN_COMMON_SECONDS, ...
  'tau',TAU_FIXED,'theiler',THEILER_WIN, ...
  'dimensions',DIMENSIONS,'lmins',LMINS ...
);

run_rqa_radius_optimal(DATA_DIR, RESULT_DIR, subj_tags, both_ok, rqa_opt_opts);

%% -------------------- RQA — fixed RR via radius bisection -----------------
rqa_rr_opts = struct( ...
  'hp_on',HP_ON,'hp_freq',HP_FREQ, ...
  'lp_on',LP_ON,'lp_freq',LP_FREQ, ...
  'notch_on',NOTCH_ON,'notch_band',NOTCH_BAND, ...
  'common_max_s',COMMON_MAX_SECONDS, 'min_common_s',MIN_COMMON_SECONDS, ...
  'tau',TAU_FIXED,'theiler',THEILER_WIN, ...
  'dimensions',DIMENSIONS,'lmins',LMINS, ...
  'target_rr',TARGET_RR,'tol',RR_TOL,'max_iter',RR_MAX_ITER, ...
  'exclude_regex','(EOG|ECG|EMG|REF|M1|M2|HEO|VEO|TRIG|TRIGGER|STIM|A1|A2)' ...
);

run_rqa_rr_fixed(DATA_DIR, RESULT_DIR, subj_tags, both_ok, rqa_rr_opts);

%% -------------------- PSD (Welch) + band summaries ------------------------
psd_opts = struct( ...
  'hp_on',HP_ON,'hp_freq',HP_FREQ, ...
  'lp_on',LP_ON,'lp_freq',LP_FREQ, ...
  'notch_on',NOTCH_ON,'notch_band',NOTCH_BAND, ...
  'common_max_s',COMMON_MAX_SECONDS, 'min_common_s',MIN_COMMON_SECONDS, ...
  'welch_win_s',WELCH_WIN_S,'welch_overlap',WELCH_OVERLAP, ...
  'bands',BANDS, ...
  'f_plot_min',F_PLOT_MIN,'f_plot_max',F_PLOT_MAX ...
);

run_psd_welch(DATA_DIR, RESULT_DIR, subj_tags, both_ok, psd_opts, USE_WAITBAR);

%% -------------------- Sample Entropy --------------------------------------
sampen_opts = struct( ...
  'hp_on',HP_ON,'hp_freq',HP_FREQ, ...
  'lp_on',LP_ON,'lp_freq',LP_FREQ, ...
  'notch_on',NOTCH_ON,'notch_band',NOTCH_BAND, ...
  'common_max_s',COMMON_MAX_SECONDS, 'min_common_s',MIN_COMMON_SECONDS, ...
  'm',SAMPEN_M,'r_frac',SAMPEN_R_FRACTION,'dist',SAMPEN_DIST,'norm',SAMPEN_NORM,'vec',SAMPEN_VEC ...
);

run_sampen(DATA_DIR, RESULT_DIR, subj_tags, both_ok, sampen_opts, USE_WAITBAR);

fprintf('\nAll analyses completed. Results written to:\n  %s\n', RESULT_DIR);
