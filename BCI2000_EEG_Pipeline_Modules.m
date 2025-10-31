% =========================================================================
%  BCI2000_EEG_Pipeline_Modules.m
%  Helper functions used by BCI2000_EEG_Pipeline_RunAll.m
%  (FieldTrip I/O, common-window extraction, label cleaning, and
%   the five analysis modules: time delay, RQA-opt, RQA-RRfixed, PSD, SampEn)
%
%  NOTE: This file assumes you have FieldTrip on your MATLAB path and
%        your custom functions: phasespace1, radsel_mai2024, RQAmac, ami, sampen.
% =========================================================================

function [subj_tags, both_ok] = bci2000_list_subjects(data_dir, smin, smax)
subj_tags = arrayfun(@(s) sprintf('S%03d', s), smin:smax, 'uni', 0);
has_R01 = false(size(subj_tags)); has_R02 = false(size(subj_tags));
for i = 1:numel(subj_tags)
    st = subj_tags{i};
    has_R01(i) = exist(fullfile(data_dir, st, sprintf('%sR01.edf', st)), 'file') ~= 0;
    has_R02(i) = exist(fullfile(data_dir, st, sprintf('%sR02.edf', st)), 'file') ~= 0;
end
both_ok = find(has_R01 & has_R02);
end

% ----------------------------- TIME DELAY ---------------------------------
function run_time_delay_ami_acf(data_dir, result_dir, subj_tags, both_ok, O)
if exist(O.pdf_file,'file'), delete(O.pdf_file); end

% collect per-electrode vectors
elec_labels = {};
EO_AMI = {}; EC_AMI = {}; EO_ACF = {}; EC_ACF = {};
init_collectors = false;

for k = 1:numel(both_ok)
    sid_num  = parse_sid(subj_tags, both_ok(k));
    subj_tag = sprintf('S%03d', sid_num);
    [raw1, raw2, Fs] = load_pair(data_dir, subj_tag);

    [dat1, dat2] = ft_filter_pair(raw1, raw2, Fs, O.hp_on,O.hp_freq,O.lp_on,O.lp_freq,O.notch_on,O.notch_band);

    % common window at the BEGINNING
    [XEO, XEC, L] = common_begin_window(dat1, dat2, Fs, O.win_seconds, 5);
    if L < O.min_pts, fprintf(2,'[%s] Too short (%d pts), skip.\n', subj_tag, L); continue; end

    if ~init_collectors
        elec_labels = dat1.label(:);
        nElec = numel(elec_labels);
        EO_AMI = cell(1,nElec); EC_AMI = cell(1,nElec);
        EO_ACF = cell(1,nElec); EC_ACF = cell(1,nElec);
        init_collectors = true;
    else
        nElec = numel(elec_labels);
    end

    for e = 1:nElec
        xEO = XEO(e,:);  xEC = XEC(e,:);

        % AMI (first minimum via peaks on -AMI)
        [ami_eo, ~] = ami(xEO, 50, 40);
        td_ami_eo = first_min_lag_from_negpeaks(ami_eo);

        [ami_ec, ~] = ami(xEC, 50, 40);
        td_ami_ec = first_min_lag_from_negpeaks(ami_ec);

        % ACF threshold 1/e
        td_acf_eo = acf_threshold_lag(xEO, O.acf_thresh);
        td_acf_ec = acf_threshold_lag(xEC, O.acf_thresh);

        if ~isnan(td_ami_eo), EO_AMI{e}(end+1) = td_ami_eo; end
        if ~isnan(td_ami_ec), EC_AMI{e}(end+1) = td_ami_ec; end
        EO_ACF{e}(end+1) = td_acf_eo;
        EC_ACF{e}(end+1) = td_acf_ec;

        % Optional 3x2 demo pages for selected subjects/electrodes
        if any(sid_num==O.demo_subjects)
            key = lower(regexprep(elec_labels{e},'[^A-Za-z0-9]',''));
            if ismember(key, lower(regexprep(O.demo_electrodes,'[^A-Za-z0-9]','')))
                plot_time_delay_demo_3x2(O.pdf_file, sid_num, elec_labels{e}, ...
                    xEO, ami_eo, td_acf_eo, O.acf_thresh, ...
                    xEC, ami_ec, td_acf_ec);
            end
        end
    end
    fprintf('[%s] OK: 20 s, Fs=%.1f, %d electrodes.\n', subj_tag, Fs, nElec);
end

% export medians per electrode
if ~isempty(elec_labels)
    nE = numel(elec_labels);
    medEO_AMI = nan(nE,1); medEC_AMI = nan(nE,1);
    medEO_ACF = nan(nE,1); medEC_ACF = nan(nE,1);
    td_median_total = nan(nE,1);

    for e = 1:nE
        if ~isempty(EO_AMI{e}), medEO_AMI(e) = median(EO_AMI{e},'omitnan'); end
        if ~isempty(EC_AMI{e}), medEC_AMI(e) = median(EC_AMI{e},'omitnan'); end
        if ~isempty(EO_ACF{e}), medEO_ACF(e) = median(EO_ACF{e},'omitnan'); end
        if ~isempty(EC_ACF{e}), medEC_ACF(e) = median(EC_ACF{e},'omitnan'); end
        combo = [EO_ACF{e}, EC_ACF{e}];
        if ~isempty(combo), td_median_total(e) = median(combo,'omitnan'); end
    end

    T = table(elec_labels, medEO_AMI, medEC_AMI, medEO_ACF, medEC_ACF, td_median_total, ...
        'VariableNames', {'Electrode','Median_AMI_EO','Median_AMI_EC','Median_ACF_EO','Median_ACF_EC','td_median_total'});
    writetable(T, O.tsv_out, 'FileType','text','Delimiter','\t');
    fprintf('>> Saved: %s (%d electrodes)\n', O.tsv_out, height(T));

    % Global histograms (optional, consistent with your script)
    plot_td_global_histograms(O.pdf_file, EO_AMI, EC_AMI, EO_ACF, EC_ACF);
end
end

% ----------------------------- RQA — optimal radius -----------------------
function run_rqa_radius_optimal(data_dir, result_dir, subj_tags, both_ok, O)
for m = O.dimensions
for lmin = O.lmins
    out_file = fullfile(result_dir, sprintf('Processing_RQA_Radjust_m%d_td%d_lmin%d.txt', m, O.tau, lmin));
    rows = {};
    hdr  = {'Subject','Electrode','Condition','Run','Duration_s', ...
            'RR','DET','Lmean','Lmax','ENTR','LAM','Vmax', ...
            'Dimension','Lmin','Radius','Tau','Fs'};

    for k = 1:numel(both_ok)
        sid_num  = parse_sid(subj_tags, both_ok(k));
        subj_tag = sprintf('S%03d', sid_num);
        try
            [raw1, raw2, Fs] = load_pair(data_dir, subj_tag);
            [dat1, dat2] = ft_filter_pair(raw1, raw2, Fs, O.hp_on,O.hp_freq,O.lp_on,O.lp_freq,O.notch_on,O.notch_band);
            [XEO, XEC, L] = common_begin_window(dat1, dat2, Fs, O.common_max_s, O.min_common_s);
            if L < O.min_common_s*Fs, continue; end

            labels = clean_labels(dat1.label);
            for e = 1:numel(labels)
                % EO
                [R, radius] = rqa_with_opt_radius(XEO(e,:), m, O.tau, lmin, O.theiler);
                rows(end+1,:) = {subj_tag, labels{e}, 'EyesOpen','R01', L/Fs, ...
                    R(1),R(2),R(3),R(4),R(5),R(6),R(7), m,lmin,radius,O.tau,Fs};

                % EC
                [R, radius] = rqa_with_opt_radius(XEC(e,:), m, O.tau, lmin, O.theiler);
                rows(end+1,:) = {subj_tag, labels{e}, 'EyesClosed','R02', L/Fs, ...
                    R(1),R(2),R(3),R(4),R(5),R(6),R(7), m,lmin,radius,O.tau,Fs};
            end
            fprintf('[%s] RQA-opt: %.1f s, %d ch, Fs=%.1f\n', subj_tag, L/Fs, numel(labels), Fs);
        catch ME
            fprintf(2,'[%s] RQA-opt ERROR: %s -> skip\n', subj_tag, ME.message);
        end
    end

    if ~isempty(rows)
        T = cell2table(rows,'VariableNames',hdr);
        writetable(T, out_file, 'Delimiter','\t');
        fprintf('>> Saved: %s (%d rows)\n', out_file, height(T));
    else
        fprintf(2,'[RQA-opt] No rows for m=%d, lmin=%d\n', m, lmin);
    end
end
end
end

% ----------------------------- RQA — fixed RR via bisection ---------------
function run_rqa_rr_fixed(data_dir, result_dir, subj_tags, both_ok, O)
for m = O.dimensions
for lmin = O.lmins
    out_file = fullfile(result_dir, sprintf('Processing_RQA_RRfixed_m%d_td%d_lmin%d.txt', m, O.tau, lmin));
    rows = {};
    hdr  = {'Subject','Electrode','Condition','Run','Duration_s', ...
            'RR','DET','Lmean','Lmax','ENTR','LAM','Vmax', ...
            'Dimension','Lmin','Radius','Tau','Fs','RR_target','RR_reached','RR_abs_err','N_iter'};

    for k = 1:numel(both_ok)
        sid_num  = parse_sid(subj_tags, both_ok(k));
        subj_tag = sprintf('S%03d', sid_num);
        try
            [raw1, raw2, Fs] = load_pair(data_dir, subj_tag);
            [dat1, dat2] = ft_filter_pair(raw1, raw2, Fs, O.hp_on,O.hp_freq,O.lp_on,O.lp_freq,O.notch_on,O.notch_band);
            [XEO, XEC, L] = common_begin_window(dat1, dat2, Fs, O.common_max_s, O.min_common_s);
            if L < O.min_common_s*Fs, continue; end

            [labels, idx_keep] = intersect_clean_eeg_labels(dat1.label, dat2.label, O.exclude_regex);

            for ei = 1:numel(idx_keep)
                e = idx_keep(ei);
                % EO
                [R_adj, r_adj, rr_adj, n_it] = rqa_radius_for_fixed_rr(XEO(e,:), m, O.tau, lmin, O.theiler, O.target_rr, O.tol, O.max_iter);
                rows(end+1,:) = {subj_tag, labels{ei}, 'EyesOpen','R01', L/Fs, ...
                    R_adj(1),R_adj(2),R_adj(3),R_adj(4),R_adj(5),R_adj(6),R_adj(7), ...
                    m,lmin,r_adj,O.tau,Fs,O.target_rr, rr_adj, abs(rr_adj-O.target_rr), n_it};

                % EC
                [R_adj, r_adj, rr_adj, n_it] = rqa_radius_for_fixed_rr(XEC(e,:), m, O.tau, lmin, O.theiler, O.target_rr, O.tol, O.max_iter);
                rows(end+1,:) = {subj_tag, labels{ei}, 'EyesClosed','R02', L/Fs, ...
                    R_adj(1),R_adj(2),R_adj(3),R_adj(4),R_adj(5),R_adj(6),R_adj(7), ...
                    m,lmin,r_adj,O.tau,Fs,O.target_rr, rr_adj, abs(rr_adj-O.target_rr), n_it};
            end
            fprintf('[%s] RQA-RR: %.1f s, %d ch, Fs=%.1f\n', subj_tag, L/Fs, numel(idx_keep), Fs);
        catch ME
            fprintf(2,'[%s] RQA-RR ERROR: %s -> skip\n', subj_tag, ME.message);
        end
    end

    if ~isempty(rows)
        T = cell2table(rows,'VariableNames',hdr);
        writetable(T, out_file, 'Delimiter','\t');
        fprintf('>> Saved: %s (%d rows)\n', out_file, height(T));
    else
        fprintf(2,'[RQA-RR] No rows for m=%d, lmin=%d\n', m, lmin);
    end
end
end
end

% ----------------------------- PSD (Welch) ---------------------------------
function run_psd_welch(data_dir, result_dir, subj_tags, both_ok, O, use_waitbar)
fig_dir = fullfile(result_dir,'PSD'); if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
rows = {};
hdr = {'Subject','Electrode','Condition','Run','Duration_s','Band', ...
       'Absolute_Total_Power','Relative_Total_Power','Mean_Frequency'};
out_txt = fullfile(result_dir,'Processing_PSD.txt');

hwb = []; if use_waitbar, hwb = waitbar(0,'Initializing…','Name','PSD Analysis'); end

have_ref = false; labels_ref = {}; f_ref = []; sumEO=[]; sumEC=[]; n_by_elec=[];

for k = 1:numel(both_ok)
    if use_waitbar, waitbar(k/numel(both_ok), hwb, sprintf('Processing %d/%d',k,numel(both_ok))); end
    sid_num  = parse_sid(subj_tags, both_ok(k));
    subj_tag = sprintf('S%03d', sid_num);
    try
        [raw1, raw2, Fs] = load_pair(data_dir, subj_tag);
        [dat1, dat2] = ft_filter_pair(raw1, raw2, Fs, O.hp_on,O.hp_freq,O.lp_on,O.lp_freq,O.notch_on,O.notch_band);
        [XEO, XEC, L] = common_begin_window(dat1, dat2, Fs, O.common_max_s, O.min_common_s);
        if L < O.min_common_s*Fs, continue; end

        labels = clean_labels(dat1.label);
        nChan  = numel(labels);

        nwin  = max(1, round(O.welch_win_s*Fs));
        nover = round(O.welch_overlap*nwin);

        y_max = 0; FF = []; PEO = cell(nChan,1); PEC = cell(nChan,1);

        for e = 1:nChan
            [Pxx_EO, f] = pwelch(XEO(e,:), nwin, nover, [], Fs);
            [Pxx_EC, ~] = pwelch(XEC(e,:), nwin, nover, [], Fs);
            if isempty(FF), FF = f; end

            msk30 = (f >= O.f_plot_min) & (f <= O.f_plot_max);
            if any(msk30), y_max = max([y_max, max(Pxx_EO(msk30)), max(Pxx_EC(msk30))]); end

            % band summaries
            for b = 1:numel(O.bands.name)
                br = O.bands.rng(b,:); br(2) = min(br(2), Fs/2);
                mb = (f >= br(1)) & (f <= br(2)); if ~any(mb), continue; end
                f_band = f(mb);

                abs_pow = sum(Pxx_EO(mb));
                rel_pow = abs_pow / sum(Pxx_EO);
                mean_f  = sum(f_band .* Pxx_EO(mb)) / max(abs_pow,eps);
                rows(end+1,:) = {subj_tag, labels{e}, 'EyesOpen','R01', L/Fs, O.bands.name{b}, abs_pow, rel_pow, mean_f};

                abs_pow = sum(Pxx_EC(mb));
                rel_pow = abs_pow / sum(Pxx_EC);
                mean_f  = sum(f_band .* Pxx_EC(mb)) / max(abs_pow,eps);
                rows(end+1,:) = {subj_tag, labels{e}, 'EyesClosed','R02', L/Fs, O.bands.name{b}, abs_pow, rel_pow, mean_f};
            end

            PEO{e} = Pxx_EO; PEC{e} = Pxx_EC;
        end

        % per-subject 64ch figure
        plot_psd_subject(fig_dir, subj_tag, FF, PEO, PEC, labels, O.f_plot_min, O.f_plot_max, y_max);

        % accumulate for group mean
        if ~have_ref
            f_ref = FF(:).'; labels_ref = labels(:);
            sumEO = zeros(nChan, numel(f_ref)); sumEC = zeros(nChan, numel(f_ref));
            n_by_elec = zeros(nChan,1); have_ref = true;
        else
            if numel(FF) ~= numel(f_ref) || any(abs(FF(:)-f_ref(:))>1e-12)
                for ee = 1:nChan
                    PEO{ee} = interp1(FF, PEO{ee}, f_ref, 'linear', 'extrap');
                    PEC{ee} = interp1(FF, PEC{ee}, f_ref, 'linear', 'extrap');
                end
            end
        end
        for ee = 1:nChan
            sumEO(ee,:) = sumEO(ee,:) + reshape(PEO{ee},1,[]);
            sumEC(ee,:) = sumEC(ee,:) + reshape(PEC{ee},1,[]);
        end
        n_by_elec = n_by_elec + 1;

    catch ME
        fprintf(2,'[%s] PSD ERROR: %s -> skip\n', subj_tag, ME.message);
    end
end

% write table
if ~isempty(rows)
    T = cell2table(rows,'VariableNames',hdr);
    writetable(T, out_txt, 'Delimiter','\t');
    fprintf('>> Saved: %s (%d rows)\n', out_txt, height(T));
end

% group mean figure
if have_ref && any(n_by_elec>0)
    meanEO = sumEO ./ max(n_by_elec,1);
    meanEC = sumEC ./ max(n_by_elec,1);
    plot_psd_groupmean(fig_dir, f_ref, meanEO, meanEC, labels_ref, O.f_plot_min, O.f_plot_max);
end

if ~isempty(hwb) && ishghandle(hwb), try close(hwb); catch, end; end
end

% ----------------------------- Sample Entropy ------------------------------
function run_sampen(data_dir, result_dir, subj_tags, both_ok, O, use_waitbar)
out_file = fullfile(result_dir, 'Processing_SampEn.txt');
rows = {};
hdr = {'Subject','Electrode','Condition','Run','Duration_s','SampEn','StdErr','Dimension','Fs'};

hwb = []; tStart = tic;
if use_waitbar
    hwb = waitbar(0,'Initializing…','Name','SampEn Analysis');
    setappdata(hwb,'Canceling',0);
end

for k = 1:numel(both_ok)
    sid_num  = parse_sid(subj_tags, both_ok(k));
    subj_tag = sprintf('S%03d', sid_num);
    if use_waitbar
        p = k/numel(both_ok);
        elapsed = toc(tStart);
        remain  = elapsed*(1-p)/max(p,eps);
        waitbar(p, hwb, sprintf('Processing %d/%d | %.1f%% | ETA %s', k, numel(both_ok), 100*p, seconds(remain)));
    end

    try
        [raw1, raw2, Fs] = load_pair(data_dir, subj_tag);
        [dat1, dat2] = ft_filter_pair(raw1, raw2, Fs, O.hp_on,O.hp_freq,O.lp_on,O.lp_freq,O.notch_on,O.notch_band);
        [XEO, XEC, L] = common_begin_window(dat1, dat2, Fs, O.common_max_s, O.min_common_s);
        if L < O.min_common_s*Fs, continue; end

        labels = clean_labels(dat1.label);
        for e = 1:numel(labels)
            % EO
            x = double(XEO(e,:));
            [E, SE] = sampen(x, O.m, O.r_frac, O.dist, O.norm, O.vec);   % check your sampen signature
            rows(end+1,:) = {subj_tag, labels{e}, 'EyesOpen','R01', L/Fs, E(O.m), SE(O.m), O.m, Fs};
            % EC
            x = double(XEC(e,:));
            [E, SE] = sampen(x, O.m, O.r_frac, O.dist, O.norm, O.vec);
            rows(end+1,:) = {subj_tag, labels{e}, 'EyesClosed','R02', L/Fs, E(O.m), SE(O.m), O.m, Fs};
        end
        fprintf('[%s] SampEn: %.1f s, %d ch, Fs=%.1f\n', subj_tag, L/Fs, numel(labels), Fs);
    catch ME
        fprintf(2,'[%s] SampEn ERROR: %s -> skip\n', subj_tag, ME.message);
    end
end

if ~isempty(rows)
    T = cell2table(rows,'VariableNames',hdr);
    writetable(T, out_file, 'Delimiter','\t');
    fprintf('>> Saved: %s (%d rows)\n', out_file, height(T));
end

if ~isempty(hwb) && ishghandle(hwb), try close(hwb); catch, end; end
end

% ============================== UTILITIES =================================
function sid_num = parse_sid(subj_tags, idx_in_bothok)
sid_num = str2double( regexprep(subj_tags{idx_in_bothok}, '[^\d]', '') );
end

function [raw1, raw2, Fs] = load_pair(data_dir, subj_tag)
cfg = []; cfg.dataset = fullfile(data_dir, subj_tag, sprintf('%sR01.edf', subj_tag));
raw1 = ft_preprocessing(cfg);
cfg = []; cfg.dataset = fullfile(data_dir, subj_tag, sprintf('%sR02.edf', subj_tag));
raw2 = ft_preprocessing(cfg);
Fs = raw1.fsample;
if abs(Fs - raw2.fsample) > 1e-6
    error('Fs mismatch: R01=%.6f, R02=%.6f', Fs, raw2.fsample);
end
end

function [dat1, dat2] = ft_filter_pair(raw1, raw2, Fs, hp_on,hp, lp_on,lp, notch_on,notch_band)
Nyq = Fs/2;
cfgf = [];
if hp_on && hp>0 && hp<Nyq, cfgf.hpfilter='yes'; cfgf.hpfreq=hp; else, cfgf.hpfilter='no'; end
if lp_on && lp>0 && lp<Nyq, cfgf.lpfilter='yes'; cfgf.lpfreq=lp; else, cfgf.lpfilter='no'; end
if notch_on && Nyq>50 && notch_band(1)>0 && notch_band(2)<Nyq
    cfgf.bsfilter='yes'; cfgf.bsfreq=notch_band; cfgf.bsinstabilityfix='reduce';
else
    cfgf.bsfilter='no';
end
dat1 = ft_preprocessing(cfgf, raw1);
dat2 = ft_preprocessing(cfgf, raw2);
end

function [XEO, XEC, L] = common_begin_window(dat1, dat2, Fs, max_s, min_s)
N1 = size(dat1.trial{1},2); N2 = size(dat2.trial{1},2);
dur = min([N1,N2]/Fs);
dur = min(dur, max_s);
if dur < min_s
    L=0; XEO=[]; XEC=[]; return;
end
L = floor(dur*Fs);
XEO = dat1.trial{1}(:,1:L);
XEC = dat2.trial{1}(:,1:L);
end

function labels = clean_labels(labels_raw)
labels = cellfun(@(s) regexprep(s,'[^A-Za-z0-9]+',''), labels_raw, 'UniformOutput', false);
end

function [labels_keep, idx_keep] = intersect_clean_eeg_labels(lbl1_raw, lbl2_raw, exclude_regex)
clean = @(C) cellfun(@(s) regexprep(s,'[^A-Za-z0-9]+',''), C, 'UniformOutput', false);
upperc = @(C) cellfun(@upper, C, 'UniformOutput', false);
L1 = upperc(clean(lbl1_raw)); L2 = upperc(clean(lbl2_raw));
common_up = intersect(L1, L2, 'stable');
is_bad = ~cellfun('isempty', regexp(common_up, exclude_regex, 'once'));
common_up = common_up(~is_bad);
[~, idx_keep] = ismember(common_up, L1);
labels_keep = clean(lbl1_raw(idx_keep));
end

function lag = first_min_lag_from_negpeaks(ami_vec)
[~, locs] = findpeaks(-ami_vec);
if ~isempty(locs), lag = locs(1)-1; else, lag = NaN; end
end

function lag = acf_threshold_lag(x, thr)
[ac, lags] = xcorr(x, 'normalized');
pos = lags >= 0; L = lags(pos); A = ac(pos);
idx = find(A <= thr, 1, 'first'); if isempty(idx), idx = 1; end
lag = L(idx);
end

function plot_time_delay_demo_3x2(pdf_file, sid, elab, xEO, ami_eo, td_acf_eo, thr, xEC, ami_ec, td_acf_ec)
colBlue = [100 149 237]/255; colRed  = [238  44  44]/255;
f = figure('Visible','off');
subplot(3,2,1); plot(xEO,'Color',colBlue); grid on; title(sprintf('S%d - %s - EO (20 s)', sid, elab));
subplot(3,2,3); [ac1,l1]=xcorr(xEO,'normalized'); p=l1>=0; plot(l1(p),ac1(p),'Color',colBlue); hold on; xline(td_acf_eo,'--k',sprintf('%g',td_acf_eo),'LineWidth',2); yline(thr,':k','1/e','LineWidth',2); xlim([0,40]); title('EO — ACF');
subplot(3,2,5); plot(0:numel(ami_eo)-1, ami_eo,'Color',colBlue); grid on; xlim([0,40]); title('EO — AMI');
subplot(3,2,2); plot(xEC,'Color',colRed);  grid on; title(sprintf('S%d - %s - EC (20 s)', sid, elab));
subplot(3,2,4); [ac2,l2]=xcorr(xEC,'normalized'); p=l2>=0; plot(l2(p),ac2(p),'Color',colRed);  hold on; xline(td_acf_ec,'--k',sprintf('%g',td_acf_ec),'LineWidth',2); yline(thr,':k','1/e','LineWidth',2); xlim([0,40]); title('EC — ACF');
subplot(3,2,6); plot(0:numel(ami_ec)-1, ami_ec,'Color',colRed); grid on; xlim([0,40]); title('EC — AMI');
exportgraphics(f, pdf_file, 'ContentType','vector','Append',true); close(f);
end

function plot_td_global_histograms(pdf_file, EO_AMI, EC_AMI, EO_ACF, EC_ACF)
colBlue = [100 149 237]/255; colRed  = [238  44  44]/255;
all_eo_ami = cell2mat(EO_AMI(~cellfun('isempty',EO_AMI)));
all_ec_ami = cell2mat(EC_AMI(~cellfun('isempty',EC_AMI)));
all_eo_acf = cell2mat(EO_ACF(~cellfun('isempty',EO_ACF)));
all_ec_acf = cell2mat(EC_ACF(~cellfun('isempty',EC_ACF)));
fGA = figure('Visible','off');
subplot(2,2,1); histogram(all_eo_ami,'BinWidth',1,'FaceColor',colBlue);  title('GLOBAL — EO (AMI)'); grid on;
subplot(2,2,2); histogram(all_eo_acf,'BinWidth',1,'FaceColor',colBlue);  title('GLOBAL — EO (ACF)'); grid on;
subplot(2,2,3); histogram(all_ec_ami,'BinWidth',1,'FaceColor',colRed );  title('GLOBAL — EC (AMI)'); grid on;
subplot(2,2,4); histogram(all_ec_acf,'BinWidth',1,'FaceColor',colRed );  title('GLOBAL — EC (ACF)'); grid on;
exportgraphics(fGA, pdf_file, 'ContentType','vector','Append',true); close(fGA);

fGB = figure('Visible','off');
histogram(all_eo_acf,'BinWidth',1,'FaceAlpha',0.5,'FaceColor',colBlue); hold on;
histogram(all_ec_acf,'BinWidth',1,'BinWidth',1,'FaceAlpha',0.5,'FaceColor',colRed);
title('GLOBAL — ACF overlay (EO vs EC)'); grid on;
exportgraphics(fGB, pdf_file, 'ContentType','vector','Append',true); close(fGB);
end

function plot_psd_subject(fig_dir, subj_tag, f, PEO, PEC, labels, fmin, fmax, y_max)
nChan = numel(labels); cols = ceil(sqrt(nChan)); rows = ceil(nChan/cols);
msk = (f>=fmin)&(f<=fmax);
fh = figure('Name',sprintf('PSD — %s (EO vs EC)',subj_tag),'Color','w','Position',[100 100 1400 900]);
tl = tiledlayout(rows, cols, 'Padding','compact','TileSpacing','compact');
for e = 1:nChan
    nexttile;
    plot(f(msk), PEO{e}(msk), '-'); hold on;
    plot(f(msk), PEC{e}(msk), '-'); hold off;
    xlim([fmin fmax]); if y_max>0, ylim([0 y_max]); end
    title(labels{e}, 'Interpreter','none');
    if e > (rows-1)*cols, xlabel('Hz'); end
    if mod(e-1, cols)==0, ylabel('PSD'); end
end
title(tl, sprintf('PSD (Welch) — %s — EO vs EC', subj_tag));
out = fullfile(fig_dir, sprintf('PSD_%s_64ch', subj_tag));
savefig(fh, [out '.fig']); exportgraphics(fh, [out '.png'], 'Resolution', 200); close(fh);
end

function plot_psd_groupmean(fig_dir, f, meanEO, meanEC, labels, fmin, fmax)
nChan = numel(labels); cols = ceil(sqrt(nChan)); rows = ceil(nChan/cols);
msk = (f>=fmin)&(f<=fmax);
ymax = max([max(meanEO(:,msk),[],'all'), max(meanEC(:,msk),[],'all')]); if ~isfinite(ymax)||ymax<=0, ymax=1; end
fh = figure('Name','PSD Mean — 64ch (EO vs EC)','Color','w','Position',[80 80 1400 900]);
tl = tiledlayout(rows, cols, 'Padding','compact','TileSpacing','compact');
for e = 1:nChan
    nexttile;
    plot(f(msk), meanEO(e,msk), '-'); hold on;
    plot(f(msk), meanEC(e,msk), '-'); hold off;
    xlim([fmin fmax]); ylim([0 ymax]);
    title(labels{e}, 'Interpreter','none');
    if e > (rows-1)*cols, xlabel('Hz'); end
    if mod(e-1, cols)==0, ylabel('PSD'); end
end
title(tl, 'PSD mean (Welch) — EO vs EC');
out = fullfile(fig_dir, 'PSD_mean_64ch');
savefig(fh, [out '.fig']); exportgraphics(fh, [out '.png'], 'Resolution', 200); close(fh);
end

% ----------------------------- RQA helpers --------------------------------
function [R, radius] = rqa_with_opt_radius(x, m, tau, lmin, theiler)
% Heuristic radius from embedded space, then compute RQA
[Y, ~] = phasespace1(x, m, tau);
radius = radsel_mai2024(Y, 2, m);
R = RQAmac(x, m, tau, radius, lmin, theiler);
end

function [R_adj, r_adj, rr_adj, n_iter] = rqa_radius_for_fixed_rr(x, m, tau, lmin, theiler, target_rr, tol, max_iter)
% Bisection on radius to reach target RR (robust bounds expansion)
[Y, ~] = phasespace1(x, m, tau);
% Start radius heuristic (slightly different per your practice)
if m==1, r0 = radsel_mai2024(Y, 2, m)/2; else, r0 = radsel_mai2024(Y, 2, m)*2; end
[R0, rr0] = rqa_rr(x, m, tau, r0, lmin, theiler);
n_iter = 1;

if ~isfinite(rr0)
    r0 = max(r0*0.5, eps);
    [R0, rr0] = rqa_rr(x, m, tau, r0, lmin, theiler);
    n_iter = n_iter + 1;
end

r_lo = r0; rr_lo = rr0; r_hi = r0; rr_hi = rr0;
grow=1.5; shrink=0.5;

while rr_hi < target_rr && isfinite(rr_hi) && n_iter < max_iter
    r_hi = r_hi*grow; [R1, rr_hi] = rqa_rr(x, m, tau, r_hi, lmin, theiler); n_iter=n_iter+1;
    if ~isfinite(rr_hi), r_hi=r_hi/grow; break; end
end
while rr_lo > target_rr && isfinite(rr_lo) && n_iter < max_iter
    r_lo = max(r_lo*shrink, eps); [R1, rr_lo] = rqa_rr(x, m, tau, r_lo, lmin, theiler); n_iter=n_iter+1;
    if ~isfinite(rr_lo), r_lo=r_lo/shrink; break; end
end

if ~(isfinite(rr_lo)&&isfinite(rr_hi)&&rr_lo<=target_rr&&rr_hi>=target_rr)
    [~,ix] = min([abs(rr_lo-target_rr), abs(rr_hi-target_rr)]);
    if ix==1, r_adj=r_lo; R_adj=rqa_rr_full(x,m,tau,r_adj,lmin,theiler); rr_adj=R_adj(1);
    else,     r_adj=r_hi; R_adj=rqa_rr_full(x,m,tau,r_adj,lmin,theiler); rr_adj=R_adj(1);
    end
    return;
end

for it = 1:max_iter
    r_mid = 0.5*(r_lo + r_hi);
    [R_mid, rr_mid] = rqa_rr(x, m, tau, r_mid, lmin, theiler); n_iter=n_iter+1;
    if ~isfinite(rr_mid)
        r_hi = r_mid; continue;
    end
    if rr_mid < target_rr, r_lo=r_mid; rr_lo=rr_mid; else, r_hi=r_mid; rr_hi=rr_mid; end
    if abs(rr_mid-target_rr)<=tol || abs(r_hi-r_lo)<=eps
        r_adj = r_mid; R_adj = rqa_rr_full(x,m,tau,r_adj,lmin,theiler); rr_adj=R_adj(1); return;
    end
end
% fallback
[~,ix] = min([abs(rr_lo-target_rr), abs(rr_hi-target_rr)]);
if ix==1, r_adj=r_lo; else, r_adj=r_hi; end
R_adj = rqa_rr_full(x,m,tau,r_adj,lmin,theiler); rr_adj=R_adj(1);
end

function [R, rr] = rqa_rr(x, m, tau, r, lmin, theiler)
R = RQAmac(x, m, tau, r, lmin, theiler);
rr = R(1);
end

function R = rqa_rr_full(x, m, tau, r, lmin, theiler)
R = RQAmac(x, m, tau, r, lmin, theiler);
end
