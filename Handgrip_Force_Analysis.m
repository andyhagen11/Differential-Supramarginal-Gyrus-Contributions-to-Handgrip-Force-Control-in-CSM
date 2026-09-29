%% Handgrip Force Metrics for fMRI Task 
% Participants performed repeated handgrip squeezes during fMRI. 4 blocks of 18 seconds with 18 seconds of rest in between 
% Analyzes handgrip force (kgf) from physio .txt files 
% 
% 1. 4th Order Zero-Lag Butterworth Lowpass Filter (10Hz)
% 2. Baseline correction based on rest blocks (zeroing)
% 3. Peak/valley detection and block segmentation
% 4. Force, timing and contraction shape metrics 
% 5. Fatigability measures (Linear/Exp fits, fatigability index, AUC fatigability index)
% 6. Participant figures and XLSX and Group-level exports

clear; 

% Configuration 
root_dir = '/Volumes/orthofbfs/SPARC/Research_Drive/Hagen/CSM_fMRI/Handgrip/Raw_Data';
output_dir = '/Volumes/orthofbfs/SPARC/Research_Drive/Hagen/CSM_fMRI/Handgrip/Analysis_Results/';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

FS = 2000; 
FILTER_CUTOFF = 10; % Hz

% Load in previously processed data (if it exists) so you can continue appending group results
if exist(fullfile(output_dir, 'GroupDataBackup.mat'), 'file')
    load(fullfile(output_dir, 'GroupDataBackup.mat'), 'DataStruct', 'group_metrics_tbl', 'group_cycle_tbl');
    disp('Loaded existing group data files');
else
    DataStruct = struct();
    group_metrics_tbl = table();
    group_cycle_tbl = table();
    disp('No backup found, starting group data fresh.');
end

%% Loop through subject files
files = dir(fullfile(root_dir, '*.txt'));

for i = 1:numel(files)
    file_name = files(i).name;
    tokens = split(file_name, '_');
    subj = tokens{1};
    subj_field = matlab.lang.makeValidName(subj);

     % Skip subjects you already processed
    if isfield(DataStruct, subj_field) 
        fprintf('Skipping %s (already processed)\n', subj);
        continue
    end
    fprintf('Processing %s...\n', subj);
    
    raw_data = readmatrix(fullfile(root_dir, file_name), 'FileType', 'text', 'NumHeaderLines', 15);
    n_samples = size(raw_data, 1);
    time_sec  = (0 : n_samples-1).' / FS;  % perfectly uniform, no rounding error like in txt data file
    force_raw = raw_data(:,3) * 9.80665;  % convert kgf to Newtons

    if strcmp(subj,"sub-023"); force_raw = force_raw(1216480:end); end % Fix for acquisition that started early

    if strcmp(subj,"sub-011")|| strcmp(subj,"sub-016"); AMP_THRESH = 20; DIST_THRESH = 400; % Participant fix for smaller amplitude peaks
    elseif strcmp(subj,"sub-019"); AMP_THRESH = 10; DIST_THRESH = 200;
    elseif strcmp(subj,"sub-021"); AMP_THRESH = 100; DIST_THRESH = 400;
    else 
        AMP_THRESH = 35; % 35 N as default
        DIST_THRESH = 400; % data points
    end

    %% 1. Filtering: 4th Order Zero-Lag Butterworth
    [b, a] = butter(4, FILTER_CUTOFF/(FS/2), 'low');
    % filtfilt provides zero-phase (zero-lag) filtering
    filtered_force = filtfilt(b, a, force_raw);

    % Also tried a Savitzky–Golay filter which could be appropriate for
    % this data
    % % Rule of thumb: window = ~50-100ms, polynomial order 2-4
    % window_ms = 100;   % ms
    % window_pts = round(window_ms/1000 * FS); % convert to samples
    % if mod(window_pts, 2) == 0
    %     window_pts = window_pts + 1; % must be odd
    % end
    % poly_order = 3;
    % filtered_force = sgolayfilt(force_raw, poly_order, window_pts);

    %% 2. Baseline Correction (Zeroing)
    % Identify rest blocks: gaps where no peaks exist
    [temp_pks, temp_locs] = findpeaks(filtered_force, 'MinPeakProminence', AMP_THRESH, 'MinPeakDistance', DIST_THRESH);
    % Find regions between peaks that are longer than 3 seconds
    diff_pks = diff(temp_locs);
    % This is a simplified approximation of rest blocks
    % We define rest as the signal in the gaps between the first and last peak that are significantly wide.
    rest_indices = [];
    for j = 1:length(diff_pks)
        if diff_pks(j) > 10000 % 5 seconds at 2000Hz
            rest_indices = [rest_indices, temp_locs(j)+1 : temp_locs(j+1)-1];
        end
    end
    baseline_from_rest = prctile(filtered_force(rest_indices), 5); % anchors 0 to the bottom 5% of the rest signal
    task_indices = setdiff(temp_locs(1):temp_locs(end), rest_indices); 
    baseline_from_task = min(filtered_force(task_indices)); % Get minimum from task (non-rest) signal as alternative floor
    baseline_val = min(baseline_from_rest, baseline_from_task); % Pick lowest between rest and task baselines
    force_signal = filtered_force - baseline_val;
    if strcmp(subj,"sub-021") % Participant fix for an outlier valley (take second lowest instead)
        baseline_from_task = mink(filtered_force(task_indices), 2);
        force_signal = filtered_force - baseline_from_task(2); 
    end 
 
    %% 3. Peak/Valley Detection and Block Segmentation
    % Re-run peak detection on zeroed/filtered signal
    % Rule: Peaks/Valleys must be separated by at least 200ms (400 pts at 2000Hz) and exceed a promineince of 35 N
    [pks, pks_locs] = findpeaks(force_signal, 'MinPeakProminence', AMP_THRESH, 'MinPeakDistance', DIST_THRESH); 
    val_locs_array = islocalmin(force_signal, 'MinSeparation', DIST_THRESH, 'MinProminence', AMP_THRESH); 
    val_locs = find(val_locs_array);

    % Find start of trial (looking for first window of 500 pts increasing)
    st_pt = 1;
    window_size = 500;
    if strcmp(subj,"sub-001"); window_size = 1000; end % participant fix
    if strcmp(subj,"sub-023"); window_size = 250; end % participant fix
    for k = 1:(length(force_signal) - window_size)
        window = force_signal(k:k+window_size);
        if all(diff(window) > 0)
            st_pt = k;
            break;
        end
    end        
    
    % Block Segmentation
    % Identify blocks by finding gaps between peaks > 5 seconds (10000 pts)
    block_sts_pks = [];
    block_ends_pks = [];
    
    diff_pks = diff(pks_locs);
    block_sts_pks = pks_locs(1);
    ct = 1;
    for m = 1:length(pks_locs)-1
        if diff_pks(m) > 10000
            ct = ct + 1;
            block_ends_pks(ct-1) = pks_locs(m);
            block_sts_pks(ct) = pks_locs(m+1);
        end
    end
    block_ends_pks(ct) = pks_locs(end);

    % Refine block ends
    % Search 5 seconds after last block peak and find the local minimum once you approach zero
    block_ends = [];
    zero_threshold = 0.05 * max(force_signal);  
    sustain_window = 500; 
    if strcmp(subj,"sub-007"); zero_threshold = 0.15 * max(force_signal); end
    if strcmp(subj,"sub-009") || strcmp(subj,"sub-017"); zero_threshold = 0.10 * max(force_signal); end
    if strcmp(subj,"sub-019"); zero_threshold = 0.30 * max(force_signal); sustain_window = 50; end
    if strcmp(subj,"sub-023"); zero_threshold = 0.10 * max(force_signal); sustain_window = 100; end
               
    for n = 1:length(block_ends_pks)
        search_start = block_ends_pks(n);
        search_end = block_ends_pks(n)+10000;
        search_sig = force_signal(search_start:search_end);
        % Step 1: find first sample that drops below threshold
        below_thresh = find(search_sig < zero_threshold, 1, 'first');
        if ~isempty(below_thresh)
            % Step 2: take 1000pt window after that point, find minimum
            win_start = below_thresh;
            win_end = below_thresh + sustain_window;
            [~, min_offset] = min(search_sig(win_start:win_end));
            end_pt_cur = search_start + win_start + min_offset - 2; % -2 corrects for two 1-based index offsets stacked
        else
            % Fallback: 0.5 seconds after last peak
            end_pt_cur = search_start + 500;
            warning('Block %d end fallback used for subject %s', n, subj);
        end
        block_ends(n) = end_pt_cur;
    end

    % Refine block starts
    % Search 5 seconds before first block peak and find first window of 500 data points that are all increasing
    block_sts = [];
    block_sts(1) = st_pt; % First block start already defined
    window_size = 500;
    if strcmp(subj,"sub-003") || strcmp(subj,"sub-005") || strcmp(subj,"sub-009") || strcmp(subj,"sub-011")|| strcmp(subj,"sub-016") || strcmp(subj,"sub-018") || strcmp(subj,"sub-019")|| strcmp(subj,"sub-022") || strcmp(subj,"sub-023") 
        window_size = 250; % participant fix for fast squeezers
    end 
    for p = 1:length(block_ends)-1
        start_pt_cur = [];
        search_range_start = force_signal(block_sts_pks(p+1)-10000 : block_sts_pks(p+1));
        for q = 1:(length(search_range_start) - window_size)
            window = search_range_start(q:q+window_size);
            if all(diff(window) > 0)
                start_pt_cur = q + (block_sts_pks(p+1)-10000);
                break;
            end
        end
        block_sts(p+1) = start_pt_cur;
    end
    
    % Number of cycles to include in early/late comparions - 3% of number of peaks rule (round up) with floor of 2
    n_cycles_pct = max(2, ceil(0.03 * numel(pks_locs)));
    fprintf('Subject %s: using %d cycles for early/late windows (3%% of %d total cycles)\n', subj, n_cycles_pct, numel(pks_locs));

    % Initialize subject-level outputs
    subj_dir = fullfile(output_dir, subj); 
    if ~exist(subj_dir, 'dir')
        mkdir(subj_dir); 
    end
    SubjFits = struct();
    Metrics = struct();
    subj_cycle_tbl = table();
    pooled_num_cycles = [];
    pooled_pk_forces = [];
    pooled_ISI = [];
    pooled_rise_time = [];
    pooled_RFD = [];
    pooled_fall_time = [];
    pooled_duty_cycle = [];
    signal_fig = figure;
    model_fig = figure;
    AUC_fig  = figure;

    % Loop through each block
    for ii = 1:length(block_sts)
        block_name = ['Block_', num2str(ii)];
        cur_block_pks = pks_locs(pks_locs > block_sts(ii) & pks_locs < block_ends(ii));
        cur_block_vals = val_locs(val_locs > block_sts(ii) & val_locs < block_ends(ii));
        
        % Delete any false valleys that occur before the first block peak
        bad_vals = cur_block_vals(cur_block_vals < cur_block_pks(1));
        
        if ~isempty(bad_vals)
            fprintf('Subject %s Block %d: removed %d valley(s) before first peak\n', ...
                subj, ii, numel(bad_vals));
            
            cur_block_vals(cur_block_vals < cur_block_pks(1)) = [];
        end

        % Delete any drop off final peaks (due to early termination at block end)
         if force_signal(cur_block_pks(end)) <= 0.65 * (force_signal(cur_block_pks(end-1))) % using 0.65 * prior peak as threshold
            % remove last peak
            cur_block_pks(end) = [];
            fprintf('Subject %s Block %d: removed last peak (termination artifact)\n', subj, ii);
            % set block end to last remaining valley, then remove last valley
            block_ends(ii) = cur_block_vals(end);
            cur_block_vals(end) = [];
         end
          
         % Fix for messy rest data to match block ends criteria above
         if strcmp(subj,'sub-019') && ii==2 
             block_ends(ii) = cur_block_vals(end-1) - 1717; 
             cur_block_vals = cur_block_vals(1:end-2);
             cur_block_pks = cur_block_pks(1:end-2);
         end 

        % Confirm that block has one more peak than valley so metrics are accurate
        assert(length(cur_block_vals) == length(cur_block_pks) - 1, 'Expected N-1 valleys for N peaks is not true');

        % Set peak number to start calculations with (always skip first peak is a ramp up and not a representative cycle)
        pk_skip = 2;  % default to skip first peak - keeping consistent for all participants

        
    %% 4. Force, Timing and Contraction Shape Metrics 
        % For counting cycles, we include all peaks here
        num_cycles_in_block = length(cur_block_pks);

        % Peak force per cycle - Skip first peak since often not a full contraction
        pk_forces = force_signal(cur_block_pks(pk_skip:end));
        pk_times = time_sec(cur_block_pks(pk_skip:end));
                
        % Inter-squeeze interval (ISI): current peak to subsequent peak 
        ISI = diff(time_sec(cur_block_pks(pk_skip:end)));
        
        % Rise Time: previous valley to current peak
        rise_time = time_sec(cur_block_pks(pk_skip:end)) - time_sec(cur_block_vals(pk_skip-1:end));

        % Rate of Force Development (RFD): Peak force minus valley force, divided by rise time (N/s)
        valley_forces = force_signal(cur_block_vals(pk_skip-1:end));
        RFD = (pk_forces - valley_forces) ./ rise_time;
        
        % Fall Time: current peak to subsequent valley - (no valley identified after last peak)
        fall_time = time_sec(cur_block_vals(pk_skip:end)) - time_sec(cur_block_pks(pk_skip:end-1));
       
        % Duty Cycle (Rise / Total Cycle) - exclude last rise time because missing fall time for last peak
        duty_cycle = rise_time(1:end-1) ./ (rise_time(1:end-1) + fall_time(1:end));

        % Store Block Metrics (Skipping 1st peak for Means/CV)
        Metrics.(block_name).NumCycles = num_cycles_in_block;
        Metrics.(block_name).PeakForce_Mean = mean(pk_forces);
        Metrics.(block_name).PeakForce_CV = std(pk_forces)/mean(pk_forces);
        Metrics.(block_name).ISI_Mean = mean(ISI);
        Metrics.(block_name).ISI_CV = std(ISI)/mean(ISI);
        Metrics.(block_name).RiseTime_Mean = mean(rise_time);
        Metrics.(block_name).RiseTime_CV = std(rise_time)/mean(rise_time);
        Metrics.(block_name).RFD_Mean = mean(RFD);
        Metrics.(block_name).RFD_CV = std(RFD) / mean(RFD);
        Metrics.(block_name).FallTime_Mean = mean(fall_time);
        Metrics.(block_name).FallTime_CV = std(fall_time)/mean(fall_time);
        Metrics.(block_name).DutyCycle_Mean = mean(duty_cycle);
        Metrics.(block_name).DutyCycle_CV = std(duty_cycle)/mean(duty_cycle);

        % Pooled metrics
        pooled_num_cycles = [pooled_num_cycles; num_cycles_in_block(:)];
        pooled_pk_forces = [pooled_pk_forces; pk_forces(:)];
        pooled_ISI = [pooled_ISI; ISI(:)];
        pooled_rise_time = [pooled_rise_time; rise_time(:)];
        pooled_RFD = [pooled_RFD; RFD(:)];
        pooled_fall_time = [pooled_fall_time; fall_time(:)];
        pooled_duty_cycle = [pooled_duty_cycle; duty_cycle(:)];

    %% 5. Fatigability Measures
        % Fatigability Model Fits (Linear & Exponential) for peak force across cycles
        % Note: by keeping x and y in each block aligned with trial timing,
        % model metrics(e.g. lm intercept and exp intital) become
        % meaningless/uncomparable
        x = pk_times(:);
        y = pk_forces(:);
        n = numel(y);

        % Linear Model
        lin_model = fitlm(x, y);
        yhat_lin = predict(lin_model, x);
        RSS_lin = sum((y - yhat_lin).^2);
        AIC_lin = n * log(RSS_lin/n) + 4;
        BIC_lin = n * log(RSS_lin/n) + 2*log(n);
        lin_intercept = lin_model.Coefficients.Estimate(1);
        lin_slope = lin_model.Coefficients.Estimate(2);
        R2_lin = lin_model.Rsquared.Ordinary;

        % Exponential Model
        % Model: y = A * exp(-k*x) + C
        exp_model_fn = @(b,x) b(1)*exp(-b(2)*x) + b(3);
        % Initial guesses: [Amplitude, Decay, Offset]
        t_range = x(end) - x(1);
        expected_pct_drop = 0.5;
        beta0 = [max(y)-min(y), -log(1 - expected_pct_drop)/t_range, min(y)];
        try
            warning('off','all') % Turn off warnings for bad fits
            exp_model = nlinfit(x, y, exp_model_fn, beta0); 
            warning('on','all')
            A_exp = exp_model(1); 
            k_exp = exp_model(2); 
            C_exp = exp_model(3);
            yhat_exp = exp_model_fn(exp_model, x);
            RSS_exp = sum((y - yhat_exp).^2);
            AIC_exp = n * log(RSS_exp/n) + 6;
            BIC_exp = n * log(RSS_exp/n) + 3*log(n);
            R2_exp = 1 - (RSS_exp/sum((y - mean(y)).^2));
        catch
            A_exp=NaN; k_exp=NaN; C_exp=NaN; AIC_exp=NaN; BIC_exp=NaN; R2_exp=NaN;
            exp_model = [NaN, NaN, NaN];
        end

        % Generate envelopes starting at 2nd peak (pk_times(1)) and ending at last peak
        t_start = pk_times(1);
        t_end = pk_times(end);
        block_t = (t_start : (1/FS) : t_end).'; % transpose to Nx1
        yhat_lin_block = predict(lin_model, block_t);
        yhat_exp_block = exp_model_fn(exp_model, block_t);

        % Store model metrics
        Metrics.(block_name).LM_Intercept = lin_intercept;
        Metrics.(block_name).LM_Slope = lin_slope;
        Metrics.(block_name).LM_R2 = R2_lin;
        Metrics.(block_name).LM_AIC = AIC_lin;
        Metrics.(block_name).LM_BIC = BIC_lin;
        Metrics.(block_name).Exp_A = A_exp;
        Metrics.(block_name).Exp_k = k_exp;
        Metrics.(block_name).Exp_C = C_exp;
        Metrics.(block_name).Exp_Initial = A_exp + C_exp;
        Metrics.(block_name).Exp_R2 = R2_exp;
        Metrics.(block_name).Exp_AIC = AIC_exp;
        Metrics.(block_name).Exp_BIC = BIC_exp;

        % Store model fits
        SubjFits.(block_name).BlockTime = block_t;
        SubjFits.(block_name).LinEnvelope = yhat_lin_block;
        SubjFits.(block_name).ExpEnvelope = yhat_exp_block;

        % Fatigability Index
        % Squeezes Early (after skip) vs Late per block - Peak Force, ISI, and Rise Time
        early_force = mean(pk_forces(1:n_cycles_pct));
        late_force = mean(pk_forces(end-n_cycles_pct+1:end));
        Metrics.(block_name).FatigIdx_Force = (early_force - late_force) / early_force;
        
        early_ISI = mean(ISI(1:n_cycles_pct));
        late_ISI = mean(ISI(end-n_cycles_pct+1:end));
        Metrics.(block_name).FatigIdx_ISI = (early_ISI - late_ISI) / early_ISI;

        early_RFD = mean(RFD(1:n_cycles_pct));
        late_RFD  = mean(RFD(end-n_cycles_pct+1:end));
        Metrics.(block_name).FatigIdx_RFD = (early_RFD - late_RFD) / early_RFD;
        
        if length(pk_forces) < 6; warning('Subject %s block %d only has %d peak(s)', subj, ii, num_cycles_in_block); end
        
        % Impulse (AUC) Fatigability Index
        % Actual AUC (from first valley to block end)
        actual_auc_idx = cur_block_vals(pk_skip-1) : block_ends(ii);
        t_actual = time_sec(actual_auc_idx);
        f_actual = force_signal(actual_auc_idx);
        actual_auc = trapz(t_actual, f_actual);

        % Calculate per cycle AUC for fatigability index
        cycle_auc = [];
        for c = (pk_skip-1):(pk_skip-1 + length(pk_forces)-1)
            v_start = cur_block_vals(c);
            if c == length(pk_forces)
                v_end = block_ends(ii); % make last cycle end at block end
            else
                v_end = cur_block_vals(c+1);
            end
            t_actual_seg = time_sec(v_start:v_end);
            f_actual_seg = force_signal(v_start:v_end);
            cycle_auc(c) = trapz(t_actual_seg, f_actual_seg);
        end
        early_auc = mean(cycle_auc(1:n_cycles_pct));
        late_auc  = mean(cycle_auc(end-n_cycles_pct+1 : end));
        cycle_auc = cycle_auc(:); % ensure column vector

        % Expected AUC - create template of early cycles and tile template from first valley to block end
        template_sum = 0;
        count_template = 0;
        n_template = round(early_ISI * FS);  % average number of samples in one cycle
        for c2 = (pk_skip-1):(pk_skip-1 + n_cycles_pct-1) % Early valley window
            v_start = cur_block_vals(c2);
            v_end = cur_block_vals(c2+1);
            t_seg = time_sec(v_start:v_end);
            f_seg = force_signal(v_start:v_end);
            t_new = linspace(t_seg(1), t_seg(end), n_template);
            resampled = interp1(t_seg, f_seg, t_new, 'linear');
            template_sum = template_sum + resampled;
            count_template = count_template + 1;
        end
        avg_template = template_sum / count_template;

        % Build expected curve tiled over block duration (time-based)
        block_duration = time_sec(block_ends(ii)) - t_actual(1);
        t_extended     = (t_actual(1) : 1/FS : t_actual(1) + block_duration).';
        expected_curve = zeros(size(t_extended));
        c2 = 1;
        while true
            t_cycle_start = t_extended(1) + (c2-1) * early_ISI;
            if t_cycle_start >= t_extended(end), break; end
            template_t = linspace(0, early_ISI, n_template);
            cycle_mask = t_extended >= t_cycle_start & t_extended < (t_cycle_start + early_ISI);
            t_local = t_extended(cycle_mask) - t_cycle_start;
            expected_curve(cycle_mask) = interp1(template_t, avg_template, t_local, 'linear', 0);
            c2 = c2 + 1;
        end
        % Use same curve for trapz calculation
        expected_auc = trapz(t_extended, expected_curve);

        % Store Block 1 for global use
        if ii == 1
            b1_early_auc = early_auc;
            b1_avg_template = avg_template;
            b1_early_ISI = early_ISI;
            b1_template_t = linspace(0, early_ISI, n_template);
        end

        % Store AUC measures
        Metrics.(block_name).AUC_Actual = actual_auc;
        Metrics.(block_name).AUC_Expected = expected_auc;
        Metrics.(block_name).FatigIdx_AUC_Actual = (early_auc - late_auc) / early_auc;
        Metrics.(block_name).FatigIdx_AUC_Expected = (expected_auc - actual_auc)/expected_auc;

        % Trail-level measures when the last block is reached 
        if ii == length(block_sts)
            % Trial Level Means and CVs
            Metrics.Trial.NumCycles = sum(pooled_num_cycles);
            Metrics.Trial.NumCycles_Early = n_cycles_pct;
            Metrics.Trial.PeakForce_Mean = mean(pooled_pk_forces);
            Metrics.Trial.PeakForce_CV = std(pooled_pk_forces) / mean(pooled_pk_forces);
            Metrics.Trial.ISI_Mean = mean(pooled_ISI);
            Metrics.Trial.ISI_CV = std(pooled_ISI) / mean(pooled_ISI);
            Metrics.Trial.RFD_Mean = mean(pooled_RFD);
            Metrics.Trial.RFD_CV = std(pooled_RFD) / mean(pooled_RFD);
            Metrics.Trial.RiseTime_Mean = mean(pooled_rise_time);
            Metrics.Trial.RiseTime_CV = std(pooled_rise_time) / mean(pooled_rise_time);
            Metrics.Trial.FallTime_Mean = mean(pooled_fall_time);
            Metrics.Trial.FallTime_CV = std(pooled_fall_time) / mean(pooled_fall_time);
            Metrics.Trial.DutyCycle_Mean = mean(pooled_duty_cycle);
            Metrics.Trial.DutyCycle_CV = std(pooled_duty_cycle) / mean(pooled_duty_cycle);Metrics.Trial.PeakForce_Mean = mean([Metrics.Block_1.PeakForce_Mean, Metrics.Block_2.PeakForce_Mean, Metrics.Block_3.PeakForce_Mean, Metrics.Block_4.PeakForce_Mean]);
            Metrics.Trial.LM_Slope_Mean = mean([Metrics.Block_1.LM_Slope, Metrics.Block_2.LM_Slope, Metrics.Block_3.LM_Slope, Metrics.Block_4.LM_Slope]);
            Metrics.Trial.LM_Slope_CV = std([Metrics.Block_1.LM_Slope, Metrics.Block_2.LM_Slope, Metrics.Block_3.LM_Slope, Metrics.Block_4.LM_Slope]) / mean([Metrics.Block_1.LM_Slope, Metrics.Block_2.LM_Slope, Metrics.Block_3.LM_Slope, Metrics.Block_4.LM_Slope]);
            Metrics.Trial.Exp_k_Mean = mean([Metrics.Block_1.Exp_k, Metrics.Block_2.Exp_k, Metrics.Block_3.Exp_k, Metrics.Block_4.Exp_k]);
            Metrics.Trial.Exp_k_CV = std([Metrics.Block_1.Exp_k, Metrics.Block_2.Exp_k, Metrics.Block_3.Exp_k, Metrics.Block_4.Exp_k]) / mean([Metrics.Block_1.Exp_k, Metrics.Block_2.Exp_k, Metrics.Block_3.Exp_k, Metrics.Block_4.Exp_k]);
            Metrics.Trial.FatigIdx_Force_Mean = mean([Metrics.Block_1.FatigIdx_Force, Metrics.Block_2.FatigIdx_Force, Metrics.Block_3.FatigIdx_Force, Metrics.Block_4.FatigIdx_Force]);
            Metrics.Trial.FatigIdx_Force_CV = std([Metrics.Block_1.FatigIdx_Force, Metrics.Block_2.FatigIdx_Force, Metrics.Block_3.FatigIdx_Force, Metrics.Block_4.FatigIdx_Force]) / mean([Metrics.Block_1.FatigIdx_Force, Metrics.Block_2.FatigIdx_Force, Metrics.Block_3.FatigIdx_Force, Metrics.Block_4.FatigIdx_Force]);
            Metrics.Trial.FatigIdx_ISI_Mean = mean([Metrics.Block_1.FatigIdx_ISI, Metrics.Block_2.FatigIdx_ISI, Metrics.Block_3.FatigIdx_ISI, Metrics.Block_4.FatigIdx_ISI]);
            Metrics.Trial.FatigIdx_ISI_CV = std([Metrics.Block_1.FatigIdx_ISI, Metrics.Block_2.FatigIdx_ISI, Metrics.Block_3.FatigIdx_ISI, Metrics.Block_4.FatigIdx_ISI]) / mean([Metrics.Block_1.FatigIdx_ISI, Metrics.Block_2.FatigIdx_ISI, Metrics.Block_3.FatigIdx_ISI, Metrics.Block_4.FatigIdx_ISI]);
            Metrics.Trial.FatigIdx_RFD_Mean = mean([Metrics.Block_1.FatigIdx_RFD, Metrics.Block_2.FatigIdx_RFD, Metrics.Block_3.FatigIdx_RFD, Metrics.Block_4.FatigIdx_RFD]);
            Metrics.Trial.FatigIdx_RFD_CV = std([Metrics.Block_1.FatigIdx_RFD, Metrics.Block_2.FatigIdx_RFD, Metrics.Block_3.FatigIdx_RFD, Metrics.Block_4.FatigIdx_RFD]) / mean([Metrics.Block_1.FatigIdx_RFD, Metrics.Block_2.FatigIdx_RFD, Metrics.Block_3.FatigIdx_RFD, Metrics.Block_4.FatigIdx_RFD]);
            Metrics.Trial.AUC_Actual_Sum = sum([Metrics.Block_1.AUC_Actual, Metrics.Block_2.AUC_Actual, Metrics.Block_3.AUC_Actual, Metrics.Block_4.AUC_Actual]);
            Metrics.Trial.AUC_Actual_Mean = mean([Metrics.Block_1.AUC_Actual, Metrics.Block_2.AUC_Actual, Metrics.Block_3.AUC_Actual, Metrics.Block_4.AUC_Actual]);
            Metrics.Trial.AUC_Actual_CV = std([Metrics.Block_1.AUC_Actual, Metrics.Block_2.AUC_Actual, Metrics.Block_3.AUC_Actual, Metrics.Block_4.AUC_Actual]) / mean([Metrics.Block_1.AUC_Actual, Metrics.Block_2.AUC_Actual, Metrics.Block_3.AUC_Actual, Metrics.Block_4.AUC_Actual]);
            Metrics.Trial.AUC_Expected_Sum = sum([Metrics.Block_1.AUC_Expected, Metrics.Block_2.AUC_Expected, Metrics.Block_3.AUC_Expected, Metrics.Block_4.AUC_Expected]);
            Metrics.Trial.AUC_Expected_Mean = mean([Metrics.Block_1.AUC_Expected, Metrics.Block_2.AUC_Expected, Metrics.Block_3.AUC_Expected, Metrics.Block_4.AUC_Expected]);
            Metrics.Trial.AUC_Expected_CV = std([Metrics.Block_1.AUC_Expected, Metrics.Block_2.AUC_Expected, Metrics.Block_3.AUC_Expected, Metrics.Block_4.AUC_Expected]) / mean([Metrics.Block_1.AUC_Expected, Metrics.Block_2.AUC_Expected, Metrics.Block_3.AUC_Expected, Metrics.Block_4.AUC_Expected]);
            Metrics.Trial.FatigIdx_AUC_Actual_Mean = mean([Metrics.Block_1.FatigIdx_AUC_Actual, Metrics.Block_2.FatigIdx_AUC_Actual, Metrics.Block_3.FatigIdx_AUC_Actual, Metrics.Block_4.FatigIdx_AUC_Actual]);
            Metrics.Trial.FatigIdx_AUC_Actual_CV = std([Metrics.Block_1.FatigIdx_AUC_Actual, Metrics.Block_2.FatigIdx_AUC_Actual, Metrics.Block_3.FatigIdx_AUC_Actual, Metrics.Block_4.FatigIdx_AUC_Actual]) / mean([Metrics.Block_1.FatigIdx_AUC_Actual, Metrics.Block_2.FatigIdx_AUC_Actual, Metrics.Block_3.FatigIdx_AUC_Actual, Metrics.Block_4.FatigIdx_AUC_Actual]);
            Metrics.Trial.FatigIdx_AUC_Expected_Mean = mean([Metrics.Block_1.FatigIdx_AUC_Expected, Metrics.Block_2.FatigIdx_AUC_Expected, Metrics.Block_3.FatigIdx_AUC_Expected, Metrics.Block_4.FatigIdx_AUC_Expected]);
            Metrics.Trial.FatigIdx_AUC_Expected_CV = std([Metrics.Block_1.FatigIdx_AUC_Expected, Metrics.Block_2.FatigIdx_AUC_Expected, Metrics.Block_3.FatigIdx_AUC_Expected, Metrics.Block_4.FatigIdx_AUC_Expected]) / mean([Metrics.Block_1.FatigIdx_AUC_Expected, Metrics.Block_2.FatigIdx_AUC_Expected, Metrics.Block_3.FatigIdx_AUC_Expected, Metrics.Block_4.FatigIdx_AUC_Expected]);
            Metrics.Trial.FatigIdx_AUC_Expected_Sum = (Metrics.Trial.AUC_Expected_Sum  - Metrics.Trial.AUC_Actual_Sum ) / Metrics.Trial.AUC_Expected_Sum ;

            % Trial Level Fatigability Indexes 
            % (Block 1 - Block 4) / Block 1
            Metrics.Trial.TrialFatigIdx_Force = (Metrics.Block_1.PeakForce_Mean - Metrics.Block_4.PeakForce_Mean) / Metrics.Block_1.PeakForce_Mean;
            Metrics.Trial.TrialFatigIdx_ISI = (Metrics.Block_1.ISI_Mean - Metrics.Block_4.ISI_Mean) / Metrics.Block_1.ISI_Mean;
            Metrics.Trial.TrialFatigIdx_RFD = (Metrics.Block_1.RFD_Mean - Metrics.Block_4.RFD_Mean) / Metrics.Block_1.RFD_Mean;
            Metrics.Trial.TrialFatigIdx_LM_Slope = (Metrics.Block_1.LM_Slope - Metrics.Block_4.LM_Slope) / Metrics.Block_1.LM_Slope; 
            Metrics.Trial.TrialFatigIdx_Exp_k = (Metrics.Block_1.Exp_k - Metrics.Block_4.Exp_k) / Metrics.Block_1.Exp_k;
            Metrics.Trial.TrialFatigIdx_AUC_Actual = (Metrics.Block_1.AUC_Actual - Metrics.Block_4.AUC_Actual) / Metrics.Block_1.AUC_Actual;
            Metrics.Trial.TrialFatigIdx_AUC_Expected = (Metrics.Block_1.AUC_Expected - Metrics.Block_4.AUC_Expected) / Metrics.Block_1.AUC_Expected;

            % Calculate global fatigability measures using compressed timeline across all blocks    
            % Goal: Remove rest gaps and 1st squeezes to create a continuous fatiguability timeline
            global_pk_times = [];
            global_pk_forces = [];
            global_val_times = [];
            global_val_forces = [];
            task_time_combined = [];
            task_sig_combined = [];
            cumulative_task_time = 0;
            
            % Once the last block is reached...
            for b_idx = 1:length(block_sts)
                % Identify peaks for this block
                b_pks = pks_locs(pks_locs > block_sts(b_idx) & pks_locs < block_ends(b_idx));
                b_vals = val_locs(val_locs > block_sts(b_idx) & val_locs < block_ends(b_idx));
                
                % Delete bad valleys that occur before block peak 1
                b_vals = b_vals(b_vals > b_pks(1));

                % Extract peaks and valleys (skip 1st peak)
                b_forces = force_signal(b_pks(pk_skip:end));
                b_valley_forces = force_signal(b_vals(pk_skip-1:end));  
                b_times = time_sec(b_pks(pk_skip:end)) - time_sec(b_vals(pk_skip-1));
                b_times_vals = time_sec(b_vals(pk_skip-1:end)) - time_sec(b_vals(pk_skip-1));
    
                % Append to global peak arrays using the cumulative offset
                global_pk_forces = [global_pk_forces; b_forces(:)];
                global_pk_times = [global_pk_times; (b_times(:) + cumulative_task_time)];
                global_val_forces = [global_val_forces; b_valley_forces(:)];
                global_val_times = [global_val_times; (b_times_vals(:) + cumulative_task_time)];
    
                % Extract force signal segment for this block
                seg_sig = force_signal(b_vals(pk_skip-1):block_ends(b_idx));
                seg_time = time_sec(b_vals(pk_skip-1):block_ends(b_idx));
                if b_idx == 1
                    seg_time = seg_time - seg_time(1);
                else
                    seg_time = seg_time - seg_time(1) + task_time_combined(end);
                end
    
                % Append to global signal arrays (exlcuding rest gaps)
                task_sig_combined = [task_sig_combined; seg_sig];
                task_time_combined = [task_time_combined; seg_time];
                
                % Update the offset for the next block
                % Offset = total time from 2nd peak to last peak of current block
                cumulative_task_time = cumulative_task_time + (time_sec(block_ends(b_idx)) - time_sec(b_vals(1)));
            end

            % Define model fit variables
            y_global = global_pk_forces(:);
            x_global = global_pk_times(:);
            n_global = length(y_global);
            
            % Global Linear Model 
            global_lin_model = fitlm(x_global, y_global);
            global_yhat_lin = predict(global_lin_model, x_global);
            global_RSS_lin = sum((y_global - global_yhat_lin).^2);
            
            global_lin_intercept = global_lin_model.Coefficients.Estimate(1);
            global_lin_slope = global_lin_model.Coefficients.Estimate(2);
            global_R2_lin = global_lin_model.Rsquared.Ordinary;
            global_AIC_lin = n_global * log(global_RSS_lin/n_global) + 4;
            global_BIC_lin = n_global * log(global_RSS_lin/n_global) + 2*log(n_global);
            
            % Global Exponential Model
            exp_model_fn = @(b,x) b(1)*exp(-b(2)*x) + b(3);
            % Initial guesses: [Amplitude, Decay, Offset]
            global_t_range = x(end) - x(1);
            global_expected_pct_drop = 0.6;
            global_beta0 = [max(y)-min(y), -log(1 - global_expected_pct_drop)/global_t_range, min(y)];
            
            try
                warning('off','all') % Turn off warnings for bad fits
                global_exp_model = nlinfit(x_global, y_global, exp_model_fn, global_beta0);
                warning('on','all') 
                global_yhat_exp = exp_model_fn(global_exp_model, x_global);
                global_RSS_exp = sum((y_global - global_yhat_exp).^2);
                global_R2_exp = 1 - (global_RSS_exp/sum((y_global - mean(y_global)).^2));
                
                global_exp_A = global_exp_model(1);
                global_exp_k = global_exp_model(2);
                global_exp_C = global_exp_model(3);
                global_exp_initial = global_exp_A + global_exp_C;
            catch
                global_exp_model = [NaN, NaN, NaN];
                global_RSS_exp = NaN; global_R2_exp = NaN;
                global_exp_A = NaN; global_exp_k = NaN; global_exp_C = NaN; global_exp_initial = NaN;
            end
            
            global_exp_AIC = n_global * log(global_RSS_exp/n_global) + 6;
            global_exp_BIC = n_global * log(global_RSS_exp/n_global) + 3*log(n_global);

            % Store global model metrics
            Metrics.Global.LM_Intercept = global_lin_intercept;
            Metrics.Global.LM_Slope = global_lin_slope;
            Metrics.Global.LM_R2 = global_R2_lin;
            Metrics.Global.LM_AIC = global_AIC_lin;
            Metrics.Global.LM_BIC = global_BIC_lin;
            Metrics.Global.Exp_A = global_exp_A;
            Metrics.Global.Exp_k = global_exp_k;
            Metrics.Global.Exp_C = global_exp_C;
            Metrics.Global.Exp_Initial = global_exp_initial;
            Metrics.Global.Exp_R2 = global_R2_exp;
            Metrics.Global.Exp_AIC = global_exp_AIC;
            Metrics.Global.Exp_BIC = global_exp_BIC;

            % Store model fits
            task_time_combined = task_time_combined(:); % ensure column
            SubjFits.Global.BlockTime = task_time_combined;
            SubjFits.Global.LinEnvelope = predict(global_lin_model, task_time_combined);
            SubjFits.Global.ExpEnvelope = exp_model_fn(global_exp_model, task_time_combined);

            % Global Fatigability Index
            % Take early cycles and late cycles from task global measures and calculate indexes
            global_early_force = mean(global_pk_forces(1:n_cycles_pct));
            global_late_force = mean(global_pk_forces(end-n_cycles_pct+1:end));
            Metrics.Global.FatigIdx_Force = (global_early_force - global_late_force) / global_early_force;
            
            global_early_ISI = mean(diff(global_pk_times(1:n_cycles_pct))); 
            global_late_ISI = mean(diff(global_pk_times(end-n_cycles_pct+1:end))); 
            Metrics.Global.FatigIdx_ISI = (global_early_ISI - global_late_ISI) / global_early_ISI;

            global_rise_times = global_pk_times - global_val_times;
            global_RFD = (global_pk_forces - global_val_forces) ./ global_rise_times;
            global_early_RFD = mean(global_RFD(1:n_cycles_pct));
            global_late_RFD  = mean(global_RFD(end-n_cycles_pct+1:end));
            Metrics.Global.FatigIdx_RFD = (global_early_RFD - global_late_RFD) / global_early_RFD;
            
            % Global Impulse (AUC) Fatigability Index
            % Golbal actual AUC
            global_actual_auc = trapz(task_time_combined, task_sig_combined); 

            % Global expected AUC - tile B1 template across full combined task duration
            global_duration = task_time_combined(end) - task_time_combined(1);
            t_global_extended = (task_time_combined(1) : 1/FS : task_time_combined(end)).';
            
            global_expected_curve = zeros(size(t_global_extended));
            c3 = 1;
            while true
                t_cycle_start = t_global_extended(1) + (c3-1) * b1_early_ISI;
                if t_cycle_start >= t_global_extended(end), break; end
                cycle_mask  = t_global_extended >= t_cycle_start & t_global_extended < (t_cycle_start + b1_early_ISI);
                t_local = t_global_extended(cycle_mask) - t_cycle_start;
                global_expected_curve(cycle_mask) = interp1(b1_template_t, b1_avg_template, t_local, 'linear', 0);
                c3 = c3 + 1;
            end
            global_expected_auc = trapz(t_global_extended, global_expected_curve);

            % Store global AUC metrics
            Metrics.Global.AUC_Actual = global_actual_auc;
            Metrics.Global.AUC_Expected = global_expected_auc;
            Metrics.Global.FatigIdx_AUC_Actual = (b1_early_auc - late_auc) / b1_early_auc;
            Metrics.Global.FatigIdx_AUC_Expected = (global_expected_auc - global_actual_auc) / global_expected_auc;
        end

    %% 6. Data and Figure Export
        % Collect per-cycle data for xlsx
        % NaN pad for equal lengths
        ISI = [ISI; NaN];
        fall_time = [fall_time; NaN];
        duty_cycle = [duty_cycle; NaN];
        n = numel(pk_forces);
        block_tbl = table(repmat({subj}, n, 1),repmat(ii, n, 1),(1:n).',pk_forces,pk_times,ISI,rise_time,RFD,fall_time,duty_cycle,cycle_auc, ...
            'VariableNames', {'ID','Block','Cycle','PeakForce','PeakForceTimes','ISI','RiseTime','RFD','FallTime','DutyCycle', 'AUC'});
        subj_cycle_tbl = [subj_cycle_tbl; block_tbl];

        % Figure 1. Peaks/valleys quality control plotting
        figure(signal_fig);
        set(gca, "FontSize", 12);
        plot(force_signal, 'k'); 
        hold on;
        plot(cur_block_vals, force_signal(cur_block_vals), 'b*');
        plot(cur_block_pks, force_signal(cur_block_pks), 'm*');
        plot(block_sts(ii), force_signal(block_sts(ii)), 'c*');
        plot(block_ends(ii), force_signal(block_ends(ii)), 'r*');
        xlabel('Frame (2000 Hz)');
        ylabel('Force (N)');
        title([subj, ' Handgrip Force'])

        % Figure 2: Plot the linear and exponential fits
        figure(model_fig);
        subplot(5, 1, ii);
        hold on;
        box off;
        plot(time_sec(cur_block_vals(1):block_ends(ii)), force_signal(cur_block_vals(1):block_ends(ii)), 'k');
        if all(isnan(SubjFits.(block_name).ExpEnvelope))
            plot(nan,nan,'w'); % placeholder so legend spacing stays consistent
        else
            plot(SubjFits.(block_name).BlockTime, SubjFits.(block_name).ExpEnvelope, 'b-', 'LineWidth', 2);
        end
        plot(SubjFits.(block_name).BlockTime, SubjFits.(block_name).LinEnvelope, 'r-.',  'LineWidth', 2);
        title([subj ' Block ' num2str(ii) ' Fatigability Models'], 'Interpreter', 'none');
        xlabel('Time (s)'); 
        xlim([(time_sec(cur_block_vals(1))-0.5),((time_sec(block_ends(ii)))+0.5)]);
        ylabel('Force (N)');
        ylim([0,(max(force_signal)*1.05)])
        legend(['Block ' num2str(ii) ' Signal'], sprintf('Exp:  k=%.3f  R²=%.2f', k_exp, R2_exp), sprintf('Linear:  m=%.3f  R²=%.2f', lin_slope, R2_lin), 'Location', 'eastoutside');

        % Add global fit to last subplot once last block is reached
        if ii == length(block_sts)
            subplot(5, 1, 5);
            plot(task_time_combined, task_sig_combined, 'k'); 
            hold on;
            box off;
            plot(task_time_combined, exp_model_fn(global_exp_model, task_time_combined), 'b-', 'LineWidth', 2);
            plot(task_time_combined, predict(global_lin_model, task_time_combined), 'r-.', 'LineWidth', 2);
            title([subj ' Global Fatigability Models'], 'Interpreter', 'none');
            xlabel('Cumulative Task Time (s)'); 
            xlim([(task_time_combined(1)-0.5),(task_time_combined(end)+0.5)]);
            ylabel('Force (N)');
            ylim([0,(max(force_signal)*1.05)])
            legend('Combined Signal',sprintf('Exp: k=%.2f R^2=%.2f', global_exp_k, global_R2_exp), sprintf('Linear: m=%.2f R^2=%.2f', global_lin_slope, global_R2_lin), 'Location', 'eastoutside');
        end

        % Figure 3. Plot the expected and actual AUC as shaded areas
        figure(AUC_fig);
        subplot(5, 1, ii);
        hold on;
        box off;
        %plot(t_extended, expected_curve, 'b-', 'LineWidth', 0.5);
        p1 = area(t_extended, expected_curve, 'FaceColor', 'b','FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot(t_actual, f_actual, 'r-', 'LineWidth', 0.5);
        p2 = area(t_actual, f_actual, 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        p3 = plot(nan,nan, 'w'); % dummy plot for legend
        legend([p1,p2,p3],... 
        sprintf('Expected Impulse: %.1f N·s', expected_auc), ...
        sprintf('Actual Impulse: %.1f N·s', actual_auc), ...
        sprintf('Fatigability Index: %.1f%%', Metrics.(block_name).FatigIdx_AUC_Expected * 100), ...
        'Location', 'eastoutside');
        title([subj ' Block ' num2str(ii) ' Expected vs Actual Impulse'], 'Interpreter', 'none');
        xlabel('Time (s)'); 
        xlim([(time_sec(cur_block_vals(1))-0.5),((time_sec(block_ends(ii)))+0.5)]);
        xlim([(t_actual(1) - 0.5), (t_actual(end) + 0.5)]);
        ylabel('Force (N)');
        ylim([0,(max(force_signal)*1.05)]);

        % Add global AUC to last subplot once last block is reached
        if ii == length(block_sts)
            subplot(5, 1, 5);
            hold on;
            box off;
            %plot(t_global_extended, global_expected_curve, 'b-', 'LineWidth', 0.5);
            p4 = area(t_global_extended, global_expected_curve, 'FaceColor', 'b','FaceAlpha', 0.3, 'EdgeColor', 'none');
            plot(task_time_combined, task_sig_combined, 'r-', 'LineWidth', 0.5);
            p5 = area(task_time_combined, task_sig_combined, 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            p6 = plot(nan,nan, 'w'); % dummy plot for legend
            legend([p4,p5,p6],... 
            sprintf('Expected Impulse: %.1f N·s', global_expected_auc), ...
            sprintf('Actual Impulse: %.1f N·s', global_actual_auc), ...
            sprintf('Fatigability Index: %.1f%%', Metrics.Global.FatigIdx_AUC_Expected * 100), ...
            'Location', 'eastoutside');
            title([subj ' Global Expected vs Actual Impulse'], 'Interpreter', 'none');
            xlabel('Cumulative Time (s)');
            xlim([(task_time_combined(1)-0.5),(task_time_combined(end)+0.5)]);
            ylabel('Force (N)');
            ylim([0,(max(force_signal)*1.05)])
        end
    end

    % Create Subject and Group Level Metrics Table
    % Define variables to pull from Metrics structure
    block_metric_names = { ...
    'NumCycles', 'PeakForce_Mean','PeakForce_CV', ...
    'ISI_Mean','ISI_CV', ...
    'RFD_Mean','RFD_CV', ...
    'RiseTime_Mean','RiseTime_CV', ...
    'FallTime_Mean','FallTime_CV', ...
    'DutyCycle_Mean','DutyCycle_CV', ...
    'LM_Slope','LM_Intercept','LM_R2','LM_AIC','LM_BIC', ...
    'Exp_k','Exp_A','Exp_C','Exp_Initial','Exp_R2','Exp_AIC','Exp_BIC', ...
    'FatigIdx_Force','FatigIdx_ISI','FatigIdx_RFD', ...
    'AUC_Actual','AUC_Expected','FatigIdx_AUC_Actual', 'FatigIdx_AUC_Expected'};
    trial_metric_names = { ...
    'NumCycles', 'NumCycles_Early', ...
    'PeakForce_Mean','PeakForce_CV', ...
    'ISI_Mean','ISI_CV', ...
    'RFD_Mean','RFD_CV', ...
    'RiseTime_Mean','RiseTime_CV', ...
    'FallTime_Mean','FallTime_CV', ...
    'DutyCycle_Mean','DutyCycle_CV', ...
    'LM_Slope_Mean','LM_Slope_CV', ...
    'Exp_k_Mean','Exp_k_CV', ...
    'FatigIdx_Force_Mean','FatigIdx_Force_CV', ...
    'FatigIdx_ISI_Mean','FatigIdx_ISI_CV', ...
    'FatigIdx_RFD_Mean','FatigIdx_RFD_CV', ...
    'AUC_Actual_Sum','AUC_Actual_Mean','AUC_Actual_CV', ...
    'AUC_Expected_Sum','AUC_Expected_Mean','AUC_Expected_CV', ...
    'FatigIdx_AUC_Actual_Mean','FatigIdx_AUC_Actual_CV', ...
    'FatigIdx_AUC_Expected_Mean','FatigIdx_AUC_Expected_CV', 'FatigIdx_AUC_Expected_Sum'...
    'TrialFatigIdx_Force','TrialFatigIdx_ISI','TrialFatigIdx_RFD', ...
    'TrialFatigIdx_LM_Slope','TrialFatigIdx_Exp_k', ...
    'TrialFatigIdx_AUC_Actual','TrialFatigIdx_AUC_Expected'};
    global_metric_names = { ...
    'LM_Slope','LM_Intercept','LM_R2','LM_AIC','LM_BIC', ...
    'Exp_k','Exp_A','Exp_C','Exp_Initial','Exp_R2','Exp_AIC','Exp_BIC', ...
    'FatigIdx_Force','FatigIdx_ISI','FatigIdx_RFD', ...
    'AUC_Actual','AUC_Expected','FatigIdx_AUC_Actual','FatigIdx_AUC_Expected'};
    
    % Build table dynamically
    vars = {};
    var_names = {};
        vars{end+1} = {subj};
    var_names{end+1} = 'ID';
    
    % Block level variables
    nBlocks = 4;
    for b1 = 1:nBlocks
        block_name = ['Block_' num2str(b1)];
        for m1 = 1:length(block_metric_names)
            metric = block_metric_names{m1};
            vars{end+1} = Metrics.(block_name).(metric);
            var_names{end+1} = ['B' num2str(b1) '_' metric];
        end
    end
    
    % Trial level variables
    for m1 = 1:length(trial_metric_names)
        metric = trial_metric_names{m1};
        vars{end+1} = Metrics.Trial.(metric);
        var_names{end+1} = ['Trial_' metric];
    end
    
    % Global level variables
    for m1 = 1:length(global_metric_names)
        metric = global_metric_names{m1};
        vars{end+1} = Metrics.Global.(metric);
        var_names{end+1} = ['Global_' metric];
    end
    
    % Convert to table and save
    subj_metrics_tbl = cell2table(vars, 'VariableNames', var_names);
    writetable(subj_metrics_tbl, fullfile(subj_dir, [subj, '_BlockMetrics.xlsx']));
    group_metrics_tbl = [group_metrics_tbl; subj_metrics_tbl];
    group_metrics_tbl = sortrows(group_metrics_tbl, 'ID');

    % Save Subject and Group Level Cycles
    writetable(subj_cycle_tbl, fullfile(subj_dir, [subj, '_CycleData.xlsx']));
    group_cycle_tbl = [group_cycle_tbl; subj_cycle_tbl];
    group_cycle_tbl = sortrows(group_cycle_tbl, 'ID');

    % Save figures as .fig and .png 
    savefig(signal_fig, fullfile(subj_dir, [subj, '_Signal_QC.fig']));
    exportgraphics(signal_fig, fullfile(subj_dir, [subj, '_Signal_QC.png']), 'Resolution', 600);
    savefig(model_fig, fullfile(subj_dir, [subj, '_Model_Viz.fig']));
    exportgraphics(model_fig, fullfile(subj_dir, [subj, '_Model_Viz.png']), 'Resolution', 600);
    savefig(AUC_fig, fullfile(subj_dir, [subj, '_AUC.fig']));
    exportgraphics(AUC_fig, fullfile(subj_dir, [subj, '_AUC.png']), 'Resolution', 600);
    
    % Save data to structure
    DataStruct.(subj_field).Metrics = Metrics;
    DataStruct.(subj_field).Cycles = subj_cycle_tbl;
    DataStruct.(subj_field).Fits = SubjFits;
    DataStruct = orderfields(DataStruct);

    % Save to data backup
    save(fullfile(output_dir, 'GroupDataBackup.mat'), 'DataStruct', 'group_metrics_tbl', 'group_cycle_tbl');

end

% Group Level Compilation
save(fullfile(output_dir, 'Handgrip_GroupResults.mat'), 'DataStruct');
writetable(group_metrics_tbl, fullfile(output_dir, 'Group_BlockMetrics.xlsx'));
writetable(group_cycle_tbl, fullfile(output_dir, 'Group_CycleData.xlsx'));
disp('Processing complete and group data exported!');