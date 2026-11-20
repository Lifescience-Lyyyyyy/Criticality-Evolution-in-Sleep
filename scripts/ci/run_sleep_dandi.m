function run_sleep_dandi
% Run the pinned estimator on audited population counts without joining bouts.
config = jsondecode(fileread('config/sleep_dandi.json'));
manifest = jsondecode(fileread('manifests/sources.json'));
receipt = jsondecode(fileread('data/fixtures/sleep_dandi_input.json'));
src = fullfile(pwd,'external','prox_crit_toolkit','src');
addpath(src);
input = load(fullfile('data','fixtures',receipt.fixture),'cohort');
ids = fieldnames(input.cohort);
assert(numel(ids) == receipt.records);
timer = tic;
controls_mode = isfield(receipt,'mode') && strcmp(receipt.mode,'controls');
pilot = numel(receipt.asset_ids) == 1 && receipt.swr_windows == 0 && ~controls_mode;
reference = strcmp(arrayfun(@(i) variant(input.cohort.(ids{i})),1:numel(ids),'UniformOutput',false),'reference');
original = nan(size(ids));
if pilot
    for i = find(reference(:))'
        item = input.cohort.(ids{i});
        [original(i),~,~,~,~,~,status] = calc_db(counts(item),config.ar_order, ...
            config.bin_seconds,config.beta,struct('fit_method',config.fit_method));
        assert(status == 0 && isfinite(original(i)));
    end
end

cleanup = trace_solver(src);
if controls_mode
    rows = controls(input.cohort,ids,config);
else
    rows = struct([]);
    for i = 1:numel(ids)
        item = input.cohort.(ids{i});
        m = jsondecode(item.metadata);
        c = counts(item);
        fit = estimate(c,config.ar_order,config);
        if pilot && reference(i)
            assert(is_valid(fit) && abs(fit.d2-original(i)) < 1e-8, ...
                'Pilot fit or original-vs-diagnostic replay failed.');
        end
        row = struct('record',m,'fit',fit,'long_reach',[]);
        if reference(i)
            row.long_reach = estimate(c,config.long_ar_order,config);
        end
        rows = [rows,row]; %#ok<AGROW>
        if reference(i) || mod(i,50)==0
            fprintf('%s %s %s %s, %d segments, flag %g, d2 %.8g\n', ...
                m.session,m.region,m.context,m.state,numel(c),fit.solver_flag,fit.d2);
        end
    end
end
clear cleanup
result = struct('schema_version',1,'experiment','sleep-dandi','input',receipt, ...
    'settings',config,'toolkit_commit',manifest.estimator.commit,'project_commit',getenv('GITHUB_SHA'), ...
    'matlab',version,'run_id',getenv('GITHUB_RUN_ID'),'elapsed_seconds',toc(timer), ...
    'timestamp_utc',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')), ...
    'pilot_original_replay_max_error',NaN, ...
    'records',rows);
if pilot
    primary = [rows(reference).fit];
    result.pilot_original_replay_max_error = max(abs([primary.d2]-original(reference)'));
else
    result.pilot_original_replay_max_error = NaN;
end
out = fullfile('results','sleep_dandi');
if ~isfolder(out), mkdir(out); end
writelines(string(jsonencode(result,'PrettyPrint',true)),fullfile(out,'matlab_fits.json'));
if any(reference)
    primary = [rows(reference).fit];
    assert(any(arrayfun(@is_valid,primary)), 'No valid primary estimates. See diagnostics.');
else
    assert(any(arrayfun(@is_valid,[rows.fit])), 'No valid controls or window estimates.');
end
end

function name = variant(row)
m = jsondecode(row.metadata);
name = m.variant;
end

function out = counts(row)
out = cellfun(@(x)reshape(double(x),1,[]),row.counts(:),'UniformOutput',false);
end

function pass = is_valid(f)
pass = isfinite(f.d2) && f.d2 >= -1e-10 && f.exit_status == 0 && f.solver_flag > 0 ...
    && f.min_ar_root > 1 && f.innovation_variance > 0;
end

function out = estimate(segments,order,config)
setappdata(0,'estimator_solver',struct('flag',NaN,'output',struct()));
[d2,~,kernel,sigma,~,~,status] = calc_db(segments,order,config.bin_seconds, ...
    config.beta,struct('fit_method',config.fit_method));
solver = getappdata(0,'estimator_solver');
rmappdata(0,'estimator_solver');
sizes = cellfun(@numel,segments)';
population = [segments{:}];
mu = mean(population);
v = var(population,1);
a = kernel(:)';
out = struct('d2',d2,'ar_order',order,'ar_coefficients',a,'innovation_variance',sigma^2, ...
    'exit_status',status,'solver_flag',solver.flag,'solver_output',solver.output, ...
    'min_ar_root',min(abs(roots([-flip(a),1]))),'segment_bins',sizes, ...
    'duration_s',sum(sizes)*config.bin_seconds,'rate_hz',mu/config.bin_seconds, ...
    'count_variance',v,'fano',v/mu,'zero_fraction',mean(population==0), ...
    'observed_acf',segment_acf(segments,20),'residual_acf',[], ...
    'model_acf',[],'acf_rmse',NaN,'residual_excess_kurtosis',NaN, ...
    'residual_abs_tail_over_3sd',NaN,'silent_run_counts',0,'max_silence_s',0,'silence_p95_s',0);
silent = cellfun(@zero_runs,segments,'UniformOutput',false);
silent = [silent{:}];
if ~isempty(silent)
    out.silent_run_counts = numel(silent);
    out.max_silence_s = max(silent)*config.bin_seconds;
    out.silence_p95_s = quantile(silent,.95)*config.bin_seconds;
end
if numel(a) ~= order || v == 0
    return
end
residuals = cellfun(@(x)filter([1,-a],1,(x-mu)/sqrt(v)),segments,'UniformOutput',false);
residuals = cellfun(@(x)x(order+1:end),residuals,'UniformOutput',false);
error_values = [residuals{:}];
out.residual_acf = segment_acf(residuals,20);
dev = error_values-mean(error_values);
out.residual_excess_kurtosis = mean(dev.^4)/mean(dev.^2)^2-3;
out.residual_abs_tail_over_3sd = mean(abs(dev)>3*sqrt(mean(dev.^2)));
if out.min_ar_root > 1
    transition = [a;eye(order-1),zeros(order-1,1)];
    shock = zeros(order);
    shock(1) = out.innovation_variance;
    stationary = reshape((eye(order^2)-kron(transition,transition))\shock(:),order,order);
    model = stationary;
    out.model_acf = zeros(1,20);
    for lag = 1:20
        model = transition*model;
        out.model_acf(lag) = model(1,1)/stationary(1,1);
    end
    out.acf_rmse = sqrt(mean((out.observed_acf-out.model_acf).^2,'omitnan'));
end
end

function out = zero_runs(x)
edges = find(diff([false,x==0,false]));
out = edges(2:2:end)-edges(1:2:end);
end

function acf = segment_acf(series,lags)
values = [series{:}];
mu = mean(values);
v = mean((values-mu).^2);
num = zeros(1,lags);
den = zeros(1,lags);
for i = 1:numel(series)
    x = series{i}-mu;
    for j = 1:min(numel(x)-1,lags)
        num(j) = num(j) + dot(x(1:end-j),x(j+1:end));
        den(j) = den(j) + numel(x)-j;
    end
end
acf = num./den/v;
end

function rows = controls(cohort,ids,config)
% Equal-exposure sensitivity cohorts, prespecified from audit metadata.
meta = cell(numel(ids),1);
spike_matrices = cell(size(meta));
means = zeros(size(meta));
groups = strings(size(meta));
for i = 1:numel(ids)
    item = cohort.(ids{i});
    meta{i} = jsondecode(item.metadata);
    r = meta{i};
    spike_matrices{i} = item.counts(:);
    c = cellfun(@(x)sum(double(x),1),spike_matrices{i},'UniformOutput',false);
    means(i) = mean([c{:}]);
    groups(i) = string(r.asset_id) + "|" + string(r.context) + "|" + string(r.region);
end
[~,~,index] = unique(groups);
rows = struct([]);
for i = 1:numel(ids)
    m = meta{i};
    cells = spike_matrices{i};
    source = cellfun(@(x)sum(double(x),1),cells,'UniformOutput',false);
    target = min(means(index == index(i)));
    if target <= 0, continue, end
    baseline = light_fit(source,config);
    r = m;
    r.variant = 'matched_controls_reference';
    r.seed = NaN;
    r.iteration = 0;
    rows = [rows,struct('record',r,'fit',baseline,'long_reach',[])]; %#ok<AGROW>
    for k = 1:config.repetitions
        seed = config.seed + 1000*index(i)+k;
        n = size(cells{1},1);
        stream = RandStream('mt19937ar','Seed',seed);
        selected = randperm(stream,n,m.target_units);
        fitted = baseline;
        if n ~= m.target_units
            sampled = cellfun(@(x)sum(double(x(selected,:)),1),cells,'UniformOutput',false);
            fitted = light_fit(sampled,config);
        end
        r = m;
        r.variant = 'unit_count';
        r.unit_count = m.target_units;
        r.iteration = k;
        r.seed = seed;
        rows = [rows,struct('record',r,'fit',fitted,'long_reach',[])]; %#ok<AGROW>

        r = m;
        r.variant = 'rate_thinning';
        r.unit_count = n;
        r.iteration = k;
        r.seed = seed + 10^6 + 100000*strcmp(m.state,'Wake');
        r.target_rate_hz = target/config.bin_seconds;
        fitted = baseline;
        if target < means(i)
            stream = RandStream('mt19937ar','Seed',r.seed);
            fitted = light_fit(thin(source,target/means(i),stream),config);
        end
        rows = [rows,struct('record',r,'fit',fitted,'long_reach',[])]; %#ok<AGROW>
    end
    fprintf('%s %s %s %s: %d repetitions and matched units %d\n', ...
            m.session,m.context,m.region,m.state,config.repetitions,m.target_units);
end
end

function f = light_fit(c,config)
setappdata(0,'estimator_solver',struct('flag',NaN,'output',struct()));
[d2,~,kernel,sigma,~,~,code] = calc_db(c,config.ar_order,config.bin_seconds,config.beta,struct('fit_method',config.fit_method));
solver = getappdata(0,'estimator_solver');
rmappdata(0,'estimator_solver');
x = [c{:}];
f = struct('d2',d2,'exit_status',code,'solver_flag',solver.flag, ...
    'min_ar_root',min(abs(roots([-flip(kernel(:)); 1]))),'innovation_variance',sigma^2, ...
    'rate_hz',mean(x)/config.bin_seconds,'count_variance',var(x,1), ...
    'observed_acf',segment_acf(c,5),'duration_s',numel(x)*config.bin_seconds);
end

function result = thin(series,p,stream)
% Invert the binomial CDF in base MATLAB. Preserve every segment boundary.
x = [series{:}];
u = rand(stream,size(x));
draws = zeros(size(x));
for k = unique(x)
    if k==0, continue, end
    idx = x==k;
    counts = 0:k-1;
    cdf = betainc(1-p,k-counts,counts+1);
    draws(idx) = sum(u(idx)' > cdf,2)';
end
result = mat2cell(draws,1,cellfun(@numel,series)');
result = result(:);
end
