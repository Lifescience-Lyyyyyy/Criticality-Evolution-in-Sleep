function run_zebra_finch
paths = jsondecode(fileread('config/paths.json'));
config = jsondecode(fileread('config/zebra_finch.json'));
manifest = jsondecode(fileread(paths.source_manifest));
src = fullfile(pwd, 'external', 'prox_crit_toolkit', 'src');
addpath(src);
input = load(fullfile(paths.fixture_dir, 'zebra_lman_population.mat'), 'cohort');
timer = tic;
states = {'ACTIVE', 'NON_ACTIVE', 'SWS'};
birds = fieldnames(input.cohort)';
bin = config.primary_bin_seconds;
n = config.reanalysis.repetitions;
order = config.primary_ar_order;

primary = struct([]);
for bird = birds
    for state = states
        row = input.cohort.(bird{1}).(state{1});
        series = pool(row.segments, true(size(row.broad_unit_ids)));
        fit = estimate(series, order, bin, false);
        primary = [primary, record(fit, bird{1}, state{1}, 'broad', numel(row.broad_unit_ids))];
    end
end

cleanup = trace_solver(src);
diagnostics = struct([]);
controls = struct([]);
windows = struct([]);
for b = 1:numel(birds)
    bird = birds{b};
    cohort = input.cohort.(bird);
    state_means = cellfun(@(s) state_mean(cohort.(s)), states);
    target = min(state_means);
    for j = 1:numel(states)
        state = states{j};
        row = cohort.(state);
        full = true(size(row.broad_unit_ids));
        series = pool(row.segments, full);
        fprintf('%s %s: estimating %d repetitions for each control\n', bird, state, n);
        common = estimate(series, order, bin, true);
        credible = model_checks(series, common, config.seed + b*100 + j);
        diagnostics = [diagnostics, record(credible, bird, state, 'broad', nnz(full))];
        long_fit = model_checks(series, estimate(series, 10, bin, true), config.seed + b*100 + j);
        diagnostics = [diagnostics, record(long_fit, bird, state, 'reach_0p40', nnz(full))];
        if ~isempty(row.strict)
            n_strict = size(row.strict{1}, 1);
            selected = pool(row.strict, true(1, n_strict));
            fit = model_checks(selected, estimate(selected, order, bin, true), config.seed + b*100 + j);
            diagnostics = [diagnostics, record(fit, bird, state, 'strict_unpaired', n_strict)];
        end
        if j < 3
            for kind = {'unmasked', 'literal_mask'}
                variant = pool(row.(kind{1}), full);
                fit = model_checks(variant, estimate(variant, order, bin, true), config.seed + b*100 + j);
                diagnostics = [diagnostics, record(fit, bird, state, kind{1}, nnz(full))];
            end
        end
        for k = 1:n
            seed = config.seed + b*10000 + k;
            selected = logical(row.unit_masks(k, :));
            counts = pool(row.segments, selected);
            fit = estimate(counts, order, bin, true);
            controls = [controls, control_row(fit, counts, bird, state, 'unit_subsampling', ...
                k, double(row.unit_seeds(k)), nnz(selected), NaN)];
            stream = RandStream('mt19937ar', 'Seed', seed + j*1000);
            counts = thin_counts(series, target / state_means(j), stream);
            fit = estimate(counts, order, bin, true);
            controls = [controls, control_row(fit, counts, bird, state, 'rate_thinning', ...
                k, seed + j*1000, nnz(full), target / bin)];
        end
        for k = 1:numel(row.windows)
            counts = pool(row.windows(k), full);
            fit = estimate(counts, order, bin, true);
            window = control_row(fit, counts, bird, state, 'window_60s', k, NaN, nnz(full), NaN);
            window.start_s = row.window_starts(k);
            windows = [windows, window];
        end
    end
end
clear cleanup

result = struct('schema_version', 2, 'experiment', 'zebra-lman-reanalysis', ...
    'strict_primary_status', 'BLOCKED_INSUFFICIENT_GOOD_UNITS', ...
    'broad_state_estimates', primary, 'diagnostics', diagnostics, 'resampling', controls, ...
    'windows', windows, 'settings', config, 'data_version', manifest.zebra_finch.commit, ...
    'toolkit_commit', manifest.estimator.commit, 'project_commit', getenv('GITHUB_SHA'), ...
    'run_id', getenv('GITHUB_RUN_ID'), 'matlab', version, 'elapsed_seconds', toc(timer), ...
    'timestamp_utc', char(datetime('now', 'TimeZone', 'UTC', ...
        'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''')));
writelines(string(jsonencode(result, 'PrettyPrint', true)), fullfile(paths.results_dir, 'zebra_finch_reanalysis.json'));

replay = diagnostics(strcmp({diagnostics.variant}, 'broad'));
assert(numel(primary) == 6 && numel(controls) == 6*2*n);
assert(all([primary.exit_status] == 0) && all(isfinite([primary.d2])) && all([primary.d2] >= -1e-10));
assert(all([replay.solver_flag] > 0) && max(abs([replay.d2] - [primary.d2])) < 1e-8);
assert(all([primary.innovation_variance] > 0) && all([primary.fitted_min_root] > 1));
fprintf('Validated %d broad state fits, %d control fits, and %d timed windows\n', ...
    numel(primary), numel(controls), numel(windows));
end

function series = pool(cells, selected)
series = cellfun(@(x) sum(double(x(selected, :)), 1), cells(:), 'UniformOutput', false);
end

function rate = state_mean(row)
series = pool(row.segments, true(size(row.broad_unit_ids)));
rate = mean([series{:}]);
end

function out = estimate(series, order, bin, debug)
if debug
    setappdata(0, 'estimator_solver', struct('flag', NaN, 'output', struct()));
end
[d2, ~, kernel, sigma, ~, ~, status] = calc_db(series, order, bin, 2, struct('fit_method', 'YuleWalker'));
solver = struct('flag', NaN, 'output', struct());
if debug
    solver = getappdata(0, 'estimator_solver');
    rmappdata(0, 'estimator_solver');
end
out = struct('d2', d2, 'ar_order', order, 'ar_coefficients', kernel, ...
    'innovation_variance', sigma^2, 'exit_status', status, 'solver_flag', solver.flag, ...
    'solver_output', solver.output, 'fitted_min_root', min(abs(roots([-flip(kernel(:)); 1]))), ...
    'segments', numel(series), 'seconds', sum(cellfun(@numel, series))*bin);
end

function out = record(fit, bird, state, variant, units)
out = fit;
out.bird = bird;
out.state = state;
out.variant = variant;
out.units = units;
end

function out = model_checks(series, fit, seed)
out = fit;
pool = [series{:}];
mean_count = mean(pool);
out.firing_rate_hz = mean_count / 0.04;
out.count_variance = var(pool, 1);
out.fano = out.count_variance / mean_count;
out.zero_fraction = mean(pool == 0);
normalized = cellfun(@(x) (x - mean_count) / std(pool), series, 'UniformOutput', false);
a = fit.ar_coefficients(:)';
resid = cellfun(@(x) filter([1, -a], 1, x), normalized, 'UniformOutput', false);
resid = cellfun(@(x) x((numel(a)+1):end), resid, 'UniformOutput', false);
r = [resid{:}];
out.residual_acf = autocorr_segments(resid, 20);
out.residual_excess_kurtosis = mean((r-mean(r)).^4) / var(r, 1)^2 - 3;
out.residual_quantiles_innovation_sd = quantile(r, [0.001, 0.01, 0.05, 0.50, 0.95, 0.99, 0.999]) / sqrt(fit.innovation_variance);
out.residual_tail_above_3sd = mean(abs(r-mean(r)) > 3*std(r));
out.observed_acf = autocorr_segments(series, 20);
stream = RandStream('mt19937ar', 'Seed', seed);
noise = sqrt(fit.innovation_variance)*randn(stream, 1, numel(pool)+2000);
sample = filter(1, [1, -a], noise);
simulation = mat2cell(sample(2001:end), 1, cellfun(@numel, series)');
out.gaussian_ar_generated_acf = autocorr_segments(simulation, 20);
out.acf_rmse = sqrt(mean((out.observed_acf - out.gaussian_ar_generated_acf).^2));
out.zero_run_lengths_bins = zeros(1, 0);
for row = series(:)'
    starts_stops = find(diff([false, row{1} == 0, false]));
    out.zero_run_lengths_bins = [out.zero_run_lengths_bins, starts_stops(2:2:end)-starts_stops(1:2:end)];
end
end

function acf = autocorr_segments(series, maxlag)
center = mean([series{:}]);
v = cellfun(@(x) x-center, series, 'UniformOutput', false);
variance = mean([v{:}].^2);
acf = zeros(1, maxlag);
for lag = 1:maxlag
    sums = 0;
    samples = 0;
    for row = v(:)'
        x = row{1};
        sums = sums + dot(x(1:end-lag), x(1+lag:end));
        samples = samples + max(numel(x)-lag, 0);
    end
    acf(lag) = sums / samples / variance;
end
end

function out = control_row(fit, series, bird, state, kind, iteration, seed, units, target)
counts = [series{:}];
out = struct('bird', bird, 'state', state, 'variant', kind, 'iteration', iteration, ...
    'seed', seed, 'units', units, 'target_rate_hz', target, 'actual_rate_hz', mean(counts)/0.04, ...
    'count_variance', var(counts,1), 'fano', var(counts,1)/mean(counts), ...
    'zero_fraction', mean(counts==0), 'acf', autocorr_segments(series, 5), ...
    'd2', fit.d2, 'exit_status', fit.exit_status, 'solver_flag', fit.solver_flag);
end

function result = thin_counts(series, p, stream)
if p >= 1
    result = series;
    return
end
x = [series{:}];
y = zeros(size(x));
u = rand(stream, size(x));
for count = unique(x)
    if count == 0
        continue
    end
    idx = x == count;
    k = 0:count-1;
    cdf = betainc(1-p, count-k, k+1);
    y(idx) = sum(u(idx)' > cdf, 2)';
end
result = mat2cell(y, 1, cellfun(@numel, series)');
result = result(:);
end
