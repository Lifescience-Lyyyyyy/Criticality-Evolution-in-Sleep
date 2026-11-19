function run_zebra_finch
paths = jsondecode(fileread('config/paths.json'));
config = jsondecode(fileread('config/zebra_finch.json'));
manifest = jsondecode(fileread(paths.source_manifest));
src = fullfile(pwd, 'external', 'prox_crit_toolkit', 'src');
addpath(src);
fixture = load(fullfile(paths.fixture_dir, 'zebra_lman_population.mat'), 'cohort');
timer = tic;

names = {'primary', 'short_windows', 'bin20', 'long_reach', 'matched_units', 'matched_rate'};
primary = evaluate_cohort(fixture.cohort, names(1), config, false);
cleanup = trace_solver(src);
diagnostics = evaluate_cohort(fixture.cohort, names, config, true);
clear cleanup

result = struct('experiment', 'zebra-finch-lman-wake-sws', 'primary', primary, ...
    'diagnostics', diagnostics, 'settings', config, 'data_version', manifest.zebra_finch.commit, ...
    'toolkit_commit', manifest.estimator.commit, 'project_commit', getenv('GITHUB_SHA'), ...
    'run_id', getenv('GITHUB_RUN_ID'), 'matlab', version, 'elapsed_seconds', toc(timer), ...
    'timestamp_utc', char(datetime('now', 'TimeZone', 'UTC', ...
        'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''')));
writelines(string(jsonencode(result, 'PrettyPrint', true)), fullfile(paths.results_dir, 'zebra_finch.json'));

primary_replay = diagnostics(strcmp({diagnostics.variant}, 'primary'));
assert(numel(primary) >= 4 && all([primary.exit_status] == 0));
assert(all(isfinite([primary.d2_bits_per_second])) && all([primary.d2_bits_per_second] >= -1e-10));
assert(all([primary.innovation_sd] > 0) && all([primary.fitted_min_root] > 1));
assert(all([primary_replay.solver_flag] > 0));
assert(max(abs([primary.d2_bits_per_second] - [primary_replay.d2_bits_per_second])) < 1e-8);
end

function rows = evaluate_cohort(cohort, variants, config, trace)
options = struct('fit_method', 'YuleWalker');
row = 0;
for bird = fieldnames(cohort)'
    for state = {'WAKE', 'SWS'}
        inputs = cohort.(bird{1}).(state{1});
        for variant = variants
            order = config.primary_ar_order;
            bin = config.primary_bin_seconds;
            segments = cellfun(@double, inputs.primary, 'UniformOutput', false);
            switch variant{1}
                case 'short_windows'
                    parts = cellfun(@(x) split_window(x, round(10/bin)), segments, 'UniformOutput', false);
                    segments = vertcat(parts{:});
                case 'bin20'
                    bin = 0.02;
                    order = 10;
                    segments = cellfun(@double, inputs.bin20, 'UniformOutput', false);
                case 'long_reach'
                    order = 10;
                case {'matched_units', 'matched_rate'}
                    segments = cellfun(@double, inputs.(variant{1}), 'UniformOutput', false);
            end
            if trace
                setappdata(0, 'estimator_solver', struct('flag', NaN, 'output', struct()));
            end
            [estimate, ~, kernel, sigma, ~, ~, status] = calc_db(segments, order, bin, 2, options);
            solver = struct('flag', NaN, 'output', struct());
            if trace
                solver = getappdata(0, 'estimator_solver');
                rmappdata(0, 'estimator_solver');
            end
            pooled = [segments{:}];
            row = row + 1;
            rows(row) = struct('bird', bird{1}, 'state', state{1}, 'variant', variant{1}, ...
                'd2_bits_per_second', estimate, 'ar_order', order, 'bin_seconds', bin, ...
                'ar_kernel', kernel, 'innovation_sd', sigma, 'exit_status', status, ...
                'solver_flag', solver.flag, 'solver_output', solver.output, ...
                'fitted_min_root', min(abs(roots([-flip(kernel(:)); 1]))), ...
                'residual_lag1_correlation', residual_lag1(segments, pooled, kernel), ...
                'count_mean', mean(pooled), 'count_variance', var(pooled), ...
                'zero_bin_fraction', mean(pooled == 0), 'segments', numel(segments), ...
                'samples', numel(pooled));
            fprintf('%s %s %s: d2=%.12g bits/s, upstream status %d, SQP flag %g\n', ...
                bird{1}, state{1}, variant{1}, estimate, status, solver.flag);
        end
    end
end
end

function chunks = split_window(series, width)
lengths = repmat(width, 1, floor(numel(series) / width));
tail = mod(numel(series), width);
if tail > 0
    lengths(end + 1) = tail;
end
chunks = mat2cell(series, 1, lengths)';
end

function lag1 = residual_lag1(segments, pooled, kernel)
center = mean(pooled);
scale = std(pooled);
order = numel(kernel);
numerator = 0;
denominator = 0;
for row = segments'
    residual = filter([1, -kernel(:)'], 1, (row{1} - center) / scale);
    tail = residual((order + 1):end);
    numerator = numerator + sum(tail(1:end-1) .* tail(2:end));
    denominator = denominator + sum(tail .^ 2);
end
lag1 = numerator / denominator;
end
