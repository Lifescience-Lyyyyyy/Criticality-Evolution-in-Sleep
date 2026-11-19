function run_pnas_mouse
paths = jsondecode(fileread('config/paths.json'));
config = jsondecode(fileread('config/pnas.json'));
manifest = jsondecode(fileread(paths.source_manifest));
src = fullfile(pwd, 'external', 'prox_crit_toolkit', 'src');
addpath(src);
data = load(fullfile(paths.fixture_dir, 'mouse_population_40ms.mat'), 'cohort');
expected = string({manifest.mouse.files.session});
assert(isequal(sort(string(fieldnames(data.cohort))'), sort(expected)), ...
    'The count fixture must contain exactly the 13 published sessions.');
timer = tic;

primary = estimate_states(data.cohort, expected, config, false);
original_y10 = data.cohort.YMV10.NREM;
zero_counts = num2cell(double(data.cohort.YMV10.NREM_zero_duration_counts(:)));
with_zero = [cellfun(@double, original_y10, 'UniformOutput', false); zero_counts];
options = struct('fit_method', config.fit_method);
[y10_with_zero, ~, ~, ~, ~, ~, y10_status] = calc_db(with_zero, ...
    config.ar_order, config.bin_seconds, config.beta, options);
assert(y10_status == 0 && isfinite(y10_with_zero));
cleanup = trace_solver(src);
diagnostics = estimate_states(data.cohort, expected, config, true);
clear cleanup

delta = [primary.d2_bits_per_second];
delta = delta(2:2:end) - delta(1:2:end);
result = struct('experiment', 'pnas-mouse-wake-nrem', 'primary', primary, ...
    'optimizer_replay', diagnostics, 'paired_nrem_minus_wake', delta, ...
    'ymv10_nrem_including_zero_duration', y10_with_zero, ...
    'mean_paired_difference', mean(delta), 'median_paired_difference', median(delta), ...
    'settings', config, 'data_version', manifest.mouse.doi, ...
    'toolkit_commit', manifest.estimator.commit, 'project_commit', getenv('GITHUB_SHA'), ...
    'run_id', getenv('GITHUB_RUN_ID'), 'matlab', version, 'elapsed_seconds', toc(timer), ...
    'timestamp_utc', char(datetime('now', 'TimeZone', 'UTC', ...
        'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''')));
result_path = fullfile(paths.results_dir, 'mouse.json');
writelines(string(jsonencode(result, 'PrettyPrint', true)), result_path);
disp(jsonencode(result));
assert(all([primary.exit_status] == 0) && all([diagnostics.solver_flag] > 0), ...
    'At least one fit or optimization failed. See the per-session results.');
assert(all(isfinite([primary.d2_bits_per_second])) && all([primary.d2_bits_per_second] >= -1e-10));
assert(all([primary.innovation_sd] > 0) && all([primary.fitted_min_root] > 1));
assert(max(abs([primary.d2_bits_per_second] - [diagnostics.d2_bits_per_second])) < 1e-8);
assert(mean(delta) > 0, 'The published mean Wake < NREM direction was not reproduced.');
end

function rows = estimate_states(cohort, sessions, config, trace)
options = struct('fit_method', config.fit_method);
row = 0;
for session = sessions
    for index = 1:numel(config.states)
        state = config.states{index};
        input = cohort.(session).(state);
        segments = cellfun(@double, input, 'UniformOutput', false);
        if trace
            setappdata(0, 'estimator_solver', struct('flag', NaN, 'output', struct()));
        end
        [estimate, ~, kernel, sigma, ~, critical_kernel, status] = calc_db(segments, ...
            config.ar_order, config.bin_seconds, config.beta, options);
        solver = struct('flag', NaN, 'output', struct());
        if trace
            solver = getappdata(0, 'estimator_solver');
            rmappdata(0, 'estimator_solver');
        end
        row = row + 1;
        rows(row) = struct('session', session, 'state', state, ...
            'd2_bits_per_second', estimate, 'ar_kernel', kernel, ...
            'innovation_sd', sigma, 'critical_kernel', critical_kernel, ...
            'exit_status', status, 'solver_flag', solver.flag, 'solver_output', solver.output, ...
            'fitted_min_root', min(abs(roots([-flip(kernel(:)); 1]))), ...
            'segments', numel(segments), 'samples', sum(cellfun(@numel, segments)));
        fprintf('%s %s d2 = %.12g bits/s, upstream status %d, SQP flag %g\n', ...
            session, state, estimate, status, solver.flag);
    end
end
end
