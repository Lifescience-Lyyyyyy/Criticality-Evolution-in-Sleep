function run_synthetic
toolkit = fullfile(pwd, 'external', 'prox_crit_toolkit');
src = fullfile(toolkit, 'src');
demo = fullfile(toolkit, 'demo', 'smoke_demo.m');
addpath(src);
seed = 2026;
timer = tic;

first = official_example(demo, seed);
repeated = official_example(demo, seed);

% A temporary diagnostic replay only requests two extra outputs from the
% same fmincon call. The primary values come from the original source.
original = fileread(fullfile(src, 'calc_db.m'));
call = '[x, db] = fmincon(objective, kernc0, [], [], [], [], [], [], nonlcon, options);';
logged = strrep(call, '[x, db]', '[x, db, solver_flag, solver_output]');
instrumented = strrep(original, call, sprintf('%s\n    %s', logged, ...
    'setappdata(0, ''synthetic_solver'', struct(''flag'', solver_flag, ''output'', solver_output));'));
assert(~strcmp(original, instrumented), 'Pinned optimizer call was not found.');
diagnostic_dir = tempname;
mkdir(diagnostic_dir);
writelines(string(instrumented), fullfile(diagnostic_dir, 'calc_db.m'));
addpath(diagnostic_dir, '-begin');
clear calc_db
replay = official_example(demo, seed);
solver = getappdata(0, 'synthetic_solver');
rmappdata(0, 'synthetic_solver');
rmpath(diagnostic_dir);
clear calc_db

% For stationary Gaussian AR(1), the distance to a unit root with optimal
% innovation variance is 0.5*log2(2/(1+phi)) bits per sample.
rng(seed, 'twister');
noise = randn(1, 12000);
rho = [0.50, 0.95];
options = struct('fit_method', 'YuleWalker');
for i = 1:numel(rho)
    x = filter(1, [1, -rho(i)], noise);
    [estimate, ~, phi, sigma, ~, ~, status] = calc_db(x(2001:end), 1, 1, 2, options);
    ar1(i) = struct('generating_phi', rho(i), 'fitted_phi', phi, ...
        'd2_bits_per_sample', estimate, 'fitted_analytic_d2', 0.5*log2(2/(1+phi)), ...
        'true_analytic_d2', 0.5*log2(2/(1+rho(i))), ...
        'innovation_sd', sigma, 'exit_status', status);
end

result = struct('experiment', 'pnas-d2-synthetic', 'official', first, ...
    'official_repeated', repeated, 'diagnostic_replay', replay, 'optimizer', solver, ...
    'ar1', ar1, 'ar1_samples', 10000, 'ar1_burnin', 2000, 'seed', seed, ...
    'elapsed_seconds', toc(timer), 'matlab', version, ...
    'project_commit', getenv('GITHUB_SHA'), 'run_id', getenv('GITHUB_RUN_ID'), ...
    'timestamp_utc', char(datetime('now', 'TimeZone', 'UTC', ...
        'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''')));
writelines(string(jsonencode(result, 'PrettyPrint', true)), 'results/synthetic.json');
disp(jsonencode(result));
estimates = [ar1.d2_bits_per_sample];
assert(first.exit_status == 0 && repeated.exit_status == 0 && replay.exit_status == 0);
assert(isfinite(first.d2_bits_per_second) && first.d2_bits_per_second >= -1e-10);
assert(max(abs(first.ar_kernel - [0.8; 0.2; -0.1])) < 0.05 && first.innovation_sd > 0);
assert(abs(first.d2_bits_per_second - repeated.d2_bits_per_second) < 1e-12);
assert(abs(first.d2_bits_per_second - replay.d2_bits_per_second) < 1e-12 && solver.flag > 0);
assert(all([ar1.exit_status] == 0) && all(isfinite(estimates)) && all(estimates >= -1e-10));
assert(estimates(2) < estimates(1), 'Near-critical AR(1) should have lower d2.');
assert(max(abs(estimates - [ar1.fitted_analytic_d2])) < 1e-8);
assert(max(abs(estimates - [ar1.true_analytic_d2])) < 0.015);
end

function result = official_example(demo, seed)
rng(seed, 'twister');
run(demo);
result = struct('d2_bits_per_second', db, 'exit_status', exit_status, ...
    'ar_kernel', kern(:), 'innovation_sd', sigma, 'critical_kernel', kernc(:), ...
    'ar_order', order, 'beta', b, 'samples', T, 'bin_seconds', deltaT, ...
    'fit_method', dbopt.fit_method);
end
