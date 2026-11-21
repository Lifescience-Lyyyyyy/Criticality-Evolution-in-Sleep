function fit_matlab
% Fit checked synthetic spike-count fixtures with the unchanged sleep estimator.
plan = jsondecode(fileread('benchmark_work/inputs/fit_plan.json'));
assert(numel(plan) == 5 || (numel(plan) == 1 && strcmp(plan.name,'empirical')), ...
    'Expected the frozen five model specifications or one empirical specification.');
for i = 1:numel(plan)
    work = plan(i);
    run_sleep_dandi(work.receipt,fullfile('benchmark_work','fits',work.name), ...
        ['sleep-model-benchmark-' work.name],work.analysis_spec);
end
end
