function run_atlas_validation
% Isolated Phase 11 spot checks using the pinned sleep fit and diagnostics.
scenario = getenv('ATLAS_SCENARIO');
assert(ismember(scenario, {'v1','qc'}));
receipt = fullfile('data','fixtures','phase11',scenario,'validation_input.json');
run_sleep_dandi(receipt, fullfile('results','validation','runner'), 'atlas-targeted-validation');
end
