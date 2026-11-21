function run_atlas_validation
% Isolated Phase 11 spot checks using the pinned sleep fit and diagnostics.
scenario = getenv('ATLAS_SCENARIO');
assert(ismember(scenario, {'v1','qc','shorter','pairwise','extension'}));
batch = 'phase11';
if ismember(scenario, {'pairwise','extension'})
    batch = 'phase11_1';
end
receipt = fullfile('data','fixtures',batch,scenario,'validation_input.json');
run_sleep_dandi(receipt, fullfile('results','validation','runner'), 'atlas-targeted-validation');
end
