function run_phase10
% Reuse the validated Phase 8 fit implementation and diagnostics.
run_sleep_dandi(fullfile('data','fixtures','phase10','phase10_input.json'), ...
    'results/phase10','phase10-sleep');
end
