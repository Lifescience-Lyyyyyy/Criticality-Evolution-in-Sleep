function run_crcns_sleep
% Reuse the segment-aware estimator and diagnostics used in Phase 8.
run_sleep_dandi(fullfile('data','fixtures','phase9','crcns_input.json'), ...
                fullfile('results','sleep_crcns'),'sleep-crcns');
end
