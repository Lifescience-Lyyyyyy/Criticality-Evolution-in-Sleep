function cleanup = trace_solver(src)
% Observe optional fmincon outputs in a temporary copy of the pinned code.
original = fileread(fullfile(src, 'calc_db.m'));
call = '[x, db] = fmincon(objective, kernc0, [], [], [], [], [], [], nonlcon, options);';
logged = strrep(call, '[x, db]', '[x, db, solver_flag, solver_output]');
instrumented = strrep(original, call, sprintf('%s\n    %s', logged, ...
    'setappdata(0, ''estimator_solver'', struct(''flag'', solver_flag, ''output'', solver_output));'));
assert(~strcmp(original, instrumented), 'Pinned optimizer call was not found.');
diagnostic_dir = tempname;
mkdir(diagnostic_dir);
writelines(string(instrumented), fullfile(diagnostic_dir, 'calc_db.m'));
addpath(diagnostic_dir, '-begin');
clear calc_db
cleanup = onCleanup(@() restore_path(diagnostic_dir));
end

function restore_path(diagnostic_dir)
rmpath(diagnostic_dir);
clear calc_db
end
