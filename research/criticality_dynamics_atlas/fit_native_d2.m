function fit_native_d2
% Use the pinned objective and optimizer with sufficient statistics of the counts.
input = 'research/criticality_dynamics_atlas/native_d2_input.json';
r = jsondecode(fileread(input));
src = fullfile(pwd,'external','prox_crit_toolkit','src');
addpath(src);
copy = tempname;
mkdir(copy);
cleanup = onCleanup(@() rmdir(copy, 's'));
original = fileread(fullfile(src, 'calc_db.m'));
anchor = "if strcmp(fit_method, 'MaxLikelihood')";
assert(numel(strfind(original, anchor))==1);
code = "if isfield(dbopt,'gamma')" + newline + ...
    "    g = dbopt.gamma(1:order); g = g(:);" + newline + ...
    "    kern = toeplitz([1; g(1:end-1)])\g;" + newline + ...
    "    sigma = sqrt(1-g'*kern); H = nan;" + newline + ...
    "elseif strcmp(fit_method, 'MaxLikelihood')";
writelines(string(strrep(original,anchor,code)),fullfile(copy,'calc_db.m'));
solver = trace_solver(copy);

% Compare the new entry path with the unmodified time-series path in the same runner.
rng(124);
x = {filter(1,[1,-.6],randn(1,541)),filter(1,[1,-.6],randn(1,435))};
pooled = [x{:}];
centered = cellfun(@(v) (v-mean(pooled))/std(pooled),x,'UniformOutput',false);
g = zeros(1,10);
for k=1:numel(g)
    g(k) = sum(cellfun(@(v) dot(v(1:end-k),v(1+k:end)),centered))/sum(cellfun(@numel,x)-k);
end
for p=[1,5,10]
    [dx,~,ax,sx] = calc_db(x,p,.04,2,struct('fit_method','YuleWalker'));
    [dg,~,ag,sg] = calc_db(0,p,.04,2,struct('gamma',g));
    assert(abs(dx-dg)<1e-6 && max(abs(ax-ag))<1e-10 && abs(sx-sg)<1e-10);
end
fits = repmat(struct('id','','fit',[],'alternate',[]),numel(r),1);
for j=1:numel(r)
    q = r(j);
    fits(j) = struct('id',q.id,'fit',fit(q,q.order),'alternate',fit(q,q.alternate_order));
end
clear solver
result = struct('input_sha256',getenv('ATLAS_INPUT_SHA256'),'matlab',version, ...
    'toolkit_commit','f26c782dca3b1c9ab18cc92c6449591c44e0d507', ...
    'run_id',getenv('GITHUB_RUN_ID'),'commit',getenv('GITHUB_SHA'),'rows',fits);
writelines(string(jsonencode(result)),fullfile('atlas_work','native-d2-fits.json'));
end

function y = fit(x,p)
setappdata(0,'estimator_solver',struct('flag',NaN));
[d2,~,a,s,~,~,status] = calc_db(0,p,x.width_ms/1000,2,struct('gamma',x.gamma));
flag = getappdata(0,'estimator_solver');
rmappdata(0,'estimator_solver');
y = struct('d2',d2,'flag',flag.flag,'exit_status',status,'min_root',min(abs(roots([-flip(a(:));1]))), ...
    'ar_coefficients',a(:)','innovation_variance',s^2,'duration',x.bins*x.width_ms/1000);
end
