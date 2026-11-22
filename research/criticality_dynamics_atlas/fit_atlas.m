function fit_atlas
plan = jsondecode(fileread('atlas_work/plan.json'));
source = fullfile(pwd,'external','prox_crit_toolkit','src');
addpath(source);
pin = 'f26c782dca3b1c9ab18cc92c6449591c44e0d507';
for j=1:numel(plan)
    r = plan(j);
    name = r.name;
    rows = struct([]);
    if strcmp(name,'null_ar1')
        % p=beta/2 uses the toolkit integral directly and has no fmincon flag.
        file = fullfile(fileparts(r.receipt),'counts.mat');
        raw = load(file,'cohort');
        ids = fieldnames(raw.cohort);
        for i=1:numel(ids)
            item = raw.cohort.(ids{i});
            m = jsondecode(item.metadata);
            [d,~,a,s,~,~,exit_flag] = calc_db(item.counts,1,.04,2,struct('fit_method','YuleWalker'));
            f = struct('d2',d,'solver_flag',NaN,'exit_status',exit_flag, ...
                'min_ar_root',1/abs(a(1)),'innovation_variance',s^2, ...
                'ar_coefficients',a,'duration_s',sum(cellfun(@numel,item.counts))*.04);
            rows = [rows,fit_row(m,f,[])]; %#ok<AGROW>
        end
    else
        folder = fullfile('atlas_work','full_fits',name);
        run_sleep_dandi(r.receipt,folder,'criticality-dynamics-atlas',r.spec);
        data = jsondecode(fileread(fullfile(folder,'matlab_fits.json')));
        assert(strcmp(data.toolkit_commit,pin) && strcmp(data.input.sha256,r.input_sha256));
        for i=1:numel(data.records)
            m = data.records(i);
            rows = [rows,fit_row(m.record,m.fit,m.long_reach)]; %#ok<AGROW>
        end
    end
    assert(numel(rows)==r.records);
    receipt = struct('plan',r,'matlab',version,'toolkit_commit',pin, ...
        'run_id',getenv('GITHUB_RUN_ID'),'commit',getenv('GITHUB_SHA'),'rows',rows);
    writelines(string(jsonencode(receipt)),fullfile('atlas_work',[name '.fits.json']));
end
end

function row = fit_row(m,f,alt)
row=struct('id',m.id,'fit',brief(f),'alternate',[]);
if ~isempty(alt)
    row.alternate=brief(alt);
end
end

function out=brief(f)
out=struct('d2',f.d2,'flag',f.solver_flag,'exit_status',f.exit_status, ...
    'min_root',f.min_ar_root,'innovation_variance',f.innovation_variance, ...
    'ar_coefficients',f.ar_coefficients,'duration',f.duration_s);
if isfield(f,'residual_acf')
    out.residual_acf_peak=max(abs(f.residual_acf));
    out.residual_excess=f.residual_excess_kurtosis;
    out.acf_rmse=f.acf_rmse;
end
end
