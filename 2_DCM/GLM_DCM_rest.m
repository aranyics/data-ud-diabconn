
set(0,'DefaultFigureVisible','off');

% !! directories and includes
% diabetes and obese cohorts are analysed independently. Comment lines to change group!
%source_path  = '/path/to/diabetes/GIFT-ICA/input/obese/';     % input for obese group
source_path  = '/path/to/diabetes/GIFT-ICA/input/diabetes/';   % input for diabetes group
project_path = '/path/to/diabetes/DCM/largescale/';            % output dir
%regions_path = '/path/to/diabetes/GIFT-ICA/obese_network/';   % adjusted regions for obese group
regions_path = '/path/to/diabetes/GIFT-ICA/diabetes_network/'; % adjusted regions for diabetes group
DCM_dir = 'conn/';                                             % output subdir
addpath('./'); % add the directory where this script is found to paths

% !! script run - it is advised to run each batches (GLM, VOI, DCM) separately to check results
argGLM = 1; % run GLM batch
argVOI = 1; % run time-series extraction batch
argDCM = 1; % run DCM batch
motcorr = 0; % turned off motion regression, as it was done during pre-processing
resting_state = 1; % select for resing-state fMRI data

SessID = 'null';
ParID = 'null';
SubjID = 'null';
SubjName = 'null';

% !! fmri specifications
VolDim4 = 256; % fMRI length in scans
frameCut = 0; % delete the first 'frameCut' volumes
TR = 2.3; % Repetition time (seconds)
hpf = 120; % high-pass filter cutoff for GLM (in seconds)
lpf = 'AR(1)'; % 'none' | 'AR(1)' | 'FAST'
globnorm = 'None';

% !! fmri filename specifications, based on previously arranged input directory
fmri_prefix = '';
fmri_postfix = 'rfmri_filtered';
fmri_ext = '.nii';
gzipped = 0;

% !! set conditions for task fMRI (not used here)
cond_unit = 'scans';
%cond_onset = [7:20:5*20];     %cond_onset = [7 27 47 67 87];
%cond_dur = repmat(10, [1, 5]);  %cond_dur = [10 10 10 10 10];

if motcorr == 1
    contVec = [0 0 0 0 0 0 1];
else
    contVec = 1;
end

% !! adjusted region coordinates file
regions_prefix = '';
regions_postfix = 'networks';
regions_ext = '.csv';

% !! logging of VOI extraction 
VOIerrCnt = 2;
clear VOIerr;
VOIerr(1).subj = 'VOI-errors';

% !! directory for time-series output (optional, when modelCnt is set to 99)
%ts_path = strcat(project_path, DCM_dir, '/ts/');
ts_path = strcat(project_path, '/ts/');



% MAIN SCRIPT
subjects = dir(strcat(source_path));
subjects = subjects(3:end);
networkdirs = dir(strcat(regions_path, 'ROI_*'));


for subjCnt = [1:43] %for 43 diabetes subjects
%for subjCnt = [1:27] %for 27 obese subjects
    subjCnt

    % !! Get subject names based on input directory
    SubjID = strsplit( subjects(subjCnt).name, '.' );
    SubjID = SubjID(1);
    SubjID = SubjID{1};
    
    % !! Get region names and adjusted coordinates
    regions_file = strcat( regions_path, networkdirs(subjCnt).name, '/', regions_prefix, regions_postfix, regions_ext );
    regions = importdata(regions_file);
    Rnames = strcat(regions.textdata(2:end,3), '_', regions.textdata(2:end,4));
    RXmm = regions.data(:,1);
    RYmm = regions.data(:,2);
    RZmm = regions.data(:,3);
        
for paradigm = [1] %no different tasks used for this project
    
    if paradigm == 1
        ParID = 'rfMRI';
    end
    
for sessionCnt = [1] %no multiple sessions used for this project
    
    if sessionCnt == 1
        SessID = 's01';
    end
    
    
    try
        
        
    %SubjName = strcat(SubjID, '-', ParID, '_', SessID);
    SubjName = strcat(SubjID);
    
    fprintf('\n\n%s \n', SubjName);

    % !! Set data and output paths
    data_filename = strcat(fmri_prefix, fmri_postfix, fmri_ext);
    %data_filename = subjects(subjCnt).name;
    data_path_pref = strcat(source_path, subjects(subjCnt).name, '/');
    %data_path_post = strcat('/rfmri/');
    %data_path = strcat(data_path_pref);
    data_DCM_path = strcat(project_path, DCM_dir);
    %data_DCM_subdir = strcat(ParID, '_', SessID);
    %data_filepath = strcat(data_path_pref);
    data_filepath = fullfile(data_path_pref);
    results_DCM_path = strcat(project_path);
    mkdir(data_DCM_path);
    %mkdir(data_filepath);
    
    
%       if argGLM || argVOI
%         if exist(strcat(data_filepath, data_filename, '.gz'), 'file') ~= 2
%             fromfmri = strcat(source_path, source_disks(disk).name, '/', subjects(subjCnt).name, '/MNINonLinear/Results/', ParID, '_', SessID, '/', ParID, '_', SessID, '.nii.gz');
%             tofmri   = strcat(data_filepath, data_filename, '.gz');
%             link_cmd = strcat({'ln -s '}, {fromfmri}, {' '}, {tofmri});
% 
%             fprintf('Linking %s to %s\n', fromfmri, tofmri);
%             system(link_cmd{:});
% 
%             if (gzipped)
%                 fprintf('Gunzip fmri to %s \n', data_filepath);
%                 gunzip(strcat(data_filepath, data_filename, '.gz'), data_filepath);
%             end
% 
%         end
%       end  
    
    % !! Listing each volumes in a format usable for SPM
    stdVolList = {};
    for i = (frameCut+1):VolDim4
        tmpVol = int2str(i);
        stdVolList{end+1} = strcat(data_filepath, data_filename, ',', tmpVol);
    end
    
    stdVolList = transpose(stdVolList); 
    
    catch
        continue
        
    end

   % try
% Initialise SPM
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');
%spm_get_defaults('cmdline',1);



clear matlabbatch
matlabbatch{1}.cfg_basicio.cfg_mkdir.parent = { data_DCM_path };
matlabbatch{1}.cfg_basicio.cfg_mkdir.name = SubjName;
spm_jobman('run',matlabbatch);
clear matlabbatch
matlabbatch{1}.cfg_basicio.cfg_cd.dir = { strcat(data_DCM_path, '/', SubjName) };
spm_jobman('run',matlabbatch);



%%
% GLM for 1 session
%--------------------------------------------------------------------------

if argGLM == 1

clear matlabbatch
matlabbatch{1}.spm.stats.fmri_spec.dir = { strcat(data_DCM_path, '/', SubjName) };
matlabbatch{1}.spm.stats.fmri_spec.timing.units = cond_unit;
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
matlabbatch{1}.spm.stats.fmri_spec.sess.scans =    stdVolList  ;
if resting_state == 1
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
else
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond.name = 'act';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond.onset = cond_onset;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond.duration = cond_dur;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond.tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond.orth = 1;
end
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
if motcorr == 1
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = { motion_regress_path };
else
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = { '' };
end
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = hpf;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = globnorm;
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.4;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = lpf;

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'F';
if motcorr == 1
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = eye(7,7);
else
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = [1];
end
matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
if resting_state == 0
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'act>';
    if motcorr == 1
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1 0 0 0 0 0 0 0];
    else
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1 0];
    end
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
end
matlabbatch{3}.spm.stats.con.delete = 0;
%matlabbatch{4}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
%matlabbatch{4}.spm.stats.con.consess{1}.fcon.name = 'F';
%matlabbatch{4}.spm.stats.con.consess{1}.fcon.weights =  [-1 1];
%matlabbatch{4}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
%matlabbatch{4}.spm.stats.con.delete = 0;


% profile clear;
%     profile on -history;
%     spm_jobman('run',matlabbatch);
%     pinfo = profile('info');
%     assignin('base', 'profileInfo', pinfo);
%     save('profile_GLM.mat', 'pinfo');
%     profile off;
%     profile clear;

spm_jobman('run',matlabbatch);

end

%%
% VOI extraction
%--------------------------------------------------------------------------
if argVOI == 1

% img origo     [-2.5 2.5 2.5]
% spm origo     [85 120 65]

for roi = [1:length(Rnames)]
    
    if exist(strcat(regions_path, networkdirs(subjCnt).name, '/', Rnames{roi}, '.nii.gz'), 'file') == 2
        gunzip(strcat(regions_path, networkdirs(subjCnt).name, '/', Rnames{roi}, '.nii.gz'), strcat(regions_path, networkdirs(subjCnt).name))
    end
    roi_mask = strcat(regions_path, networkdirs(subjCnt).name, '/', Rnames{roi}, '.nii');

    clear matlabbatch
    try
    matlabbatch{1}.spm.util.voi.spmmat(1) = cellstr(fullfile(data_DCM_path, SubjName,'SPM.mat'));
    matlabbatch{1}.spm.util.voi.adjust = NaN;
    matlabbatch{1}.spm.util.voi.session = 1;
    matlabbatch{1}.spm.util.voi.name = Rnames{roi};
    %matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [RXmm(roi), RYmm(roi), RZmm(roi)];
    %matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 8;
    %matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
    matlabbatch{1}.spm.util.voi.roi{1}.mask.image = cellstr(roi_mask);
    matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
    matlabbatch{1}.spm.util.voi.roi{2}.mask.image = cellstr('mask.nii');
    matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
    matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';
    spm_jobman('run',matlabbatch);
    catch
       VOIerr(VOIerrCnt).subj = SubjName; 
       VOIerr(VOIerrCnt).voi = matlabbatch{1}.spm.util.voi.name; 
       VOIerrCnt = VOIerrCnt + 1;
    end
end

end




% !! Select models to estimate (currently 1 stands for full model and 99 is used for saving VOI time-series in CSV
for modelCnt = [1]
    
    if modelCnt == 99
        mkdir(ts_path);
    end

    


% !! DCM specification and estimation starts here

if argDCM == 1
    
    
clear DCM

% !! load SPM.mat
load(fullfile(data_DCM_path, SubjName,'SPM.mat'));

% !! load VOIs
k = 1;
%for i = [2 3 4 5 7 8 9 23:25 27 28] % you can use a subset of regions if you wish
for i = [1:length(Rnames)]
    load( fullfile(data_DCM_path, SubjName, strcat('VOI_', Rnames{i}, '_1.mat')), 'xY' );
    DCM.xY(k) = xY;
    k = k+1;
end

% !! make DCM structure
DCM.n = length(DCM.xY);      % number of regions
DCM.v = length(DCM.xY(1).u); % number of time points
if (resting_state == 0)
    uN = size(SPM.Sess.U, 2);    % number of inputs
end

DCM.Y.dt  = SPM.xY.RT;
DCM.Y.X0  = DCM.xY(1).X0;

for i = 1:DCM.n
    DCM.Y.y(:,i)  = DCM.xY(i).u;
    DCM.Y.name{i} = DCM.xY(i).name;
end

DCM.Y.Q    = spm_Ce(ones(1,DCM.n)*DCM.v);

if (resting_state == 0)
    DCM.U.dt   =  SPM.Sess.U(1).dt;
    DCM.U.name = [SPM.Sess.U.name];
    DCM.U.u    = [SPM.Sess.U(1).u(1:end,1)];
else
    % DCM.U.dt   = 0;
    DCM.U.name = cellstr('null');
    DCM.U.u    = zeros(VolDim4-frameCut,1);
end

DCM.delays = repmat(SPM.xY.RT,DCM.n,1);
DCM.TE     = 0.04;

DCM.a = zeros(DCM.n, DCM.n, 1);     % intrinsic
DCM.b = zeros(DCM.n, DCM.n, 0);     % condition effects
DCM.c = zeros(DCM.n, 1);            % model inputs
DCM.d = zeros(DCM.n, DCM.n, 0);     % nonlinear effects

DCM.options.nonlinear  = 0;
DCM.options.two_state  = 0;
DCM.options.stochastic = 0;
DCM.options.endogenous = 0;
DCM.options.induced    = 0;
DCM.options.centre     = 1;
DCM.options.nograph    = 1;

if (resting_state == 0)
    DCM.c = zeros(DCM.n, uN);
    DCM.b = zeros(DCM.n, DCM.n, uN);
end
if (DCM.options.nonlinear == 1)
    DCM.d = zeros(DCM.n, DCM.n, DCM.n);
end

if modelCnt == 1 % model number one is specified as the Mfull model for now
    
    % ----model specification-------------------------------------------|-
    % ----edit this part only-------------------------------------------*-
    
    %modelName = 'M2NW-test'; fprintf('\n%s\n', modelName);
    modelName = 'Mfull'; fprintf('\n%s\n', modelName);

    DCM.a = ones(DCM.n, DCM.n, 1);
    
    DCM_mode = 'csd';

    % DCM options are set based on estimation method
    if strcmp(DCM_mode, 'sto') % stochastic DCM
      DCM.options.endogenous = 1;
      DCM.options.stochastic = 1;
    end
    if strcmp(DCM_mode, 'csd') % cross-spectra DMC
      DCM.options.induced    = 1;
    end
    
    DCM_filename = strcat('DCM_', SubjName, '_', modelName,'_', DCM_mode, '.mat');
       
    % ------------------------------------------------------------------^-
    
    % !! saving DCM structure to MAT file
    save( fullfile(data_DCM_path,SubjName,DCM_filename),'DCM');
    continue
    
    % !! estimate DCM parameters - this overwrites the original DCM.mat file
    if strcmp(DCM_mode, 'sto')
        DCM = spm_dcm_estimate( fullfile(data_DCM_path,SubjName,DCM_filename) );
    end
    if strcmp(DCM_mode, 'csd')
        DCM = spm_dcm_fmri_csd( fullfile(data_DCM_path,SubjName,DCM_filename) );
    end

    % !! saving results in text format as well
    saveDCMResults(DCM, data_DCM_path, results_DCM_path, modelName, SubjName, resting_state, DCM_filename, DCM_mode);

end


% When 99 is selected for the model to estimate, it instead writes time-series data to CSV
if modelCnt == 99
    

try

    mkdir(ts_path);
    
    header = DCM.Y.name(1);
    for v = 2:size(DCM.Y.name, 2)
        header = strcat( header, ',', DCM.Y.name(v));
    end
    
    csvname = strcat(ts_path, SubjName, '_ts.csv');
    
    dlmwrite(csvname, header, 'delimiter', '');
    dlmwrite(csvname, DCM.Y.y, '-append', 'delimiter', ',');
    
catch
    
    print('Warning: write timeseries\n')
    
end
    
end


end



end

%     catch
% 
%         try
%             rm_cmd = strcat({'rm '}, {data_filepath}, {data_filename});
%             %system(rm_cmd{:});
%         catch
% 
%         end
%         
%         continue
%     end
% 
%     
%     try
%         rm_cmd = strcat({'rm '}, {data_filepath}, {data_filename});
%         %system(rm_cmd{:});
%     catch
% 
%     end


end


end



assignin('base', 'voierr', VOIerr);

end

fprintf('\nfinished \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINISH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
