
msg = 'building PEB model'
% all estimated DCM files copied to a single directory (RSN_36/DCM)
DCM_dir = '/path/to/diabetes/DCM/largescale/RSN_36/DCM';
cd(DCM_dir);
files = dir(fullfile(DCM_dir, 'DCM*'));
files = struct2cell(files);
files = files(1,1:70)'; % 1-72 obese diabetes, 1-42 50 diab, 43-49 51-70 obese
nregions = 36;
nparams = nregions*nregions;

% BMR_dir = '/mnt/raid6_data/user/aranyics/diabetes/hu/process/DCM/largescale/estimated/csd/BMR';
% cd(BMR_dir);
% BMR_files = dir(fullfile(BMR_dir, 'BMR_DCM*'));
% BMR_files = struct2cell(BMR_files);
% BMR_files = BMR_files(1,:)';

PEB_dir = '/path/to/diabetes/DCM/largescale/RSN_36/PEB/';

field = {'A'};



L1d.GCM = files( [1:42 51] );
L1d.X = [ ones(length(L1d.GCM), 1) ];
L1d.Xnames = {'diab'};
L1d.comparison = 'hierarch_diab_obese_avg_diff_missproba';

L1o.GCM = files( [45:50 52:72] );
L1o.X = [ ones(length(L1o.GCM), 1) ];
L1o.Xnames = {'obes'};
L1o.comparison = 'hierarch_diab_obese_avg_diff_missproba';

L2.X = [1 1; 1 -1];
L2.Xnames = {'mean', 'diab-obes'};
L2.comparison = 'hierarch_diab_obese_avg_diff_missproba';

L1do.GCM = [ L1d.GCM; L1o.GCM ];
L1do.X = [ ones(length(L1d.GCM), 1) ones(length(L1d.GCM), 1); ...
           ones(length(L1o.GCM), 1) -1*ones(length(L1o.GCM), 1)];
L1do.X(:,2) = L1do.X(:,2) - mean(L1do.X(:,2));
L1do.Xnames = {'mean', 'diab-obes'};
L1do.comparison = 'avg_diff_missproba';


cd(DCM_dir);



% Group-level models of diabetes and obesity with PEB and BMR
msg='estimating PEB model and saving PEB and DCM structures'
[PEB1 DCM1] = spm_dcm_peb( L1d.GCM, L1d, field );
[PEB2 DCM2] = spm_dcm_peb( L1o.GCM, L1o, field );
mkdir( fullfile(PEB_dir, L2.comparison) );
save( fullfile(PEB_dir, L2.comparison, strcat('PEB_DCM_', L2.comparison, '_L1_diab')), 'PEB1' );
save( fullfile(PEB_dir, L2.comparison, strcat('PEB_DCM_', L2.comparison, '_L1_obes')), 'PEB2' );
save( fullfile(PEB_dir, L2.comparison, strcat('DCM_', L2.comparison, '_L1_diab')), 'DCM1' );
save( fullfile(PEB_dir, L2.comparison, strcat('DCM_', L2.comparison, '_L1_obes')), 'DCM2' );

msg='saving PEB connectivity matrix'
A1 = reshape( PEB1.Ep, nregions, nregions );
A2 = reshape( PEB2.Ep, nregions, nregions );
dlmwrite( fullfile(PEB_dir, L2.comparison, strcat('peb_', L2.comparison, '_L1_', 'diab', '_A.csv')), A1, 'delimiter', ',' );
dlmwrite( fullfile(PEB_dir, L2.comparison, strcat('peb_', L2.comparison, '_L2_', 'obese', '_A.csv')), A2, 'delimiter', ',' );

msg='computing BMR on PEB results and saving BMR structures'
% load( fullfile(PEB_dir, L2.comparison, strcat('PEB_DCM_', L2.comparison)), 'PEB' );
[BMAPEB1] = spm_dcm_peb_bmc(PEB1);
save( fullfile(PEB_dir, L2.comparison, strcat('BMAPEB_DCM_', L2.comparison, '_diab')), 'BMAPEB1' );
[BMAPEB2] = spm_dcm_peb_bmc(PEB2);
save( fullfile(PEB_dir, L2.comparison, strcat('BMAPEB_DCM_', L2.comparison, '_obese')), 'BMAPEB2' );

msg='saving BMR connectivity and probability matrix'
A1 = reshape( BMAPEB1.Ep, nregions, nregions );
A2 = reshape( BMAPEB2.Ep, nregions, nregions );
pA1 = reshape( BMAPEB1.Pp, nregions, nregions );
pA2 = reshape( BMAPEB2.Pp, nregions, nregions );
dlmwrite( fullfile(PEB_dir, L2.comparison, strcat('bmra_', L2.comparison, '_L1_', 'diab', '_A.csv')), A1, 'delimiter', ',' );
dlmwrite( fullfile(PEB_dir, L2.comparison, strcat('bmra_', L2.comparison, '_L2_', 'obese', '_A.csv')), A2, 'delimiter', ',' );
dlmwrite( fullfile(PEB_dir, L2.comparison, strcat('bmra_', L2.comparison, '_L1_', 'diab', '_pA.csv')), pA1, '-append', 'delimiter', ',' );
dlmwrite( fullfile(PEB_dir, L2.comparison, strcat('bmra_', L2.comparison, '_L2_', 'obese', '_pA.csv')), pA2, '-append', 'delimiter', ',' );



%% Population mean and group-level differences PEB and BMR
msg='estimating PEB model and saving PEB and DCM structures'
[PEB DCM] = spm_dcm_peb( L1do.GCM, L1do, field);
mkdir( fullfile(PEB_dir, L1do.comparison) );
save( fullfile(PEB_dir, L1do.comparison, strcat('PEB_DCM_', L1do.comparison, '_L2_diff')), 'PEB' );
%save( fullfile(PEB_dir, L1do.comparison, strcat('DCM_', L1do.comparison, '_L2_diff')), 'DCM' );

msg='saving improved DCM connectivity and probability matrices'
mkdir( fullfile(PEB_dir, L1do.comparison, 'DCM') );
for subj = 1:length(DCM)
   [p b e] = fileparts(L1do.GCM{subj});
	dlmwrite( fullfile(PEB_dir, L1do.comparison, 'DCM', strcat('econn_dcm_peb_est_', b, '.csv')), DCM{subj}.Ep.A, '-append', 'delimiter', ',' );
	dlmwrite( fullfile(PEB_dir, L1do.comparison, 'DCM', strcat('peconn_dcm_peb_est_', b, '.csv')), DCM{subj}.Pp.A, '-append', 'delimiter', ',' );
	dlmwrite( fullfile(PEB_dir, L1do.comparison, 'DCM', strcat('ceconn_dcm_peb_est_', b, '_diag.csv')), reshape(diag(DCM{subj}.Cp(1:nparams, 1:nparams)), nregions, nregions), '-append', 'delimiter', ',' );
end

msg='saving PEB connectivity matrix'
A11 = reshape( PEB.Ep(1:nparams), nregions, nregions );
A12 = reshape( PEB.Ep((nparams+1):2*nparams), nregions, nregions );
dlmwrite( fullfile(PEB_dir, L1do.comparison, strcat('peb_', L1do.comparison, '_L11_', 'mean', '_A.csv')), A11, '-append', 'delimiter', ',' );
dlmwrite( fullfile(PEB_dir, L1do.comparison, strcat('peb_', L1do.comparison, '_L12_', 'diff', '_A.csv')), A12, '-append', 'delimiter', ',' );

msg='computing BMR on PEB results and saving BMR structures'
% load( fullfile(PEB_dir, L2.comparison, strcat('PEB_DCM_', L2.comparison)), 'PEB' );
[BMAPEB] = spm_dcm_peb_bmc(PEB);
save( fullfile(PEB_dir, L1do.comparison, strcat('BMAPEB_DCM_', L1do.comparison, '_mean_diff')), 'BMAPEB' );

msg='saving BMR connectivity and probability matrix'
A11 = reshape( BMAPEB.Ep(1:nparams), nregions, nregions );
A12 = reshape( BMAPEB.Ep((nparams+1):2*nparams), nregions, nregions );
pA11 = reshape( BMAPEB.Pp(1:nparams), nregions, nregions );
pA12 = reshape( BMAPEB.Pp((nparams+1):2*nparams), nregions, nregions );
dlmwrite( fullfile(PEB_dir, L1do.comparison, strcat('bmra_', L1do.comparison, '_L11_', 'mean', '_A.csv')), A11, '-append', 'delimiter', ',' );
dlmwrite( fullfile(PEB_dir, L1do.comparison, strcat('bmra_', L1do.comparison, '_L12_', 'diff', '_A.csv')), A12, '-append', 'delimiter', ',' );
dlmwrite( fullfile(PEB_dir, L1do.comparison, strcat('bmra_', L1do.comparison, '_L11_', 'mean', '_pA.csv')), pA11, '-append', 'delimiter', ',' );
dlmwrite( fullfile(PEB_dir, L1do.comparison, strcat('bmra_', L1do.comparison, '_L12_', 'diff', '_pA.csv')), pA12, '-append', 'delimiter', ',' );



% matlab -nodisplay -nosplash -nodesktop -r "run('/path/to/diabetes/DCM/batch/script_save_PEB_model.m'); exit;"
