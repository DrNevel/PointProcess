% ONLY RR
% f='p029164-2153-02-10-01-19_2153-02-11 16_00_00_1_hour_annot';
% pth='C:\Users\maxim\Dropbox\Dottorato\TesistiMagistrali\Stefanie Messner\WF';
% fpth=fullfile(pth,f);
% 
% Features=get_features_RR(fpth,'fs',125,'PPorderMono',7,'get_PP',1,'get_spectra',1,'get_TimeDomain',1,'get_Freqs',1,'get_Compl',1);


f='p002332-2151-04-11-20-58_annotations.mat';
pth='C:\Users\maxim\Dropbox\Dottorato\Codes\Matlab\Functions\Feature Package\Test\annotations';
fpth=fullfile(pth,f);

% Features_RR=get_features_RR(fpth,'sel_EKGR',1,'fs',125,'PPorderMono',9,'get_PP',1,'get_spectra',1,'get_TimeDomain',1,'get_Freqs',1,'get_Compl',1);
% Features_COV=get_features(fpth,'sel_EKGR',1,'fs',125,'PPorderBiva',7,'get_PP',1,'get_biva',1,'UndSampl',20,'get_spectra',1,'get_TimeDomain',1,'get_Freqs',1,'get_Compl',1);


f='p002332-2151-04-11-20-58_annotations2.mat';
pth='C:\Users\maxim\OneDrive - Politecnico di Milano\DOTTORATO\Studies Data\Sepsis_Identification_1stHR\Annotations2';
fpth=fullfile(pth,f);
% Features_RR=get_features_RR(fpth,'sel_EKGR',1,'fs',125,'PPorderMono',9,'get_PP',1,'UndSampl',20,'get_spectra',1,'get_TimeDomain',1,'get_Freqs',1,'get_Compl',1);
Features_COV=get_features(fpth,'sel_EKGR',1,'fs',125,'PPorderBiva',7,'get_PP',1,'get_biva',1,'UndSampl',20,'get_spectra',1,'get_TimeDomain',1,'get_Freqs',1,'get_Compl',1);


% 