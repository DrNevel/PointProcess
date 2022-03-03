clear; close all; clc;

% se usi questi, imposta 'fs' a 125 Hz
% f='p002332-2151-04-11-20-58_annotations2';
% pth='Test\annotations';
% fpth=fullfile(pth,f);

% f='p002332-2151-04-11-20-58_annotations2_short';
% pth='Sepsis_Data\Annotations2';
% fpth=fullfile(pth,f);

% se usi questo, imposta 'fs' a 250 Hz
f='12726_data.mat';
pth='Sepsis_Data\Annotations2';
fpth=fullfile(pth,f);

% about 17 seconds
% tic
% Features_RR=get_features_RR(fpth,'sel_EKGR',1,'fs',250,'PPorderMono',9,'get_PP',1,'UndSampl',20,'get_spectra',1,'get_TimeDomain',1,'get_Freqs',1,'get_Compl',1);
Features_RR=get_features_multi(fpth,'sel_EKGR',1,'fs',250,'PPorderMono',8,'get_PP',1,'UndSampl',20,'get_spectra',1);
% Features_RR=get_features_copy(fpth,'sel_EKGR',1,'fs',250,'PPorderMono',8,'get_PP',1,'get_mono',1,'UndSampl',20,'get_spectra',1);
% toc

% about 280 seconds
% tic
% Features_COV=get_features(fpth,'sel_EKGR',1,'fs',125,'PPorderBiva',7,'get_PP',1,'get_biva',1,'UndSampl',20,'get_spectra',1,'get_TimeDomain',1,'get_Freqs',1,'get_Compl',1);
Features_COV=get_features_multi(fpth,'sel_EKGR',1,'fs',250,'PPorderBiva',8,'get_PP',1,'get_biva',1,'UndSampl',20,'get_spectra',1);
% Features_COV=get_features_copy(fpth,'sel_EKGR',1,'fs',250,'PPorderBiva',8,'get_PP',1,'get_biva',1,'UndSampl',20,'get_spectra',1);
% toc

% about 610 seconds
tic
Features_TRI=get_features_multi(fpth,'sel_EKGR',1,'fs',250,'PPorderMulti',8,'get_PP',1,'get_multi',1,'UndSampl',20,'get_spectra',1,'cov_type','pat-syst','correct_diast',1);
%Features_TRI=get_features_copy(fpth,'sel_EKGR',1,'fs',250,'PPorderMulti',8,'get_PP',1,'get_multi',1,'UndSampl',20,'get_spectra',1,'cov_type','pat-syst');
toc

% 