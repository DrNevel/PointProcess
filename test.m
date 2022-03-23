f='p002332-2151-04-11-20-58_annotations2';
pth='Test\annotations';
fpth=fullfile(pth,f);

f='12726_data';
pth='Sepsis_Data\Annotations4';
fpth=fullfile(pth,f);

%Features_RR=get_features_RR(fpth,'sel_EKGR',1,'fs',125,'PPorderMono',9,'get_PP',1,'UndSampl',20,'get_spectra',1,'get_TimeDomain',1,'get_Freqs',1,'get_Compl',1);
Features_PP_Spectral=get_features_RR(fpth,'sel_EKGR',1,'fs',125,'PPorderMono',9,'get_PP',1,'UndSampl',20,'get_spectra',1,'get_TimeDomain',0,'get_Freqs',0,'get_Compl',0);

figure
subplot(2,2,1)
plot(Features_PP_Spectral.PP_U.Mono.RR.powLF)
hold on
plot(Features_PP_Spectral.PP_new.Mono.RR.powLF)
legend('Old Uncensored LF','New LF')
subplot(2,2,2)
plot(Features_PP_Spectral.PP_U.Mono.RR.powHF)
hold on
plot(Features_PP_Spectral.PP_new.Mono.RR.powHF)
legend('Old Uncensored HF ','New HF')
subplot(2,2,3)
plot(Features_PP_Spectral.PP_U.Mono.RR.powVLF)
hold on
plot(Features_PP_Spectral.PP_new.Mono.RR.powVLF)
legend('Old Uncensored VLF','New VLF')
subplot(2,2,4)
plot(Features_PP_Spectral.PP_U.Mono.RR.powTot)
hold on
plot(Features_PP_Spectral.PP_new.Mono.RR.powTot)
legend('Old Uncensored TOT','New TOT')

Features_COV=get_features(fpth,'sel_EKGR',1,'fs',125,'PPorderBiva',7,'get_PP',1,'get_biva',1,'UndSampl',20,'get_spectra',1,'get_TimeDomain',1,'get_Freqs',1,'get_Compl',1);


% prova