function testplot(Features,EKGR_old, COV_V)

figure, ax = subplot(7,1,1);plot(EKGR_old-EKGR_old(1):0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.Mu,...
    EKGR_old-EKGR_old(1):0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.meanRR),legend('Mu','meanRR'),ylabel('sec')
ax = [ax,subplot(7,1,2)]; plot(EKGR_old-EKGR_old(1),COV_V),legend('SAP'),ylabel('mmHg')
ax = [ax,subplot(7,1,3)]; imagesc(0:0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.Baro.FS(:,end),Features.PP.Cov.Baro.GAIN_2to1,[0 0.003]),ylabel('S->RR s/mmHg')
ax = [ax,subplot(7,1,4)]; imagesc(0:0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.Baro.FS(:,end),Features.PP.Cov.Baro.PHASE_2to1,[-3.14 3.14]),ylabel('S->RR rad')
ax = [ax,subplot(7,1,5)]; imagesc(0:0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.Baro.FS(:,end),Features.PP.Cov.Baro.GAIN_1to2,[0 500]),ylabel('RR->S mmHg/s')
ax = [ax,subplot(7,1,6)]; imagesc(0:0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.Baro.FS(:,end),Features.PP.Cov.Baro.PHASE_1to2,[-3.14 3.14]),ylabel('RR->S rad')
ax = [ax,subplot(7,1,7)]; imagesc(0:0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.Baro.FS(:,end),Features.PP.Cov.Baro.COH,[0 1]),ylabel('COH')
linkaxes(ax,'x')
figure, ax = subplot(7,1,1);plot(EKGR_old-EKGR_old(1):0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.Mu,...
    EKGR_old-EKGR_old(1):0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.meanRR),legend('Mu','meanRR'),ylabel('sec')
ax = [ax,subplot(7,1,2)]; plot(EKGR_old-EKGR_old(1),COV_V),legend('SAP'),ylabel('mmHg')
ax = [ax,subplot(7,1,3)]; plot(0:0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.Baro.GAIN_21_LF,...
    0:0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.Baro.GAIN_21_HF),legend('LF','HF'),ylabel('S->RR s/mmHg')
ax = [ax,subplot(7,1,5)]; plot(0:0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.Baro.GAIN_12_LF,...
    0:0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.Baro.GAIN_12_HF),legend('LF','HF'),ylabel('RR->S mmHg/s')
ax = [ax,subplot(7,1,4)]; plot(0:0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.Baro.PHASE_21_LF,...
    0:0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.Baro.PHASE_21_HF),legend('LF','HF'),ylabel('S->RR rad')
ax = [ax,subplot(7,1,6)]; plot(0:0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.Baro.PHASE_12_LF,...
    0:0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.Baro.PHASE_12_HF),legend('LF','HF'),ylabel('RR->S rad')
ax = [ax,subplot(7,1,7)]; plot(0:0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.Baro.COH_LF,...
    0:0.05:EKGR_old(end)-EKGR_old(1),Features.PP.Cov.Baro.COH_HF),legend('LF','HF'),ylabel('COH')
linkaxes(ax,'x')
