clear,clc
% pth = 'C:\Users\maxim\OneDrive - Politecnico di Milano\DOTTORATO\Studies Data\Sepsis_Identification_1stHR\Annotations';
% save_pth_PP = 'C:\Users\maxim\OneDrive\Desktop\PointProcessed_13';
% save_pth_PP_biva = 'E:\Sepsis1HRStudy\PPData\PointProcessed_9_biva_S';
% save_pth_hour = 'C:\Users\maxim\OneDrive\Desktop\Featured_1Hour';
% save_pth_15mins = 'C:\Users\maxim\OneDrive\Desktop\Featured_15Mins';
% % mkdir(save_pth_PP);
% mkdir(save_pth_PP_biva);
% mkdir(save_pth_hour);
% mkdir(save_pth_15mins);

% MAC has to swtich \ with / !!!!! N.B: PP won't work in MAC!!
pth = strcat(cd,'\Test\annotations');
save_pth_PP_biva =  strcat(cd,'\Test\PP');
cd(pth)
mkdir(save_pth_PP_biva);

%%

list = what(pth);
list = cell2mat(list.mat);

% list2 = cell2mat(what(lst_pth).mat);
% list = setdiff(string(list(:,1:end-16)),string(list2(:,1:end-4)));
% list = strcat(list,'_annotations.mat')
% list = str2mat(list)

% for i = 1:size(REC,1)
%     s = num2str(REC{i,:});
%     wfdb2mat(s);
% end

%% Extract PP and Features
cov_type = 'syst'; % bivariate PP Mode
to_correct = [];
get_PPs = 1; % put 1 To Run PointProcess

for i = 1%:size(list,1)
    
    fprintf('%d / %d \n',i,size(list,1))
    clear Features RR_correction
    selected_wf = list(i,:);
    
    [Feat,RR_correction] = get_features(selected_wf,'rsmpl_ts',0,...
        'get_TimeDomain',1,'get_Freqs',1,'get_Compl',1,...
        'get_PP',get_PPs,'get_mono',0,'get_biva',1,'get_spectra',1,'fres',512,...
        'PPorderMono',13,'PPorderBiva',4,'UndSampl',10,'cov_type',cov_type);
    if isfield(Feat,'error')
        to_correct = [to_correct,i];
        continue
    end
    % Saving in existent files
    %     if exist(strcat(save_pth_PP_biva,'\',selected_wf(1:end-16),'.mat'),'file')
    %         load(strcat(save_pth_PP_biva,'\',selected_wf(1:end-16)))
    %         if strcmp(cov_type,'syst') && get_PPs == 1
    %             Features.PP = Feat.PP;
    %         elseif strcmp(cov_type,'diast') && get_PPs == 1
    %             Features.PP_DIAST = Feat.PP;
    %         elseif get_PPs == 0 && get_class == 1
    %             Features.CLASS = Feat.CLASS;
    %             Features.SIGS = Feat.SIGS;
    %         end
    %         if get_class == 1 && get_PPs == 1
    %             Features.CLASS = Feat.CLASS;
    %             Features.SIGS = Feat.SIGS;
    %         end
    %     else
    %         warning('New File is going to be written, press a button to continue...')
    % %         pause
    Features = Feat;
    %     end
    save(strcat(save_pth_PP_biva,'\',selected_wf(1:end-16)),'Features','RR_correction','-v7.3');
    
end
to_correct
%