% clear all
% load('OutputData.mat')
% CN = (CN - CN(1,:))./CN(1,:);
% CN1 = CN;
% itr1 = timestep;
% save('CNp25.mat','CN1','itr1')
% % 
% clear all
% load('OutputData.mat')
% CN = (CN - CN(1,:))./CN(1,:);
% CN2 = CN;
% itr2 = timestep;
% save('CNp25.mat','CN2','itr2','-append')
% % 
% % 
% clear all
% load('OutputData.mat')
% CN = (CN - CN(1,:))./CN(1,:);
% CN3 = CN;
% itr3 = timestep;
% save('CNp25.mat','CN3','itr3','-append')
% 
% i = 600;
% CNm = (CN1(1:i,:) + CN2(1:i,:) + CN3(1:i,:))/3;
% itrm = (itr1(1:i) + itr2(1:i) + itr3(1:i))/3;
% save('CNp25.mat','CNm','itrm','-append')


clear all
load('OutputData.mat')
CN = (CN - CN(1,:));
CN1 = CN;
itr1 = timestep;
save('CNp25_TN.mat','CN1','itr1')
% 
clear all
load('OutputData.mat')
CN = (CN - CN(1,:));
CN2 = CN;
itr2 = timestep;
save('CNp25_TN.mat','CN2','itr2','-append')
% 
% 
clear all
load('OutputData.mat')
CN = (CN - CN(1,:));
CN3 = CN;
itr3 = timestep;
save('CNp25_TN.mat','CN3','itr3','-append')

i = 600;
CNm = (CN1(1:i,:) + CN2(1:i,:) + CN3(1:i,:))/3;
itrm = (itr1(1:i) + itr2(1:i) + itr3(1:i))/3;
save('CNp25_TN.mat','CNm','itrm','-append')


% 

