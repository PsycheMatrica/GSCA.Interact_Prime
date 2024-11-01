%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Illustration for InteractGSCA Prime package                             %
%   Author: Heungsun Hwang & Gyeongcheol Cho                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                            %
%   - This code aims to illustrate how to use InteractGSCA Prime package. %
%   - The dataset is generated from the model used in Shen, Cho and Hwang %
%     (under review).                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References                                                              %
%     * Hwang, H., Ho, M.-H. R., & Lee, J. (2010). Generalized structured %
%         component analysis with latent interactions. Psychometrika,     %
%         75(2), 228–242. https://doi.org/10.1007/s11336-010-9157-5       %
%     * Shen, Z., Cho, G., & Hwang, H. (under review)  Comparison of      %
%         component-based structural equation modeling methods in testing %
%         component interaction effects.                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
load('Example_data.mat')
z0=Dataset;
P = 3;
P_int = 1;
P_total = 4;
Jp = 5;
J = Jp * P;

wp  = ones(Jp,1);
W0 = blkdiag(wp,wp,wp)*99;
C0 = W0';
nnlv_index=[1,2];
B0 = zeros(P_total,P);
B0([1,2,4],3)=99;
ind_sign=[1,6,11];
N_Boot = 0;
Flag_Parallel=false;

Max_iter = 100; 
Min_limit = .00001;

[INI, TABLE, ETC] = InteractGSCA(z0, W0, C0, B0, nnlv_index, ind_sign, N_Boot, Max_iter, Min_limit,Flag_Parallel);

INI.W 
INI.C
INI.B 
INI.b0

TABLE.W
TABLE.C
TABLE.B
TABLE.b0


