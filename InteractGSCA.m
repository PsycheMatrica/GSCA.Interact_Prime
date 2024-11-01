function [INI, TABLE, ETC] =  InteractGSCA(z0, W0, C0, B0, nnlv_index, ind_sign, N_Boot, Max_iter, Min_limit,Flag_Parallel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% InteractGSCA() - MATLAB function to perform Generalized Structured      %
%                  Component Analysis (GSCA) with Component Interactions  %
% Author: Heungsun Hwang & Gyeongcheol Cho                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:                                                        %
%   Data = an N by J matrix of scores for N individuals on J indicators   %
%   W0 = a J by P matrix of weight parameters                             %
%   C0 = a P by J matrix of loading parameters                            %
%   B0 = a (P + P_int) by P matrix of path coefficients, where P_int      %
%          represents the number of component interaction terms           %     
%   nnlv_index = a P_int by 2 matrix of indices for component interactions%
%   ind_sign = a P by 1 vector whose p-th element represents the number   %
%               of the sign-fixing indicator for the p-th component       % 
%   N_Boot = Integer representing the number of bootstrap samples for     %
%            calculating standard errors (SE) and 95% confidence          %
%            intervals (CI)                                               %
%   Max_iter = Maximum number of iterations for the Alternating Least     % 
%              Squares (ALS) algorithm                                    %
%   Min_limit = Tolerance level for ALS algorithm                         %
%   Flag_Parallel = Logical value to determine whether to use parallel    %
%                   computing for bootstrapping                           %
% Output arguments:                                                       %
%   INI: Structure array containing goodness-of-fit values, R-squared     % 
%        values, and matrices parameter estimates                         %
%     .W: a J by P matrix of weight estimates                             %
%     .C: a P by J matrix of loading estimates                            %
%     .b0: a 1 by P matrix of path coefficient estimates                  %
%     .B: a (P + P_int) by P matrix of path coefficient estimates         %
%  TABLE: Structure array containing tables of parameter estimates, their %
%         SEs, 95% CIs,and other statistics                               %
%     .W: Table for weight estimates                                      %
%     .C: Table for loading estimates                                     %
%     .b0: Table for intercept estimates                                  %
%     .B: Table for path coefficients estimates                           %
%  ETC: Structure array including bootstrapped parameter estmates         %
%     .W_Boot: Matrix of bootstrapped weight estimates                    %
%     .C_Boot: Matrix of bootstrapped loading estimates                   %
%     .b0_Boot: Matrix of bootstrapped intercept estimates                %
%     .B_Boot: Matrix of bootstrapped path coefficient estimates          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z0 = zscore(z0,1);
N = size(z0,1);
ones_N=ones(N,1);
Z_p1 = [ones_N,z0];

W0=(W0~=0); Nw=sum(sum(W0,2),1);
C0=(C0~=0); Nc=sum(sum(C0,2),1);
B0=(B0~=0); Nb=sum(sum(B0,2),1);

W0_p1 = blkdiag(true,W0);
C0_p1 = blkdiag(true,C0);

b0 = true(1,size(B0,2));
b0(1,sum(B0,1)==0) = false; Nb0=sum(b0,2);
B0_ext = [false(size(B0,1)+1,1),[b0;B0]];
nnlv_index=nnlv_index+1;

[J_p1,P_p1] = size(W0_p1);
P_int = size(nnlv_index,1);       
C0_ext = [C0_p1;false(P_int,J_p1)];  
A0_ext = [C0_ext,B0_ext];

V001 = eye(J_p1);
V00 = [V001,W0_p1];

V00 = V00~=0;

W_p1 = W0_p1; A_ext = A0_ext;V = V00;
W_p1=double(W_p1);
A_ext=double(A_ext);
V=double(V);

WE = we(nnlv_index, P_p1);

ind_sign_ext=[0,ind_sign+1];

[est_W, est_C, est_B, est_b0]  = als_lint(Z_p1, W0_p1, A0_ext, W_p1, A_ext, V, nnlv_index,ind_sign_ext, WE, J_p1, P_p1, Max_iter, Min_limit);

INI.W=est_W;
INI.C=est_C;
INI.B=est_B;
INI.b0=est_b0;

if N_Boot<100
    TABLE.W=[est_W(W0),NaN(Nw,5)];
    TABLE.C=[est_C(C0),NaN(Nc,5)];
    TABLE.B=[est_B(B0),NaN(Nb,5)];
    TABLE.b0=[est_b0(b0),NaN(Nb0,5)];
    ETC.W_Boot=[];
    ETC.C_Boot=[];
    ETC.B_Boot=[];  
    ETC.b0_Boot=[];  
else
    W_Boot=zeros(Nw,N_Boot);
    C_Boot=zeros(Nc,N_Boot);
    B_Boot=zeros(Nb,N_Boot);  
    b0_Boot=zeros(Nb0,N_Boot);
    if Flag_Parallel
        parfor b=1:N_Boot
            [Z_ib,~]=GC_Boot(z0); 
            Z_ib_ext = [ones_N,zscore(Z_ib,1)];
            [W_b, C_b, B_b, b0_b] = als_lint(Z_ib_ext, W0_p1, A0_ext, W_p1, A_ext, V, nnlv_index,ind_sign_ext, WE, J_p1, P_p1, Max_iter, Min_limit);
            W_Boot(:,b)=W_b(W0);
            C_Boot(:,b)=C_b(C0);
            B_Boot(:,b)=B_b(B0);
            b0_Boot(:,b)=b0_b(b0);
        end
    else
        for b=1:N_Boot
            [Z_ib,~]=GC_Boot(z0); 
            Z_ib_ext = [ones_N,zscore(Z_ib,1)];
            [W_b, C_b, B_b, b0_b] = als_lint(Z_ib_ext, W0_p1, A0_ext, W_p1, A_ext, V, nnlv_index,ind_sign_ext, WE, J_p1, P_p1, Max_iter, Min_limit);
            W_Boot(:,b)=W_b(W0);
            C_Boot(:,b)=C_b(C0);
            B_Boot(:,b)=B_b(B0);
            b0_Boot(:,b)=b0_b(b0);
        end
    end
    alpha=.05;
    CI=[alpha/2,alpha,1-alpha,1-(alpha/2)];
    loc_CI=round(CI*(N_Boot-1))+1; % .025 .05 .95 .975
    
    TABLE.W=para_stat(est_W(W0),W_Boot,loc_CI);
    TABLE.C=para_stat(est_C(C0),C_Boot,loc_CI);
    TABLE.B=para_stat(est_B(B0),B_Boot,loc_CI); 
    TABLE.b0=para_stat(est_b0(b0),b0_Boot,loc_CI);

    ETC.W_Boot=W_Boot;
    ETC.C_Boot=C_Boot;
    ETC.B_Boot=B_Boot;        
    ETC.b0_Boot=b0_Boot;
end
end
function Table=para_stat(est_mt,boot_mt,CI_mp)
    boot_mt=sort(boot_mt,2);
    SE=std(boot_mt,0,2);
    Table=[est_mt,SE,boot_mt(:,CI_mp(1,1)),boot_mt(:,CI_mp(1,4))]; 
end
function [in_sample,out_sample,index,N_oob]=GC_Boot(Data)
    N=size(Data,1); 
    index=ceil(N*rand(N,1));
    in_sample=Data(index,:); 
    index_oob=(1:N)'; index_oob(index)=[];
    out_sample=Data(index_oob,:);
    N_oob=length(index_oob);
end

