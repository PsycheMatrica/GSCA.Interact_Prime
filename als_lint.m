function [W, C, B, b0, Gamma] = als_lint(Z, W0, A0, W, A, V, nnlv_index, ind_sign, WE, nvar, nlv, itmax, ceps)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ALS_Basic() - MATLAB function to implement the ALS algorithm for        %  
    %               Generalized Structured Component Analysis (GSCA) with     % 
    %               component interactions.                                   %
    % Author: Heungsun Hwang & Gyeongcheol Cho                                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nobs = size(Z,1);
    ntv = nvar + nlv;
    nnlv = size(nnlv_index,1);       
    ntlv = nlv + nnlv;                
    
    Gamma = Z*W;
    Gg_sd=sqrt(diag(Gamma'*Gamma)'./nobs);
    Gamma=Gamma./repmat(Gg_sd,[nobs,1]);
    NLV = lint_score(nnlv_index,Gamma);
    Gamma=[Gamma, NLV];
    
    W=W./repmat(Gg_sd,[nvar,1]);
    V(:,(nvar+1):end)=W;
    
    Psi = Z*V; 
    it = 0;                 
    imp = 100000;         
    f0 = 10^10;            
    while it <= itmax && imp > ceps 
          it = it+1;                    
          for t = 1:ntv  
              H1 = eye(ntv); 
              H1(t,t) = 0;
              aindex = A0(:,t);
              na=sum(aindex,1);
              if na>0 
                  a = A(:,t);
                  a(aindex) = 0;           
                  e = zeros(1,ntv); 
                  e(t) = 1;
                  Y = Psi(:,t) - Gamma*(A*H1 + a*e); 
                  X = Gamma(:,aindex);
                  A(aindex,t) = (X'*X)\(X'*Y*e');
              end
          end
          for p = 1:nlv
              windex = W0(:,p);
              nw = sum(windex,1);
              if nw>1
                  Zp = Z(:,windex);
                  we_index = find(WE(p,:));
                  nnzwe = length(we_index);              
                  H1 = ones(1,ntlv);
                  H2 = ones(1,ntv);
                  if nnzwe > 1                          
                     M1 = zeros(ntv*nobs,nw);
                     for i = 1:nnzwe
                         e = zeros(1,ntv); 
                         if we_index(i) > nlv           
                            num_lv = WE(p,we_index(i));  
                            LV = Gamma(:,num_lv);
                            DIAGLV = LV*ones(1,nw);
                            DZj = DIAGLV.*Zp;              
                            a = A(we_index(i),:);
                            beta = -a;
                            M1 = M1 + kron(beta',DZj);
                         else                 
                            e(nvar+p) = 1;
                            a = A(p,:);
                            beta = e-a;
                            M1 = M1 + kron(beta',Zp);
                         end
                         H1(we_index(i)) = 0;
                         if (nvar + we_index(i)) <= ntv
                             H2(nvar + we_index(i)) = 0;
                         end
                     end
                     H1 = diag(H1);
                     H2 = diag(H2);
                     Y = Gamma*H1*A - Psi*H2;  
                     M2 = reshape(Y,nobs*ntv,1);
                 else      
                     e = zeros(1,ntv); 
                     e(nvar+p) = 1;
                     a = A(p,:);
                     beta = e - a;
                     q = beta*beta';
                     H1(we_index) = 0;
                     H2(nvar + we_index) = 0;
                     H1 = diag(H1);
                     H2 = diag(H2);
                     Y = Gamma*H1*A - Psi*H2; 
                     M1 = q*(Zp'*Zp);
                     M2 = Zp'*Y*beta';
                 end          
                 theta = M1\M2;     
                 zw = Zp*theta; 
                 theta = theta/sqrt((zw'*zw)/nobs);             
                 W(windex,p) = theta;
                 V(windex,nvar+p) = theta;
              end
          end
          Gamma1 = Z*W;
          NLV = lint_score(nnlv_index,Gamma1);
          Gamma = [Gamma1, NLV];
          Psi = Z*V;       
          dif = Psi-Gamma*A;          
          f = sum(sum(dif.*dif)); 
          imp = f0-f;
          f0 = f;
    end

    for p=1:nlv
       if ind_sign(p)~=0
           if Z(:,ind_sign(p))'*Gamma1(:,p) < 0
              Gamma1(:,p)=-Gamma1(:,p);
              W(:,p)=-W(:,p);
           end
       end
    end    
    NLV = lint_score(nnlv_index,Gamma1);
    Gamma = [Gamma1, NLV];

    for t = 1:ntv 
       H1 = eye(ntv); 
       H1(t,t) = 0;
       aindex = A0(:,t);       
       a = A(:,t);
       a(aindex) = 0;                     
       e = zeros(1,ntv); 
       e(t) = 1;
       Y = Psi - Gamma*(A*H1 + a*e);      
       X = Gamma(:,aindex);
       A(aindex,t) = (X'*X)\(X'*Y*e');
    end

    W=W(2:end,2:end);
    C=A(2:(end-nnlv),2:nvar);
    b0=A(1,(nvar+2):end);
    B=A(2:end,(nvar+2):end);
