% This script file test the gradient descent maxent (GDME) algorithm
%


clear all;
showflag = 1;

% function [E_rmse, E_aad, E_aid, E_sad, E_sid, E_time] = gdme_test(SNR,showflag)

%% reflectance from USGS library
load A;
load BANDS; % BANDS: selected 188 band index from original 224 bands
type = 5;
c = 4; estc = c;
A = A(BANDS,1:c);

% Generate simulated data
[mixed, abf] = getSynData(A, type, 7, 1, c-1, 1);
[M,N,Band] = size(mixed);


% Add Gaussian noise
SNR = 10; 
variance = sum(mixed(:).^2)/10^(SNR/10)/M/N/Band;
n = sqrt(variance)*randn([M,N,Band]);
mixed = mixed+n;          
mixed = reshape(mixed,M*N,Band)';  % column:bands, row:samples          

% Test different algorithms
[A_gdme, s_gdme, t_gdme] = gdme(mixed, SNR, 1, estc, A);
[A_varnt, s_varnt, t_varnt] = gdme(mixed, SNR, 2, estc, A);
    
E_time = [t_gdme, t_varnt];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for method = 1:2
    if method == 1
        Aest = A_gdme; sest = s_gdme;
    else method == 2
        Aest = A_varnt; sest = s_varnt;
    end
    
    
    % Permute Results
    CRD = corrcoef([A Aest]);
    DD = abs(CRD(c+1:2*c,1:c));  
    perm_mtx = zeros(c,c);
    aux=zeros(c,1);
    for i=1:c
        [ld cd]=find(max(DD(:))==DD); 
        ld=ld(1);cd=cd(1); % in the case of more than one maximum
        perm_mtx(ld,cd)=1; 
        DD(:,cd)=aux; DD(ld,:)=aux';
    end
    Aest = Aest*perm_mtx;
    sest = sest'*perm_mtx;
    Sest = reshape(sest,[M,N,c]);
    sest = sest';
    
    % Postprocessing abundance maps
    %         if SNR < 20
    %             H = ones(3,3)/9;
    %             S_new = zeros(M,N,c);
    %             for i=1:4
    %                 S_new(:,:,i) = conv2(Sest(:,:,i), H, 'same');
    %             end
    %             Sest = Snew;
    %             sest = reshape(Sest, M*N,c);
    %             sest = sest./repmat(sum(sest,2), [1 4]);
    %         end
    
    % Visualize the estimation results
    if showflag,
        figure, 
        for i=1:c
            subplot(c,4,4*i-3),
            plot(A(:,i),'r'); axis([0 Band 0 1])
            if i==1 title('True end-members'); end
            subplot(c,4,4*i-2),
            plot(Aest(:,i),'g');axis([0 Band 0 1])
            if i==1 title('Estimated end-members'); end
            subplot(c,4,4*i-1),
            imagesc(reshape(abf(i,:),M,N));
            if i==1 title('True abundance'); end
            subplot(c,4,4*i),
            imagesc(Sest(:,:,i));
            if i==1 title('Estimated abundance'); end
        end
    end
    
    % Quantitative evaluation of spectral signature and abundance
    % Rmse error of abundances
    E_rmse(method) = sqrt(sum(sum(((abf-sest).*(abf-sest)).^2))/(M*N*c));
    
    % The angle between abundances (AAD)
    nabf = diag(abf*abf'); 
    nsest = diag(sest*sest');
    ang_beta = 180/pi*acos( diag(abf*sest')./sqrt(nabf.*nsest));
    E_aad(method) = mean(ang_beta.^2)^.5;
    
    % Cross entropy between abundance (AID)
    E_entropy = sum(abf.*log((abf+1e-9)./(sest+1e-9))) + sum(sest.*log((sest+1e-9)./(abf+1e-9)));
    E_aid(method) = mean(E_entropy.^2)^.5;
    
    % The angle between material signatures
    nA = diag(A'*A);
    nAest = diag(Aest'*Aest);
    ang_theta = 180/pi*acos( diag(A'*Aest)./sqrt(nA.*nAest) );
    E_sad(method) = mean(ang_theta.^2)^.5;
    
    % The spectral information divergence
    pA = A./(repmat(sum(A),[length(A(:,1)) 1]));
    qA = Aest./(repmat(sum(Aest),[length(A(:,1)) 1])); 
    qA = abs(qA); 
    SID = sum(pA.*log(pA./qA)) + sum(qA.*log(qA./pA));
    E_sid(method) = mean(SID.^2)^.5;
    
end

E_rmse
E_aad
E_aid
E_sad
E_sid
E_time