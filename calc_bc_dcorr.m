function [Y] = calc_bc_dcorr(X)
% This function can be used to calculate a bias-corrected distance correlation matrix for 
% fMRI time series data. X: time*feature
% Cf. Kinsey, S. et al. (2026). Explicitly nonlinear fMRI networks reveal hidden trajectories of infant brain development. bioRxiv: the preprint server for biology, 2026.04.07.716703. https://doi.org/10.64898/2026.04.07.716703
% skinsey8@gsu.edu
% This code file is licensed under the MIT License (see LICENSE).
% Implementation based on Székely, G., & Rizzo, M. L. (2014). Partial distance correlation with methods for dissimilarities. The Annals of Statistics, 42(6), 2382-2412. https://doi.org/10.1214/009053607000000505

[nT,nF] = size(X);
Y = zeros(nT*nT,nF);

for ii = 1:nF
    x = X(:,ii);
    a = pdist2(x,x);
    A = a - (sum(a)/(nT-2)) - (sum(a)/(nT-2))' + (sum(a(:))/((nT-1)*(nT-2)));
    A(1:nT+1:end) = 0;
    Y(:,ii) = A(:); % Y = U-centered, modified distance matrices. Y: nT^2*nF
end
clear X A

Y = single(Y);
I = magic(nT);
I = eye(size(I));
diagidx = logical(I);

Y = (Y(~diagidx,:)'*Y(~diagidx,:))/(nT*(nT-3)); % Y = Unbiased squared distance covariance
u = sqrt(diag(Y));
u = u*u';

Y = (Y./u);
Y(isnan(Y)) = 0; % Y = Squared bias-corrected distance correlation (bcdCor)
clear u;
Y(Y<0) = 0;
Y = sqrt(Y); % Y = Bias-corrected distance correlation (truncated)

end