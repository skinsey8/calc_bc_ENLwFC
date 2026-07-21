function [ENLwFC,NULL_ENLwFC,NLwFC,LINwFC] = calc_bc_ENLwFC(data)
% This function can be used to calculate an explicitly nonlinear functional
% connectivity matrix for fMRI time series data. data: time*feature
% Cf. Kinsey, S. et al. (2026). Explicitly nonlinear fMRI networks reveal hidden trajectories of infant brain development. bioRxiv: the preprint server for biology, 2026.04.07.716703. https://doi.org/10.64898/2026.04.07.716703
% skinsey8@gsu.edu
% This code file is licensed under the MIT License (see LICENSE).

[nT,nF] = size(data);
data = zscore(data);

LINwFC = single(corr(data));
LINwFC(isnan(LINwFC)) = 0;
NLwFC = calc_bc_dcorr(data);
NLwFC = reshape(NLwFC,[nF*nF 1]);
LINwFC = reshape(LINwFC,[nF*nF 1]);
NULL_ENLwFC = sqrt(((LINwFC.*asin(LINwFC))+sqrt(1-LINwFC.^2)-LINwFC.*asin(LINwFC/2)-sqrt(4-LINwFC.^2)+1)/(1+pi/3-sqrt(3))); % Transformation based on Edelmann, D., Móri, T. F., & . Székely, G. J. (2021). On relationships between the Pearson and the distance correlation coefficients. Statistics & Probability Letters, 169. https://doi.org/10.1016/j.spl.2020.108960
ENLwFC = NLwFC - NULL_ENLwFC;

ENLwFC = reshape(ENLwFC,[nF nF]);
NULL_ENLwFC = reshape(NULL_ENLwFC,[nF nF]);
NLwFC = reshape(NLwFC,[nF nF]);
LINwFC = reshape(LINwFC,[nF nF]);

end