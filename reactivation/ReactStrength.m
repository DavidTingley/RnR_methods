function [R,phi] = ReactStrength(Qref,Qtar,varargin)

% Computes reactivation strength with various methods.
% 
% USAGE
%
% [R,phi] = reactStrength(Qref,Qtar,optionName,optionValue)
%
% INPUTS
%   Qref            binned spike train during the reference epoch (time x
%                   cells)
%   Qtar            binned spike train during the target epoch
%   <options>       optional list of property-value pairs (see table below)
%
% OUTPUTS
%   R               reactivation strength (time x significant components)
%   phi             encoding strength (normalized eigenvalues)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'method'     'pca' or 'ica' (default = 'pca')
%    =========================================================================

% Copyright (C) 2019 Adrien Peyrache & David Tingley.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

nCells = size(Qref,2);
if nCells ~= size(Qtar,2)
    error('''Qred'' and ''Qtar'' should have the same number of columns. (type ''help ReactStrength'' for details).') 
end

%Defaults
method =  'pca';

if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help ReactStrength'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i+3) ' is not a property (type ''help ReactStrength'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'method',
      method = varargin{i+1};
      if strcmp(method,'pca') && strcmp(method,'ica')
        error('Incorrect value for property ''method'' (type ''help ReactStrength'' for details).');
      end
  end
end

Cref = corrcoef(Qref);
[PCs,lambdas] = pcacov(Cref);

%Marcenko-Pastur threshold
lMax = (1 + sqrt(nCells / size(Qref,1))).^2;

nPCs    = sum(lambdas>lMax);
phi     = lambdas(1:nPCs)./lMax;

if strcmp(method,'ica')
    PCs         = fast_ica(zscore(Qref),nPCs);
end

scoreTar    = zscore(Qtar)*PCs(:,1:nPCs);

% This trick is used to get rid of the diagonal term in react. strength
tmp     = (zscore(Qtar).^2)*(PCs(:,1:nPCs).^2);
R       = scoreTar.^2 - tmp;
