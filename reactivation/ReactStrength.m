function [R,phi] = ReactStrength(Qref,Qtar,varargin)

% Computes reactivation strength with various methods
% 
% USAGE
%
% reactStrength(Qref,Qtar,optionName,optionValue)
%
% INPUTS
%   Qref            binned spike train during the reference epoch (time x
%                   cells)
%   Qtar            binned spike train during the target epoch
%   <options>       optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'methods'     'pca' or 'ica' (default = 'pca')
%    =========================================================================

% Copyright (C) 2019 Adrien Peyrache & David Tingley.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

%Defaults
methods =  'pca';

if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help ReactStrength'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i+3) ' is not a property (type ''help ReactStrength'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'methods',
      duration = varargin{i+1};
      if strcmp(methods,'pca') && strcmp(methods,'ica')
        error('Incorrect value for property ''methods'' (type ''help ReactStrength'' for details).');
      end
  end
end

Cref = corrcoef(Qref);
[PCs,lambdas] = pcacov(Cref);

scorePRE = zscore(Qpre)*PCweights(:,1);
scorePOST = zscore(Qpost)*PCweights(:,1);

singleNpre = (zscore(Qpre).^2)*(PCweights(:,1).^2);
singleNpost = (zscore(Qpost).^2)*(PCweights(:,1).^2);

reactPRE = scorePRE.^2 - singleNpre;
reactPOST = scorePOST.^2 - singleNpost;
