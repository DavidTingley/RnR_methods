function out=rebin(in,N_bins,dim)
% Copyright 2004 T.L. van Vuure
% This M-file allows the user to rebin a histogram or matrix or histograms.
% The matlab functions hist or histc take distributed data as input and
% divide it over bins. However they assume the data is distributed as
% floats. When they are integers (e.g. from a Multi Channel Analyzer
% spectrum) this leads to strange artifacts in the spectrum when hist or
% histc are used, except if the bins exactly match the bins the data were
% in already. 
% This M-file implements weighted distribution of the contents of a bin in
% the original spectrum over more than one target bins, or more than one
% original bin in one target bin. It can be considered linear
% interpolation, assuming that the original data are samples from a
% continuum (floats).
% N.B. While an elegant solution (the entire operation comes down to a single
% matrix multiplication) this is not optimized for speed. It is not
% recommended for more than a couple of hundred spectra at a time.
if nargin==0 || nargin>3
    error(['Rebins a matrix of histograms.\n' ...
    'matrix out = reb_gen(matrix in , number of bins , dimension containing the histograms).\n' ...
    'If the dimension containing the histograms is not specified, the last dimension is assumed.\n' ...
    'In the case of a single vector, both orientations are accepted.\n'],1);
end

dims=ndims(in);
if nargin == 2 % if dim is not indicated, take the last one
    if isvector(in) && size(in,2)==1 % is it a column vector?
        dim=1; % then the algorithm finds one histogram rather than n very short ones.
    else
        dim=dims;
    end
end
if dim>dims || round(dim)~=dim
    error('Invalid dimension specified.');
end
if dim==dims
    reorder_flag=false;     % histograms in last dimension
else
    in=shiftdim(in,dim);    % move the dimension containing the histograms to the end
                            % that way we can use the reshape function efficiently below.
    reorder_flag=true;
end

% N-Dimensional matrix is collapsed to a 2D one.
original_size=size(in);
working_size=[ prod(original_size(1:dims-1))   size(in,dims) ];
in = reshape(in,working_size);
% To simplify treating all histograms one by one, we rearrange them in a 2D
% matrix
target_size=[ original_size(1:dims-1) N_bins ]; % for the output matrix


if round(N_bins) ~= N_bins
    error('Non-integer number of histograms. Check your math.');
end
if N_bins < size(in,2)
    % This is the normal case, rebinning to fewer bins. This
    % will produce a matrix with weights that indicate which fraction of
    % the original bin should be transferred to the destination bin.
    binning_larger_flag=false;
    few_bins=N_bins;
    many_bins=size(in,2);
else
    % As a shortcut, the rebinning matrix for the opposite case, rebinning
    % to more bins, will be produced from the normal case by producing the
    % same matrix as above, and then transposing and normalizing it.
    % This is of course not a recommended procedure as it introduces even
    % more errors than when binning to fewer bins. However, it is the best
    % possible method without a priori knowledge of the contents of the
    % histograms.
    binning_larger_flag=true;
    few_bins=size(in,2);
    many_bins=N_bins;
end

% Assembly of the rebinning matrix with weights that indicate which fraction of
% the original bin should be transferred to the destination bin.
M=spalloc(few_bins,many_bins,few_bins+many_bins); % pre-allocate space
cnt2=1;
for cnt1 = 1:many_bins
    if cnt1 <= cnt2*many_bins/few_bins
        M(cnt2,cnt1) = 1;
    else
        M(cnt2,cnt1) = cnt2*many_bins/few_bins - (cnt1 - 1);
        M(cnt2+1,cnt1) = cnt1 - cnt2*many_bins/few_bins;
        cnt2 = cnt2 + 1;
    end
end
if binning_larger_flag
    norm_factor=sum(M(1,:),2);
    M=M'./norm_factor;
end
clear cnt1 cnt2 % Clean up counters for safe future use of variable names

in=double(in); % Making sure input histograms are in double float format
for cnt1=1:size(in,1)
    out( cnt1 , 1:N_bins ) = ( M * in(cnt1,:)' )'; 
end
% Linear algebra defines the matrix multiplication
% with column vectors, while we use row vectors, hence the double
% transposition above. Matlab is very quick at this.

out = reshape(out,target_size); % restoring the original number of dimensions
if reorder_flag
    out=shiftdim(out,dims-dim); % move the dimension containing the histograms
                              % back to its original position
end
out=full(out);
end