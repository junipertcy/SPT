% histStat.m
%
% Calculates central moments of arbitrary order from histogram
% data, including matrices. Note the binning size will not greatly affect
% the calculated moments. For an example sample of size N=15000, choosing
% binning size of 100 is adequate.
% 
% Syntax:
% [allMeans allVars ...]=histStat(bincenters,counts)
% 
% INPUT:
%
% OUTPUT:
%
%
%
% Tzu-Chi Yen
% Apt. 29, 2014
%
% Modified from MATLAB File Exchange,
% code originally written by Kaare.


function varargout=histStat(bincenters, counts)
% Assumes each column in "counts" is a histogram. NaN-values are set to
% zero. The first output is the mean.


if numel(bincenters)~=size(counts,1)
    error('Error: each column in counts must be as long as bincenters.')
end

%removing NaNs:
counts(isnan(counts))=0;

bincenters=bincenters(:);

EXn=zeros(nargout,size(counts,2)); % non-central moments


for n=1:nargout
    % Calculating regular moment:
    EXn(n,:)=moments(bincenters,n);
    
    % Calculating central moment:
    varargout{n}=EXn(1,:).^(n); % k=0 term in binomial theorem
    if n>1 % not centering the mean
        for k=1:n
            varargout{n}=varargout{n}+nchoosek(n,k)*EXn(k,:).*EXn(1,:).^(n-k); % nchooseek = binomial coefficient
        end
    end
end

    function output=moments(bincenters,n)
        %returns E[x^n]
        output=sum(counts.*((bincenters(:).^n*ones(1,size(counts,2)))),1)./sum(counts,1);
        
    end
end