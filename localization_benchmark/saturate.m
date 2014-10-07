% saturate.m
% 
% function called by detect_particles_IAMS_centroid.m
%====================================================================== 
%
% SATURATE: Saturates a vector of values to lower/upper bounds
%
% SYNTAX:  res = saturate(in, min, max)
%
% INPUTS:  in       a vector (or scalar) of numbers to be saturated
%          min      lower bound. 
%          max      upper bound
%
% The return value of a vector of the same size as in with all
% values >=min and <=max. Input values between min and max are
% left unchanged.
%
% Ivo Sbalzarini, 12.2.2003
% Institute of Computational Science, Swiss Federal
% Institute of Technology (ETH) Zurich. 
% E-mail: sbalzarini@inf.ethz.ch
%
%====================================================================== 

function res=saturate(in,min,max)

res = in;
res(find(in>max))=max;
res(find(in<min))=min;

return
