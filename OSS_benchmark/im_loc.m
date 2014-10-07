% im_loc.m
% 
% Immobile locator that finds the location of possibly blinking
% immobile reporters in a typical tracking experiment.
%
% input:
%           loc: locations of the center of the trajectories
%           cut1: cutoff for specifying distinct detections (unit: pixel)
%                 Example: say if two detected locations is less
%                 than cut1=1 pixels, then we calculate the mean
%                 location of the two points, representative of
%                 only one immobile particle.
%
% output:
%           tart: locations of the particles
%
% usual usage: 
%           tart = im_loc(loc,1);
%
% Tzu-Chi Yen
% Jul. 25, 2014

function [ tart ] = im_loc(loc, cut1)
 tart =[];
 [d1 d2] = kNearestNeighbors(loc,loc,2);
 loc(d1(d2(:,2)>cut1),:)=[];

 [d1 d2] = kNearestNeighbors(loc,loc,length(loc));
 criterion = max(diff(d2(1,:))>1);
 
 while criterion == 1
 tart = [tart; mean(loc(d1(1,1:find(diff(d2(1,:))>cut1,1)),:))];
 loc(d1(1,1:find(diff(d2(1,:))>1,1)),:)=[];
 [d1 d2] = kNearestNeighbors(loc,loc,length(loc));
 criterion = max(diff(d2(1,:))>1);
 end
 tart = [tart; mean(loc)];
end
