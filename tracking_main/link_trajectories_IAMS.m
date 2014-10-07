% link_trajectories_IAMS.m
%
% Original IAMS code from Hsieh's lab
%
%


%====================================================================== 
%
% LINK_TRAJECTORIES: links particle positions into trajectories
%
% SYNTAX:  peaks = link_trajectories(peaks, L, viz, nfig)
%
% INPUTS:  peaks     peaks cell list of particle positions in all
%                    frames and zero and second order intensity
%                    moments as produced by detect_particles
%          L         maximum allowed particle displacement between
%                    two subsequent frames (cutoff)
%          viz       =1 if intermediate visualization is needed
%          nfig      figure number for first viz image
%
% The particle matching is done such that the sum of squared
% displacements between two frames is minimized and also the
% quadratic differences in the zero and second order intensity
% moments of the particles are minimized (more precise: between
% the particle properties peaks{i}(:,3) and peaks{i}(:,4)).
%
% After linkage, the same cell list "peaks" is returned but with
% the fields peaks{i}(:,6) now filled. All other fields are left
% untouched. peaks{i}(j,6) contains the index of the one
% particle in frame i+1 corresponding to particle j in frame i.
% If no correspondence was found (i.e. the trajectory terminates)
% -1 is given.
%
% Ivo Sbalzarini, 12.2.2003
% Institute of Computational Science, Swiss Federal
% Institute of Technology (ETH) Zurich. 
% E-mail: sbalzarini@inf.ethz.ch
%
% Based on the logistic transportation algorithm described in:
%     http://www.damtp.cam.ac.uk/user/fdl/people/sd/digimage/...
%        document/track2d/matching.htm#L_2_PARTICLE_MATCHING
%====================================================================== 

%====================================================================== 
% STEP 5: Linking locations into trajectories
%====================================================================== 

function peaks = link_trajectories_IAMS(peaks, L, viz, nfig)

nframe = length(peaks);

for iframe = 2:nframe,
    % initialize all linked lists to -1
    peaks{iframe-1}(:,6) = -1;
    disp(sprintf('Linking paths in frame %d of %d',iframe,nframe))
    m = length(peaks{iframe-1}(:,1)); % previous frame: p_i, i=1,..,m
    n = length(peaks{iframe}(:,1));   % this frame: q_j, j=1,..,n
    % empty association matrix
    A = zeros(m+1,n+1);               % a_{ij}=1 <=> p_i = q_j 
                                      % dummy particles a_{(m+1)j},
				      % a_{i(n+1)}
    % delta(i,j): quadratic distance between p_i and q_j
    xrep = repmat(peaks{iframe}(:,1)',m,1);
    yrep = repmat(peaks{iframe}(:,2)',m,1);
    xdiff = xrep - repmat(peaks{iframe-1}(:,1),1,n);
    ydiff = yrep - repmat(peaks{iframe-1}(:,2),1,n);
    delta = (xdiff.^2)+(ydiff.^2);

    % --- Cost for peak intensity and shape closed --- 
%     % dm0(i,j): quadratic difference between m0 moments of p_i and q_j
%     xrep = repmat(peaks{iframe}(:,3)',m,1);
%     xdiff = xrep - repmat(peaks{iframe-1}(:,3),1,n);
%     dm0 = xdiff.^2;
% 
%     % dm2(i,j): quadratic difference between m2 moments of p_i and q_j
%     xrep = repmat(peaks{iframe}(:,4)',m,1);
%     xdiff = xrep - repmat(peaks{iframe-1}(:,4),1,n);
%     dm2 = xdiff.^2;

    % C(i,j): cost function for link p_i, q_j
    C = L*L*ones(m+1,n+1);   % set broken dummy links to L^2
%    C(1:m,1:n) = (delta+dm0+dm2);
    C(1:m,1:n) = delta;

    % set cost of matchings that will never occur to Inf
    C1 = C(1:m,1:n) - repmat(C(m+1,1:n),m,1);
    C2 = C(1:m,1:n) - repmat(C(1:m,n+1),1,n);
    set1 = find(C1 > 0);
    set2 = find(C2 > 0);
    s = union(set1,set2);
    [i,j] = ind2sub([m n],s);
    C(sub2ind(size(C),i,j)) = Inf;

    % initialize link matrix A
    for i=1:m,
	% sort costs of real particles
	[srtcst,srtidx] = sort(C(i,:));
	% append index of dummy particle
	iidx = 1;
	dumidx = find(srtidx==(n+1));
	% search for available particle of smallest cost or dummy
	while and(sum(A(:,srtidx(iidx)))~=0, iidx<dumidx),
	    iidx = iidx + 1;
	end;
	A(i,srtidx(iidx)) = 1;
    end;
    % set dummy particle for columns with no entry
    s = sum(A,1);
    A(m+1,find(s < 1)) = 1;
    % dummy always corresponds to dummy
    A(m+1,n+1) = 1;

    % consistency check for matrix A
    s = sum(A(:,1:n),1);
    if find(s(1:n)~=1),
	disp('Inconsistent initial matrix A. Columns: ');
        find(s(1:n)~=1)
    end;
    s = sum(A(1:m,:),2);
    if find(s(1:m)~=1),
	disp('Inconsistent initial matrix A. Rows:');
	find(s(1:m)~=1)
    end;

    if viz,
        figure(nfig)
	clf
	hist(reshape(delta+dm0+dm2,m*n,1),100)
	xlabel('inter-particle cost distance')
    end;

    % iteration loop for the logistic transportation algorithm
    finished = 0;
    mincost = [];
    while ~finished,
	% non-set links of finite cost
	todo = intersect(find(A(1:m,1:n)==0),find(C(1:m,1:n)<Inf));

	% determine induced changes and reduced cost Cred for each
	% candidate link insertion
	[Icand,Jcand] = ind2sub([m n],todo);
	Cred = zeros(size(Icand));
	Xcand = zeros(size(Icand));
	Ycand = zeros(size(Icand));
	for ic=1:length(Icand),
	    Cred(ic) = C(Icand(ic),Jcand(ic));
	    Xcand(ic) = find(A(Icand(ic),:)==1);
	    Ycand(ic) = find(A(:,Jcand(ic))==1);
	    Cred(ic) = Cred(ic)-C(Icand(ic),Xcand(ic))-C(Ycand(ic),Jcand(ic));
	    Cred(ic) = Cred(ic)+C(Ycand(ic),Xcand(ic));
	end;

	% find minimum cost and corresponding action
	[minc,mini] = min(Cred);
	mincost = [mincost, minc];

	% if minimum is < 0, link addition is favorable
        if minc < 0,
	    % add link and update dependencies to preserve topology
	    A(Icand(mini),Jcand(mini)) = 1;
	    A(Ycand(mini),Jcand(mini)) = 0;
	    A(Icand(mini),Xcand(mini)) = 0;
	    A(Ycand(mini),Xcand(mini)) = 1;
	else
	    % done if best change is no more an improvement
	    finished = 1;
	end;

	% consistency check for matrix A
	s = sum(A(:,1:n),1);
	if find(s(1:n)~=1),
	    disp('Inconsistent matrix A during optimization. Columns: ');
	    find(s(1:n)~=1)
	end;
	s = sum(A(1:m,:),2);
	if find(s(1:m)~=1),
	    disp('Inconsistent matrix A during optimization. Rows:');
	    find(s(1:m)~=1)
	end;
    end;

    if viz,
	if (length(mincost) > 0),
	    figure(nfig+1)
	    clf
	    plot(mincost)
	    xlabel('iteration')
	    ylabel('reduced cost for best link insertion')
	end;
    end;

    % Convert link matrix to linked list representation
    [Ilink,Jlink] = find(A(1:m,:));
    % if link is to dummy particle, set index to -1
    Jlink(find(Jlink==(n+1))) = -1;
    % set linked list indices
    peaks{iframe-1}(Ilink,6) = Jlink;
end; 

% terminate all linked lists at the very end
peaks{nframe}(:,6) = -1;

return;
