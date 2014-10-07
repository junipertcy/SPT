% sort_trajectories_IAMS.m
%
% Original IAMS code from Hsieh's lab
%
%

function [trajectories, all_part] = sort_trajectories_IAMS(peaks)

    N_img = length(peaks);
    
    trajectories    = cell(1,1);
   % trajectories{1} = peaks{1}(:,1:2);

    n_curr_part     =  size(peaks{1},1);
    
    all_part   = [(1:n_curr_part)' ones(size(peaks{1},1),1) N_img*ones(size(peaks{1},1),1)];
                            % row-index: particle ID
                            % 1st column: current pos in "peaks"
                            % 2nd column: frame# that the particle appeared
                            % 3rd column: frame# that the particle left
    
        
    
    for i_frame = 1:N_img
    
        n_curr_part    =  size(peaks{i_frame},1);               % total number of particles in the frame
        last_part_IDs  =  find(all_part(:,1)~=-1);              % IDs of the ones that made it from the last frame
        n_new_part     =  n_curr_part - length(last_part_IDs);  % --> number of new particles
                        
        if n_new_part > 0
            for i_new_pos = 1:n_curr_part                       % check which position in "peaks" 
                if sum( find(all_part(:,1)==i_new_pos)) == 0    %    the new particle has
                    all_part = [all_part; i_new_pos i_frame 0]; % append new particle info at the end of "all_part"
                end
            end
        end
        
        curr_pos = all_part(:,1); % needed for next loop, since all_part(:,1) will be reassigned
        
        for i_part = 1:n_curr_part % in the order they are in "peaks"

            ID = find(curr_pos == i_part); % returns row-index, i.e. ID
            
            trajectories{ID}(i_frame+1-all_part(ID,2),1) = peaks{i_frame}(i_part,1); % x_position
            trajectories{ID}(i_frame+1-all_part(ID,2),2) = peaks{i_frame}(i_part,2); % y_position
            trajectories{ID}(i_frame+1-all_part(ID,2),3) = peaks{i_frame}(i_part,3); % empty so far (SD-x???)
            trajectories{ID}(i_frame+1-all_part(ID,2),4) = peaks{i_frame}(i_part,4); % empty so far (SD-y???)
            trajectories{ID}(i_frame+1-all_part(ID,2),5) = peaks{i_frame}(i_part,5); % contrast

            all_part(ID,1) = peaks{i_frame}(i_part,6); % position in the _NEXT_ frame

            if all_part(ID,1) == -1 % if this is the last frame for this particle
                all_part(ID,3) = i_frame; 
            end
        end

    end

    all_part = all_part(:,2:3);
    
return
