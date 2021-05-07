function [poses_ls, landmarks_ls] = slamLS(landmarks, pose0, transitions, observations, N_ITERATIONS)

    % initialize outputs
    poses_ls = struct();
    landmarks_ls = landmarks;

    % initialize first pose
    curr_pose = [pose0.x; pose0.y; pose0.theta];
    curr_pose_id = pose0.id;

    % H and b dimensions parameters
    pose_dim = length(curr_pose);

    for i=1:length(transitions)
        % reset H and b
        H = zeros(3, 3);
        b = zeros(3, 1);

        % observations in current pose
        curr_obs = searchObservationsById(observations, curr_pose_id);

        for j=1:N_ITERATIONS
            for k=1:length(curr_obs)
                curr_land = searchById(landmarks, curr_obs(k).id);
                if !isnumeric(curr_land)  % if it is a number such landmark is not available
                    land_pose = [curr_land.x_pose; curr_land.y_pose];
                    [e, J] = errorAndJacobian(curr_pose, curr_obs(k).range, land_pose);
                    H += J.'*J;
                    b += J.'*e;
                endif
            endfor

            % update estimate of the current pose
            curr_pose -= H\b;
        endfor

        % append pose after LS optimization to output
        poses_ls(end+1) = pose(
            curr_pose_id,
            curr_pose(1),
            curr_pose(2),
            curr_pose(3)
        );

        curr_pose = transition_function(curr_pose, transitions(i));
        curr_pose_id++;
    endfor

    % remove first, empty and meaningless, cell of struct
    poses_ls = poses_ls(2:end);

end

function [e, J] = errorAndJacobian(curr_pose, range_obs, land_pose)
    z_hat = measurement_function(curr_pose, land_pose);
    e = range_obs - z_hat;
    J = [curr_pose(1)-land_pose(1)/z_hat curr_pose(2)-land_pose(2)/z_hat 0];
end

function new_pose = transition_function(curr_pose, transition)
    % Gather pose and transition data
    curr_x = curr_pose(1);
    curr_y = curr_pose(2);
    curr_theta = curr_pose(3);
    v = transition.v(1);
    w = transition.v(3);

    % Compute new pose
    new_pose = [
        curr_x + v*cos(curr_theta);
        curr_y + v*sin(curr_theta);
        curr_theta + w
    ];
end

function z_hat = measurement_function(curr_pose, land_pose)
    z_hat = norm([curr_pose(1:2)-land_pose]);
end

function out = searchObservationsById(observations, pose_id)
    index = find([observations.pose_id] == pose_id);
    if (index)
        out = observations(index).observation;
    else
        out = [];
    end
end
