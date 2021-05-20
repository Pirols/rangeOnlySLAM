function [poses_ls, landmarks_ls, chi_stats, num_inliers] = slamLS(landmarks_ig,
                poses_ig, observations, n_iterations, damping, kernel_threshold)

    % initialize outputs
    poses_ls = struct();
    landmarks_ls = struct();
    chi_stats = zeros(1, n_iterations);
    num_inliers = zeros(1, n_iterations);

    % no poses provided
    n_poses = length(poses_ig);
    if !n_poses
        return
    endif

    % build initial state
    pose0 = poses_ig(1);
    curr_pose_id = pose0.id;
    state = buildState(pose0, landmarks_ig);
    state_dim = length(state);

    % array of relevant landmarks' ids
    land_ids = unique([landmarks_ig.id]);

    for i=2:n_poses

        % update state
        curr_pose = poses_ig(i);
        curr_pose_id = curr_pose.id;
        state(1:3) = [curr_pose.x; curr_pose.y; curr_pose.theta];

        % gather observations in current state
        curr_pose_obs = searchObservationsById(observations, curr_pose_id);

        for iter=1:n_iterations
            % reset H and b
            H = zeros(state_dim, state_dim);
            b = zeros(state_dim, 1);
            for k=1:length(curr_pose_obs)
                % Is the observation regarding a relevant landmarks?
                if any(land_ids == curr_pose_obs(k).id)
                    [e, J] = errorAndJacobian(state, curr_pose_obs(k));
                    chi = e.'*e;
                    if chi > kernel_threshold
                        e *= sqrt(kernel_threshold/chi);
                        chi = kernel_threshold;
                    else
                        num_inliers(iter)++;
                    endif
                    chi_stats(iter) += chi;
                    H += J.'*J;
                    b += J.'*e;
                endif
            endfor

            % apply damping
            H += eye(state_dim)*damping;

            % update estimate of the current pose
            state -= H\b;

            % normalize pose angle
            %theta = state(3);
            %state(3) = atan2(sin(theta), cos(theta));
        endfor

        % append pose after LS optimization to output
        poses_ls(end+1) = pose(
            curr_pose_id,
            state(1),
            state(2),
            state(3)
        );
    endfor

    % for in state - fill landmarks_ls
    for i=1:length(land_ids)
        land_id = land_ids(i);
        land_x = state(4+2*land_id-2);
        land_y = state(4+2*land_id-1);
        landmarks_ls(end+1) = landmark(land_id, [land_x, land_y]);
    endfor

    % remove first, empty field of struct
    poses_ls(1) = [];
    landmarks_ls(1) = [];

end

function [e, J] = errorAndJacobian(state, obs)
    % find landmark
    curr_land_id = obs.id;
    curr_land_idx_x = (4+2*curr_land_id-2);
    curr_land_idx_y = (4+2*curr_land_id-1);
    curr_land_pose = [
        state(curr_land_idx_x);
        state(curr_land_idx_y)
    ];

    % compute error and build e
    z_hat = measurement_function(state(1:2), curr_land_pose);
    e = z_hat - obs.range;

    % compute jacobian
    J = zeros(length(e), length(state));
    J(1:2) = [(state(1)-curr_land_pose(1))/z_hat (state(2)-curr_land_pose(2))/z_hat];
    J(curr_land_idx_x) = (curr_land_pose(1)-state(1))/z_hat;
    J(curr_land_idx_y) = (curr_land_pose(2)-state(2))/z_hat;
end

function new_state = transition_function(state, transition)
    % Gather pose and transition data
    curr_theta = state(3);
    v = transition.v(1);
    w = transition.v(3);

    % landmarks stay fixed
    new_state = state;

    % Compute new state
    new_state(1:3) += [
        v*cos(curr_theta);
        v*sin(curr_theta);
        w
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

function [state] = buildState(pose0, landmarks)
    state = [pose0.x; pose0.y; pose0.theta];

    for i=1:length(landmarks)
        curr_land_id = landmarks(i).id;
        state(4+2*curr_land_id-2) = landmarks(i).x_pose;
        state(4+2*curr_land_id-1) = landmarks(i).y_pose;
    endfor
end
