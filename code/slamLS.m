function [poses_ls, landmarks_ls] = slamLS(landmarks, pose0, transitions, observations, N_ITERATIONS)

    % initialize poses struct
    poses_ls = struct();

    % build state
    state = buildState(pose0, landmarks);
    state_dim = length(state);
    curr_pose_id = pose0.id;

    % array of relevant landmarks' ids
    land_ids = unique([landmarks.id]);

    for i=1:length(transitions)

        % transition to next state (first state have no observations)
        state = transition_function(state, transitions(i));
        ++curr_pose_id;

        % gather observations in current state
        curr_obs = searchObservationsById(observations, curr_pose_id);

        for j=1:N_ITERATIONS
            % reset H and b
            H = zeros(state_dim, state_dim);
            b = zeros(state_dim, 1);
            for k=1:length(curr_obs)
                % Is the observation regarding a relevant landmarks?
                if any(land_ids == curr_obs(k).id)
                    [e, J] = errorAndJacobian(state, curr_obs(k));
                    H += J.'*J;
                    b += J.'*e;
                endif
            endfor

            % update estimate of the current pose
            state -= H\b;
        endfor

        % append pose after LS optimization to output
        poses_ls(end+1) = pose(
            curr_pose_id,
            state(1),
            state(2),
            state(3)
        );
    endfor

    landmarks_ls = struct();
    % for in state - build landmarks_ls
    for i=1:length(land_ids)
        land_id = land_ids(i);
        land_x = state(4+2*land_id-2);
        land_y = state(4+2*land_id-1);
        landmarks_ls(end+1) = landmark(land_id, [land_x, land_y]);
    endfor

    % remove first, empty and meaningless, cell of struct
    poses_ls = poses_ls(2:end);
    landmarks_ls = landmarks_ls(2:end);

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
