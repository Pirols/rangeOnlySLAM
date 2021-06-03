function [XR_ls, XL_ls, chi_stats_l, num_inliers_l, ...
          chi_stats_r, num_inliers_r] = slamLS(
        XR_ig, XL_ig,
        Zl, land_associations,
        Zr, pose_associations,
        num_iterations, damping, kernel_threshold)

    % retrieve cardinalities
    global POSE_DIM;
    global LAND_DIM;
    global POSE_NUM;
    global LAND_NUM;
    SYSTEM_SIZE = POSE_DIM*POSE_NUM + LAND_DIM*LAND_NUM;

    % initialize outputs
    XR_ls = XR_ig;
    XL_ls = XL_ig;
    chi_stats_l   = zeros(1, num_iterations);
    num_inliers_l = zeros(1, num_iterations);
    chi_stats_r   = zeros(1, num_iterations);
    num_inliers_r = zeros(1, num_iterations);

    % variables used to avoid repeating the same computations
    DAMPING_MATRIX = eye(SYSTEM_SIZE)*damping;

    for i=1:num_iterations
        % reset H and b (dampening H)
        H = DAMPING_MATRIX;
        b = zeros(SYSTEM_SIZE, 1);

        % pose-landmark measurement correction
        for j=1:size(Zl, 2)
            % retrieve measurement information
            pose_index = land_associations(1, j);
            land_index = land_associations(2, j);

            % relevant pose and landmark
            Xr = XR_ls(:, pose_index);
            Xl = XL_ls(:, land_index);
            z  = Zl(:, j);

            % compute error and jacobian(s)
            [e, Jr, Jl] = errorAndJacobianLand(Xr, Xl, z);

            % check if outlier
            chi = e.'*e;
            if chi > kernel_threshold
                e *= sqrt(kernel_threshold/chi);
                chi = kernel_threshold;
            else
                num_inliers_l(i)++;
            endif
            chi_stats_l(i) += chi;

            % Update H and b
            Hrl = Jr.'*Jl;

            pose_matrix_index = poseMatrixIndex(pose_index);
            land_matrix_index = landMatrixIndex(land_index);

            H(pose_matrix_index:pose_matrix_index+POSE_DIM-1,
              pose_matrix_index:pose_matrix_index+POSE_DIM-1) += Jr.'*Jr;

            H(pose_matrix_index:pose_matrix_index+POSE_DIM-1,
              land_matrix_index:land_matrix_index+LAND_DIM-1) += Hrl;

            H(land_matrix_index:land_matrix_index+LAND_DIM-1,
              pose_matrix_index:pose_matrix_index+POSE_DIM-1) += Hrl.';

            H(land_matrix_index:land_matrix_index+LAND_DIM-1,
              land_matrix_index:land_matrix_index+LAND_DIM-1) += Jl.'*Jl;

            b(pose_matrix_index:pose_matrix_index+POSE_DIM-1) += Jr.'*e;
            b(land_matrix_index:land_matrix_index+LAND_DIM-1) += Jl.'*e;
        endfor

        % pose-pose odometry correction
        for j=1:size(Zr, 2)
            % retrieve transition information
            pose_i_index = pose_associations(1, j);
            pose_j_index = pose_associations(2, j);

            % relevant pose and landmark
            Xi = XR_ls(:, pose_i_index);
            Xj = XR_ls(:, pose_j_index);
            z  = Zr(:, j);

            % compute error and jacobian(s)
            [e, Ji, Jj] = errorAndJacobianPose(Xi, Xj, z);

            % check if outlier
            chi = e.'*e;
            if chi > kernel_threshold
                e *= sqrt(kernel_threshold/chi);
                chi = kernel_threshold;
            else
                num_inliers_r(i)++;
            endif
            chi_stats_r(i) += chi;

            % Update H and b
            pose_i_matrix_index = poseMatrixIndex(pose_i_index);
            pose_j_matrix_index = poseMatrixIndex(pose_j_index);

            H(pose_i_matrix_index:pose_i_matrix_index+POSE_DIM-1,
              pose_i_matrix_index:pose_i_matrix_index+POSE_DIM-1) += Ji.'*Ji;

            H(pose_i_matrix_index:pose_i_matrix_index+POSE_DIM-1,
              pose_j_matrix_index:pose_j_matrix_index+POSE_DIM-1) += Ji.'*Jj;

            H(pose_j_matrix_index:pose_j_matrix_index+POSE_DIM-1,
              pose_i_matrix_index:pose_i_matrix_index+POSE_DIM-1) += Jj.'*Ji;

            H(pose_j_matrix_index:pose_j_matrix_index+POSE_DIM-1,
              pose_j_matrix_index:pose_j_matrix_index+POSE_DIM-1) += Jj.'*Jj;

            b(pose_i_matrix_index:pose_i_matrix_index+POSE_DIM-1) += Ji.'*e;
            b(pose_j_matrix_index:pose_j_matrix_index+POSE_DIM-1) += Jj.'*e;
        endfor

        % compute correction and update state
        [XR_ls, XL_ls] = boxPlus(XR_ls, XL_ls, -H\b);
    endfor
endfunction

function [e, Ji, Jj] = errorAndJacobianPose(Xi, Xj, z)
    global POSE_DIM;

    Xj_hat = v2t(trans_fun(Xi, z));
    Xi_mat = v2t(Xi);
    Xj_mat = v2t(Xj);

    % compute error
    e = inv(Xi_mat)*Xj_mat - inv(Xi_mat)*Xj_hat;
    e = reshape(e(1:2, :), 6, 1);

    % get values
    x_i     = Xi(1);
    y_i     = Xi(2);
    theta_i = Xi(3);
    x_j     = Xj(1);
    y_j     = Xj(2);
    theta_j = Xj(3);
    % compute trigonometric values
    c_i  = cos(theta_i);
    s_i  = sin(theta_i);
    c_j  = cos(theta_j);
    s_j  = sin(theta_j);
    c_ij = cos(theta_i-theta_j);
    s_ij = sin(theta_i-theta_j);

    %% compute jacobians
    % w.r.t. x_i, y_i, theta_i
    Ji = [
        0,    0,    -s_ij;
        0,    0,    -c_ij;
        0,    0,    c_ij;
        0,    0,    -s_ij;
        -c_i, -s_i, y_j*c_i - y_i*c_i + x_i*s_i - x_j*s_i;
        s_i,  -c_i, x_i*c_i - x_j*c_i + y_i*s_i - y_j*s_i;
    ];
    % w.r.t. x_j, y_j, theta_j
    Jj = [
        0,    0,   s_ij;
        0,    0,   c_ij;
        0,    0,   -c_ij;
        0,    0,   s_ij;
        c_i,  s_i, 0;
        -s_i, c_i, 0;
    ];
endfunction

function [e, Jr, Jl] = errorAndJacobianLand(Xr, Xl, z)
    global POSE_DIM;
    global LAND_DIM;

    % compute error
    sub = [Xr(1) - Xl(1), Xr(2) - Xl(2)];
    z_hat = norm(sub);
    e = z_hat - z;

    % compute jacobians
    Jr = zeros(length(z), POSE_DIM);
    Jl = zeros(length(z), LAND_DIM);
    Jr(1:2) = sub/z_hat;
    Jl(1:2) = -sub/z_hat;
endfunction

function [XR, XL] = boxPlus(XR, XL, dx)
    global POSE_DIM;
    global LAND_DIM;
    global POSE_NUM;
    global LAND_NUM;

    for i=1:POSE_NUM
        pose_matrix_index = poseMatrixIndex(i);
        XR(:, i) += dx(pose_matrix_index:pose_matrix_index+POSE_DIM-1);
    endfor

    for i=1:LAND_NUM
        land_matrix_index = landMatrixIndex(i);
        XL(:, i) += dx(land_matrix_index:land_matrix_index+LAND_DIM-1);
    endfor
endfunction

function pose_matrix_index = poseMatrixIndex(pose_index)
    global POSE_NUM;

    if pose_index > POSE_NUM
        pose_matrix_index = -1;
    else
        global POSE_DIM;
        pose_matrix_index = 1 + (pose_index-1)*POSE_DIM;
    endif
endfunction

function land_matrix_index = landMatrixIndex(land_index)
    global LAND_NUM;

    if land_index > LAND_NUM
        land_matrix_index = -1;
    else
        global POSE_DIM;
        global LAND_DIM;
        global POSE_NUM;
        land_matrix_index = 1 + POSE_NUM*POSE_DIM + (land_index-1)*LAND_DIM;
    endif
endfunction

function A = v2t(v)
    c = cos(v(3));
    s = sin(v(3));
    A = [c, -s, v(1);
         s,  c, v(2);
         0,  0, 1];
endfunction

function new_pose = trans_fun(pose, trans)
    v = trans(1);
    theta = pose(3);
    new_pose = reshape(
        [pose(1) + v*cos(theta);
         pose(2) + v*sin(theta);
         theta + trans(2)],
        size(pose));
endfunction
