function [XR_ls, XL_ls, chi_stats, num_inliers] = slamLS(XR_ig, XL_ig, Z,
                 associations, num_iterations, damping, kernel_threshold)

    % retrieve cardinalities
    global POSE_DIM;
    global LAND_DIM;
    global POSE_NUM;
    global LAND_NUM;
    SYSTEM_SIZE = POSE_DIM*POSE_NUM + LAND_DIM*LAND_NUM;

    % initialize outputs
    chi_stats = zeros(1, num_iterations);
    num_inliers = zeros(1, num_iterations);
    XR_ls = XR_ig;
    XL_ls = XL_ig;

    % convenient variable
    DAMPING_MATRIX = eye(SYSTEM_SIZE)*damping;

    for i=1:num_iterations
        % reset H and b (dampening H)
        H = DAMPING_MATRIX;
        b = zeros(SYSTEM_SIZE, 1);
        for j=1:size(Z, 2)
            % retrieve measurement information
            pose_index = associations(1, j);
            land_index = associations(2, j);
            z = Z(:, j);

            % relevant pose and landmark
            Xr = XR_ls(:, :, pose_index);
            Xl = XL_ls(:, land_index);

            % compute error and jacobian(s)
            [e, Jr, Jl] = errorAndJacobian(Xr, Xl, z);

            % check if outlier
            chi = e.'*e;
            if chi > kernel_threshold
                e *= sqrt(kernel_threshold/chi);
                chi = kernel_threshold;
            else
                num_inliers(i)++;
            endif
            chi_stats(i) += chi;

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

        % compute correction and update state
        [XR_ls, XL_ls] = boxPlus(XR_ls, XL_ls, -H\b);
    endfor
endfunction

function [e, Jr, Jl] = errorAndJacobian(Xr, Xl, z)

    global POSE_DIM;
    global LAND_DIM;

    % compute error
    sub = [Xr(1, 3) - Xl(1), Xr(2, 3) - Xl(2)];
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
        dxr = dx(pose_matrix_index:pose_matrix_index+POSE_DIM-1);
        XR(:, :, i) = v2t(dxr)*XR(:, :, i);
    endfor

    for i=1:LAND_NUM
        land_matrix_index = landMatrixIndex(i);
        XL(:, i) += dx(land_matrix_index:land_matrix_index+LAND_DIM-1, :);
    endfor

endfunction

function A = v2t(v)
    c = cos(v(3));
    s = sin(v(3));
    A = [c, -s, v(1);
         s,  c, v(2);
         0,  0,   1];
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
