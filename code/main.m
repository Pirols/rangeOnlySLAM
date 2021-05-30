close all
clear
clc

addpath 'g2o_wrapper'
source 'mse.m'

%%% load datasets
%% initial guess
[~, poses_ig, ~, observations] = loadG2o('../datasets/slam2d_range_only_initial_guess.g2o');
% remove first, empty field of structs
poses_ig(1) = [];
observations(1) = [];
% trilaterate landmarks
landmarks_ig = initLandmarks(poses_ig, observations);
%% ground truth
[landmarks_gt, poses_gt, ~, ~] = loadG2o('../datasets/slam2d_range_only_ground_truth.g2o');
% remove first, empty field of structs
landmarks_gt(1) = [];
poses_gt(1) = [];

%%% parameters
% cardinalities
global POSE_DIM = 3;
global LAND_DIM = 2;
global POSE_NUM = length(poses_ig);
global LAND_NUM = length(landmarks_ig);

%%% build XR, XL, Z and associations to perform ICP
% Initialize XR and XL
XR_ig = zeros(3, 3, POSE_NUM);
XR_gt = zeros(3, 3, POSE_NUM);
XL_ig = zeros(LAND_DIM, LAND_NUM);
XL_gt = zeros(LAND_DIM, LAND_NUM);

% XR
for i=1:POSE_NUM
    %%% initial guess
    theta = poses_ig(i).theta;
    c = cos(theta);
    s = sin(theta);
    XR_ig(1, 1, i) = c;
    XR_ig(1, 2, i) = -s;
    XR_ig(1, 3, i) = poses_ig(i).x;
    XR_ig(2, 1, i) = s;
    XR_ig(2, 2, i) = c;
    XR_ig(2, 3, i) = poses_ig(i).y;
    XR_ig(3, 3, i) = 1;

    %%% ground truth
    theta = poses_gt(i).theta;
    c = cos(theta);
    s = sin(theta);
    XR_gt(1, 1, i) = c;
    XR_gt(1, 2, i) = -s;
    XR_gt(1, 3, i) = poses_gt(i).x;
    XR_gt(2, 1, i) = s;
    XR_gt(2, 2, i) = c;
    XR_gt(2, 3, i) = poses_gt(i).y;
    XR_gt(3, 3, i) = 1;
endfor

% XL
for i=1:LAND_NUM
    XL_ig(1:2, i) = [landmarks_ig(i).x_pose; landmarks_ig(i).y_pose];
    land_id = landmarks_ig(i).id;
    land_gt = searchById(landmarks_gt, land_id);
    if isnumeric(land_gt)
        disp("Error: landmark found in landmarks_ig but not in landmarks_gt");
        return
    endif
    XL_gt(1:2, i) = [land_gt.x_pose; land_gt.y_pose];
endfor

% Z and associations
Z = [];
associations = [];
for i=1:length(observations)
    pose_num = find([poses_ig.id] == observations(i).pose_id);
    for j=1:length(observations(i).observation)
        land_num = find([landmarks_ig.id] == observations(i).observation(j).id);
        if !length(land_num)
            % observation not relative to a relevant landmark
            continue
        endif
        Z(:, end+1) = observations(i).observation(j).range;
        associations(:, end+1) = [pose_num; land_num];
    endfor
endfor

%%% least squares optimization
% LS parameters
LS_ITERATIONS = 100;
DAMPING = 0.01;
KERNEL_THRESHOLD = 1.0;

% LS
[XR_ls, XL_ls, chi_stats, num_inliers] = slamLS(
    XR_ig,
    XL_ig,
    Z,
    associations,
    LS_ITERATIONS,
    DAMPING,
    KERNEL_THRESHOLD
);

%%% Evaluation
% Map
mse_lands_ig = mse_land(XL_ig, landmarks_ig, landmarks_gt);
mse_lands_ls = mse_land(XL_ls, landmarks_ig, landmarks_gt);
printf("MSE over landmarks went from %d to %d (after LS optimization)\n",
        mse_lands_ig, mse_lands_ls);
% Trajectory (only XY)
mse_poses_xy_ig = mse_pose(XR_ig, poses_ig, poses_gt);
mse_poses_xy_ls = mse_pose(XR_ls, poses_ig, poses_gt);
printf("MSE over poses (only XY) went from %d to %d (after LS optimization)\n",
        mse_poses_xy_ig, mse_poses_xy_ls);

%%% Plot results
% Map
figure();
title("Map");
hold on;
plot(XL_ig(1, :), XL_ig(2, :), 'r*', "linewidth", 2);
plot(XL_ls(1, :), XL_ls(2, :), 'bx', "linewidth", 2);
plot(XL_gt(1, :), XL_gt(2, :), 'go', "linewidth", 2);
legend("Initial guess", "Optimized", "Ground truth");

% Trajectory
figure();
title("Trajectory");
hold on;
plot(reshape(XR_ig(1, 3, :), 1, POSE_NUM), reshape(XR_ig(2, 3, :), 1, POSE_NUM),
     'r-', 'linewidth', 3);
plot(reshape(XR_ls(1, 3, :), 1, POSE_NUM), reshape(XR_ls(2, 3, :), 1, POSE_NUM),
     'b-', 'linewidth', 3);
plot(reshape(XR_gt(1, 3, :), 1, POSE_NUM), reshape(XR_gt(2, 3, :), 1, POSE_NUM),
     'g-', 'linewidth', 3);
legend("Initial guess", "Optimized", "Ground truth");

% Chi evolution
figure();
title("Chi evolution");
hold on;
plot(chi_stats, 'b-', 'linewidth', 2)
