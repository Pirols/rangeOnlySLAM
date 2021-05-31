close all
clear
clc

addpath 'g2o_wrapper'

%%% load datasets
disp("----- Loading datasets and initializing landmarks...")
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
disp("----- Parsing loaded data...");
% Initialize XR and XL
XR_ig = zeros(POSE_DIM, POSE_NUM);
XR_gt = zeros(POSE_DIM, POSE_NUM);
XL_ig = zeros(LAND_DIM, LAND_NUM);
XL_gt = zeros(LAND_DIM, LAND_NUM);

% XR
disp("---------- XR");
for i=1:POSE_NUM
    %%% initial guess
    XR_ig(:, i) = [poses_ig(i).x; poses_ig(i).y; poses_ig(i).theta];

    %%% ground truth
    XR_gt(:, i) = [poses_gt(i).x; poses_gt(i).y; poses_gt(i).theta];
endfor

% XL
disp("---------- XL");
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
disp("---------- Z and associations");
Z = [];
associations = [];

% check if there is noise in the observations
mean_err_obs = 0;

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

        % add current observation's noise
        land = searchById(landmarks_gt, observations(i).observation(j).id);
        pose = searchById(poses_gt, observations(i).pose_id);
        val = norm([pose.x - land.x_pose; pose.y - land.y_pose]);
        mean_err_obs += abs(Z(:, end) - val);
    endfor
endfor

% get mean observations' noise
mean_err_obs /= size(Z, 2);
printf("--------------- Mean error in the observations is %f\n", mean_err_obs);

%%% least squares optimization
disp("----- Starting LS optimization...");
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

%%% Plot results
disp("----- Plotting results...");
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
plot(XR_ig(1, :), XR_ig(2, :), 'r-', 'linewidth', 3);
plot(XR_ls(1, :), XR_ls(2, :), 'b-', 'linewidth', 3);
plot(XR_gt(1, :), XR_gt(2, :), 'g-', 'linewidth', 3);
legend("Initial guess", "Optimized", "Ground truth");

% Chi evolution
figure();
title("Chi evolution");
hold on;
plot(chi_stats, 'b-', 'linewidth', 2)
