close all
clear
clc

addpath 'g2o_wrapper'

%%% load datasets
disp("----- Loading datasets and initializing landmarks...")
%% initial guess
[~, poses_ig, transitions, observations] = loadG2o('../datasets/slam2d_range_only_initial_guess.g2o');
% remove first, empty field of structs
poses_ig(1)     = [];
observations(1) = [];
transitions(1)  = [];
% trilaterate landmarks
landmarks_ig = initLandmarks(poses_ig, observations);
%% ground truth
[landmarks_gt, poses_gt, ~, ~] = loadG2o('../datasets/slam2d_range_only_ground_truth.g2o');
% remove first, empty field of structs
landmarks_gt(1) = [];
poses_gt(1)     = [];

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

% Zl and land_associations
disp("---------- Zl and land_associations");
Zl = [];
land_associations = [];

% needed to check the noise in the observations
land_obs_mean_err = 0;

for i=1:length(observations)
    pose_num = find([poses_ig.id] == observations(i).pose_id);
    for j=1:length(observations(i).observation)
        land_num = find([landmarks_ig.id] == observations(i).observation(j).id);
        if !length(land_num)
            % observation not relative to a relevant landmark
            continue
        endif
        Zl(:, end+1) = observations(i).observation(j).range;
        land_associations(:, end+1) = [pose_num; land_num];

        % add current observation's noise
        land_gt = searchById(landmarks_gt, observations(i).observation(j).id);
        pose_gt = searchById(poses_gt, observations(i).pose_id);
        land_obs_mean_err += abs(Zl(:, end) - norm([
            pose_gt.x - land_gt.x_pose;
            pose_gt.y - land_gt.y_pose
        ]));
    endfor
endfor
% get mean observations' noise
land_obs_mean_err /= size(Zl, 2);
printf("--------------- Mean error in the range measurements is %f\n", land_obs_mean_err);

% Zl and land_associations
disp("---------- Zr and pose_associations");
Zr = [];
pose_associations = [];

for i=1:length(transitions)
    pose_num_i = find([poses_ig.id] == transitions(i).id_from);
    pose_num_j = find([poses_ig.id] == transitions(i).id_to);

    pose_associations(:, end+1) = [pose_num_i; pose_num_j];
    Zr(:, end+1) = [transitions(i).v(1); transitions(i).v(3)];
endfor

%%% least squares optimization
disp("----- Starting LS optimization...");
% LS parameters
LS_ITERATIONS    = 100;
DAMPING          = 0.01;
KERNEL_THRESHOLD = 1.0;

% LS
[XR_ls, XL_ls, chi_stats_l, num_inliers_l, chi_stats_r, num_inliers_r] = slamLS(
    XR_ig,
    XL_ig,
    Zl, land_associations,
    Zr, pose_associations,
    LS_ITERATIONS,
    DAMPING,
    KERNEL_THRESHOLD
);

%%% Plot results
disp("----- Plotting results...");

% Map
disp("---------- Map");
figure();
title("Map");
hold on;
plot(XL_ig(1, :), XL_ig(2, :), 'r*', "linewidth", 2);
plot(XL_ls(1, :), XL_ls(2, :), 'bx', "linewidth", 2);
plot(XL_gt(1, :), XL_gt(2, :), 'go', "linewidth", 2);
legend("Initial guess", "Optimized", "Ground truth");

% Trajectory (XY)
disp("---------- Trajectory (XY)");
figure();
title("Trajectory (XY)");
hold on;
plot(XR_ig(1, :), XR_ig(2, :), 'r-', 'linewidth', 3);
plot(XR_ls(1, :), XR_ls(2, :), 'b-', 'linewidth', 3);
plot(XR_gt(1, :), XR_gt(2, :), 'g-', 'linewidth', 3);
legend("Initial guess", "Optimized", "Ground truth");

% % Uncomment paragraph for theta plot
% % Trajectory (theta)
% disp("---------- Trajectory (theta)");
% figure();
% title("Trajectory (theta)");
% hold on;
% plot(XR_ig(3, :), 'r-', 'linewidth', 3);
% plot(XR_ls(3, :), 'b-', 'linewidth', 3);
% plot(XR_gt(3, :), 'g-', 'linewidth', 3);
% legend("Initial guess", "Optimized", "Ground truth");

% Chi evolution (measurements)
disp("---------- Chi evolution (measurements)");
figure();
title("Chi evolution (measurements)");
hold on;
plot(chi_stats_l, 'b-', 'linewidth', 2)

% Chi evolution (odometry)
disp("---------- Chi evolution (odometry)");
figure();
title("Chi evolution (odometry)");
hold on;
plot(chi_stats_r, 'b-', 'linewidth', 2)
