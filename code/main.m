close all
clear
clc

addpath 'g2o_wrapper'
addpath 'visualization'
source 'mse.m'

% disable singularity warnings (couldn't find a way to disable just that kind of warnings)
warning("off", "all")

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
% LS
LS_ITERATIONS = 10;
DAMPING = 0.01;
KERNEL_THRESHOLD = 1.0;

%%% least squares optimization
[poses_ls, landmarks_ls, chi_stats, num_inliers] = slamLS(
    landmarks_ig,
    poses_ig,
    observations,
    LS_ITERATIONS,
    DAMPING,
    KERNEL_THRESHOLD
);

%%% printing some evaluation metrics
mse_lands = mse_land(landmarks_gt, landmarks_ig);
mse_lands_ls = mse_land(landmarks_gt, landmarks_ls);
printf("MSE over landmarks went from %d to %d (after LS optimization)\n",
        mse_lands, mse_lands_ls);
mse_poses_xy = mse_pose(poses_gt, poses_ig);
mse_poses_xy_ls = mse_pose(poses_gt, poses_ls);
printf("MSE over poses(only xy) went from %d to %d (after LS optimization)\n",
        mse_poses_xy, mse_poses_xy_ls);

%%% Build following plots
%       * Trajectory (both ground truth and guess)
%       * Landmarks map (both ground truth and guess)
%       * Chi evolution
