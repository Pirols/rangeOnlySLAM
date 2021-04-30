close all
clear
clc

addpath 'g2o_wrapper'
addpath 'visualization'
% source "utilities/geometry_helpers_2d.m"

%%% load datasets
% initial guess
[~, poses, transitions, observations] = loadG2o('../datasets/slam2d_range_only_initial_guess.g2o');
% ground truth
[landmarks_gt, poses_gt, transitions_gt, observations_gt] = loadG2o('../datasets/slam2d_range_only_ground_truth.g2o');

%%% remove first elements of structs which are empty and meaningless
% initial guess
poses = poses(2:end);
transitions = transitions(2:end);
observations = observations(2:end);
landmarks = initLandmarks(poses, observations); % obtained using odometry edges
% ground truth
landmarks_gt = landmarks_gt(2:end);
poses_gt = poses_gt(2:end);
transitions_gt = transitions_gt(2:end);
observations_gt = observations_gt(2:end);

% %%% init stuff
% % initial pose
% mu = [0; 0; 0];

% % init trajectory
% trajectory = [mu(1), mu(2)];

% % least squares optimization

% % render graphics
% figure(1); title("rangeonly-slam");
% plotState(landmarks, mu);
