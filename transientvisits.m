% transientvisits.m computes the expected number of visits to transient states

function [W] = transientvisits(Q)

W = inv(eye(size(Q))-Q);
