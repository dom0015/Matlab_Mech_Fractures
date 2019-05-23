function [ y ] = sample_fracture( sample_params,points_normalized )
%SAMPLE_FRACTURE Summary of this function goes here
%   Detailed explanation goes here

sigma=1;
lambda=0.5;
n=100;
[ eig_value,eig_func ] = KL_1D( n,sigma,lambda );

[ y ] = frac_sample( sample_params,points_normalized,eig_value,eig_func );
y=exp(y);




end

