function [ fracture_width ] = generate_frac_fin( params, data_generator )
%GENERATE_FRAC_FIN Summary of this function goes here
%   Detailed explanation goes here

fracture_width=data_generator.prior_mean.*exp(data_generator.gen_basis*params);
end

