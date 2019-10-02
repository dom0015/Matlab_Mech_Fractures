function [ fracture_matrice] = fracture2cells_parameters( fracture_matrice,ALFA,MAT_FRAC )
%FRACTURE2CELLS Summary of this function goes here
%   Detailed explanation goes here
n=length(fracture_matrice);
for i=1:n
    fracture_matrice{i}.above_material=ALFA{i};
    fracture_matrice{i}.under_material=ALFA{i};
    fracture_matrice{i}.material=MAT_FRAC{i};
    fracture_matrice{i}.left_boundary=2; 
    fracture_matrice{i}.right_boundary=2;
end

end
