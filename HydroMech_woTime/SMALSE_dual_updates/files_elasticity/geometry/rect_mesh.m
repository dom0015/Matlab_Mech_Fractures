function [ node,elem,bdFlag ] = rect_mesh( height, width, h_elem, w_elem )
%% Creation mesh (nodes and elements) from rectangle positioned in (0,0)
% Input height and width of rectangle, number of elements in row - w_elem,
% number elements in column - h_elem.
% Output node matrix = position of verticles of mash
%        elem matrix = indexes of nodes creating elements
%        bdFlag matrix = for each element contain border-flag for each edge 

%% nodes matrix creation
x1 = linspace(0,width,w_elem+1);
n2 = h_elem+1;
x2 = linspace(0,height,h_elem+1);
[X1,X2] = meshgrid(x1,x2);
node = [X1(:) X2(:)];
%% elem matrix creation
m = 2*w_elem*h_elem;
elem = zeros(m,3,'int32');
idx = 1;
temp=reshape( repmat( 1:h_elem, 2,1 ), 1, [] );
for i=1:w_elem
    first=temp+(i-1)*n2;    % vector of first verticles of elements
    second=[temp(2:end) (h_elem+1)]+i*n2; % second verticles
    third=temp+1+(i-1)*n2;  % third verticles 
    third(1:2:end)=third(1:2:end)+n2;
    % together make elements of column in mesh
    elem(idx:(idx+h_elem*2-1),:)=[first'  second'  third']; 
    idx=idx+h_elem*2;
end
%% bdFlag matrix creation
% bdFlag - in each row=element edges of triangle (bc,ac,ab) - (0-no border)
%           (1-bottom,2-right,3-top,4-left)
bdFlag=zeros(m,3,'int8');
bdFlag(1:(h_elem*2):end,3)=1;
bdFlag((end-h_elem*2+1):2:end,1)=2;
bdFlag(h_elem*2:(h_elem*2):end,1)=3;
bdFlag(2:2:(h_elem*2),2)=4;
end