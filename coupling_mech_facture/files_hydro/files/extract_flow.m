function [ Q ] = extract_flow( A, u, node, elem, h_elem, Dirichlet_windows, downEdge, rightEdge, upEdge, leftEdge )
%EXTRACT_FLOW Summary of this function goes here
%   Detailed explanation goes here
%% EXTRACT FLOW FROM THE SOLUTION (ze slabe formulace)
t=A*u;
% t=reshape(t,h_elem+1,h_elem+1);
for i=1:size(Dirichlet_windows,1)
    a_=Dirichlet_windows(i,2);
    b_=Dirichlet_windows(i,3);
    okna_stred=(node(2:(h_elem+1),2)+node(1:(h_elem),2))/2;
    temp_=(1:h_elem)';
    okno=temp_((okna_stred>a_)&(okna_stred<b_));
    switch Dirichlet_windows(i,1)
        case 1	% okno dole
            o=downEdge(okno,:);
        case 2  % okno vpravo
            o=rightEdge(okno,:);
        case 3  % okno nahore
            o=upEdge(okno,:);
        case 4  % okno vlevo
            o=leftEdge(okno,:);
    end
    Q(i)=sum((t(unique(o))));
end

%% EXTRACT FLOW FROM THE SOLUTION (using diferences)
% for i=1:size(Dirichlet_windows,1)
%     a_=Dirichlet_windows(i,2);
%     b_=Dirichlet_windows(i,3);
%     okna_stred=(node(2:(h_elem+1),2)+node(1:(h_elem),2))/2;
%     temp_=(1:h_elem)';
%     okno=temp_((okna_stred>a_)&(okna_stred<b_));
%     switch Dirichlet_windows(i,1)
%         case 1	% okno dole
%             o=downEdge(okno,:);
%         case 2  % okno vpravo
%             o=rightEdge(okno,:);
%         case 3  % okno nahore
%             o=upEdge(okno,:);
%         case 4  % okno vlevo
%             o=leftEdge(okno,:);
%     end
%     temp = sum(ismember(elem,o),2)==2;
%     temp = elem(temp,:);
%     pressure = u(temp);
%     temp_x = node(:,1); coord_x = temp_x(temp);
%     temp_y = node(:,2); coord_y = temp_y(temp);
%     l=size(pressure,1); grad=zeros(l,2);
%     for j=1:l
%         M = [coord_x(j,:)' coord_y(j,:)' ones(3,1)];
%         b = pressure(j,:)';
%         sol = M\b;
%         grad(j,:) = sol(1:2)';
%     end
%     switch Dirichlet_windows(i,1)
%         case 1	% okno dole
%             Q(i) = -sum(grad(:,2));
%         case 2  % okno vpravo
%             Q(i) = sum(grad(:,1));
%         case 3  % okno nahore
%             Q(i) = sum(grad(:,2));
%         case 4  % okno vlevo
%             Q(i) = -sum(grad(:,1));
%     end
% end
end

