function Btot = LineBcompute(id_pr,np,x_line,X1,X2,X3,J)
num_elements = size(x_line,2);
num_tria = size(X1,2);
Btot = zeros(3,num_elements/4);
for k = 1:num_elements/4;
    if mod(k,5) == 0
        fprintf('=')
    endif
    kk = k+id_pr*num_elements/4;
    for ii = 1:num_tria
        [gps,gws] = tri_gauss_points(1,X1(:,ii),X2(:,ii),X3(:,ii));
        Btot(:,k) += gws.*cross(J(:,ii),x_line(:,kk)-gps)./(norm(x_line(:,kk)-gps)^3);
    end
end
Btot = {Btot'};
