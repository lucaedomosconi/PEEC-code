function [Btot] = parallel_losses(id_pr,np,prec_in_small_tria,...
            distance_near,prec_in_near_tria,prec_in_far_tria,...
            midpoint,num_elements,J,X1,X2,X3)
Btot = zeros(3,num_elements/np);
for iii = 1:num_elements/np
    ii = id_pr*(num_elements/np) + iii;
  if mod(iii,20) == 0
    fprintf("=", id_pr,iii,num_elements/4)
  endif
  for jj = 1:num_elements
    if jj == ii
    %{
      fprintf("c")
      x1tr = X1(:,ii);
      x2tr = x1tr + (X2(:,ii) - x1tr)/3;
      x3tr = x1tr + (X3(:,ii) - x1tr)/3;
      [gpst,gwst] = tri_gauss_points(prec_in_small_tria,x1tr,x2tr,x3tr);
      gps = squeeze(gpst);
      gws = squeeze(gwst);
      for kk = 1:size(gps,2)
        Btot(:,iii) += gws(kk).*cross(J(:,ii),midpoint(:,ii)-gps(:,kk))./(norm(midpoint(:,ii)-gps(:,kk))^3);
      end
      x1tr = X2(:,ii);
      x2tr = x1tr + (X3(:,ii) - x1tr)/3;
      x3tr = x1tr + (X1(:,ii) - x1tr)/3;
      [gpst,gwst] = tri_gauss_points(prec_in_small_tria,x1tr,x2tr,x3tr);
      gps = squeeze(gpst);
      gws = squeeze(gwst);
      for kk = 1:size(gps,2)
        Btot(:,iii) += gws(kk).*cross(J(:,ii),midpoint(:,ii)-gps(:,kk))./(norm(midpoint(:,ii)-gps(:,kk))^3);
      end
      x1tr = X3(:,ii);
      x2tr = x1tr + (X1(:,ii) - x1tr)/3;
      x3tr = x1tr + (X2(:,ii) - x1tr)/3;
      [gpst,gwst] = tri_gauss_points(prec_in_small_tria,x1tr,x2tr,x3tr);
      gps = squeeze(gpst);
      gws = squeeze(gwst);
      for kk = 1:size(gps,2)
        Btot(:,iii) += gws(kk).*cross(J(:,ii),midpoint(:,ii)-gps(:,kk))./(norm(midpoint(:,ii)-gps(:,kk))^3);
      end
      %}
    elseif norm(midpoint(ii) - midpoint(jj)) < distance_near
      [gpst,gwst] = tri_gauss_points(prec_in_near_tria,X1(:,jj),X2(:,jj),X3(:,jj));
      gps = squeeze(gpst);
      gws = squeeze(gwst);
      for kk = 1:size(gps,2)
        Btot(:,iii) += gws(kk).*cross(J(:,ii),midpoint(:,ii)-gps(:,kk))./(norm(midpoint(:,ii)-gps(:,kk))^3);
      end
    else
      [gpst,gwst] = tri_gauss_points(prec_in_far_tria,X1(:,jj),X2(:,jj),X3(:,jj));
      gps = squeeze(gpst);
      gws = squeeze(gwst);
      for kk = 1:size(gps,2)
        Btot(:,iii) += gws(kk).*cross(J(:,ii),midpoint(:,ii)-gps(:,kk))./(norm(midpoint(:,ii)-gps(:,kk))^3);
      end
    endif
  end
end
Btot = {Btot'};
end