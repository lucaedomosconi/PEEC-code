function Btot = parallel_losses()

for ii = 1:num_elements
  if mod(ii,50) == 0
    disp(ii)
  endif
  for jj = 1:num_elements
    if jj == ii
      x1tr = X1(:,ii);
      x2tr = x1tr + (X2(:,ii) - x1tr)/3;
      x3tr = x1tr + (X3(:,ii) - x1tr)/3;
      [gpst,gwst] = tri_gauss_points(prec_in_small_tria,x1tr,x2tr,x3tr);
      gps = squeeze(gpst);
      gws = squeeze(gwst);
      for kk = 1:size(gps,2)
        Btot(:,ii) += gws(kk).*cross(J(:,ii),gps(:,kk)-midpoint(:,ii))./(norm(gps(:,kk)-midpoint(:,ii))^3);
      end
      x1tr = X2(:,ii);
      x2tr = x1tr + (X3(:,ii) - x1tr)/3;
      x3tr = x1tr + (X1(:,ii) - x1tr)/3;
      [gpst,gwst] = tri_gauss_points(prec_in_small_tria,x1tr,x2tr,x3tr);
      gps = squeeze(gpst);
      gws = squeeze(gwst);
      for kk = 1:size(gps,2)
        Btot(:,ii) += gws(kk).*cross(J(:,ii),gps(:,kk)-midpoint(:,ii))./(norm(gps(:,kk)-midpoint(:,ii))^3);
      end
      x1tr = X3(:,ii);
      x2tr = x1tr + (X1(:,ii) - x1tr)/3;
      x3tr = x1tr + (X2(:,ii) - x1tr)/3;
      [gpst,gwst] = tri_gauss_points(prec_in_small_tria,x1tr,x2tr,x3tr);
      gps = squeeze(gpst);
      gws = squeeze(gwst);
      for kk = 1:size(gps,2)
        Btot(:,ii) += gws(kk).*cross(J(:,ii),gps(:,kk)-midpoint(:,ii))./(norm(gps(:,kk)-midpoint(:,ii))^3);
      end
    elseif norm(midpoint(ii) - midpoint(jj)) < distance_near
      [gpst,gwst] = tri_gauss_points(prec_in_near_tria,X1(:,jj),X2(:,jj),X3(:,jj));
      gps = squeeze(gpst);
      gws = squeeze(gwst);
      for kk = 1:size(gps,2)
        Btot(:,ii) += gws(kk).*cross(J(:,ii),gps(:,kk)-midpoint(:,ii))./(norm(gps(:,kk)-midpoint(:,ii))^3);
      end
    else
      [gpst,gwst] = tri_gauss_points(prec_in_far_tria,X1(:,jj),X2(:,jj),X3(:,jj));
      gps = squeeze(gpst);
      gws = squeeze(gwst);
      for kk = 1:size(gps,2)
        Btot(:,ii) += gws(kk).*cross(J(:,ii),gps(:,kk)-midpoint(:,ii))./(norm(gps(:,kk)-midpoint(:,ii))^3);
      end
    endif
  end
end