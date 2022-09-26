function [curlw] = curlw(k, X1, X2, X3)
  % Numbering of the nodes
  %
  %  eta2
  %   ^
  %   |
  %   3
  %   |`\
  %   |  `\
  %   1----2 --> eta1
  %
  E1 = X3 - X2;
  E2 = X1 - X3;
  E3 = X2 - X1;
  A2 = norm(cross(E2, E3, 1), 2, 'cols');
  switch (k)
    case 1
      curlw = E1 ./ A2;
    case 2
      curlw = E2 ./ A2;
    case 3
      curlw = E3 ./ A2;
    otherwise
      error('Bad k')
  end
end
