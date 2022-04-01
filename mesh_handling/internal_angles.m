function angles = internal_angles(V, F)
%angles = internal_angles(V, F)
%
% Parameters
% ----------
% V : #vertices x dim float
%   mesh vertices vertices
% F : #faces x 3 int
%   mesh connectivity list, indices into V of connected vertices
%
% Returns
% -------
% angles : #tri x 3 float
%   internal angles for angles across from vertices 23, 13, and 12
%
% NPMitchell 2022

if size(V, 2) == 3
    % 3d mesh vertices
    s12 = vecnorm(V(F(:,2),:) - V(F(:,1),:), 2, 2);
    s31 = vecnorm(V(F(:,3),:) - V(F(:,1),:), 2, 2);
    s23 = vecnorm(V(F(:,3),:) - V(F(:,2),:), 2, 2);
    a23 = acos((s12.^2 + s31.^2 - s23.^2)./(2.*s12.*s31));
    a31 = acos((s23.^2 + s12.^2 - s31.^2)./(2.*s23.*s12));
    a12 = acos((s31.^2 + s23.^2 - s12.^2)./(2.*s31.*s23));
    angles = [a23 a31 a12];
elseif size(V, 2) == 2
    error('handle 2d case here')
else
    error('handle N-dim case here')
end

