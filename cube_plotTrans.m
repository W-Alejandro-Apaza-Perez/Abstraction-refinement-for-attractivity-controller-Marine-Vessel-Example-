function cube_plotTrans(origin,X,Y,Z,color)
% CUBE_PLOT plots a cube with dimension of X, Y, Z.
%

%         List of colors
%         b blue
%         g green
%         r red
%         c cyan
%         m magenta
%         y yellow
%         k black
%         w white
% OUPUTS:
% Plot a figure in the form of cubics.
%
% EXAMPLES
% cube_plot([0,0,0],2,3,4,'red')
%
% ------------------------------Code Starts Here------------------------------ %
% Define the vertexes of the unit cubic
ver = [1 1 0;
    0 1 0;
    0 1 1;
    1 1 1;
    0 0 1;
    1 0 1;
    1 0 0;
    0 0 0];
%  Define the faces of the unit cubic
fac = [1 2 3 4;
    4 3 5 6;
    6 7 8 5;
    1 2 8 7;
    6 7 1 4;
    2 3 5 8];
cube = [ver(:,1)*X+origin(1),ver(:,2)*Y+origin(2),ver(:,3)*Z+origin(3)];
patch('Faces',fac,'Vertices',cube,'FaceColor',color,'FaceAlpha',.2,'EdgeAlpha',.2);
view(3)
end