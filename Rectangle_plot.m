function Rectangle_plot(origin,X,Y,color)
% RECTANGLE_PLOT plots a cube with dimension of X, Y.
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
% Plot a figure in the form of rectangles.
%
% EXAMPLES
% Rectangle_plot([0 0],2,3,'red')
%
% ------------------------------Code Starts Here------------------------------ %
% Define the vertexes of the unit cubic
ver = [0 1;
       1 1;
       1 0;
       0 0];
%  Define the faces of the unit cubic
fac = [1 2 3 4];
rectangle = [ver(:,1)*X+origin(1),ver(:,2)*Y+origin(2)];
patch('Faces',fac,'Vertices',rectangle,'FaceColor',color);
view(2)
end