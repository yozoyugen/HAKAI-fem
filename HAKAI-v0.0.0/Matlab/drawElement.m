function drawElement(coordmat, elementmat, flag_text)

% HAKAI - A 3-dimensional finite element program
% Copyright (c) 2024 Yozo Yugen
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see http://www.gnu.org/licenses/.


    figure

    nNode = size(coordmat,1);
    plot3(coordmat(:,1), coordmat(:,2), coordmat(:,3), 'bo' );
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis equal;
    hold on;

    if flag_text == 1
        for i = 1 : nNode
            text(coordmat(i,1), coordmat(i,2), coordmat(i,3), sprintf('#%d', i) );
        end
    end


    nElement = size(elementmat,1);
    for i = 1 : nElement
        nd = elementmat(i,:);
        ii = nd([1 2 3 4 1]);
        line( coordmat( ii,1),  coordmat( ii,2), coordmat( ii,3), 'color', 'b' )
        ii = nd([5 6 7 8 5]);
        line( coordmat( ii,1),  coordmat( ii,2), coordmat( ii,3), 'color','b' )
        ii = nd([1 5]);
        line( coordmat( ii,1),  coordmat( ii,2), coordmat( ii,3), 'color','b' )
        ii = nd([2 6]);
        line( coordmat( ii,1),  coordmat( ii,2), coordmat( ii,3), 'color','b' )
        ii = nd([3 7]);
        line( coordmat( ii,1),  coordmat( ii,2), coordmat( ii,3), 'color','b' )
        ii = nd([4 8]);
        line( coordmat( ii,1),  coordmat( ii,2), coordmat( ii,3), 'color','b' )
    end

end
