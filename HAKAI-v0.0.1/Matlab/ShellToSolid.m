function [newPART] = ShellToSolid(PART)

% ShellToSolid - A three-dimensional finite element re-mesh generator
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

    newPART = [];

    coordmat = PART.coordmat'; % row major
    elementmat = PART.elementmat'; % row major

    nNode = PART.nNode;
    nElement = PART.nElement;

    nodeNormal = zeros(nNode, 3);

    figure;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis equal;
    hold on;

        for i = 1 : nNode
            text(coordmat(i,1), coordmat(i,2), coordmat(i,3), sprintf('#%d', i) );
        end


    mag = 0;
    for i = 1 : nElement
        nd = elementmat(i,:);
        ii = nd([1 2 3 4 1]);
        line( coordmat( ii,1),  coordmat( ii,2), coordmat( ii,3), 'color', 'b' )

        v1 = coordmat(ii(2), :) - coordmat(ii(1), :);
        v4 = coordmat(ii(4), :) - coordmat(ii(1), :);
        %mag = vecnorm(v1);
        n = cross(v1,v4);
        n = n / vecnorm(n);

        %ctr = sum( coordmat(ii(1:4), :)) /4;
        %quiver3( ctr(1), ctr(2), ctr(3), mag*n(1), mag*n(2), mag*n(3)  );

        nodeNormal(ii(1), :) = nodeNormal(ii(1), :) + n;
        nodeNormal(ii(2), :) = nodeNormal(ii(2), :) + n;
        nodeNormal(ii(3), :) = nodeNormal(ii(3), :) + n;
        nodeNormal(ii(4), :) = nodeNormal(ii(4), :) + n;

        mag = mag + (vecnorm(v1) + vecnorm(v4))/2;
    end
    mag = mag / nElement;


    for i = 1 : nNode
        %mag = 1;
        quiver3( coordmat(i,1), coordmat(i,2), coordmat(i,3),...
                 mag*nodeNormal(i,1), mag*nodeNormal(i,2), mag*nodeNormal(i,3)  );
    end



    %      8       7
    %   5       6
    %
    %      4       3   
    %   1       2

    cd_temp = zeros(nElement*8, 3);
    el_temp = zeros(nElement, 8);  
    th = PART.shell_thickness;

    for i = 1 : nElement

        el_temp(i, :) = [1 2 3 4 5 6 7 8] + (i-1)*8;

        el = elementmat(i,:);
        n1 = nodeNormal(el(1), :);
        n2 = nodeNormal(el(2), :);
        n3 = nodeNormal(el(3), :);
        n4 = nodeNormal(el(4), :);
        n1 = n1 / vecnorm(n1);
        n2 = n2 / vecnorm(n2);
        n3 = n3 / vecnorm(n3);
        n4 = n4 / vecnorm(n4);

        cd_temp(1+(i-1)*8,:) = coordmat(el(1),:) - th * 0.5 * n1;
        cd_temp(2+(i-1)*8,:) = coordmat(el(2),:) - th * 0.5 * n2;
        cd_temp(3+(i-1)*8,:) = coordmat(el(3),:) - th * 0.5 * n3;
        cd_temp(4+(i-1)*8,:) = coordmat(el(4),:) - th * 0.5 * n4;
        cd_temp(5+(i-1)*8,:) = coordmat(el(1),:) + th * 0.5 * n1;
        cd_temp(6+(i-1)*8,:) = coordmat(el(2),:) + th * 0.5 * n2;
        cd_temp(7+(i-1)*8,:) = coordmat(el(3),:) + th * 0.5 * n3;
        cd_temp(8+(i-1)*8,:) = coordmat(el(4),:) + th * 0.5 * n4;
        
    end    

    % figure;
    % xlabel('X');
    % ylabel('Y');
    % zlabel('Z');
    % axis equal;
    % hold on;
    % 
    % plot3(cd_temp(1:8,1), cd_temp(1:8,2), cd_temp(1:8,3), 'bo' );

    node_modify = [1:8*nElement];
    final_coord = cd_temp(1:8,:);
    
    tol = 1E-10;

    for i = 8+1 : nElement * 8

        % if i > size(coord,1)
        %     break;
        % end

        ci = cd_temp(i,:);
        for j = 1 : size(final_coord,1)
            cj = final_coord(j,:);
            v = vecnorm(ci-cj);
            if v < tol
                node_modify(i) = j;
                break
            end

            if j == size(final_coord,1)
                final_coord = [final_coord; ci];
                node_modify(i) = size(final_coord,1);
            end
        end
    end

    size(final_coord)

    for i = 1 : size(el_temp,1)
        for j = 1 : 8
            el_temp(i,j) = node_modify( el_temp(i,j) );
        end
    end


    newPART.coordmat = final_coord;  % row major
    newPART.elementmat = el_temp;  % row major
    
    drawElement(newPART.coordmat', newPART.elementmat', 1)


    fname = sprintf("shellsolid_temp.txt", i);
    out = fopen(fname,"w");
    fprintf(out,"*Node\n");
    
    for i = 1 : size(final_coord,1)
        fprintf(out, sprintf("%d,   %.6e,   %.6e,   %.6e\n", i, final_coord(i,1), final_coord(i,2), final_coord(i,3)) );
    end

    fprintf(out,"*Element, type=C3D8R\n");
    
    for i = 1 : size(el_temp,1)
        fprintf(out, sprintf("%d, %d, %d, %d, %d, %d, %d, %d, %d\n", i, ...
                             el_temp(i,1), el_temp(i,2), el_temp(i,3), el_temp(i,4), ...
                             el_temp(i,5), el_temp(i,6), el_temp(i,7), el_temp(i,8) ) );
    end

    fclose(out);
    

end