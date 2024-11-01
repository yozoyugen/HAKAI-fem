function [newPART] = Gilgamsh(PART)

% Gilgamsh - A three-dimensional finite element re-mesh generator
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

    final_coord = [];
    coord = [];
    element = [];
    node_modify = [];

    type = PART.element_type;
    nUnit = 27;
    nEN = 8; 

    cm = PART.coordmat'; % row major
    em = PART.elementmat'; % row major
    nElement = PART.nElement;


    if contains(type, 'C3D8')

        % PART.nNode = 8;
        % PART.coordmat = ...
        %     [0  0  0
        %      1  0  0
        %      1  1  0
        %      0  1  0
        %      0  0  1
        %      1  0  1
        %      1  1  1
        %      0  1  1];
        % 
        % PART.nElement = 1;
        % PART.elementmat = ...
        %     [1 2 3 4 5 6 7 8];
    
        coord = zeros(nElement*27, 3);
        element = zeros(nElement*8, 8);
    
        %          8       25       7
        %      26       27      24
        %   5       23      6
        %
        %          17      20       16 
        %      21      22       19
        %   14      18      15
        %
        %          4       11       3
        %      12       13      10
        %   1       9       2    
    
    
        for i = 1 : nElement
    
            element([1:8]+(i-1)*8, :) = ...
                [1 9 13 12 14 18 22 21
                 9 2 10 13 18 15 19 22
                 13 10 3 11 22 19 16 20
                 12 13 11 4 21 22 20 17
                 14 18 22 21 5 23 27 26
                 18 15 19 22 23 6 24 27
                 22 19 16 20 27 24 7 25
                 21 22 20 17 26 27 25 8] + (i-1)*27;
    
            %offset = nNode + (27-8)*(i-1);
            c9 = (cm(em(i,1),:) + cm(em(i,2),:)) * 0.5;
            c10 = (cm(em(i,2),:) + cm(em(i,3),:)) * 0.5;
            c11 = (cm(em(i,3),:) + cm(em(i,4),:)) * 0.5;
            c12 = (cm(em(i,1),:) + cm(em(i,4),:)) * 0.5;
            c13 = (c10 + c12) * 0.5;
    
            c14 = (cm(em(i,1),:) + cm(em(i,5),:)) * 0.5;
            c15 = (cm(em(i,2),:) + cm(em(i,6),:)) * 0.5;
            c16 = (cm(em(i,3),:) + cm(em(i,7),:)) * 0.5;
            c17 = (cm(em(i,4),:) + cm(em(i,8),:)) * 0.5;
            c18 = (c14 + c15) * 0.5;
            c19 = (c15 + c16) * 0.5;
            c20 = (c16 + c17) * 0.5;
            c21 = (c14 + c17) * 0.5;
            c22 = (c19 + c21) * 0.5;
    
            c23 = (cm(em(i,5),:) + cm(em(i,6),:)) * 0.5;
            c24 = (cm(em(i,6),:) + cm(em(i,7),:)) * 0.5;
            c25 = (cm(em(i,7),:) + cm(em(i,8),:)) * 0.5;
            c26 = (cm(em(i,5),:) + cm(em(i,8),:)) * 0.5;
            c27 = (c24 + c26) * 0.5;
    
            coord([1:8]+(i-1)*27,:) = cm(em(i,:),:);
            coord([9:27]+(i-1)*27,:) = [c9;c10;c11;c12;c13;c14;c15;c16;c17;c18;c19;c20;
                                      c21;c22;c23;c24;c25;c26;c27];
    
        end
    
        node_modify = [1:27*nElement];
        final_coord = coord(1:27,:);

    elseif contains(type, 'S4')

        %   4       7       3
        %
        %   8       9       6 
        %   
        %   1       5       2    

        coord = zeros(nElement*9, 3);
        element = zeros(nElement*4, 4);

        for i = 1 : nElement

            element((1:4)+(i-1)*4, :) = ...
                [1 5 9 8
                 5 2 6 9
                 9 6 3 7
                 8 9 7 4] + (i-1)*9;
    
            c5 = (cm(em(i,1),:) + cm(em(i,2),:)) * 0.5;
            c6 = (cm(em(i,2),:) + cm(em(i,3),:)) * 0.5;
            c7 = (cm(em(i,3),:) + cm(em(i,4),:)) * 0.5;
            c8 = (cm(em(i,1),:) + cm(em(i,4),:)) * 0.5;
            c9 = (c6 + c8) * 0.5;

            coord((1:4)+(i-1)*9,:) = cm(em(i,:),:);
            coord((5:9)+(i-1)*9,:) = [c5;c6;c7;c8;c9];
            
        end

        node_modify = [1:9*nElement];
        final_coord = coord(1:9,:);        

        nUnit = 9;
        nEN = 4; 
    
    end

    tol = 1E-10;

    for i = nUnit+1 : nElement * nUnit
        if i > size(coord,1)
            break;
        end

        ci = coord(i,:);
        for j = 1 : size(final_coord,1)
            cj = final_coord(j,:);
            v = norm(ci-cj);
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

    for i = 1 : size(element,1)
        for j = 1 : nEN
            element(i,j) = node_modify( element(i,j) );
        end
    end


    newPART.coordmat = final_coord'; % column major
    newPART.elementmat = element'; % column major
    newPART.nNode = size(final_coord,1);
    newPART.nElement = size(element,1);
    newPART.shell_thickness = [];
    if contains(type, 'S4')
        newPART.shell_thickness = PART.shell_thickness;
    end


    % str = "";
    % for i = 1 : size(final_coord)
    %     str = str + ...
    %     sprintf("%d,   %.6e,   %.6e,   %.6e\n", i, final_coord(i,1), final_coord(i,2), final_coord(i,3));
    % end
    % str
    % 
    % str = "";
    % for i = 1 : size(element)
    %     str = str + ...
    %     sprintf("%d, %d, %d, %d, %d, %d, %d, %d, %d\n", i, ...
    %             element(i,1), element(i,2), element(i,3), element(i,4), ...
    %             element(i,5), element(i,6), element(i,7), element(i,8) );
    % end
    % str

    fname = sprintf("mesh_temp.txt", i);
    out = fopen(fname,"w");
    fprintf(out,"*Node\n");
    
    for i = 1 : size(final_coord,1)
        fprintf(out, sprintf("%d,   %.6e,   %.6e,   %.6e\n", i, final_coord(i,1), final_coord(i,2), final_coord(i,3)) );
    end

    if contains(type, 'C3D8')

        fprintf(out,"*Element, type=C3D8R\n");
        
        for i = 1 : size(element,1)
            fprintf(out, sprintf("%d, %d, %d, %d, %d, %d, %d, %d, %d\n", i, ...
                                 element(i,1), element(i,2), element(i,3), element(i,4), ...
                                 element(i,5), element(i,6), element(i,7), element(i,8) ) );
        end

    elseif contains(type, 'S4')



    end

    fclose(out);

end