%This function is useful to plot meshing of any region  within domain in 3d
%It reads data from binary files created by Hercules. One can plot domain in
%3d and color the elements according to Vp, Vs or rho values. It is also
%useful when it is desired to see which processor posseses which element.
%All these options are available through parameters_for_matlab.in. Usage of
%parameters file is added below.

%A sample parameters_for_matlab.in file

%1-where to start plotting (x coordinate) : 0 ** These are coordinates of
%2-where to end plotting (x coordinate) : 800    domain where you want to
%3-where to start plotting (y coordinate) : 0    plot.
%4-where to end plotting (y coordinate) : 800
%5-where to start plotting (z coordinate) : 0
%6-where to end plotting (z coordinate) : 800
%7-which of these to plot as 4th dimension Vs(1) Vp(2) or Rho(3): 1 ** Put 1, 2, or 3 only
%8-number of processors used : 8  **This is essential
%9-path to directory where coordinates are stored :/Users/testrun/outputfiles/For_Matlab
%10-path to directory where data are stored : /Users/testrun/outputfiles/For_Matlab
%11-plot processor(p) or data(d) : d
%** If you put p, colors will show which processor posseses which element
%and parameter 7 is disregarded.
%If you put d, colors will show Vs, Vp or Rho according to parameter 7.

%The ordering of parameters are important and should not be changed. You are
%free to change the names of parameters as long as followed by a colon(:).
%Finally, each parameter should appear in a single line and there should not
%be any blank lines in between parameters. When calling fuction, one must
%include the name of path to parameters between quotation marks i.e.
%plot3d_Hercules('/Users/Desktop/parameters_for_matlab.in')

function  plotmesh( path_to_parameters )
    %Parsing parameters from parameters_for_matlab.in
    fid = fopen(path_to_parameters, 'r');
    i=0;
    temp_numbers = zeros(8,1);
    temp_paths = cell(3,1);
    while 1
        i=i+1;
        tline = fgetl(fid);
        if ~ischar(tline)
            break
        end
        if (i ~= 9 && i ~= 10 && i~=11)
            C = textscan(tline, '%s %d', 'delimiter',':');
            temp_numbers(i) = C{2};
        else
            C = textscan(tline, '%s %s', 'delimiter',':');
            temp_paths(i-8) = C{2};
        end
        
    end
    fclose(fid);
    
    %Initiliazing of variables
    %One can plot any region within domain with these inputs
    x_min = temp_numbers(1);
    x_max = temp_numbers(2);
    y_min = temp_numbers(3);
    y_max = temp_numbers(4);
    z_min = temp_numbers(5);
    z_max = temp_numbers(6);
    Forth_dim = temp_numbers(7);
    
    %This is done to grab data from different text files written by different
    %processors.Need to know how many processors are used.
    
    number_processors = temp_numbers(8);
    
    directory_coord = char(temp_paths(1));
    directory_data  = char(temp_paths(2));
    
    %plot processors or data
    plot_processor_p_or_data_d = char(temp_paths(3));
    
    if (strcmp(plot_processor_p_or_data_d,'p'))
        plot_processors = 1;
    else
        plot_processors = 0;
    end
    
    %-----------------------Parsing and Initiliazing Ends Here-----------------
    %There are 24 coordinates and 3 data values. A is coordinate matrix.(number
    %of points by 24). B is data matrix(number of points by 3)

    % geid is the global element id at the beginning of each element record.
    geid_byte_size = 8; % 8 byte integer (int64)
    coordinate_byte_size = 8; % 8 byte float (float64)
    data_byte_size = 4; % 4 byte float (float32)
    total_bytes_per_coordinate = geid_byte_size + 3*coordinate_byte_size;
    skip_bytes_per_coordinate = total_bytes_per_coordinate - coordinate_byte_size;
    total_bytes_per_data = geid_byte_size + 3*data_byte_size;
    skip_bytes_per_data = total_bytes_per_data - data_byte_size;

    j = 0;    
    for i = 0:number_processors-1
        mesh_coordinate_file = [directory_coord '/mesh_coordinates.' num2str(i)];
        if (exist(mesh_coordinate_file ,'file'))
            fid1 = fopen(mesh_coordinate_file);
            file_byte_size = dir(mesh_coordinate_file).bytes;
            %Coord_Matrix is a auxillary matrix for holding coordinates and data.
            Coord_Matrix = zeros(3, file_byte_size/total_bytes_per_coordinate);
            for direction = 1:3
                % skip the first geid
                fseek(fid1, geid_byte_size + (direction-1)*coordinate_byte_size, 'bof');
                % Read the rest of the file
                Coord_Matrix(direction, :) = fread(fid1, [1, inf], 'float64', skip_bytes_per_coordinate);
            end
            % Reshape the matrix to 24xinf
            Coord_Matrix = reshape(Coord_Matrix, 24, []);
            
            %Converting Coord_Matrix to a (number of points) by (24) matrix
            auxilary_coord = Coord_Matrix';
            if ( j == 0 )
                A = auxilary_coord;
            else
                A = [A ; auxilary_coord];
            end
            fclose(fid1);
            [row_size, column_size] = size(auxilary_coord);
            %Repeat same procedure for data matrix if necessaery

            if (plot_processors == 0)
                mesh_data_file = [directory_data '/mesh_data.' num2str(i)];
                fid2 = fopen(mesh_data_file);
                file_byte_size = dir(mesh_data_file).bytes;
                Data_Matrix = zeros(3, file_byte_size/total_bytes_per_data);
                for direction = 1:3
                    % skip the first geid
                    fseek(fid2, geid_byte_size + (direction-1)*data_byte_size, 'bof');
                    % Read the rest of the file
                    Data_Matrix(direction, :) =  fread(fid2, [1, inf], 'float32', skip_bytes_per_data);
                end
                
                fclose(fid2);
                auxilary_data  = Data_Matrix';
                if ( j == 0 )
                    B = auxilary_data;
                else
                    B = [B ; auxilary_data];
                end
            end
            
            %Note that B is row_size by 3
            if (plot_processors == 1)
                if ( j == 0 )
                    B = i*ones(row_size,3);
                else
                    B = [B ; i*ones(row_size,3)];
                end
                
            end
            j=j+1 ;
        end
    end
    
    [r, c]=size(A);
    
    faces_matrix=[1 3 4 2;
    5 7 8 6;
    7 8 4 3;
    5 6 2 1;
    5 7 3 1;
    6 8 4 2;];
    
    %If the coordinates are in the max and min limits.Comparasion with left corner coordinates
    for i=1:r
        if( A(i,1) >= x_min && A(i,1) < x_max &&...
            A(i,2) >= y_min && A(i,2) < y_max &&...
            A(i,3) >= z_min && A(i,3) < z_max )
            vertex_matrix=zeros(8,3);
            
            for j=0:7
                
                vertex_matrix(j+1,1)=A(i,j*3+1);
                vertex_matrix(j+1,2)=A(i,j*3+2);
                vertex_matrix(j+1,3)=-A(i,j*3+3);%z_downwards
            end
            
            patch('Vertices',vertex_matrix,'Faces',faces_matrix,'FaceColor','flat','FaceVertexCData',B(i,Forth_dim));
            
            patch('Vertices',vertex_matrix,'Faces',faces_matrix,'FaceColor','flat','FaceVertexCData',B(i,Forth_dim));
            hold on;
        end
    end
    %set(gcf,'Menubar','none','Name','Cube', ...
    %'NumberTitle','off','Position',[10,350,1000,1000], ...
    %'Color',[0.5 0.8 0.3]);
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1])
    colorbar('EastOutside');
    %light('Position',[100 100 0],'Style','local');
    light('Position',[0 0 -100]);
    material shiny;
    alpha(0.4);
    alphamap('rampdown');
    camlight(45,45);
    lighting phong
    view(30,30);
    %zoom(2);
    axis([x_min x_max y_min y_max -z_max -z_min]); %z downwards in the Hercules
    xlabel('X, N-S (m)');
    ylabel('Y, E-W (m)');
    zlabel('Z (m)');
    
    if(plot_processors==1)
        colors = 'id of processors';
    else
        if(Forth_dim==1)
            colors = 'Vs';
        elseif(Forth_dim==2)
            colors = 'Vp';
        elseif(Forth_dim==3)
            colors = 'Rho';
        end
    end
    title(['Colors represent ' colors]);
    
end