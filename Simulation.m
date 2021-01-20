classdef Simulation
    
    properties
        N;
        g_min;      % Minimum of graph drawing
        g_max;      % Max of graph for drawing
        r_c;        % robot centers
        t_c;        % task centers
        r_colors;   % colors of the robots
        radii;      % radius of the circle shape of the robots
        width,      % width of the square shape of the tasks
        height      % height of the square shape of the tasks
    end
    
    methods
        function obj = Simulation(N_, g_min_, g_max_)
            % This is constructor of the class
            obj.N = N_;
            obj.g_min = g_min_;
            obj.g_max = g_max_;
            obj.width = 1;
            obj.height = 1;
            obj.r_c = zeros(N_,2);
            obj.t_c = zeros(N_,2);
            obj.r_colors = zeros(N_,3);
            obj.radii = ones(N_,1)*0.5;
        end
        
        function [] = draw_lines(obj, x)
            % Draw lines from robots to the tasks
            b = 1;
            e = obj.N;
            for ii = 1:obj.N
                loc = find(x(1, b:e) == 1);
                plot([obj.r_c(ii,1),obj.t_c(loc,1),],[obj.r_c(ii,2),obj.t_c(loc,2)],'--','color',obj.r_colors(ii,:));
                b = b + obj.N;
                e = e + obj.N;
            end
        end
        
        function [obj, c, A_eq, b_eq, LB, UB] = get_constraint_matrices(obj)
            % Get robot, task coordinates. Get constraint matrices
            r_min = obj.g_min+5;
            r_max = obj.g_max-5;
            
            obj.r_colors = rand(obj.N,3);
            
            for ii = 1:obj.N
                p_x = r_min+rand(1,1)*(r_max-r_min);
                p_y = r_min+rand(1,1)*(r_max-r_min);

                if ii ~= 1
                    while ~is_valid(p_x, p_y, ii, obj.r_c)
                        p_x = r_min+rand(1,1)*(r_max-r_min);
                        p_y = r_min+rand(1,1)*(r_max-r_min);
                    end
                    obj.r_c(ii,1) = p_x;
                    obj.r_c(ii,2) = p_y;
                else
                    obj.r_c(ii,1) = p_x;
                    obj.r_c(ii,2) = p_y;
                end
            end

            for ii = 1:obj.N
                p_x = r_min+rand(1,1)*(r_max-r_min);
                p_y = r_min+rand(1,1)*(r_max-r_min);

                if ii ~= 1
                    while ~(is_valid(p_x, p_y, ii, obj.t_c) && is_valid(p_x, p_y, obj.N+1, obj.r_c))
                        p_x = r_min+rand(1,1)*(r_max-r_min);
                        p_y = r_min+rand(1,1)*(r_max-r_min);
                    end
                    obj.t_c(ii,1) = p_x;
                    obj.t_c(ii,2) = p_y;
                else
                    obj.t_c(ii,1) = p_x;
                    obj.t_c(ii,2) = p_y;
                end
            end
            
            c = get_c_vector(obj);
            LB = zeros((obj.N*obj.N),1);
            UB = ones((obj.N*obj.N),1);
            
            A_eq = zeros(2*obj.N,obj.N*obj.N);
            ind = 0;
            st = size(A_eq, 2) / obj.N;
            for ii=1:obj.N
                A_eq(ii,ii+ind:ii*st) = ones(1, obj.N);
                A_eq(obj.N+1:end, ii+ind:ii*st) = eye(obj.N);
                ind = ind + obj.N-1;
            end
            b_eq = ones(2*obj.N, 1);
            
        end
        
        function [] = create_environment(obj)
            % Create simulation environment
            for ii = 1:obj.N
                xlim([obj.g_min obj.g_max])
                ylim([obj.g_min obj.g_max])

                axis square

                viscircles(obj.r_c(ii,:),obj.radii(ii),'Color', obj.r_colors(ii,:));
                draw_square(obj,obj.t_c(ii,1), obj.t_c(ii,2));
            end
        end
        
    end
end

function res = get_c_vector(obj)
    res = zeros(1, obj.N*obj.N);
    ind = 1;
    for ii = 1:obj.N
        for jj = 1:obj.N
            res(ind) = round(pdist([obj.r_c(ii,1),obj.r_c(ii,2);obj.t_c(jj,1),obj.t_c(jj,2)], 'euclidean'));
            ind = ind + 1;
        end
    end
end

function val = is_valid(p_x, p_y, ind, set)
    % Check whether there is enough distance between points in order not to
    % collide
    flg = 1;
    for ii = 1:ind-1
        dist_ = pdist([p_x,p_y;set(ii,1),set(ii,2)], 'euclidean');
        if dist_ <= 5.0
            flg = 0;
            break
        end
    end

    if flg == 1
        val = true;
    else
        val = false;
    end
end

function [] = draw_square(obj, c_x, c_y)
    xCenter = c_x;
    yCenter = c_y;
    xLeft = xCenter - obj.width/2;
    yBottom = yCenter - obj.height/2;

    rectangle('Position', [xLeft, yBottom, obj.width, obj.height], 'EdgeColor', 'k', 'LineWidth', 0.5);
end
