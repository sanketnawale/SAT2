function [Distances, Elevation_Angles, Ground_Distances, Visibility, Num_Visible_Sats] = ...
    Satellite_Geometry(H, Node_Coordinates, oev, Earth_Radius, Time_Vector)

    Num_Satellites = size(oev, 1); 
    num_steps = length(Time_Vector);
    num_nodes = size(Node_Coordinates, 1);

    % 🌍 Convert Ground Node Positions to Cartesian (ECEF)
    Node_Positions = zeros(num_nodes, 3);
    for n = 1:num_nodes
        lat = Node_Coordinates(n, 1);
        lon = Node_Coordinates(n, 2);
        Node_Positions(n, :) = [
            Earth_Radius * cosd(lat) * cosd(lon),
            Earth_Radius * cosd(lat) * sind(lon),
            Earth_Radius * sind(lat)
        ];
    end

    % ✅ Initialize Output Arrays
    Distances = zeros(num_nodes, Num_Satellites, num_steps);
    Elevation_Angles = zeros(num_nodes, Num_Satellites, num_steps);
    Ground_Distances = zeros(num_nodes, Num_Satellites, num_steps);
    Visibility = zeros(num_nodes, Num_Satellites, num_steps);
    Num_Visible_Sats = zeros(num_nodes, num_steps);

    % 🔄 Loop Over Time Steps
    for t = 1:num_steps
        
       fprintf('\n⏳ Time %.2f min:\n', Time_Vector(t) / 60);  % Convert sec to min

        
        for s = 1:Num_Satellites
            % ✅ Compute satellite position dynamically at time t
            mean_motion = 2 * pi / 5400;  % Assume 5400s orbital period
            true_anomaly = oev(s, 6) + mean_motion * Time_Vector(t);

            [Sat_X, Sat_Y, Sat_Z] = orbital_to_cartesian(...
                oev(s, 1), oev(s, 2), oev(s, 3), ...
                oev(s, 5), oev(s, 4), true_anomaly);
            Satellite_Position = [Sat_X, Sat_Y, Sat_Z];

            % 🛰️ Print Satellite Position
            fprintf('🛰️ Sat %d Pos -> X: %.2f km, Y: %.2f km, Z: %.2f km\n', ...
                    s, Sat_X / 1e3, Sat_Y / 1e3, Sat_Z / 1e3);

            for n = 1:num_nodes
                % ✅ Compute Distance Between Node and Satellite
                Vector = Satellite_Position - Node_Positions(n, :);
                Distances(n, s, t) = norm(Vector);

                % ✅ Compute Ground Distance
                Ground_Distances(n, s, t) = sqrt(Distances(n, s, t)^2 - H^2);

                % ✅ Compute Elevation Angle
                Elevation_Angles(n, s, t) = asind(dot(Vector, Node_Positions(n, :)) / ...
                                               (norm(Vector) * norm(Node_Positions(n, :))));

                % ✅ Visibility Condition
                if Elevation_Angles(n, s, t) > 10  % Satellites above 10° elevation are visible
                    Visibility(n, s, t) = 1;
                end
            end
        end

        % ✅ Count Visible Satellites per Node
        Num_Visible_Sats(:, t) = sum(Visibility(:, :, t), 2);

        % ✅ Print visibility results
        fprintf('🔍 Rome sees %d satellites, Milan sees %d satellites\n', ...
                Num_Visible_Sats(1, t), Num_Visible_Sats(2, t));
    end
end
