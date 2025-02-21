clear all; close all; clc;

%% 🌍 Simulation Parameters
tic
Simulation_T = 110 * 60;  % Total simulation time (110 minutes → seconds)
Time_Step = 60;           % Time step (1 min = 60 sec)
MonteCarlo = 5000;        % Monte Carlo simulation (number of packets per node)
Nodes = 3;                % Number of ground nodes (Rome, Milan, NodeRM)
Pkct_Per_Hour = 100;      % Packets per hour for each node

%% 🌍 Ground Nodes (Rome, Milan, NodeRM)
Node_Coordinates = [
    41.9028, 12.4964;  % Rome (Node 1)
    45.4642, 9.1900;   % Milan (Node 2)
    41.9, 12.5         % NodeRM (Node 3, near Rome)
];

%% 🛰️ Satellite Constellation (Walker)
Num_Satellites = 2;
Num_Planes = 6;
Orbital_Inclination = deg2rad(87); % Inclination in Radians
H = 1200e3; % Satellite altitude (meters)
Earth_Radius = 6378e3; % Earth radius (meters)
Time_Vector = 0:Time_Step:Simulation_T; % Time simulation vector

% 🛰️ Generate Walker Delta Constellation
oev = walker_delta(Num_Satellites, Num_Planes, 1, pi, Earth_Radius + H, Orbital_Inclination);

%% 📡 Call Updated Satellite Geometry Function
[Distances, Elevation_Angles, Ground_Distances, Visibility, Num_Visible_Sats, Sat_IDs] = ...
    Satellite_Geometry(H, Node_Coordinates, oev, Earth_Radius, Time_Vector);

%% 📊 Initialize Satellite Visibility Matrix
Visible_Sat_Matrix = zeros(length(Time_Vector), 3);  % Columns: [Time (min), Rome Visible Sats, Milan Visible Sats]

%% 📡 Random Access Logic for Packet Transmission, Reception & Collision Detection

% Initialize Success Rates and Collisions
SuccessRate = zeros(Nodes, length(Time_Vector));  % Success rate per node per time step
Collisions = zeros(Nodes, length(Time_Vector));   % Collision counts per node per time step

% Initialize Received Packets for NodeRM
Received_Packets_NodeRM = zeros(1, length(Time_Vector));

% Simulate Random Access and Reception Logic
for t = 1:length(Time_Vector)
    current_time_min = Time_Vector(t) / 60;  % Convert time to minutes
    fprintf('\n⏳ Time %.2f min: \n', current_time_min);

    % 🌍 Store Visibility Data in Matrix
    Visible_Sat_Matrix(t, :) = [current_time_min, Num_Visible_Sats(1, t), Num_Visible_Sats(2, t)];

    % Print visibility status
    fprintf('📡 Rome sees %d satellites: %s\n', Num_Visible_Sats(1, t), mat2str(Sat_IDs{1, t}));
    fprintf('📡 Milan sees %d satellites: %s\n', Num_Visible_Sats(2, t), mat2str(Sat_IDs{2, t}));
    fprintf('📡 NodeRM sees %d satellites: %s\n', Num_Visible_Sats(3, t), mat2str(Sat_IDs{3, t}));


    % 🔍 Identify Common Satellites
    Common_Sats = intersect(Sat_IDs{1, t}, Sat_IDs{2, t});
    if ~isempty(Common_Sats)
        fprintf('🔄 Common Satellites seen by both: %s\n', mat2str(Common_Sats));
    end

    % Tracking packets received by NodeRM
    NodeRM_Packet_Times = [];

    for n = 1:Nodes-1  % Only Rome (Node 1) and Milan (Node 2) transmit
        if Num_Visible_Sats(n, t) == 0  
            % 🛰️ No satellites visible, no packet transmission
            fprintf('🚫 No satellites visible for Node %d at %.2f min, no packets sent.\n', n, current_time_min);
            continue;
        end

        % Rome sends multiple identical packets (10)
        Num_Packets = 100; 
        Transmission_Times = rand(1, Num_Packets) * Time_Step; % Packets randomly sent within time step
        sorted_times = sort(Transmission_Times);

        % Select visible satellites
        Visible_Sats = Sat_IDs{n, t};

        % 🛑 **Check if any satellites are available**
        if isempty(Visible_Sats)
            fprintf('🚫 No satellites visible for Node %d at %.2f min, skipping transmission.\n', n, current_time_min);
            continue;
        end

        % Initialize Satellite Reception Time Storage
        Sat_Receive_Times = cell(Num_Satellites, 1);
        for s = 1:Num_Satellites
            Sat_Receive_Times{s} = [];
        end

        % Assign packets to all visible satellites
        for pkt = 1:Num_Packets
            for chosen_sat = Visible_Sats  % Send the same packet to all visible satellites
                if chosen_sat <= Num_Satellites
                    Sat_Receive_Times{chosen_sat} = [Sat_Receive_Times{chosen_sat}, sorted_times(pkt)];
                end
            end
        end

        % 🔥 Check for Collisions at Each Satellite
% 🔥 Check for Collisions at Each Satellite
        for s = 1:Num_Satellites
            if ~isempty(Sat_Receive_Times{s})
                sat_tx_times = sort(Sat_Receive_Times{s});

        % Detect collisions (packets arriving too close)
            collisions = sum(diff(sat_tx_times) < 0.01);  % Collision if packets arrive within 10ms
            total_packets = length(sat_tx_times);
    
        % Store results
            Collisions(n, t) = collisions;  % Store number of collisions
            SuccessRate(n, t) = total_packets - collisions;  % Only successful packets

             % ✅ **Fix the Relay Logic for NodeRM**
                if n == 1  
                % 🚀 Ensure we only take packets that were successfully received
                    if SuccessRate(n, t) > 0
                    successful_packets = sat_tx_times((collisions+1):end);  % Ensure valid range
                    NodeRM_Packet_Times = [NodeRM_Packet_Times, successful_packets]; 
                    end
                end
            end
        end

        % Print Debug Info
        fprintf('📊 Node %d transmitted %d packets, %d collisions\n', n, Num_Packets, Collisions(n, t));
    end

    % 🚀 **NodeRM Packet Reception Logic**
   % 🚀 **NodeRM Packet Reception Logic (Fixed)**
if ~isempty(NodeRM_Packet_Times)
    NodeRM_Packet_Times = sort(NodeRM_Packet_Times);

    % ✅ **Only accept packets from satellites visible to NodeRM**
    NodeRM_Visible_Sats = Sat_IDs{3, t};  % Get visible satellites for NodeRM
    if ~isempty(intersect(NodeRM_Visible_Sats, Sat_IDs{1, t}))  % Rome's satellites
        Received_Packets_NodeRM(t) = 1;  % Mark first packet as successful
        first_packet_time = NodeRM_Packet_Times(1);

        % 🚨 Collision Warning for Extra Packets
        if length(NodeRM_Packet_Times) > 1
            fprintf('🚨 COLLISION at NodeRM: %d packets received, but only 1 is accepted (Time: %.2f min)\n', ...
                length(NodeRM_Packet_Times), first_packet_time / 60);
        else
            fprintf('📡 NodeRM successfully received a packet at %.2f min\n', first_packet_time / 60);
        end
    else
        fprintf('📡 NodeRM Received No Packets at %.2f min (No successful relay, no visible satellites).\n', current_time_min);
    end
else
    fprintf('📡 NodeRM Received No Packets at %.2f min (No successful relay).\n', current_time_min);
end

end

%% 📊 Print Satellite Visibility Matrix
fprintf('\n==== 📊 Satellite Visibility Data ====\n');
disp(array2table(Visible_Sat_Matrix, 'VariableNames', {'Time_Min', 'Rome_Sats', 'Milan_Sats'}));

%% 📊 Final Summary
fprintf('\n==== 📊 Final Results ====\n');
fprintf('✅ Overall Success Rate for Rome: %.2f%%\n', mean(SuccessRate(1, :)) / Num_Packets * 100);
fprintf('✅ Overall Success Rate for Milan: %.2f%%\n', mean(SuccessRate(2, :)) / Num_Packets * 100);
fprintf('📡 Total Packets Successfully Received by NodeRM: %d\n', sum(Received_Packets_NodeRM));

toc;
