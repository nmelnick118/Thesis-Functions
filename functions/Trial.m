% data holds the pole values
data = [];
% time holds the timestamps
time = [];
% strain holds the strain info
strain = [];
%%
% Start reading in files
cd /Users/Nitai/Desktop/'Final Data Sets'/'TimeCourse'/CB15N;
files = dir('*pill_MESH.mat');
for i = 1:length(files)
    cells = load(files(i).name);
    C = strsplit(files(i).name, {'t', '_'});
    timestamp = char(C(3));
    timestamp = str2double(timestamp);
    all_cells = cells.frame.object;
    % max_min = TwotoCenter(cells, 0);
    % max_min = TwoPointWidth(cells);
    max_min = AreaPerimRatio(cells); % (larger value is LESS pointy!)
    max_min = transpose(max_min);
    
    temptime = [];
    for i = 1:length(max_min)
        temptime(i) = timestamp;
        strain = [strain, C(1)];
    end
    time = [time, temptime];
    data = vertcat(data, max_min);
end

%%

cd /Users/Nitai/Desktop/'Final Data Sets'/'TimeCourse'/43;
files = dir('*pill_MESH.mat');

for i = 1:length(files)
    cells = load(files(i).name);
    C = strsplit(files(i).name, {'t', '_'});
    timestamp = char(C(3));
    timestamp = str2double(timestamp);
    % max_min = TwotoCenter(cells, 0);
    
    max_min = AreaPerimRatio(cells); % (larger value is LESS pointy!)
    max_min = transpose(max_min);
    
    temptime = [];
    for i = 1:length(max_min)
        temptime(i) = timestamp;
        strain = [strain, C(1)];
    end
    time = [time, temptime];
    data = vertcat(data, max_min);
end

%%
cd /Users/Nitai/Desktop/'Final Data Sets'/'TimeCourse'/1218;
files = dir('*pill_MESH.mat');

for i = 1:length(files)
    cells = load(files(i).name);
    C = strsplit(files(i).name, {'t', '_'});
    timestamp = char(C(3));
    timestamp = str2double(timestamp);
    % max_min = TwotoCenter(cells, 0);
    
    max_min = AreaPerimRatio(cells); % (larger value is LESS pointy!)
    max_min = transpose(max_min);
    
    temptime = [];
    for i = 1:length(max_min)
        temptime(i) = timestamp;
        strain = [strain, C(1)];
    end
    time = [time, temptime];
    data = vertcat(data, max_min);
end

%%
cd /Users/Nitai/Desktop/'Final Data Sets'/'TimeCourse'/3205;
files = dir('*pill_MESH.mat');

for i = 1:length(files)
    cells = load(files(i).name);
    C = strsplit(files(i).name, {'t', '_'});
    timestamp = char(C(3));
    timestamp = str2double(timestamp);
    % max_min = TwotoCenter(cells, 0);
    
    max_min = AreaPerimRatio(cells); % (larger value is LESS pointy!)
    max_min = transpose(max_min);
    
    temptime = [];
    for i = 1:length(max_min)
        temptime(i) = timestamp;
        strain = [strain, C(1)];
    end
    time = [time, temptime];
    data = vertcat(data, max_min);
end

%%
cd /Users/Nitai/Desktop/'Final Data Sets'/'TimeCourse'/801;
files = dir('*pill_MESH.mat');

for i = 1:length(files)
    cells = load(files(i).name);
    C = strsplit(files(i).name, {'t', '_'});
    timestamp = char(C(3));
    timestamp = str2double(timestamp);
    % max_min = TwotoCenter(cells, 0);
    
    max_min = AreaPerimRatio(cells); % (larger value is LESS pointy!)
    max_min = transpose(max_min);
    
    temptime = [];
    for i = 1:length(max_min)
        temptime(i) = timestamp;
        strain = [strain, C(1)];
    end
    time = [time, temptime];
    data = vertcat(data, max_min);
end

%%
cd /Users/Nitai/Desktop/'Final Data Sets'/'TimeCourse'/858;
files = dir('*pill_MESH.mat');

for i = 1:length(files)
    cells = load(files(i).name);
    C = strsplit(files(i).name, {'t', '_'});
    timestamp = char(C(3));
    timestamp = str2double(timestamp);
    % max_min = TwotoCenter(cells, 0);
    
    max_min = AreaPerimRatio(cells); % (larger value is LESS pointy!)
    max_min = transpose(max_min);
    
    temptime = [];
    for i = 1:length(max_min)
        temptime(i) = timestamp;
        strain = [strain, C(1)];
    end
    time = [time, temptime];
    data = vertcat(data, max_min);
end

%%

%{
% cd /Users/Nitai/Desktop/TimeCourses/2012;
files = dir('*pill_MESH.mat');

for i = 1:length(files)
    cells = load(files(i).name);
    C = strsplit(files(i).name, {'t', '_'});
    timestamp = char(C(3));
    timestamp = str2double(timestamp);
    % max_min = TwotoCenter(cells, 0);
    
    max_min = AreaPerimRatio(cells); % (larger value is LESS pointy!)
    max_min = transpose(max_min);
    
    temptime = [];
    for i = 1:length(max_min)
        temptime(i) = timestamp;
        strain = [strain, C(1)];
    end
    time = [time, temptime];
    data = vertcat(data, max_min);
end

%%
% cd /Users/Nitai/Desktop/TimeCourses/801;
files = dir('*pill_MESH.mat');

for i = 1:length(files)
    cells = load(files(i).name);
    C = strsplit(files(i).name, {'t', '_'});
    timestamp = char(C(3));
    timestamp = str2double(timestamp);
    % max_min = TwotoCenter(cells, 0);
    
    max_min = AreaPerimRatio(cells); % (larger value is LESS pointy!)
    max_min = transpose(max_min);
    
    temptime = [];
    for i = 1:length(max_min)
        temptime(i) = timestamp;
        strain = [strain, C(1)];
    end
    time = [time, temptime];
    data = vertcat(data, max_min);
end

%%
% cd /Users/Nitai/Desktop/TimeCourses/858;
files = dir('*pill_MESH.mat');

for i = 1:length(files)
    cells = load(files(i).name);
    C = strsplit(files(i).name, {'t', '_'});
    timestamp = char(C(3));
    timestamp = str2double(timestamp);
    % max_min = TwotoCenter(cells, 0);
    
    max_min = AreaPerimRatio(cells); % (larger value is LESS pointy!)
    max_min = transpose(max_min);
    
    temptime = [];
    for i = 1:length(max_min)
        temptime(i) = timestamp;
        strain = [strain, C(1)];
    end
    time = [time, temptime];
    data = vertcat(data, max_min);
end
%}

time = transpose(time);
strain = transpose(strain);
data = horzcat(data, time);
table = array2table(data);
table.FourthColumn = strain;
table.Properties.VariableNames = {'SmallPole', 'LargePole', 'Time', 'Strain'};

% name = sprintf('%s.txt', char(C(1)));
name = 'timeCourseData.txt';
cd /Users/Nitai/Desktop/'Final Data Sets'/'TimeCourse'/tables;
writetable(table, name);




