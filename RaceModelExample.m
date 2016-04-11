% Race model script test/example
% http://matlaboratory.blogspot.co.uk/2015/03/race-model-inequality-script-r-ulrich-j.html

%% Generate data
% 3 columns, in Seconds

% A, V, AV
n = 50;
RTs = rand(n,3);
RTs(:,1) = RTs(:,1)+0.5;
RTs(:,2) = RTs(:,2)+0.5;
RTs(:,3) = RTs(:,3);

disp(mean(RTs))

%% Plot generated data

figure
x = (1:n)'; % Xaxis
scatter(x,RTs(:,1))
hold on
scatter(x,RTs(:,2))
scatter(x,RTs(:,3))
xlabel('n')
ylabel('RT, s')

%% Test race model

% Round and convert RT data to ms
RTs = round(RTs*1000);

% Probabilites to test at
P = 0.02:0.02:1;

% Transpose data to rows
X = RTs(:,1)';
Y = RTs(:,2)';
Z = RTs(:,3)';

% Run
figure
[Xp, Yp, Zp, Bp] = RaceModel(X, Y, Z, P, true);
