%% Draw the Map and Stops
% figure;
%
% load('usborder.mat','x','y','xx','yy');
% rng(3,'twister') % makes a plot with stops in Maine & Florida, and is reproducible
m=5; %number of salemen/vehicles

DATA= readmatrix('Dual Litter Bins_Tempe_LatLong_Distance_Time.xlsx','Sheet','Distance','Range','A1:ID238');
[Longitude, Latitude] = readvars('Dual Litter Bins_Tempe_LatLong_Distance_Time.xlsx','Sheet','Sheet2','Range','B6:C243');
nStops = length(DATA); % you can use any number, but the problem size scales as N^2
% stopsLon = zeros(nStops,1); % allocate x-coordinates of nStops
% stopsLat = stopsLon; % allocate y-coordinates
% n = 1;
% while (n <= nStops)
%     xp = rand*1.5;
%     yp = rand;
%     if inpolygon(xp,yp,x,y) % test if inside the border
%         stopsLon(n) = xp;
%         stopsLat(n) = yp;
%         n = n+1;
%     end
% end
% plot(x,y,'Color','red'); % draw the outside border
% hold on
% % Add the stops to the map
% plot(stopsLon,stopsLat,'*b')
% hold off

%%% Calculate Distances Between Points
idxs = nchoosek(1:nStops,2);


% dist = hypot(stopsLat(idxs(:,1)) - stopsLat(idxs(:,2)), ...
%              stopsLon(idxs(:,1)) - stopsLon(idxs(:,2)));
dist=zeros(length(idxs),1);
for i=1:length(idxs)
    dist(i)=DATA(idxs(i,1),idxs(i,2));
end

lendist = length(dist);

%%% Create Variables and Problem
tsp = optimproblem;
trips = optimvar('trips',lendist,1,'Type','integer','LowerBound',0,'UpperBound',1);

tsp.Objective = dist'*trips;

%%% Equality Constraints
constrips = sum(trips) == nStops+m-1; % restricts the total number of trips needed for m round trip subtours
tsp.Constraints.constrips = constrips;

constrtrips = optimconstr(nStops,1);
for stops = 1:nStops
    if stops==1 % m vehicles leave and return to the first "depot" node
        whichIdxs = (idxs == stops);
        whichIdxs = any(whichIdxs,2); 
        constrtrips(stops) = sum(trips(whichIdxs)) == 2*m;
        
        
        
    else
        whichIdxs = (idxs == stops);
        whichIdxs = any(whichIdxs,2); % one trip to and from each DL node
        constrtrips(stops) = sum(trips(whichIdxs)) == 2;
    end
end
tsp.Constraints.constrtrips = constrtrips;


consorder=zeros(length(sum(trips)),2);


routes=zeros(length(sum(trips)),2); 
j=1;
for numpaths=1:m
    
    %while 
        
        for order=1:length(trips)
            if trips(order)==1
                routes(j,:)=[idxs(order,1),idxs(order,2)];
                j=j+1;
            end
        end
        
        %numnodes
        while routeX(count,2)~=1
            %for numnodes=1:length(DATA)
            if count==1
                [row,col,val]=find(routes==1);
                routeX(1,:)=routes(row(1),:);
            else
                [row,col,val]=find(routes==routeX(count-1,2);
                routeX(count,:)=routes(row(1),:);
                if routeX(count,1)~=routeX(count-1,2)
                    routeX(count,:)=[routeX(count,2),routeX(count,1)];
                end
                
            end
            count=count+1;
            
        end
    %end
 
   if numpaths==1
    route1=routeX;
    elseif numpaths==2
    route2=routeX;
    elseif numpaths==3
    route3=routeX;
    elseif numpaths==4
    route4=routeX;
    elseif numpaths==5
    route5=routeX;
   end
   
   [row,col,val]=find(routes==1);
   routes(row(1),:)=[];
   [row,col,val]=find(routes==routeX(length(routeX),1));
    if routes(row(1),1)==1
        routes(row(1),:)=[];
    else
        routes(row(2),:)=[];
    end
   numpaths=numpaths+1; 
end
tsp.Constraints.consorder = consorder;


%%% Solve the Initial Problem
opts = optimoptions('intlinprog','Display','off');
tspsol = solve(tsp,'options',opts)

% tspsol = struct with fields:
%     trips: [19900×1 double]



%% Visualize the Solution

segments = find(tspsol.trips); % Get indices of lines on optimal path

truetrips=idxs(segments',:);
% xplot=zeros(length(truetrips),1);
% yplot=zeros(length(truetrips),1);
figure;hold on;
for i=1:length(truetrips)
    xplot=[Latitude(truetrips(i,1)),Latitude(truetrips(i,2))];
    yplot=[Longitude(truetrips(i,1)),Longitude(truetrips(i,2))];
    plot(xplot,yplot,'b')
end
%%
hold on
plot(Latitude,Longitude,'r*')
hold off
%%
lh = zeros(nStops,1); % Use to store handles to lines on plot
lh = updateSalesmanPlot(lh,tspsol.trips,idxs,Longitude,Latitude);
title('Solution with Subtours');

%% Subtour Constraints
tours = detectSubtours(tspsol.trips,idxs);
numtours = length(tours); % number of subtours
fprintf('# of subtours: %d\n',numtours);

% Index of added constraints for subtours
k = 1;
while numtours > 1 % repeat until there is just one subtour
    % Add the subtour constraints
    for ii = 1:numtours
        subTourIdx = tours{ii}; % Extract the current subtour
%         The next lines find all of the variables associated with the
%         particular subtour, then add an inequality constraint to prohibit
%         that subtour and all subtours that use those stops.
        variations = nchoosek(1:length(subTourIdx),2);
        a = false(length(idxs),1);
        for jj = 1:length(variations)
            whichVar = (sum(idxs==subTourIdx(variations(jj,1)),2)) & ...
                       (sum(idxs==subTourIdx(variations(jj,2)),2));
            a = a | whichVar;
        end
        tsp.Constraints.(sprintf('subtourconstr%i',k)) = sum(trips(a)) <= length(subTourIdx)-1;
        k = k + 1;
    end
    % Try to optimize again
    [tspsol,fval,exitflag,output] = solve(tsp,'options',opts);

    % Visualize result
    lh = updateSalesmanPlot(lh,tspsol.trips,idxs,stopsLon,stopsLat);

    % How many subtours this time?
    tours = detectSubtours(tspsol.trips,idxs);
    numtours = length(tours); % number of subtours
    fprintf('# of subtours: %d\n',numtours);
end
