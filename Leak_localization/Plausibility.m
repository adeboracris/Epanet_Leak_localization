function [Pl_n,di] = Plausibility(elec)
%{
 Criar novo PL
%}
load('Graph_s.mat')

%select sensors
elec_sensors{1}= elec;

%Graph
weights=(1.21216*10^10*length1)./(120^(1.852)*Diameter.^(4.87)); %length/diameter^5 [m][m]
G = graph(S,T,weights);

%1 Nodes between the resevoir to the sensor node
P=[];
for i=1:length(elec_sensors{1,1})
    P{i}=shortestpath(G,32,elec_sensors{1,1}(i));
end

%2 Analyze for each leak node hypothesis
node_l=[];
for i=1:31
    node_l{i}=shortestpath(G,32,i);
end

%3 Create plausability
for i=1:31
    maxi=[];
    for j=1:length(elec_sensors{1,1})
        maxi(j)= sum(ismember(P{1,j},node_l{1,i}));
        [~,IB]= ismember(P{1,j},node_l{1,i});
        Cp=node_l{1,i}(IB(2:maxi(j)));
        Pl(j,i) = distances(G,32,Cp(maxi(j)-1));
    end
end

% weights=length1; %length/diameter^5 [m][m]
% G = graph(S,T,weights);
di=[];
for j=1:31
    for i=1:length(elec_sensors{1,1})
        di(i,j)=1/distances(G,j,elec_sensors{1,1}(i));
    end
end
di(di==inf)=0;
di=di./sum(di);
di(di==0)=1;

% Pl_n= Pl;
Pl_n= Pl.*di;
% Pl_n=1./Pl_n;
Pl_n=Pl_n./sum(Pl_n);
end

