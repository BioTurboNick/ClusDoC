function [ClusterSmooth1, fig] = FunDBSCAN4DoC_GUIV2( Data_DoC,p,q,r,display1)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
%   Data = Zen format
% Routine for DBSCAN apply on Degree of colocalisation data 1
% 
        
% Calculate Lr for cumulated channels ch1
% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
% DBSCAN, plot of the contour, cluster centre identification
 % Parameters 1
%figure1=[];
fig=[];


% Data_Doc=table [X Y Lr Ch Density DofC D1_D2];
x=[Data_DoC.X Data_DoC.Y];
Density=Data_DoC.D1_D2;

if 1==display1
    fig=figure();
ax1=axes('parent',fig);
 hold on;
 plot( ax1,x(:,1), x(:,2),'Marker','.','MarkerSize',4,'LineStyle','none','color','blue');
 
 axis equal
 axis tight
end

%% Threshold for the DBSCAN_Nb_Neighbor on density
xsize=ceil (max(x(:,1))-min(x(:,1)));
ysize=ceil (max(x(:,2))-min(x(:,2)));
SizeROI=max(xsize,ysize);
AvDensity=length(x)/(xsize*ysize);
Nrandom=AvDensity*pi*r^2;

% if Nrandom<3
%     Nrandom=3;
% end

%% L(r) thresholding

if 1==0
    %% Calculate Lr for cumulated channels ch1 
        %  Threshold the Lr at Lr_Threshold
        
        Lr_Threshold=((SizeROI)^2*Nrandom/(length(x)-1)/pi).^0.5;
        [data, idx, Dis, Density ] = Lr_fun(x(:,1),x(:,2),x(:,1),x(:,2),r,SizeROI); % data=[X2 Y2 Lr Kfuncans];
        [dataT, idxT, DisT, Density_20 ] = Lr_fun(x(:,1),x(:,2),x(:,1),x(:,2),20,SizeROI); % data=[X2 Y2 Lr Kfuncans];
        % 
        data(:,5)=Data(:,12); % channel index
        data(:,6)=Density;
        data(:,7)=Density_20;
        idxthr=find(data(:,3)>Lr_Threshold); % Index corresponding to the Threshold the data with Lr
        datathr=data(idxthr,:);
        
%       Include the points count in the search radius r for a point of
%       interst

        if 1==1    
        idxthr_friends=idx(idxthr);
        idxthr_friends2=cell2mat(idxthr_friends');
        Index_thr=unique(idxthr_friends2);
        datathr=data(Index_thr,:);
        end
        
        xt(:,1)=datathr(:,1); xt(:,2)=datathr(:,2); % data threshold both channel
x=xt;
Density=datathr(:,6);
Density20=datathr(:,7);
size(Density)
size(Density20)
if display1==1
plot(ax1, x(:,1), x(:,2),'Marker','.','MarkerSize',4,'LineStyle','none','color','red');
end

end


%%
Data4dbscan=x(:,1:2);
DBSCAN_Radius=20;
DBSCAN_Nb_Neighbor=3;%ceil(Nrandom);
threads = 2;

% FAST DBSCAN CALL

class=t_dbscan(Data4dbscan(:,1), Data4dbscan(:, 2), DBSCAN_Nb_Neighbor, DBSCAN_Radius, threads);
Data4dbscan(:,3)=class;

SumofContourBig=[];
SumofContourSmall=[];
ClusterSmooth1=cell(max(class),1);

for i=1:max(class)
    % Smoothing 
    xin=Data4dbscan(Data4dbscan(:,3)==i,:); % Positions contained in the cluster i
    Data_DoCi=Data_DoC(Data4dbscan(:,3)==i,:); % Data_DoC for the cluster i
    
    if display1==1
    plot( ax1,xin(:,1), xin(:,2),'Marker','.','MarkerSize',4,'LineStyle','none','color','green');
    end

    [dataT, idxT, DisT, Density20 ] = Lr_fun(xin(:,1), xin(:,2),xin(:,1), xin(:,2),20,SizeROI); % data=[X2 Y2 Lr Kfuncans];
    %Cluster{i,1}.K=K(Data4dbscan(:,3)==i);
    % Smoothing function
    [ClusImage,  Area, Circularity, Nb, contour, edges, Cutoff_point]=Smoothing_fun4clusterV3_2(xin(:,1:2),0,0); % 0.1*max intensity
    
    %% Collect parameter
    ClusterSmooth1{i,1}.Data_DoCi=Data_DoCi;
    %ClusterSmooth1{i,1}.Image=ClusImage;
    ClusterSmooth1{i,1}.edges=edges;
    %ClusterSmooth1{i,1}.Cutoff_point=Cutoff_point;
    ClusterSmooth1{i,1}.Contour=contour;
    
    ClusterSmooth1{i,1}.TotalAreaDensity=AvDensity;
    
    % Cluster Density
    ClusterSmooth1{i,1}.Area=Area;
    ClusterSmooth1{i,1}.Circularity=Circularity;
    ClusterSmooth1{i,1}.Nb=Nb;
    ClusterSmooth1{i}.Density_Nb_A=Nb/Area;
    ClusterSmooth1{i}.RelativeDensity_Nb_A=Nb/Area/AvDensity;
    
    
    % Local Density
    ClusterSmooth1{i,1}.Density=Density(Data4dbscan(:,3)==i);
    ClusterSmooth1{i,1}.Density20=Density20;
    ClusterSmooth1{i,1}.RelativeDensity=Density(Data4dbscan(:,3)==i)/AvDensity;
    ClusterSmooth1{i,1}.RelativeDensity20=Density20/AvDensity;
    %ClusterSmooth1{i,1}.Mean_Density=mean(Density(Data4dbscan(:,3)==i));
    
    % Average density
    ClusterSmooth1{i,1}.AvRelativeDensity=mean(Density(Data4dbscan(:,3)==i)/AvDensity);
    ClusterSmooth1{i,1}.AvRelativeDensity20=mean(Density20/AvDensity);
    

    plot(ax1,contour(:,1),contour(:,2),'k');
    
    %% Clusters with DofC >0.4 
    
    
    Point_In=[Data_DoCi.X(Data_DoCi.DofC>=0.4) Data_DoCi.Y(Data_DoCi.DofC>=0.4)];
    Nb_In=size(Point_In,1);
    
    if Nb_In>1
    Density20_In=Density20(Data_DoCi.DofC>=0.4);
    [Contour_In]=Smoothing_fun4clusterV3_3(Point_In, 0,0);
    Area_In=polyarea(Contour_In(:,1),Contour_In(:,2));
    
    
    ClusterSmooth1{i,1}.Nb_In=Nb_In;
    ClusterSmooth1{i,1}.Area_In=Area_In;
    ClusterSmooth1{i,1}.AvRelativeDensity20_In=mean(Density20_In/AvDensity);
    
    plot(ax1,Contour_In(:,1),Contour_In(:,2),'r');
    
    DofCOut=Data_DoCi(Data_DoCi.DofC<0.4,:);
    Density20_Out=Density20(Data_DoCi.DofC<0.4);
    
    Nb_Out=length(Density20_Out);
    ClusterSmooth1{i,1}.Nb_Out=Nb_Out;
    ClusterSmooth1{i,1}.AvRelativeDensity20_Out=mean(Density20_Out/AvDensity);
    
    ClusterSmooth1{i,1}.DensityRatio=mean(Density20_In/AvDensity)/mean(Density20_Out/AvDensity);
    ClusterSmooth1{i,1}.Contour_In=Contour_In;
    
    else
        ClusterSmooth1{i,1}.Nb_In=Nb_In;
        ClusterSmooth1{i,1}.Area_In=0;
        ClusterSmooth1{i,1}.AvRelativeDensity20_In=0;
    
        ClusterSmooth1{i,1}.Nb_In=0;
        ClusterSmooth1{i,1}.Area_In=0;
        ClusterSmooth1{i,1}.AvRelativeDensity20_In=0;
    
        ClusterSmooth1{i,1}.DensityRatio=0;
        ClusterSmooth1{i,1}.Contour_In=0;
    end
    

end




