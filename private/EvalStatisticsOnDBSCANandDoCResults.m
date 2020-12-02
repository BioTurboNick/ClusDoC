function EvalStatisticsOnDBSCANandDoCResults(ClusterSmoothTableCh, Ch, outputFolder)

% Tabulate stats and output excel file for Clus-DoC analysis
%% Pull NbThresh from handles

handles = guidata(findobj('Tag', 'PALM GUI'));


%% Density Comparison

[row, column]=size(ClusterSmoothTableCh);

           MeanDensityDofC1=cell(row, column);
           MeanAreaDofC1=cell(row, column);
           MeanCircularityDofC1=cell(row, column);
           
           MeanDensityDofC2=cell(row, column);
           MeanAreaDofC2=cell(row, column);
           MeanCircularityDofC2=cell(row, column);
           
           MeanDensityDofC12=cell(row, column);
           MeanAreaDofC12=cell(row, column);
           MeanCircularityDofC12=cell(row, column);
           
           MeanDensity2=cell(row, column);
           MeanArea2=cell(row, column);
           MeanCircularity2=cell(row, column);
           
           MeanDensity3=cell(row, column);
           MeanArea3=cell(row, column);
           MeanCircularity3=cell(row, column);
           
           MeanNumMolsPerColoc1Cluster = cell(row, column);
           NumColoc1ClustersPerROI = cell(row, column);
           MeanNumMolsPerColoc2Cluster = cell(row, column);
           NumColoc2ClustersPerROI = cell(row, column);
           MeanNumMolsPerColoc12Cluster = cell(row, column);
           NumColoc12ClustersPerROI = cell(row, column);
           MeanNumMolsPerNonColocCluster = cell(row, column);
           NumNonColocClustersPerROI = cell(row, column);
           
           
           
for i=1:column
    
    for j=1:row
        A=ClusterSmoothTableCh{j,i};

        if ~isempty(A)
            
           % Population of cluster with Nb > NbThresh
           AA=cellfun(@(x) x(x.Nb > handles.DoC.NbThresh), A,'UniformOUtput',0);
           A=A(~cellfun('isempty', AA));
           
           % Cluster with Nb(Dof>0.4) > NbThresh      
           Cluster_DofC1=cellfun(@(x) x(x.Nb_In1 > handles.DoC.NbThresh && x.Nb_In2 <= handles.DoC.NbThresh), A,'UniformOUtput',0);
           Cluster_DofC1=Cluster_DofC1(~cellfun('isempty', Cluster_DofC1));
           
           Cluster_DofC2=cellfun(@(x) x(x.Nb_In2 > handles.DoC.NbThresh && x.Nb_In1 <= handles.DoC.NbThresh), A,'UniformOUtput',0);
           Cluster_DofC2=Cluster_DofC2(~cellfun('isempty', Cluster_DofC2));
           
           Cluster_DofC12=cellfun(@(x) x(x.Nb_In1 > handles.DoC.NbThresh && x.Nb_In2 > handles.DoC.NbThresh), A,'UniformOUtput',0);
           Cluster_DofC12=Cluster_DofC12(~cellfun('isempty', Cluster_DofC12));
           
           
           % Area and realtive density for Nb(Dof>0.4) >10
           DensityDofC1=cellfun(@(x) x.AvRelativeDensity20, Cluster_DofC1);
           AreaDofC1=cellfun(@(x) x.Area, Cluster_DofC1);
           CircularityDofC1=cellfun(@(x) x.Circularity, Cluster_DofC1);
           MeanDensityDofC1{j,i}=mean(DensityDofC1);
           MeanAreaDofC1{j,i}=mean(AreaDofC1);
           MeanCircularityDofC1{j,i}=mean(CircularityDofC1);
           
           NumMoleculesPerColoc1Cluster = cellfun(@(x) x.Nb, Cluster_DofC1);
           
           MeanNumMolsPerColoc1Cluster{j,i} = mean(NumMoleculesPerColoc1Cluster);
           NumColoc1ClustersPerROI{j,i} = numel(NumMoleculesPerColoc1Cluster);
           
           
           DensityDofC2=cellfun(@(x) x.AvRelativeDensity20, Cluster_DofC2);
           AreaDofC2=cellfun(@(x) x.Area, Cluster_DofC2);
           CircularityDofC2=cellfun(@(x) x.Circularity, Cluster_DofC2);
           MeanDensityDofC2{j,i}=mean(DensityDofC2);
           MeanAreaDofC2{j,i}=mean(AreaDofC2);
           MeanCircularityDofC2{j,i}=mean(CircularityDofC2);
           
           NumMoleculesPerColoc2Cluster = cellfun(@(x) x.Nb, Cluster_DofC2);
           
           MeanNumMolsPerColoc2Cluster{j,i} = mean(NumMoleculesPerColoc2Cluster);
           NumColoc2ClustersPerROI{j,i} = numel(NumMoleculesPerColoc2Cluster);
           
           
           DensityDofC12=cellfun(@(x) x.AvRelativeDensity20, Cluster_DofC12);
           AreaDofC12=cellfun(@(x) x.Area, Cluster_DofC12);
           CircularityDofC12=cellfun(@(x) x.Circularity, Cluster_DofC12);
           MeanDensityDofC12{j,i}=mean(DensityDofC12);
           MeanAreaDofC12{j,i}=mean(AreaDofC12);
           MeanCircularityDofC12{j,i}=mean(CircularityDofC12);
           
           NumMoleculesPerColoc12Cluster = cellfun(@(x) x.Nb, Cluster_DofC12);
           
           MeanNumMolsPerColoc12Cluster{j,i} = mean(NumMoleculesPerColoc12Cluster);
           NumColoc12ClustersPerROI{j,i} = numel(NumMoleculesPerColoc12Cluster);
           
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
           % Cluster with Nb(Dof>0.4) < NbThresh
           Cluster_Other=cellfun(@(x) x(x.Nb_In1 <= handles.DoC.NbThresh && x.Nb_In2 <= handles.DoC.NbThresh), A,'UniformOUtput',0);
           Cluster_Other=Cluster_Other(~cellfun('isempty', Cluster_Other));
                      
           Density2=cellfun(@(x) x.AvRelativeDensity20, Cluster_Other);
           Area2=cellfun(@(x) x.Area, Cluster_Other);
           Circularity2=cellfun(@(x) x.Circularity, Cluster_Other);
           MeanDensity2{j,i}=mean(Density2);
           MeanArea2{j,i}=mean(Area2);
           MeanCircularity2{j,i}=mean(Circularity2);
           
           NumMoleculesPerNonColocCluster = cellfun(@(x) x.Nb, Cluster_Other);
           
           MeanNumMolsPerNonColocCluster{j,i} = mean(NumMoleculesPerNonColocCluster);
           NumNonColocClustersPerROI{j,i} = numel(NumMoleculesPerNonColocCluster);
           
           %MeanDensity22{j,i}=mean(Density2.*WNb2);
           %MeanArea22{j,i}=mean(Area2.*WNb2);
           
           % Population of cluster with Nb<10  
           A=ClusterSmoothTableCh{j,i};
           AA=cellfun(@(x) x(x.Nb <= handles.DoC.NbThresh), A,'UniformOUtput',0);
           A=A(~cellfun('isempty', AA));
           
           Density3=cellfun(@(x) x.AvRelativeDensity20, A);
           Area3=cellfun(@(x) x.Area, A);
           Circularity3=cellfun(@(x) x.Circularity, A);

           MeanArea3{j,i}=mean(Area3);
           MeanDensity3{j,i}=mean(Density3);
           MeanCircularity3{j,i}=mean(Circularity3);
                      
        else
           MeanDensityDofC1{j,i}=[];
           MeanAreaDofC1{j,i}=[];
           MeanCircularityDofC1{j,i}=[];
           MeanDensityDofC2{j,i}=[];
           MeanAreaDofC2{j,i}=[];
           MeanCircularityDofC2{j,i}=[];
           MeanDensityDofC12{j,i}=[];
           MeanAreaDofC12{j,i}=[];
           MeanCircularityDofC12{j,i}=[];
           
           MeanDensity2{j,i}=[];
           MeanArea2{j,i}=[];
           MeanCircularity2{j,i}=[];
                
        end        
    end
    
end

Result2.DensityDofC1=MeanDensityDofC1;
Result2.AreaDofC1=MeanAreaDofC1;
Result2.CircularityDofC1=MeanCircularityDofC1;
Result2.DensityDofC2=MeanDensityDofC2;
Result2.AreaDofC2=MeanAreaDofC2;
Result2.CircularityDofC2=MeanCircularityDofC2;
Result2.DensityDofC12=MeanDensityDofC12;
Result2.AreaDofC12=MeanAreaDofC12;
Result2.CircularityDofC12=MeanCircularityDofC12;

Result2.Density2=MeanDensity2;
Result2.Area2=MeanArea2;
Result2.Circularity2=MeanCircularity2;

Result2.Density3=MeanDensity3;
Result2.Area3=MeanArea3;
Result2.Circularity3=MeanCircularity3;

switch Ch
    case 1
        ResultCh1=Result2;    
        save(fullfile(outputFolder, 'ResultCh1.mat'),'ResultCh1')
    case 2
        ResultCh2=Result2;    
        save(fullfile(outputFolder, 'ResultCh2.mat'),'ResultCh2')
    case 3
        ResultCh3=Result2;    
        save(fullfile(outputFolder, 'ResultCh3.mat'),'ResultCh3')
end
    


%% Convert to excel file
    
    % Density Area Circularity for cluster with DofC>0.4
    DensityDofC1 = cell2mat(MeanDensityDofC1(:));
    AreaDofC1 = cell2mat(MeanAreaDofC1(:));
    CircularityDofC1 = cell2mat(MeanCircularityDofC1(:));
    DensityDofC2 = cell2mat(MeanDensityDofC2(:));
    AreaDofC2 = cell2mat(MeanAreaDofC2(:));
    CircularityDofC2 = cell2mat(MeanCircularityDofC2(:));
    DensityDofC12 = cell2mat(MeanDensityDofC12(:));
    AreaDofC12 = cell2mat(MeanAreaDofC12(:));
    CircularityDofC12 = cell2mat(MeanCircularityDofC12(:));
    
    % Density Area Circularity for cluster with DofC<0.4
    Density2 = cell2mat(MeanDensity2(:));
    Area2 = cell2mat(MeanArea2(:));
    Circularity2 = cell2mat(MeanCircularity2(:));

    
    % Density Area Circularity for cluster with Nb<NbThresh 
%     Density3 = cell2mat(MeanDensity3(:));
%     Area3 = cell2mat(MeanArea3(:));
%     Circularity3=cell2mat(MeanCircularity3(:));
    
     Array=[{'Rel density in colocalised 1 clusters'}, ...
         {'Rel density in colocalised 2 clusters'}, ...
         {'Rel density in colocalised 1-2 clusters'}, ...
         {'Rel density in non-colocalised clusters'},...
         {'Average area of colicalised 1 clusters (nm^2)'}, ...	
         {'Average area of colicalised 2 clusters (nm^2)'}, ...	
         {'Average area of colicalised 1-2 clusters (nm^2)'}, ...	
         {'Average area of non-colicalised clusters (nm^2)'}, ...
         {'Circularity of colocalised 1 clusters'}, ...
         {'Circularity of colocalised 2 clusters'}, ...
         {'Circularity of colocalised 1-2 clusters'}, ...
         {'Circularity of non-colicalised clusters'}, ...
         {'Mean number of molecules per colocalised 1 cluster'}, ...
         {'Mean number of molecules per colocalised 2 cluster'}, ...
         {'Mean number of molecules per colocalised 1-2 cluster'}, ...
         {'Mean number of molecules per non-colocalised cluster'}, ...
         {'Mean number of colocalised 1 clusters per ROI'}, ...
         {'Mean number of colocalised 2 clusters per ROI'}, ...
         {'Mean number of colocalised 1-2 clusters per ROI'}, ...
         {'Mean number of non-colocalised clusters per ROI'}];
    
    Matrix_Result=[DensityDofC1...
                    DensityDofC2...
                    DensityDofC12...
                    Density2...
                    AreaDofC1...
                    AreaDofC2...
                    AreaDofC12...
                    Area2...
                    CircularityDofC1...
                    CircularityDofC2...
                    CircularityDofC12...
                    Circularity2, ...
                    cell2mat(MeanNumMolsPerColoc1Cluster(:)), ...
                    cell2mat(MeanNumMolsPerColoc2Cluster(:)), ...
                    cell2mat(MeanNumMolsPerColoc12Cluster(:)), ...
                    cell2mat(MeanNumMolsPerNonColocCluster(:)), ...
                    cell2mat(NumColoc1ClustersPerROI(:)), ...
                    cell2mat(NumColoc2ClustersPerROI(:)), ...
                    cell2mat(NumColoc12ClustersPerROI(:)), ...
                    cell2mat(NumNonColocClustersPerROI(:))];
    
    RegionName = strcat('Clus-DoC results');
    
 switch Ch
    case 1
        xlswrite(fullfile(outputFolder, 'Clus-DoC Ch1.xls'), Array, RegionName, 'A1');
        xlswrite(fullfile(outputFolder, 'Clus-DoC Ch1.xls'), Matrix_Result, RegionName, 'A2');
    case 2
        xlswrite(fullfile(outputFolder, 'Clus-DoC Ch2.xls'), Array, RegionName, 'A1');
        xlswrite(fullfile(outputFolder, 'Clus-DoC Ch2.xls'), Matrix_Result, RegionName, 'A2');
    case 3
        xlswrite(fullfile(outputFolder, 'Clus-DoC Ch3.xls'), Array, RegionName, 'A1');
        xlswrite(fullfile(outputFolder, 'Clus-DoC Ch3.xls'), Matrix_Result, RegionName, 'A2');
 end


    
end
    
    

    

   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    