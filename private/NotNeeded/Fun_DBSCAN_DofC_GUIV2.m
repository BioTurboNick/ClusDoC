function [ClusterSmoothTableCh1, ClusterSmoothTableCh2]=Fun_DBSCAN_DofC_GUIV2(ROIData, Data_DegColoc, Path_name)
% Routine to apply DBSCAN on the Degree of Colocalisation Result for


% Channel 1
% Ch1 hhhhhhh

        %load(fullfile(Path_name,'Data_for_Cluster_Analysis.mat'))
        
        if ~exist(strcat(Path_name, 'DBSCAN'),'dir')
            mkdir(fullfile(Path_name, 'DBSCAN'));
            mkdir(fullfile(Path_name, 'DBSCAN', 'Clus-DoC cluster maps', 'Ch1'));
            mkdir(fullfile(Path_name, 'DBSCAN', 'Clus-DoC cluster maps', 'Ch2'));
        end
        
        
        % Parameters to change
        %Cutoff=10; % cutoff of number of molecules per cluster
        Display1=1; % Display and save Image from DBSCAN 
        Display2=0; % Display and save Image Cluster Density Map
        r=20;       % r= raduis for the Lr_fun which calculates the local density

        ClusterSmoothTableCh1=cell(size(ROIData));
        ClusterSmoothTableCh2=cell(size(ROIData));
        
for p=1:size(ROIData,2) % index for the cell
         for q=1:size(ROIData,1) % index for the region
             
             Data=ROIData{q,p};
             if ~isempty(Data)
                 
             Data_DoC=Data_DegColoc{q,p};
                 for Ch=1:2
                     
                    Data_DoC1=Data_DoC(Data_DoC.Ch==Ch,:);            
                    % FunDBSCAN4ZEN( Data,p,q,2,r,Cutoff,Display1, Display2 )
                    % Input :
                    % -Data = Data Zen table (12 or 13 columns). 
                    % -p = index for cell
                    % -q = index for Region
                    % -Cutoff= cutoff for the min number of molecules per cluster 
                    % -Display1=1; % Display and save Image from DBSCAN 
                    % -Display2=0; % Display and save Image Cluster Density Map

                    %[ClusterSmooth2, fig,fig2,fig3] = FunDBSCAN4ZEN_V3( Data_DoC1,p,q,r,Display1,Display2);
                    [ClusterCh, fig] = FunDBSCAN4DofC_GUIV2( Data_DoC1,p,q,r,Display1);

                    % Output :
                    % -Datathr : Data after thresholding with lr_Fun +
                    % randommess Criterion
                    % -ClusterSmooth : cell/structure with individual cluster
                    % pareameter (Points, area, Nb of position Contour....)
                    % -SumofContour : cell with big(>Cutoff) and small(<Cutoff)
                    % contours. use to draw quickly all the contours at once
                    % Fig1, fig2 fig3 : handle for figures plot in the
                    % function.

                    % Save the plot and data
                    switch Ch
                        case 1
                        
                            ClusterSmoothTableCh1{q,p}=ClusterCh;
                        
                        case 2

                            ClusterSmoothTableCh2{q,p}=ClusterCh;
                    end  
                    
                    Name1 = fullfile(Path_name, 'DBSCAN', sprintf('Table_%d_Region_%d_', p, q));
                    Name2 = fullfile(Path_name, 'DBSCAN', 'Clus-DoC cluster maps', ...
                        sprintf('Ch%d', Ch), sprintf('%sClusters_Ch%d.tif', Name1, Ch));
                    
                    set(gca, 'box', 'on', 'XTickLabel', [], 'XTick', [], 'YTickLabel', [], 'YTick', [])
                    set(fig, 'Color', [1 1 1])
                    tt = getframe(fig);
                    imwrite(tt.cdata,Name2 )
                    close(fig)

                 end
            end 
        end
    end
    
     save('DBSCAN for DoC Ch1 Ch2','ClusterSmoothTableCh1','ClusterSmoothTableCh2');
end

