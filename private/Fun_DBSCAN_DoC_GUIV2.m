function [ClusterSmoothTableCh1, ClusterSmoothTableCh2]=Fun_DBSCAN_DoC_GUIV2(CellData, ROIPos, Data_DegColoc, Path_name, Chan1Color, Chan2Color)
% Routine to apply DBSCAN on the Degree of Colocalisation Result for


% Channel 1
% Ch1 hhhhhhh

%load(fullfile(Path_name,'Data_for_Cluster_Analysis.mat'))
%
%         if ~exist(strcat(Path_name, 'DBSCAN'),'dir')
%             mkdir(fullfile(Path_name, 'DBSCAN'));
%             mkdir(fullfile(Path_name, 'DBSCAN', 'Clus-DoC cluster maps', 'Ch1'));
%             mkdir(fullfile(Path_name, 'DBSCAN', 'Clus-DoC cluster maps', 'Ch2'));
%         end
%

% Parameters to change
%Cutoff=10; % cutoff of number of molecules per cluster
Display1=1; % Display and save Image from DBSCAN
Display2=0; % Display and save Image Cluster Density Map
r=20;       % r= raduis for the Lr_fun which calculates the local density

nCells = numel(unique(ROIPos(:,1)));
nROIs = numel(unique(ROIPos(:,2)));

ClusterSmoothTableCh1 = cell(nROIs, nCells);
ClusterSmoothTableCh2 = cell(nROIs, nCells);

ResultCell = cell(nROIs, nCells);


for p = 1:nCells % index for the cell
    for q = 1:nROIs % index for the region
        
        CurrentROI = ROIPos((ROIPos(:,1) == p) & (ROIPos(:,2) == q), :);
        [~, Index_In] = Cropping_Fun(CellData{p}(:,5:6), CurrentROI(3:6));
        Data = CellData{p}(Index_In,:);
        
        
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
                
                dbscanParams.epsilon = 20;
                dbscanParams.minPts = 3;
                dbscanParams.threads = 2;
                dbscanParams.UseLr_rThresh = 0;
                dbscanParams.Lr_rThreshRad = 20;
                dbscanParams.SmoothingRad = 15;
                dbscanParams.Cutoff = 10;
                dbscanParams.DoStats = 1;
                dbscanParams.CurrentChannel = Ch;
                dbscanParams.Outputfolder = Path_name;
                dbscanParams.DoCThreshold = 0.4;
                
                if Ch == 1
                    clusterColor = Chan1Color;
                elseif Ch == 2
                    clusterColor = Chan2Color;
                end
                
                
                % DBSCAN_Radius=20 - epsilon
                % DBSCAN_Nb_Neighbor=3; - minPts ;
                % threads = 2
                
                [~, ClusterCh, ~, fig, ~, ~, ResultCell{p,q}] = DBSCANHandler([Data_DoC1.X Data_DoC1.Y], ...
                    dbscanParams, q, p, true, false, clusterColor, Data_DoC1.D1_D2, Data_DoC.DoC);
                
                %                     [ClusterCh, fig] = FunDBSCAN4DoC_GUIV2(Data_DoC1, p, q, r, Display1);
                
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
                
                %                         Name1 = sprintf('_Table_%d_Region_%d_', p, q);
                %                         Name2 = fullfile(Path_name, 'DBSCAN Results', 'Clus-DoC cluster maps', ...
                %                             sprintf('Ch%d', Ch), sprintf('%sClusters_Ch%d.tif', Name1, Ch));
                %
                %                         set(gca, 'box', 'on', 'XTickLabel', [], 'XTick', [], 'YTickLabel', [], 'YTick', [])
                %                         set(fig, 'Color', [1 1 1])
                %                         tt = getframe(fig);
                %                         imwrite(tt.cdata, Name2)
                %                         close(gcf)
                
            end
        end
    end
end

ExportDBSCANDataToExcelFiles(ROIPos, ResultCell, strcat(Path_name, '\DBSCAN Results'));

save(fullfile(Path_name, 'DBSCAN Clus-DoC Results.mat'),'ClusterSmoothTableCh1','ClusterSmoothTableCh2');
end

