function [CellData, DensityROI] = DoCHandler(ROICoordinates, CellData, Lr_rad, Rmax, Step, Chan1Color, Chan2Color, Chan3Color, Outputfolder, NDatacolumns)

    % Handler function for data arrangement and plotting of DoC calculations
    % Formerly Main_Fun_DoC_GUIV2.m
    % Follows idea of:
    % Colocalization Malkusch (Histochem Cell Biol)
    % http://www.ncbi.nlm.nih.gov/pubmed/22086768
    % Moved here to avoid conflicts with private functions and data
    % handling.
    % Away from button callback for clarity.

    Data_DoC = cell(max(cellfun(@length, ROICoordinates)), length(CellData));
    DensityROI = cell(max(cellfun(@length, ROICoordinates)), length(CellData));

    for cellIter = 1:length(CellData) % cell number

        for roiIter = 1:numel(ROICoordinates{cellIter}) % ROI number

            roiHere = ROICoordinates{cellIter}{roiIter};

            % Since which ROI a point falls in is encoded in binary, decode here
            whichPointsInROI = fliplr(dec2bin(CellData{cellIter}(:,NDatacolumns + 1)));
            whichPointsInROI = whichPointsInROI(:, roiIter) == '1';
    
            dataCropped = CellData{cellIter}(whichPointsInROI, :);

            if ~isempty(dataCropped)
                
                dataCropped(isnan(dataCropped(:,12)),:)=[];

                % Send single ROI of data to DoC calculation function
                [ Data_DegColoc1, SizeROI1 ] = DoCCalc( dataCropped, Lr_rad, Rmax, Step, roiHere );

                CA12 = Data_DegColoc1.DoC1((Data_DegColoc1.Ch == 1) & (Data_DegColoc1.Lr_rAboveThresh == 1)); % Ch1 -> Ch2
                CA13 = Data_DegColoc1.DoC2((Data_DegColoc1.Ch == 1) & (Data_DegColoc1.Lr_rAboveThresh == 1)); % Ch1 -> Ch3
                CA21 = Data_DegColoc1.DoC1((Data_DegColoc1.Ch == 2) & (Data_DegColoc1.Lr_rAboveThresh == 1)); % Ch2 -> Ch1
                CA23 = Data_DegColoc1.DoC2((Data_DegColoc1.Ch == 2) & (Data_DegColoc1.Lr_rAboveThresh == 1)); % Ch2 -> Ch3
                CA31 = Data_DegColoc1.DoC1((Data_DegColoc1.Ch == 3) & (Data_DegColoc1.Lr_rAboveThresh == 1)); % Ch3 -> Ch1
                CA32 = Data_DegColoc1.DoC2((Data_DegColoc1.Ch == 3) & (Data_DegColoc1.Lr_rAboveThresh == 1)); % Ch3 -> Ch2

                handles.handles.DoCFigPerROI = figure('color', [1 1 1], 'inverthardcopy', 'off');
                
                

                handles.handles.DoCAxPerROI(1) = subplot(6,1,1);
                histHand = histogram(handles.handles.DoCAxPerROI(1), CA12, 100);
                histHand.FaceColor = Chan1Color;
                histHand.EdgeColor = rgb(52, 73, 94);
                set(handles.handles.DoCAxPerROI(1), 'XLim', [-1 1]);
                xlabel(handles.handles.DoCAxPerROI(1), 'DoC Score Ch1 -> Ch2', 'Fontsize', 12);
                ylabel(handles.handles.DoCAxPerROI(1), 'Frequency','FontSize',12);
                set(handles.handles.DoCAxPerROI(1),'FontSize',12)
                
                handles.handles.DoCAxPerROI(1) = subplot(6,1,2);
                histHand = histogram(handles.handles.DoCAxPerROI(1), CA13, 100);
                histHand.FaceColor = Chan1Color;
                histHand.EdgeColor = rgb(52, 73, 94);
                set(handles.handles.DoCAxPerROI(1), 'XLim', [-1 1]);
                xlabel(handles.handles.DoCAxPerROI(1), 'DoC Score Ch1 -> Ch3', 'Fontsize', 12);
                ylabel(handles.handles.DoCAxPerROI(1), 'Frequency','FontSize',12);
                set(handles.handles.DoCAxPerROI(1),'FontSize',12)

                handles.handles.DoCAxPerROI(2) = subplot(6,1,3);
                histHand = histogram(handles.handles.DoCAxPerROI(2), CA21, 100);
                histHand.FaceColor = Chan2Color;
                histHand.EdgeColor = rgb(52, 73, 94);
                set(handles.handles.DoCAxPerROI(2), 'XLim', [-1 1]);
                xlabel(handles.handles.DoCAxPerROI(2), 'DoC Score Ch2 -> Ch1', 'Fontsize', 12);
                ylabel(handles.handles.DoCAxPerROI(2), 'Frequency','FontSize',12);
                set(handles.handles.DoCAxPerROI(2),'FontSize',12)
                
                handles.handles.DoCAxPerROI(2) = subplot(6,1,4);
                histHand = histogram(handles.handles.DoCAxPerROI(2), CA23, 100);
                histHand.FaceColor = Chan2Color;
                histHand.EdgeColor = rgb(52, 73, 94);
                set(handles.handles.DoCAxPerROI(2), 'XLim', [-1 1]);
                xlabel(handles.handles.DoCAxPerROI(2), 'DoC Score Ch2 -> Ch3', 'Fontsize', 12);
                ylabel(handles.handles.DoCAxPerROI(2), 'Frequency','FontSize',12);
                set(handles.handles.DoCAxPerROI(2),'FontSize',12)
                
                handles.handles.DoCAxPerROI(3) = subplot(6,1,5);
                histHand = histogram(handles.handles.DoCAxPerROI(3), CA31, 100);
                histHand.FaceColor = Chan3Color;
                histHand.EdgeColor = rgb(52, 73, 94);
                set(handles.handles.DoCAxPerROI(3), 'XLim', [-1 1]);
                xlabel(handles.handles.DoCAxPerROI(3), 'DoC Score Ch3 -> Ch1', 'Fontsize', 12);
                ylabel(handles.handles.DoCAxPerROI(3), 'Frequency','FontSize',12);
                set(handles.handles.DoCAxPerROI(3),'FontSize',12)
                
                handles.handles.DoCAxPerROI(3) = subplot(6,1,6);
                histHand = histogram(handles.handles.DoCAxPerROI(3), CA32, 100);
                histHand.FaceColor = Chan3Color;
                histHand.EdgeColor = rgb(52, 73, 94);
                set(handles.handles.DoCAxPerROI(3), 'XLim', [-1 1]);
                xlabel(handles.handles.DoCAxPerROI(3), 'DoC Score Ch3 -> Ch2', 'Fontsize', 12);
                ylabel(handles.handles.DoCAxPerROI(3), 'Frequency','FontSize',12);
                set(handles.handles.DoCAxPerROI(3),'FontSize',12)
                
                set(handles.handles.DoCFigPerROI,'Position',[100 100 500 1000])
                
                drawnow;

                % Save the figure
                Name = sprintf('Table_%d_Region_%d_Hist', cellIter, roiIter);
                print(fullfile(Outputfolder, 'Clus-DoC Results', 'DoC histograms', Name), ...
                    handles.handles.DoCFigPerROI, '-dtiff');

                close gcf
                
                CellData{cellIter}(whichPointsInROI, NDatacolumns + 4) = Data_DegColoc1.DoC1;
                CellData{cellIter}(whichPointsInROI, NDatacolumns + 5) = Data_DegColoc1.DoC2;
                CellData{cellIter}(whichPointsInROI, NDatacolumns + 6) = Data_DegColoc1.Lr;
                CellData{cellIter}(whichPointsInROI, NDatacolumns + 7) = Data_DegColoc1.D1_D2;
                CellData{cellIter}(whichPointsInROI, NDatacolumns + 8) = Data_DegColoc1.Lr_rAboveThresh;
                CellData{cellIter}(whichPointsInROI, NDatacolumns + 9) = Data_DegColoc1.Density;
                
                Data_DoC{roiIter, cellIter} = Data_DegColoc1;


                % Average density for each region
                DensityROI{roiIter, cellIter} = [size([CA12;CA13;CA21;CA23;CA31;CA32],1)/SizeROI1^2, ...
                    size(CA12,1)/SizeROI1^2, ...
                    size(CA13,1)/SizeROI1^2, ...
                    size(CA21,1)/SizeROI1^2, ...
                    size(CA23,1)/SizeROI1^2, ...
                    size(CA31,1)/SizeROI1^2, ...
                    size(CA32,1)/SizeROI1^2];

            end
            %AvDensityCell(nt,:)=mean (DensityROI,1);
        end
    end

%     A = Data_DoC(~cellfun('isempty', Data_DoC(:)));
%     DoC1 = cell2mat(cellfun(@(x) x.DoC(x.Ch == 1), A(:), 'UniformOutput', false));
%     DoC2 = cell2mat(cellfun(@(x) x.DoC(x.Ch == 2), A(:), 'UniformOutput', false));
%     DoC1(DoC1==0) = [];
%     DoC2(DoC2==0) = [];
    
    % DoC1 and DoC2 are DoC scores for chan1->2 and chan2->1 for ALL
    % evaluated points above Lr_r threshold
    %  CumHist not updated for 3rd channel; I don't think it's used?
    DoC12CumHist = zeros(100, 1);
    DoC13CumHist = zeros(100, 1);
    DoC21CumHist = zeros(100, 1);
    DoC23CumHist = zeros(100, 1);
    DoC31CumHist = zeros(100, 1);
    DoC32CumHist = zeros(100, 1);
    
    assignin('base', 'CellData', CellData);
    
    for k = 1:numel(CellData)
        DoC12CumHist = DoC12CumHist + histc(CellData{k}((CellData{k}(:,NDatacolumns + 1) > 0) & (CellData{k}(:,12) == 1) & (CellData{k}(:, NDatacolumns + 8) == 1), NDatacolumns + 4), ...
            linspace(-1, 1, 100));
        DoC13CumHist = DoC13CumHist + histc(CellData{k}((CellData{k}(:,NDatacolumns + 1) > 0) & (CellData{k}(:,12) == 1) & (CellData{k}(:, NDatacolumns + 8) == 1), NDatacolumns + 5), ...
            linspace(-1, 1, 100));
        DoC21CumHist = DoC21CumHist + histc(CellData{k}((CellData{k}(:,NDatacolumns + 1) > 0) & (CellData{k}(:,12) == 2) & (CellData{k}(:, NDatacolumns + 8) == 1), NDatacolumns + 4), ...
            linspace(-1, 1, 100));
        DoC23CumHist = DoC23CumHist + histc(CellData{k}((CellData{k}(:,NDatacolumns + 1) > 0) & (CellData{k}(:,12) == 2) & (CellData{k}(:, NDatacolumns + 8) == 1), NDatacolumns + 5), ...
            linspace(-1, 1, 100));
        DoC31CumHist = DoC31CumHist + histc(CellData{k}((CellData{k}(:,NDatacolumns + 1) > 0) & (CellData{k}(:,12) == 3) & (CellData{k}(:, NDatacolumns + 8) == 1), NDatacolumns + 4), ...
            linspace(-1, 1, 100));
        DoC32CumHist = DoC32CumHist + histc(CellData{k}((CellData{k}(:,NDatacolumns + 1) > 0) & (CellData{k}(:,12) == 3) & (CellData{k}(:, NDatacolumns + 8) == 1), NDatacolumns + 5), ...
            linspace(-1, 1, 100));
    end
    
    DoC12 = DoC12CumHist/sum(DoC12CumHist(:));
    DoC13 = DoC13CumHist/sum(DoC13CumHist(:));
    DoC21 = DoC21CumHist/sum(DoC21CumHist(:));
    DoC23 = DoC23CumHist/sum(DoC23CumHist(:));
    DoC31 = DoC31CumHist/sum(DoC31CumHist(:));
    DoC32 = DoC32CumHist/sum(DoC32CumHist(:));

    % Plot summary data
    handles.handles.DoCFig = figure('color', [1 1 1], 'inverthardcopy', 'off');

    handles.handles.DoCAx(1) = subplot(6,1,1);
%     histHand = histogram(handles.handles.DoCAx(1), DoC1, 100);
    histHand = bar(handles.handles.DoCAx(1), linspace(-1, 1, 100), DoC12);
    histHand.FaceColor = Chan1Color;
    histHand.EdgeColor = rgb(52, 73, 94);
    set(handles.handles.DoCAx(1), 'XLim', [-1 1]);
    xlabel(handles.handles.DoCAx(1), 'DoC Score Ch1 -> Ch2', 'Fontsize', 12);
    ylabel(handles.handles.DoCAx(1), 'Frequency','FontSize',12);
    set(handles.handles.DoCAx(1),'FontSize',12)

    handles.handles.DoCAx(2) = subplot(6,1,2);
%     histHand = histogram(handles.handles.DoCAx(2), DoC2, 100);
    histHand = bar(handles.handles.DoCAx(2), linspace(-1, 1, 100), DoC13);
	histHand.FaceColor = Chan1Color;
    histHand.EdgeColor = rgb(52, 73, 94);
    set(handles.handles.DoCAx(2), 'XLim', [-1 1]);
    xlabel(handles.handles.DoCAx(2), 'DoC Score Ch1 -> Ch3', 'Fontsize', 12);
    ylabel(handles.handles.DoCAx(2), 'Frequency','FontSize',12);
    set(handles.handles.DoCAx(2),'FontSize',12)
    
    handles.handles.DoCAx(3) = subplot(6,1,3);
%     histHand = histogram(handles.handles.DoCAx(2), DoC2, 100);
    histHand = bar(handles.handles.DoCAx(3), linspace(-1, 1, 100), DoC21);
	histHand.FaceColor = Chan2Color;
    histHand.EdgeColor = rgb(52, 73, 94);
    set(handles.handles.DoCAx(3), 'XLim', [-1 1]);
    xlabel(handles.handles.DoCAx(3), 'DoC Score Ch2 -> Ch1', 'Fontsize', 12);
    ylabel(handles.handles.DoCAx(3), 'Frequency','FontSize',12);
    set(handles.handles.DoCAx(3),'FontSize',12)
    
    handles.handles.DoCAx(4) = subplot(6,1,4);
%     histHand = histogram(handles.handles.DoCAx(2), DoC2, 100);
    histHand = bar(handles.handles.DoCAx(4), linspace(-1, 1, 100), DoC23);
	histHand.FaceColor = Chan2Color;
    histHand.EdgeColor = rgb(52, 73, 94);
    set(handles.handles.DoCAx(4), 'XLim', [-1 1]);
    xlabel(handles.handles.DoCAx(4), 'DoC Score Ch2 -> Ch3', 'Fontsize', 12);
    ylabel(handles.handles.DoCAx(4), 'Frequency','FontSize',12);
    set(handles.handles.DoCAx(4),'FontSize',12)

    handles.handles.DoCAx(5) = subplot(6,1,5);
%     histHand = histogram(handles.handles.DoCAx(2), DoC2, 100);
    histHand = bar(handles.handles.DoCAx(5), linspace(-1, 1, 100), DoC31);
	histHand.FaceColor = Chan3Color;
    histHand.EdgeColor = rgb(52, 73, 94);
    set(handles.handles.DoCAx(5), 'XLim', [-1 1]);
    xlabel(handles.handles.DoCAx(5), 'DoC Score Ch3 -> Ch1', 'Fontsize', 12);
    ylabel(handles.handles.DoCAx(5), 'Frequency','FontSize',12);
    set(handles.handles.DoCAx(5),'FontSize',12)
    
    handles.handles.DoCAx(6) = subplot(6,1,6);
%     histHand = histogram(handles.handles.DoCAx(2), DoC2, 100);
    histHand = bar(handles.handles.DoCAx(6), linspace(-1, 1, 100), DoC32);
	histHand.FaceColor = Chan3Color;
    histHand.EdgeColor = rgb(52, 73, 94);
    set(handles.handles.DoCAx(6), 'XLim', [-1 1]);
    xlabel(handles.handles.DoCAx(6), 'DoC Score Ch3 -> Ch2', 'Fontsize', 12);
    ylabel(handles.handles.DoCAx(6), 'Frequency','FontSize',12);
    set(handles.handles.DoCAx(6),'FontSize',12)
    
    set(handles.handles.DoCFig,'Position',[100 100 500 1000])
    
    drawnow;

    % Save the figure
    try 
        print(fullfile(Outputfolder, 'Clus-DoC Results', 'DoC histograms', 'Pooled DoC histogram.tif'), ...
            handles.handles.DoCFig, '-dtiff');
    catch
        currFig = getframe(handles.handles.DoCFig);
        imwrite(currFig.cdata, fullfile(Outputfolder, 'Clus-DoC Results', 'DoC histograms', 'Pooled DoC histogram.tif'));
    end

    close gcf;
    
    save(fullfile(Outputfolder, 'Clus-DoC Results', 'Data_for_Cluster_Analysis.mat'),'Data_DoC','DensityROI'); %ROIData removed!


end
