% lasso_coord_desc_test_runner

%{

TODO :: 

- really needs an intercept.
- look more closely at behavioral model - what is justified?
- which FR inputs to put in?  Something that makes sense when
mean-subtracted?

%}


rng(100)

% Standard parameters
lasso_params = struct();
lasso_params.min_lambda     = 1E-6;
lasso_params.max_lambda     = 5;
lasso_params.n_lambdas      = 200;
lasso_params.lambda_linlog  = 'log'; % 'Linear' or 'log' defines spacing of lambda.
lasso_params.standardize    = 1;  % z-scores X and mean-centers Y.
lasso_params.verbose        = 0;    % Acts on every second � :: { 0 [do nothing] , 1 [shrink by 5x], 2 [set to zero] }

% lasso_params.n_predictors   = 8;
% lasso_params.n_dimensions   = 21;
% lasso_params.n_samples      = 50;
% lasso_params.SNR            = NaN;
% lasso_params.regularization = NaN;    % Acts on every second � :: { 0 [do nothing] , 1 [shrink by 5x], 2 [set to zero] }

% Auc right-side edge detection
% 'full_range'            :: Samples the full from Lambda 0 .. max.
% 'double_baseline_error' :: Samples the range from error(OLS) to 2x error(OLS)
% 'midpoint_RMSE'         :: Samples min(RMSE) to mean( [ min(RMSE), max(RMSE) ])
% 'midpoint_L1'           :: Samples min(RMSE) to RMSE @ mean( [ min(L1), max(L1) ])
% 'peak_inflection_point' :: Finds inflection by 2nd derivative, takes biggest peak
% 'last_inflection_point' :: Finds inflection by 2nd derivative, takes last peak
% 'cross-validation'      :: NOT IMPLEMENTED :: Find lambda with highest cross-val score.
lasso_params.auc_edge_detector = 'midpoint_RMSE';

% Misc

do_clusters          = 1;
n_random_rotations   = 50;

if do_clusters
    
    X = load('/Users/avaughan/Documents/MATLAB/_projects_CSHL/OFC_ANALYSIS_AGV/Clustering/allRegressionInputs.mat');
    X = X.allRegressionInputs_mat;
    %Strip X matrix down and make irrelevant variables NaNs.
    X = X(:, 1:9); % Toss conditional updating
    X(:, 4) = sum(X(:, 4:5),2); %Joint confidence estimate.
    X(:, 5) = [];
    for i = [1,2,4,5,6,8,]
        X(31:42,i) = NaN;
    end
    for i = [3,7,]
        X(1:30,i) = NaN;
    end
    X(:,2) = 0;
    X(31:42,2) = 1;
    X(31:42,4) =  X(31:42,3);
    X(:,3) = [];
    X(:,4) = X(:,4) + X(:,5);
    X(:,5) = [];
    
    X = nanzscore(X);
    namedFigure('Predictor matrix.')
    subplot(1,2,1); imagesc(X)
    subplot(1,2,2); imagesc(isnan(X));
    
    cluster_data = load('/Users/avaughan/Documents/MATLAB/_projects_CSHL/OFC_ANALYSIS_AGV/AGV output/cellids490_ID20_RP1/cellIDs490_ID20_RP1 -- XX_agv_raw_634cells_pC.pE.evidence_2014.2.21 - 2014-05-19 20-11-10/After fixing ppca 7.15.14/21PCs/figs_12_10_2014/workingData_2014_12_14.mat');
    cluster_normalization = 'z_score_reconstructed';
    switch cluster_normalization
        case 'minmax'
            Y_all = cluster_data.XX_norm_subsampled(cell2mat(IDsByCluster),:)';    % Range 0..1
        case 'z_score_reconstructed'
            Y_all = cluster_data.behaviorResponse_reconstructed_sortedByCluster'; % Z-scored
            %Y_all = cluster_data.XX_reconstructed(cell2mat(IDsByCluster),:);
    end
    
    %    .clusN       : [49 61 35 21 23 97 57 100 42]
    %    .cumclusN    : [49 110 145 166 189 286 343 443 485]
    cluster_indices = [0 cluster_data.cumclusN];
end

for do_average_clusters = 1
    for cluster_to_test = 2
        
        fprintf('\n------------RESULTS-------------------\n')
        fprintf('Testing control for cluster %.0f.\n',cluster_to_test)
        
        % Setup
        if do_clusters
            Y = Y_all(:,(cluster_indices(cluster_to_test)+1):cluster_indices(cluster_to_test+1));
            
            if do_average_clusters
                Y = mean(Y,2);
            end
            lasso_params.X = X;
            lasso_params.Y = Y;
            lasso_params.fig_suffix = sprintf('Cluster %.0f (Averaged: %.0f)',cluster_to_test,do_average_clusters);
        end
        
        % 1st lasso
        lasso_params.do_rotation    = 0;    % Rotate basis vectors?
        lasso_params.reset_figures  = 1;
        
        
        [X,Y,intercept_true,beta_true,RMSE,L1,auc_struct,line_handles] = lasso_coord_desc_test( lasso_params );
        
        lasso_data{do_average_clusters+1,cluster_to_test}.cluster_normalization = cluster_normalization;
        lasso_data{do_average_clusters+1,cluster_to_test}.lasso_params = lasso_params;
        lasso_data{do_average_clusters+1,cluster_to_test}.X = X;
        lasso_data{do_average_clusters+1,cluster_to_test}.Y = Y;
        lasso_data{do_average_clusters+1,cluster_to_test}.intercept_true = intercept_true;
        lasso_data{do_average_clusters+1,cluster_to_test}.beta_true = beta_true;
        lasso_data{do_average_clusters+1,cluster_to_test}.RMSE = RMSE;
        lasso_data{do_average_clusters+1,cluster_to_test}.L1 = L1;
        lasso_data{do_average_clusters+1,cluster_to_test}.auc_struct = auc_struct;
        lasso_data{do_average_clusters+1,cluster_to_test}.line_handles = line_handles;
        
        
        auc_cx = auc_struct.auc;
        lasso_params.X = X;
        lasso_params.Y = Y;
        lasso_params.intercept_true = intercept_true;
        lasso_params.beta_true = beta_true;
        lasso_params.auc_edge_RMSE = auc_struct.auc_edge_RMSE;
        lasso_params.auc_edge_L1 = auc_struct.auc_edge_L1;
        
        lasso_params.do_rotation    = 1;    % Rotate basis vectors?
        lasso_params.reset_figures  = 0;
        
        rot_stats = struct('i',num2cell(1:n_random_rotations));
        
        auc_rot = [];
        fprintf('Testing rotated basis')

        for i = 1:n_random_rotations
            fprintf('.')
            if mod(i,20) == 0
                fprintf('\n                     ')
            end
            
            [   rot_stats(i).X,...
                rot_stats(i).Y,...
                rot_stats(i).intercept_true,...
                rot_stats(i).beta_true,...
                rot_stats(i).RMSE,...
                rot_stats(i).L1,...
                rot_stats(i).auc_struct] = lasso_coord_desc_test( lasso_params );
            
            auc_rot(i) = rot_stats(i).auc_struct.auc;
        end
        
        lasso_data{do_average_clusters+1,cluster_to_test}.rot_stats = rot_stats;
        
        %%
        
        fprintf('%.0f cells\n',diff(cluster_indices(cluster_to_test + [0 1])));
        fprintf('Cluster :: %.0f\n',cluster_to_test')
        fprintf('AUC_CX \t\t:: %.3f\n',auc_cx)
        if ~isempty(auc_rot)
            %fprintf('auc_rot(%.0f) \t:: %.3f\n',length(auc_rot),auc_rot(end))
            fprintf('AUC_ROT (mean) \t:: %.3f\n',mean(auc_rot));
            fprintf('AUC_ROT (stdev)\t:: %.3f\n',std(auc_rot));
            p_sr = signrank( auc_rot - auc_cx );
            fprintf('\tp <= %.3f (Sign Rank Test)\n',p_sr);
            p_boot = sum( auc_rot <= auc_cx)/length(auc_rot);
            fprintf('\tp <= %.3f (Bootstrap)\n',p_boot);
        end
        
        %%
        % Set zdata for the control curve to the top
        
        for i = 1:length(line_handles)
            for j = 1:length(line_handles{i})
                try, % One
                    set(line_handles{i}(j),'ZData',100 * ones(size(get(line_handles{i}(j),'XData'))))
                end
            end
        end
        
        %% Plot significance curve.
        
        i = do_average_clusters +1;
        j = cluster_to_test;

        spline_length = 100;
        spline_inds = 1:(find(diff(lasso_data{i,j}.RMSE)==0,1)-1);
        lasso_data{i,j}.RMSE_spline = linspace(min(lasso_data{i,j}.RMSE),max(lasso_data{i,j}.RMSE),spline_length);
        lasso_data{i,j}.L1_spline  = spline(...
            lasso_data{i,j}.RMSE(spline_inds),...
            lasso_data{i,j}.L1(spline_inds),...
            lasso_data{i,j}.RMSE_spline);
        
        lasso_data{i,j}.rot_spline = zeros(n_random_rotations,spline_length);
        for k = 1:n_random_rotations,
            spline_inds = 1:(find(diff(lasso_data{i,j}.rot_stats(k).RMSE)==0,1,'first')-1);
            lasso_data{i,j}.rot_spline(k,:) = spline( lasso_data{i,j}.rot_stats(k).RMSE(spline_inds) ,lasso_data{i,j}.rot_stats(k).L1(spline_inds), lasso_data{i,j}.RMSE_spline);
        end
        
        namedFigure(sprintf('Significance Analysis :: Cluster %.0f\n',cluster_to_test'))

        subplot(5,1,1)
        imagesc(lasso_data{i,j}.L1_spline)
        clims = caxis;
        title('L1 Penalty','Units', 'normalized', ...
            'Position', [0.05 1], 'HorizontalAlignment', 'left')
        set(gca,'XTick',[])
        
        subplot(5,1,2)
        imagesc(lasso_data{i,j}.rot_spline)
        title('L1 Penalty (rotated','Units', 'normalized', ...
            'Position', [0.05 1], 'HorizontalAlignment', 'left')
        caxis(clims)
        set(gca,'XTick',[])
        
        colorbar
        subplot(5,1,3)
        imagesc(bsxfun(@minus,lasso_data{i,j}.rot_spline,lasso_data{i,j}.L1_spline))
        caxis([-1 1])
        colorbar
        title(' L1 (rotated) - L1 (original','Units', 'normalized', ...
            'Position', [0.05 1], 'HorizontalAlignment', 'left')
        set(gca,'XTick',[])
        
        subplot(5,1,4)
        imagesc(bsxfun(@lt,lasso_data{i,j}.rot_spline,lasso_data{i,j}.L1_spline))
        colorbar('Location','EastOutside')
        title(' L1(rotated) < L1_spline','Units', 'normalized', ...
            'Position', [0.05 1], 'HorizontalAlignment', 'left')
        set(gca,'XTick',[])
        
        subplot(5,1,5); cla; hold on
        %for significance we want real L1 to be less than the rotated.
        plot(mean(bsxfun(@lt,lasso_data{i,j}.rot_spline,lasso_data{i,j}.L1_spline)),'r')
        plot(mean(...
            bsxfun(@lt,...
            cumsum(lasso_data{i,j}.rot_spline,2),...
            cumsum(lasso_data{i,j}.L1_spline))...
            ),'k')
        hline(0.05)
        title('mean(L1(rotated) < L1_spline','Units', 'normalized', ...
            'Position', [0.05 1], 'HorizontalAlignment', 'left')
        colormap red_white_blue
        ylabel({'p( L1_{true} >=  L1_{rot} )',' [bootstrap'})
        legend('Pointwise L1','Cumulative L1')
        suptitle(sprintf(' Significance :: Averaged %.0f :: Cluster %.0f',do_average_clusters,j))
        ylim([0 1])

        %%
    end
end

lasso_plot_L_curve(lasso_data)

fprintf('\n\nFinal variables')
list_variables(lasso_data{1,1})





%{ 

:: Notes :: 

switch regularization_form with SNR = 5
case 0
    % (Do nothing.)
    - No particular effect if beta_true is not regularized.
case 1
    % (Divide some by 5.)
    - Great - 2 vs. 2.2 AUC under 50% of range
    - p ~ 0.
case 2
    % (Set some to zero.)
    - Looks good - very kinked.
    - Often significant, but not always.

%}


%%


% Debugging code for finding L1 inflection point.
% figure(3);
% subplot(3,1,1)
% plot(RMSE_spline,L1_spline)
% subplot(3,1,2)
% plot(RMSE_spline(2:end),diff(L1_spline))
% subplot(3,1,3); hold on
% plot(RMSE_spline(3:end),diff(diff(L1_spline)))
% 
% spline_ind = find(diff(RMSE)==0,1);
% RMSE_spline = linspace(min(RMSE),RMSE(spline_ind),500);
% L1_spline = spline( RMSE(1:spline_ind), L1(1:spline_ind), RMSE_spline)
% [peaks, loc] = findpeaks(  diff(diff(L1_spline)) )
% peak = peaks(find(peaks == max(peaks)))
% loc = loc(find(peaks == max(peaks)));
% 
% plot(RMSE_spline(loc+2), peak,'ro')
