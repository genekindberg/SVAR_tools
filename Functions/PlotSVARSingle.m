function [] = PlotSVARSingle(IRF_draws, varnames, varuse, diffs, irfperiods, wholetitle, shock)

%% Description
% IRF_draws from SVAR estimation, 4d matrix ordered as draws, variable,
% time, shock.
% varnames - names of variables to print on plot
% varuse - subset of variables in var to plot
% diffs - 0,1 vector, where 1 signifies that the impulse response should be
% cumulated (i.e. to go from impact on growth to impact on the level if log-differences)
% irfperiods - number of periods over which to plot the IRF

%% Prelims

FillColor05   =[.5 .5 .5];
FillColor16   =[.7 .7 .7];
FontSizeset = 8;
nvars = size(varuse,2);

%% Make shape of subplot based on size

if nvars<3
    wid = nvars;
    len = 1;
else
    wid = 3;
    len = ceil(nvars/3);
end


%% Plot variables
figure
for jj = 1:nvars
    subplot(len,wid,jj)
    if diffs(jj)==0
        fill([(1:irfperiods) fliplr((1:irfperiods))],...
        [ quantile(squeeze(IRF_draws(:,varuse(jj),:,shock)),[0.05])  fliplr(quantile(squeeze(IRF_draws(:,varuse(jj),:,shock)),[0.95]))],...
        FillColor05,'EdgeColor','none');
        hold on
        fill([(1:irfperiods) fliplr((1:irfperiods))],...
        [ quantile(squeeze(IRF_draws(:,varuse(jj),:,shock)),[0.16])  fliplr(quantile(squeeze(IRF_draws(:,varuse(jj),:,shock)),[0.84]))],...
        FillColor16,'EdgeColor','none');
        hold on
        plot(median(squeeze(IRF_draws(:,varuse(jj),:,shock)),1,'omitmissing')','-k','LineWidth',2)
        hold on
    elseif diffs(jj)==1
        fill([(1:irfperiods) fliplr((1:irfperiods))],...
        [ cumsum(quantile(squeeze(IRF_draws(:,varuse(jj),:,shock)),[0.05]))  fliplr(cumsum(quantile(squeeze(IRF_draws(:,varuse(jj),:,shock)),[0.95])))],...
        FillColor05,'EdgeColor','none');
        hold on
        fill([(1:irfperiods) fliplr((1:irfperiods))],...
        [ cumsum(quantile(squeeze(IRF_draws(:,varuse(jj),:,shock)),[0.16]))  fliplr(cumsum(quantile(squeeze(IRF_draws(:,varuse(jj),:,shock)),[0.84])))],...
        FillColor16,'EdgeColor','none');
        hold on
        plot(cumsum(median(squeeze(IRF_draws(:,varuse(jj),:,shock)),1, 'omitmissing')')','-k','LineWidth',2)
    end
    plot(zeros(irfperiods,1),'-k','LineWidth',1)
    hold off
    box on
    grid on
    set(gca,'linewidth',2)
    title([varnames{jj}])
    set(gca,'FontSize', FontSizeset)
end

sgtitle(wholetitle)

end
