function opinion_plots_multi_combined(t, convergence_metrics, ...
    convergence_metrics_max, convergence_metrics_min,  names)
    %   Plots comparison of MULTIPLE RUNS (each a batch of monte carlo) of 
    %   opinion diffusion scenario as per plot dictionary selection
    %   combines time & rate convergence metric into a single plot
    % INPUTS:
    %   node_val                    exact node values
    %   t                           time elapsed
    %   convergence metrics         array of convergence metric combining time & rate=
    %   convergence_metrics_max     max value of convergence metrics (time & rate)
    %   convergence_metrics_min     min value of convergence metrics (time & rate)
    %   names                       name plots
    % OUTPUTS:
    %  node trajectory plot
    %  convergence metric plot 
    %  percentage of values filtered for each run
    


    nb_plots = size(convergence_metrics, 1);

    colours = parula(nb_plots+1);

    f = figure();
    hold on
    box on 
    grid on

    % plot error bars
    t_err = repmat(t,nb_plots,1);
    error_bars = 20;
    for i=1:nb_plots
        semilogy(t.',convergence_metrics(i, :).', "-", 'color', colours(i, :), 'LineWidth', 2) 
        semilogy(t_err(:,i:error_bars:end),[convergence_metrics_max(i,i:error_bars:end); convergence_metrics_min(i,i:error_bars:end)],'-_', 'color', colours(i, :)) 
        semilogy(t(1,i:error_bars:end),convergence_metrics(i,i:error_bars:end),'square', 'color', colours(i, :), 'MarkerFaceColor',colours(i, :)) 
    end
    xlim([0 t(end)])
    ax = gca;
    ax.YAxis.Scale ="log";
    xlabel("Time, [s]")
    ylabel("Convergence metric")

    legend('','Location', 'northeast')
    for ii =1:length(names)            
         plot(NaN, NaN, '-', 'Color', colours(ii,:), 'DisplayName', names{ii},'LineWidth', 2)
    end

    fontsize(14,"points")

end
