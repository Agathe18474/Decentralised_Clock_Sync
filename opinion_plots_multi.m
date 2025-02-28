function  opinion_plots_multi(t, convergence_metric, ...
    convergence_metric_max, convergence_metric_min, boxplot_data, ...
    percentage_filtered, percentage_filtered_max, percentage_filtered_min,...
    plot_dict, plot_titles, line_names, zero_val, type)

    %   Plots comparison of MULTIPLE RUNS (each a batch of monte carlo) of 
    %   opinion diffusion scenario as per plot dictionary selection
    % INPUTS:
    %   node_val            exact node values
    %   t                   time elapsed
    %   plot_dict           dictionary of plot selections [boolean]. Allows: 
    %                           Scatter Plot - plots avg convergence, and
    %                               at constant intervals, boxplots and all 
    %                               convergence values 
    %                           Error Bars - plots avg convergence, and at 
    %                               constant intervals, error bars 
    %                               indicating max and min convergence
    %                               metric of the monte carlo batch                             
    %                           Filtered Values Percentage - plots avg 
    %                               percentage of values filtered out by
    %                               nodes for given runs
    %  convergence metrics   array of convergence metric for MULTIPLE runs
    %   type                is plot for meaningless values, time, or rate values?  
    % OUTPUTS:
    %  node trajectory plot
    %  convergence metric plot 
    %  percentage of values filtered for each run

    plot_length = size(convergence_metric, 1);
    max_time_step = t(end)+1;

    colours = parula(plot_length+1);
    colours_boxplot = repmat(colours(1:end-1,:), ceil((max_time_step+2)/(plot_length)),1);

    if plot_dict(1)
        % SCATTER PLOT
        f = figure();
        hold on
        for i=1:plot_length
            pp(i) = plot((t+1).',convergence_metric(i,:).', "-", 'LineWidth', 0.5, 'color', colours(i,:),'DisplayName',line_names(i));
        end
        
        boxplot(boxplot_data, 'color', colours_boxplot, 'Notch','on') 
        c = 1;
        for i=1:t(end)+1
            scatter(i, boxplot_data (:,i),[], colours_boxplot(c,:), 'filled', 'MarkerFaceAlpha', 0.25, 'jitter', 'on', 'jitterAmount', 0.15)
            if rem(i,(actual_length+1)) == 0
                c = 0;
            end
            c = c +1;
        end
        ax = gca;
        ax.YAxis.Scale ="log";
        xlabel("Time steps")
        if type == "Neutral"
            ylabel("Convergence Metric")  
        elseif type == "Time"
            ylabel("Time Convergence Metric")
        elseif type == "Rate"
            ylabel("Rate Convergence Metric")
        end
        
        
        grid on
        box on 
        fontsize(f,13, "points")
        title(plot_titles)
        legend(line_names, 'Location', 'northeastoutside')
        %legend(pp(1:actual_length+1), 'Location', 'northeastoutside')
        %legend(names, 'Location', 'best')
        hold off
    end
    
    
    
    if plot_dict(2)
        % ERROR BAR PLOT
        f = figure();
        t_cm = repmat(t, plot_length,1);
        semilogy(t_cm.',convergence_metric.', "-", 'LineWidth', 2) 
        hold on
        ax = gca; 
        ax.ColorOrder = colours;
        ax.LineStyleCyclingMethod = 'withcolor'; 
        
        % plot error bars
        t_err = repmat(t,2,1);
        for i=1:plot_length
            semilogy(t_err(:,i:plot_length:end),[convergence_metric_max(i,i:plot_length:end); convergence_metric_min(i,i:plot_length:end)],'-_', 'color', colours(i, :)) 
            semilogy(t_err(1,i:plot_length:end),convergence_metric(i,i:plot_length:end),'square', 'color', colours(i, :), 'MarkerFaceColor',colours(i, :)) 
        end
        xlim([0 t(end)])
        %legend(names, 'Location', 'northeastoutside')
        %or 
        legend(line_names, 'Location', 'best')
        %ylim([zero_val 1]) 
        grid on
        xlabel("Time steps")
        if type == "Neutral"
            ylabel("Convergence Metric")  
        elseif type == "Time"
            ylabel("Time Convergence Metric")
        elseif type == "Rate"
            ylabel("Rate Convergence Metric")
        end
        grid on
        box on 
        fontsize(f,13, "points")
        title(plot_titles)
        
    end
    if plot_dict(3)
        % FILTERED % VALUE PLOTS
        f = figure();
        t_cm = repmat(t(1:end-1), plot_length,1);
        plot(t_cm.',  percentage_filtered.', "-", 'LineWidth', 2)
        % plot error bars
        t_err = repmat(t,2,1);
        hold on 
        for i=1:plot_length
            semilogy(t_err(:,i:plot_length:end-1),[percentage_filtered_max(i,i:plot_length:end); percentage_filtered_min(i,i:plot_length:end)],'-_', 'color', colours(i, :)) 
            semilogy(t_err(1,i:plot_length:end-1),percentage_filtered(i,i:plot_length:end),'square', 'color', colours(i, :), 'MarkerFaceColor',colours(i, :)) 
        end
        %plot(parameter_tested(1):test_step_size:parameter_tested(end), percentage_filtered, "-", 'LineWidth', 2)
        ax = gca; 
        ax.ColorOrder = colours;
        ax.LineStyleCyclingMethod = 'withcolor'; 
        ylim([0 1]) 
        xlabel("Number of Disruptive Agents")
        ylabel("% of values filtered out")
        legend(line_names, 'Location', 'best')
        grid on
        title(plot_titles)
    end
end