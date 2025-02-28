function  opinion_plots_single(plotting_val, t, convergence_metric, ...
    max_difference, i_dis, plot_dict,zero_val, plot_colour, type)
    % Plots SINGLE RUN of opinion diffusion scenario scenario as per plot 
    % dictionary selection
    % INPUTS:
    %   plotting_val        node values or node rates
    %   t                   time elapsed
    %   convergence metric    array of convergence metric for a single run 
    %   i_dis               id of disruptive nodes
    %   starting counter
    %   plot_dict           dictionary of plot selections [boolean]. Allows: 
    %                           node trajectory
    %                           convergence metric 
    %   plot_colour         plot colour - monochromatic vs multichromatic
    %   type                is plot for meaningless values, time, or rate values?    
    % OUTPUTS:
    %  node trajectory plot
    %  convergence metric plot
   
    n = size(plotting_val,1);

    i_compliant = true(1,n).';
    if ~isempty(i_dis)
        i_compliant(i_dis) = false;
    end

    if ~exist('node_names','var')
        % if node_names do not exist, initialise manually 
        % initialise node names
        node_names = strings(1,n);
        node_node_spec = "Node %d";
        
        for i=1:n
            node_names(1,i) = sprintf(node_node_spec,i);
        end
        node_names(i_dis) = "Disruptive Node";
    end
    
    if plot_dict(1)
        % NODE TRAJECTORY
        f = figure();
        clf
        hold on
        fontsize(f,13, "points")

        if plot_colour == "mono"
            length_t = length(t);
            tt_normal = repmat(t,(n-length(i_dis)),1);
            plot(repmat(t,[length(i_dis),1]).', plotting_val(~i_compliant,1:length_t).', '--', 'Color','#21918c')
            h = plot(tt_normal.', plotting_val(i_compliant,1:length_t).', 'Color', '#440154');
            %legend(node_names)
            line_names = ["Compliant Nodes", "Disruptive Nodes"];
            line_widths = [1.5,1];
            line_type = ["-", "--"];
            cols =  ['#440154'; '#21918c'];
            legend('','Location', 'northeast')
            for j =1:length(line_names)            
              plot(NaN, NaN, line_type(j), 'Color', cols(j,:), 'DisplayName', line_names(j),'LineWidth', line_widths(j))
            end
            

        elseif plot_colour == "multi"
            cmap = hsv(length(i_compliant)+1); 
            node_val_compliant = plotting_val(i_compliant,:);
            for i = 1:length(i_compliant)     
                plot(t,node_val_compliant(i,:),'Color',cmap(i,:)); 
            end
            hold on 
            
            m = mean(node_val_compliant, 1);
            plot(t, m, '*k')
            rate = m./t;

            txt = ['\leftarrow rate: ' num2str(rate(round(t(end)/2)))];
            text(round(t(end)/2), m(round(t(end)/2)),txt)
            txt = ['\leftarrow rate: ' num2str(rate(3))];
            text(t(3), m(3),txt)
            txt = ['rate: ' num2str(rate(end-5)) '\rightarrow'];
            text(t(end-5),m(end-5),txt,'HorizontalAlignment','right')

            if n<=20 
                lgd=legend(node_names,'Location', 'eastoutside');
                fontsize(lgd,10,'points')   
            end
                
        end

        grid on
        xlabel("Time steps")
        if type == "Neutral"
            ylabel("Node Values")  
        elseif type == "Time"
            ylabel("Node V-times")
        elseif type == "Rate"
            ylabel("Node rates")
        end
        box on     
        xlim([0 t(end)])

        % PLOT FOR PAPERS
        %
        %legend('','Location', 'northeast')
        %cols =  ['#440154'; '#21918c'];
        %line_type = ["-", "--"];
        %for j =1:length(line_names)            
        %    plot(NaN, NaN, line_type(j), 'Color', cols(j,:), 'DisplayName', line_names(j),'LineWidth', line_widths(j))
        %end

    
    end

    if plot_dict(2)
        % RELATIVE CONVERGENCE METRIC
        f = figure();
        clf
        semilogy(t,convergence_metric, 'k') 
        
        xlim([0 t(end)])
        ylim([zero_val inf]) 
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
    end

    if plot_dict(3)
        f = figure();
        clf
        semilogy(t, max_difference)

        xlim([0 t(end)])
        ylim([zero_val inf]) 
        xlabel("Time steps")
        grid on
        box on 
        fontsize(f,13, "points")

        if type == "Neutral"
            ylabel("Maximum Node Difference")  
        elseif type == "Time"
            ylabel("Maximum Time Difference")
        elseif type == "Rate"
            ylabel("Maximum Rate Difference")
        end


    end
end