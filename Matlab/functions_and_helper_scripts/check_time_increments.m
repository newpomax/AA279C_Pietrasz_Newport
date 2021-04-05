function plot_format = check_time_increments(plot_format)
    % Based on time_increments: 's', 'min', 'hours', 'days', or 'years'
    switch plot_format.time_increments
        case 's'
            plot_format.seconds_to_increment = 1;
        case 'min'
            plot_format.seconds_to_increment = 60;
        case 'hours'
            plot_format.seconds_to_increment = 60*60;
        case 'days'
            plot_format.seconds_to_increment = 24*60*60;
        case 'years'
            plot_format.seconds_to_increment = 365.25*24*60*60;
        otherwise
            warning('Invalid time_increment "%s": defaulting to day.', ...
                plot_format.time_increments);
            plot_format.time_increments = 'day';
            plot_format.seconds_to_increment = 60;
    end
    
end