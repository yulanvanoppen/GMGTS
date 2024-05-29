function lhandle_new = add_legendentry(leg_handle, plot_handles, plot_names)
    h_list = [];
    for idx = 1:length(leg_handle.String)
        h_list = [h_list findobj(gca, 'DisplayName', leg_handle.String{idx})]; %#ok<AGROW>
    end
    h_list = [h_list plot_handles];
    if ~iscell(plot_names), plot_names = {plot_names}; end
    str_list = [leg_handle.String plot_names];
    lhandle_new = legend(h_list, str_list{:});
end