# Reproduction-number-in-heterogeneous-populations-with-structured-mixing
Code associated to with the paper "The reproduction number of infectious diseases in heterogeneous populations with structured mixing patterns" 

## How to run this code
The figures and results from the paper can be reproduced by running the scripts in this repo. Some scripts depend on temporary (pickle) files created by other scripts, as follows...

*simulation_on_network.py* requires *adjacency_list_generator.py* and creates a pickle file *R0_analytical_and_simulated.p* as output. This script takes an hour or two to run on a normal laptop computer.\

*plot_figure1.py* requires *R0_analytical_and_simulated.p* and creates *figure1.png'* as output\

*combined_method_with_vs_without_outliers.py*, *sensitivity_to_z_threshold.py*, *process_data_for_contact_matrix.py*,*process_data_for_contact_distribution.py*, and *four_methods_applied_to_survey_data* all require data files *CoMix_BE_contact_common.csv*, *CoMix_BE_participant_common.csv*, *CoMix_BE_sday.csv* from https://zenodo.org/communities/social_contact_data/, and *Euro_pops.csv* from https://population.un.org/wpp/Download/Standard/MostUsed/. \

Running *combined_method_with_vs_without_outliers.py* with *filtered=True* creates *../pickles/data_for_all_wave_plots_filtered.p*, and with *filtered=False* creates *../pickles/data_for_all_wave_plots_unfiltered.p*. These each take ten or so minutes.\

*sensitivity_to_z_threshold.py* creates *data_for_R0_vs_cutoff_plots.p*. This takes a few minutes.\

*plot_figure2.py* requires *data_for_all_wave_plots_unfiltered.p*, *data_for_all_wave_plots_filtered.p* and *data_for_R0_vs_cutoff_plots.p*, and it creates *figure2.png*. Takes a few seconds.\

*process_data_for_contact_matrix.py* creates *data_for_matrix_plots.p*. Takes a few seconds.\

*process_data_for_contact_distribution* creates *data_for_contact_distribution_plots.p*. Takes a few seconds.\

*four_methods_applied_to_survey_data* creates *data_for_R0_bootsrap_plots.p*. Takes around five minutes. \

*plot_figure2.py* requires *data_for_matrix_plots.p*, *data_for_contact_distribution_plots.p*, *data_for_R0_bootsrap_plots.p*. Takes a few seconds.\

