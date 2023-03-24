#
# Overall parameters

## Sample dataframe

A dataframe containing one line per sample on which we wish to run. Must contain at least the following columns:

- sample_id: a unique string identifying the sample. Usually a (relative) folder structure.
- sample_path: the relative path to the sample
- fs: the sampling frequency of the sample

__Warning__: fs must not change, the current dependencies do not check for an fs change.
*Note*: tests runs are achieved by giving a small number of samples in the sample dataframe

## Execution parameters

- for artefact removal:
    - version: the version we wish to execute
- for figure generation:
    - figure.version: the version of the figure generation function
    - figure.sampling: the downsampling before plotting signal
    - figure.size: the size of the output figure
- for viewing:
    - view.sampling: the downsampling before plotting signal
- for folder_structure:
    - folders.results: base folder for storing results. Usually in the form of {path}/Clean

## Variables for those parameters

- The dataframe must be stored in a variable called *samples_df*
- The execution parameters must be stored in a dictionary

# Folder structure for signals

Where each file is stored. Variables in brackets indicate that they depend on overall parameters and {date} means current date/time and {clean_params} extends to v={version}_ap={artefact_params}

- raw:             {sample_path}
- cleaned:         {result_folder}/{clean_params}/signal/{sample}_signal.npy
- summary:         {result_folder}/{clean_params}/summary/{sample}_summary.tsv
- png:             {result_folder}/{clean_params}/figures/{sample}_fp={figure_params}.png
- overallsummary : {result_folder}/{clean_params}/run_info_{date}/overallsummary.tsv
- samples :        {result_folder}/{clean_params}/run_info_{date}/samples.tsv

# Rules

## Final rules:

- all: computes the block for all signals but without generating figures
- complete: computes the block for all signals with figures
- view: does not compute anything but allows to view current state

## Rules sketch

- rule make_sample_file:
    - input:
        - None
    - output:
        - samples
    - run:
        dump sample_id column in file

- rule rm_artefacts:
    - input:
        - raw
    - output:
        - cleaned
        - summary
    - params:
        - version
        - fs
    - run:
        toolbox.rm_artefacts

- rule create_fig: (optional)
    - input:
        - raw,
        - cleaned,
        - summary
    - params:
        - figure_sampling,
        - figure_size, 
        - figure_version,
        - version,
        - fs
    - output:
        - png
    - run:
        scripts/artefact_figure

- rule overallsummary:
    - input:
        - samples
        - summary[ samples ]
    - ouput:
        - overallsummary
    - params: 
        - version
    - run:
        compute_features_from_individual_summaries

- rule all:
    - input:
        - samples
        - overallsummary
        - summary[ samples ]
        - cleaned[ samples ]
    - params:
        - version
    - output: 
        - None

- rule complete:
    - input:
        - samples
        - overallsummary
        - summary[ samples ]
        - cleaned[ samples ]
        - png[ samples ]
    - params:
        - version
        - figure.params
    - output: 
        - None

- rule view:
    - input:
        - samples
    - output: 
        - None
    - params:
        - fs
        - version
    - run:
        show_iter_with_refresh with refresh doing the following:
        get_new_samples_with_summary_and_clean_computed that are within samples
        compute_features_from_individual_summaries
        create_visualization_order (adding the new ones at end)
        update_figure