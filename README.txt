This directory contains four folders corresponding to the main components of the analysis presented in the paper.

Folder descriptions
-------------------

AdaptiveCUSUM_SimulationCode
    Contains scripts used to obtain the stable distributions, compute control limits,
    and run simulations to estimate out-of-control (OC) run lengths for the proposed
    adaptive CUSUM statistic.

EWMA_Xie
    Contains scripts used to compute control limits and run simulations for the
    EWMA charting statistic proposed by Xie et al. for the Gumbel distribution.

RealDataAnalysis
    Contains scripts used for the real data analysis on the NASA turbofan dataset.
    These scripts:
        - compute the time-varying distribution,
        - estimate the corresponding control limits, and
        - run the monitoring procedure on the data.

    The processed dataset (obtained from the raw data using the thresholds described
    in the paper) and the corresponding in-control parameter estimates are stored in:
        DataAndParamsList.RDS

RunShewhartChart
    Contains scripts used to run simulations to obtain OC run lengths for the
    Shewhart chart proposed by Zwetsloot et al.


Pre-computed objects
--------------------

For convenience, the folders also include the in-control distributions and control
limits required to reproduce the OC run length simulations. These objects can also
be recomputed using the provided scripts if desired.


Running the scripts
-------------------

The file names are intended to be self-explanatory and indicate the purpose of each
script.

For each script, the main task for the user is to set the working directory to the
appropriate folder. After this is done, all other paths in the code are specified as
relative paths and will resolve correctly.

In each script, the working directory is defined immediately under the section:

    ########## Main ##########

For example, in the folder:

    AdaptiveCUSUM_SimulationCode

the file:

    GetOC_Runlengths.R

contains the following lines:

    ############# Main #############
    setwd("/pathToFolder/AdaptiveCUSUM_SimulationCode")  # Path to this folder

When running the code on a local machine, the user should modify this line to point
to the correct location. For example:

    setwd("/home/Documents/Scripts/AdaptiveCUSUM_SimulationCode")

Once the working directory is set, all required files (such as the control limits
and stable distributions) will be loaded automatically using relative paths.
