# Covid-19 winter scenarios for Norway

This repository includes the code to procude covid-19 scenarios for the autumn/winter of 2022 in Norway. A written report describing the scenarios is published here. The scenarios have been developed by the covid-19 modelling team at the Norwegian Institute for Public Health.

The code in this repository depends on the [metapop](https://github.com/Gulfa/metapop) and [metapopnorge](https://github.com/Gulfa/metapopnorge) packages. These and the other requirements can be installed by running install_reqs.R. 

To run the scenarios edit the number of threads and simulations in the mclapply function call in run_analysis.R and run the file. This creates a all_results.RDS file that can then be plotted using plot_results.R

Any questions or comments about the code or the scenarios are welcome.
