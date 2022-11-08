# Covid-19 winter scenarios for Norway

This repository includes the code to procude covid-19 scenarios for the autumn/winter of 2022 in Norway. A written report describing the scenarios is published [here](https://www.fhi.no/contentassets/e6b5660fc35740c8bb2a32bfe0cc45d1/vedlegg/nasjonale-og-regionale-rapporter/winter_scenarios_2022_2023.pdf) with additional context in the [risk evalutation from 8th of November 2022](https://www.fhi.no/contentassets/c9e459cd7cc24991810a0d28d7803bd0/vedlegg/risikovurdering-2022-11-08.pdf). The scenarios have been developed by the covid-19 modelling team at the Norwegian Institute for Public Health.

The code in this repository depends on the [metapop](https://github.com/Gulfa/metapop) and [metapopnorge](https://github.com/Gulfa/metapopnorge) packages. These and the other requirements can be installed by running install_reqs.R. 

To run the scenarios edit the number of threads and simulations in the mclapply function call in run_analysis.R and run the file. This creates a all_results.RDS file that can then be plotted using plot_results.R

Any questions or comments about the code or the scenarios are welcome.
