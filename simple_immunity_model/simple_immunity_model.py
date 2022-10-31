#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm
from datetime import datetime
import os
pd.options.mode.chained_assignment = None  # default='warn'

def simple_immunity_model(scenario_name, age_file = "age_FHI_2022.txt", A = 0.9, w = 0.005, prob_infected_by_age_group = None, frac_pop_infected = 0.8, t_omicron_wave1_mean = 90, sd_omicron_wave1 = 30, t_omicron_wave2_mean = 212, sd_omicron_wave2 = 30, weight_wave2 = 0.35, t_out = 304, prints = False, plots=False):
    """
    Input arguments:
        scenario_name, # Name will be prepended to output files
        age_file = "age_FHI_2022.txt", # Population file
        A = 0.9, # Max value of protection (VE or infection-induced)
        w = 0.005, # Decay rate of protection
        prob_infected_by_age_group = None,
        frac_pop_infected = 0.8,
        t_omicron_wave1_mean = 90,  # days -- 1st of March 2022 counting from 1st december 2021
        sd_omicron_wave1 = 30,  # days
        t_omicron_wave2_mean = 212,  # days -- 1st of July 2022 counting from 1st december 2021
        sd_omicron_wave2 = 30,  # days
        t_out = 304, # Day to evaluate (Default October 1st counting from 2022-12-01)
        prints = False, # Whether to print info as you go
        plots=False  # Whether to plot
    """

    if prints:
        print("\n\nNow running scenario", scenario_name)


    # Make folders
    if not os.path.isdir("scenarios"):
        os.makedirs("scenarios")
    if not os.path.isdir("plots"):
        os.makedirs("plots")


    rng = np.random.default_rng() # Instantiate random number generator 

    # Load population file
    df_ages = pd.read_csv(age_file, header=None)
    df_ages.columns = ['pop']
    df_ages['age'] = range(len(df_ages))
    df_ages['age_group'] = np.minimum(8, (df_ages['age'] // 10).astype(int))
    if prints:
        print(df_ages)


    N_pop = sum(df_ages['pop'])
    if prints:
        print(N_pop)


    # Age group population
    df_ag = df_ages.groupby('age_group', as_index = False).agg({'pop': 'sum'})
    df_ag



    #%%
    # Make dataframe with one row per individual
    ages_array = np.zeros(N_pop)
    i_age = 0
    i = 0
    counter_age = 0
    ages_array_list = []
    # Make a lookup list of the starting index in ages_array corresponding to each age
    ages_index_lookup = [0]
    for i_age in range(len(df_ages)):
        ages_array_list.append(np.ones(df_ages['pop'][i_age])*i_age)
        if i_age < len(df_ages): # not last age
            ages_index_lookup.append(ages_index_lookup[-1] + df_ages['pop'][i_age]) # Add element to list which has the index of the first individual in each age
    ages_array = np.concatenate(ages_array_list)


    if prints:
        print(len(ages_array))
        print(ages_array)

    # Array -> dataframe with some more properties
    df = pd.DataFrame({
            "age": ages_array,
            "infected": np.zeros(N_pop),
            "t_infected": np.ones(N_pop)*9999,
            # "vaccinated": np.zeros(N_pop),
            "t_vaccinated": np.ones(N_pop)*9999
        })
    #%%
    def VE_individual(t, t_vax, t_inf, A, w):
        # A model where effective VE is max(VE_vax, VE_inf) * hybrid_factor, hybrid_factor = 1.1 if both vax and inf, 0 else?
        VE_vax = np.where(t > t_vax, A*np.exp(-w*(t-t_vax)), 0)
        VE_inf = np.where(t > t_inf, A*np.exp(-w*(t-t_inf)), 0)
        return np.minimum(1, np.maximum(VE_vax, VE_inf))# * np.where(((VE_vax > 0) & (VE_inf > 0)), 1.1, 1))

    if plots:
        # Test plot VE for an individual with hybrid immunity
        tt = np.linspace(0, 300, 301)
        tt_vax = np.ones(len(tt))*50
        tt_inf = np.ones(len(tt))*120
        f, ax = plt.subplots(1)
        ax.plot(tt, VE_individual(tt, tt_vax, tt_inf, A, w))
    #%%

    ## Infection-induced immunity:
    # Sample time of infection in two steps
    # Step 1. age dependent probability of having been infected during past year
    if prob_infected_by_age_group is None:
        prob_infected_by_age_group = { # key is 1/10 of lower bound of interval, e.g. 0 means 0-9 years
            0: 0.55,  # 0-9   # Assuming same as 16-29 which is the lowest age category in Symptometer 
            1: 0.55,  # 10-19 # Assuming same as 16-29 which is the lowest age category in Symptometer 
            2: 0.55,  # 20-29 # Assuming same as 16-29 which is the lowest age category in Symptometer 
            3: 0.725, # 30-39
            4: 0.7,   # 40-49  
            5: 0.625, # 50-59
            6: 0.55,  # 60-69
            7: 0.46,  # 70-79 # Same value for all 70+
            8: 0.46,  # 80+ # Same value for all 70+
            # 9: 0.46   # 90+   # Same value for all 70+
        }
    # # Renormalise overall probability to allow fraction frac_pop_infected of population to be infected in total
    raw_probs = np.array(list(prob_infected_by_age_group.values()))
    rescaled_probs = raw_probs * frac_pop_infected * np.sum(df_ag['pop']) / np.sum(raw_probs * df_ag['pop'])
    rescaled_probs = np.where(rescaled_probs < 1, rescaled_probs, 1) # Ensure no larger-than-1 probabilities
    prob_infected_by_age_group = {i: rescaled_probs[i] for i in range(9)}
    
    # Sample random people who will get infected
    df['age_group'] = (df['age'] // 10).astype(int)
    df.loc[df['age_group'] > 8, 'age_group'] = 8
    u,inv = np.unique(df['age_group'],return_inverse = True)
    df['prob_infected'] = np.array([prob_infected_by_age_group[x] for x in u])[inv]
    df['infected'] = np.random.binomial(1, df['prob_infected'])

    #%%
    # Step 2. sum of gaussians to determine time of infection
    # t_omicron_wave1_mean = 90 # days -- 1st of March 2022 counting from 1st december 2021
    # sd_omicron_wave1 = 30 # days
    # t_omicron_wave2_mean = 212 # days -- 1st of July 2022 counting from 1st december 2021
    # sd_omicron_wave2 = 30 # days

    if plots:
        # Plot gaussian
        tt = np.linspace(0, 300, 301)
        f, ax = plt.subplots(1)
        ax.plot(tt, norm.pdf(tt, loc = t_omicron_wave1_mean, scale = sd_omicron_wave1))
        plt.show()


    def time_of_infection(t_omicron_wave1_mean, sd_omicron_wave1, t_omicron_wave2_mean, sd_omicron_wave2):
        # Set t0 = 2021-12-01
        # Assume we had a wave with peak infections 1st of March, standard deviation 1 month (so width on either side of peak is ~2 months)
        # Draw timepoint from this dist.
        choose_wave_1_or_2 = rng.choice(
            np.array([1, 2]),
            size = N_pop, 
            p = np.array([1 - weight_wave2, weight_wave2]) # Relative size here is determined by sum of hosps in each of the two omicron waves. Update: Seems from Symptometer that more of the weight should be shifted towards the summer peak. Could be due to reinfections and also immunity reducing IHR compared with winter wave.
        )
        return np.where(
            choose_wave_1_or_2 == 1,
            np.random.normal(loc = t_omicron_wave1_mean, scale = sd_omicron_wave1, size = N_pop),
            np.random.normal(loc = t_omicron_wave2_mean, scale = sd_omicron_wave2, size = N_pop)
        )
    df['t_infected'] = np.where(df['infected'] == 1, time_of_infection(t_omicron_wave1_mean, sd_omicron_wave1, t_omicron_wave2_mean, sd_omicron_wave2), 9999)

    if plots:
        # Check resulting distribution of infection times
        plt.figure()
        plt.hist(df['t_infected'][df['t_infected']<9999], bins=100)
        plt.savefig("plots/"+scenario_name+"-dist_infection_times.jpg")

    #%%
    ## Sample time of vaccination. Assume random people independent of vaccination status.
    # 1. Have used SYSVAK to make an aggregated list of vaccine doses by day and by age
    # 2. For each dose, sample random individuals who get it. Set their time of vaccination to that timepoint.
    # NB! Due to privacy protection, we are not sharing the true file containing number of vaccine doses. See the readme file located in this folder.
    dfvax = pd.read_csv("sysvak_number_of_doses_by_date_age_dosenumber-template.csv")
    # Convert date to time measured from 2021-12-01:
    dfvax['time'] = (pd.to_datetime(dfvax['date_vax']) - datetime(2022, 1, 1)).astype('timedelta64[D]')
    # Select only dose 3 and 4
    dfvax_d3 = dfvax[dfvax['dose'] == 3]
    dfvax_d4 = dfvax[dfvax['dose'] == 4]
    if prints:
        print(len(dfvax))
        print("len(dfvax_d3) =", len(dfvax_d3))
        print("len(dfvax_d4) =", len(dfvax_d4))

    # As a simplistic assumption, remove as many counts from d3 as there are d4, per age 
    list_df_d3 = []
    for age in range(len(df_ages)):
        # Make long df with one row per shot, for current age, dose 3, to subsample
        dfvax_d3_curr = dfvax_d3[dfvax_d3['age'] == age]
        dfd3l_curr = dfvax_d3_curr.loc[dfvax_d3_curr.index.repeat(dfvax_d3_curr['n'])].assign(n_1=1).reset_index(drop=True).drop(columns = ['n']).rename(columns={'n_1': 'n'})
        n_subsample = sum(dfvax_d3_curr['n']) - sum(dfvax_d4['n'][dfvax_d4['age'] == age])

        if (n_subsample <= 0): # No 3rd doses left
            continue

        list_df_d3.append(dfd3l_curr.sample(n_subsample)) # Sample n_subsample rows from the one-row-per-shot dose3 dataframe for current age

    # Gather all remaining dose3 into single dataframe
    dfvax_d3_shortened = pd.concat(list_df_d3)
    # Concat it with the dose4 dataframe to give the complete list of all doses
    dfvax = pd.concat([dfvax_d3_shortened, dfvax_d4])


    #%%
    # Assign doses to people
    for age in range(len(df_ages)):
        dfvax_currage = dfvax[dfvax['age'] == age]
        n_doses = sum(dfvax_currage['n'])
        if n_doses < 1: # Skip age if no doses
            continue
        n_doses = min(df_ages['pop'][age], n_doses)
        # Draw random individuals who get this vaccine
        # First draw a random subset of n_doses individuals within range (ages_index_lookup[age], ages_index_lookup[age+1]-1)
        # Then sample the vaccination time vector with weights according to n_dose, to give each individual a vaccine time
        indices = rng.choice(
            np.arange(start = ages_index_lookup[age], stop = ages_index_lookup[age+1], step = 1),
            size = n_doses,
            replace = False
        )
        # df.loc[indices, 'vaccinated'] = 1 # Set vaccination status to 1 for chosen individuals. 
                                        # (Not strictly needed I think, could just set 
                                        # t_vaccinated = 9999 for everybody by default)
        df.loc[indices, 't_vaccinated'] = rng.choice( # Sample time of vaccination. This is easier than assigning each time to exactly n_doses persons, and should give excellent approximation.
            dfvax_currage['time'], 
            size = n_doses, 
            p = dfvax_currage['n']/np.sum(dfvax_currage['n'])
        )


    #%%
    if prints:
        # Print some vax and inf statistics
        print("Share of infected =", len(df[df['t_infected'] < 9999]) / N_pop)
        print("Share of vaccinated =", len(df[df['t_vaccinated'] < 9999]) / N_pop)
        print("Share of infected AND vaccinated =", len(df[((df['t_vaccinated'] < 9999) & (df['t_infected'] < 9999))]) / N_pop)
        print("Share of infected OR  vaccinated =", len(df[((df['t_vaccinated'] < 9999) | (df['t_infected'] < 9999))]) / N_pop)
        print("Share of infected since 2021-12-01 =", len(df[(df['t_infected'] < 9999) & (df['t_infected'] > 0)]) / N_pop)
        print("Share of vaccinated since 2021-12-01 =", len(df[(df['t_vaccinated'] < 9999) & (df['t_vaccinated'] > 0)]) / N_pop)
        print("Share of infected AND vaccinated since 2021-12-01 =", len(df[((df['t_vaccinated'] < 9999) & (df['t_vaccinated'] > 0) & (df['t_infected'] < 9999) & (df['t_infected'] > 0))]) / N_pop)
        print("Share of infected OR vaccinated since 2021-12-01 =", len(df[(((df['t_vaccinated'] < 9999) & (df['t_vaccinated'] > 0)) | ((df['t_infected'] < 9999) & (df['t_infected'] > 0)))]) / N_pop)
        print("Share of infected since 2021-06-01 =", len(df[(df['t_infected'] < 9999) & (df['t_infected'] > 182)]) / N_pop)
        print("Share of vaccinated since 2021-06-01 =", len(df[(df['t_vaccinated'] < 9999) & (df['t_vaccinated'] > 182)]) / N_pop)
        print("Share of infected AND vaccinated since 2021-06-01 =", len(df[((df['t_vaccinated'] < 9999) & (df['t_vaccinated'] > 182) & (df['t_infected'] < 9999) & (df['t_infected'] > 182))]) / N_pop)
        print("Share of infected OR vaccinated since 2021-06-01 =", len(df[(((df['t_vaccinated'] < 9999) & (df['t_vaccinated'] > 182)) | ((df['t_infected'] < 9999) & (df['t_infected'] > 182)))]) / N_pop)
        


    if plots:
        # Check resulting distribution of vaccination times (last dose)
        plt.figure()
        plt.hist(df['t_vaccinated'][df['t_vaccinated']<9999], bins=100)
        plt.savefig("plots/"+scenario_name+"-dist_vaccination_times.jpg")

    #%%
    # Calculate VE over time for each individual; aggregate immunity to age groups
    # Set Tmax to be 1st of july 2023. We define t=0 to be 1st of december 2021, so then Ndays is 577
    Ndays = 577
    t_analysis = np.linspace(0, Ndays, Ndays+1) 
    df_time = pd.DataFrame({
        "time": t_analysis
    })

    # Subsample population to save memory
    df_s = df.sample(frac=0.01)

    df_big = df_s.merge(df_time, how='cross')

    df_big['protection'] = VE_individual(
        df_big['time'], 
        df_big['t_vaccinated'], 
        df_big['t_infected'],
        A,
        w
    )

    # Plot mean VE over time
    df_mean = df_big.groupby(['time'], as_index=False).agg({'protection': 'mean'})
    df_mean['date'] = pd.date_range(start='2021-12-01', periods=len(df_mean), freq='D')

    if plots:
        f,ax = plt.subplots(1)
        ax.plot(df_mean['time'], df_mean['protection'])
        plt.show()


    if prints:
        print("Mean value of VE among all age groups on day 304 (= Oct 1st) is", df_mean['protection'][303])

    #%%
    # Quantiles
    # # Define quantile functions. Has to be done this manual, stupid way due to pandas agg function
    # def q25(x):
    #     return x.quantile(0.25)
    # def q50(x):
    #     return x.quantile(0.50)
    # def q75(x):
    #     return x.quantile(0.75)
    # #df_avg_over_time = df_big.groupby(['time'], as_index=False).agg({'protection': 'mean'})
    # df_quantiles_over_time = df_big[['time', 'protection']].groupby(['time'], as_index=False).agg({'protection': [q25, q50, q75]})
    # df_quantiles_over_time['date'] = pd.date_range(start='2021-12-01', periods=len(df_quantiles_over_time), freq='D')

    # f, ax = plt.subplots(1)
    # ax.plot(df_quantiles_over_time['date'], df_quantiles_over_time['protection']['q50'])
    # ax.fill_between(
    #     df_quantiles_over_time['date'],
    #     y1 = df_quantiles_over_time['protection']['q25'],
    #     y2 = df_quantiles_over_time['protection']['q75'],
    #     alpha = 0.2
    # )
    # ax.tick_params(axis = 'x', labelrotation = 45)
    # plt.show()

    #%%
    # Plot mean VE over time by age grup
    df_avg_over_time_by_age_group = df_big.groupby(['time', 'age_group'], as_index=False).agg({'protection': ['mean']})#,'std']})
    # df_avg_over_time_by_age_group = df_big[['time', 'age_group', 'protection']].groupby(['time', 'age_group'], as_index=False).agg({'protection': [q25, q50, q75]})
    f, ax = plt.subplots(1)
    for age_group in range(9):
        dfp = df_avg_over_time_by_age_group[df_avg_over_time_by_age_group['age_group'] == age_group]
        dfp['date'] = pd.date_range(start='2021-12-01', periods=len(dfp), freq='D')
        # ax.plot(dfp['time'], dfp['protection']['q50'], label='age_group {:d}'.format(age_group))
        ax.plot(dfp['date'], dfp['protection']['mean'], label='age_group {:d}'.format(age_group))
    ax.legend(loc='lower right', prop={'size': 6})
    ax.tick_params(axis = 'x', labelrotation = 45)
    plt.tight_layout()
    plt.savefig("plots/"+scenario_name+"-avg_protection_by_age_group.jpg", dpi=300, transparent=False)
    if plots:
        plt.show()



    # %%
    # Take out VE values on 1st of october 2022 = day 304 counting from 2021-12-01.
    # t_out = 304 # October 1st counting from 2022-12-01
    dfout = df_big[df_big['time'] == t_out]
    dfout['protection_rounded'] = np.minimum(0.9, np.round(dfout['protection']-0.05, decimals=1)) # Rounding protection into categories, in such a way that all in interval [0.8, 0.9) -> 0.8, etc., with max interval 0.9 since it's given by lower bin edge.

    # Calculate proportion in each protection category by age group
    dfout_age = dfout.groupby(['age_group'], as_index=False).size().rename(columns={'size': 'size_age'})
    dfout_age_and_protection = dfout.groupby(['age_group', 'protection_rounded'], as_index=False).size() # Number of people
    dfout = dfout_age_and_protection.merge(dfout_age, how='left', on='age_group')
    dfout['proportion'] = dfout['size'] / dfout['size_age']


    # Plot histogram of each age group's immunity distribution
    f, ax = plt.subplots(3,3)
    for age_group in range(9):
        ix = age_group // 3
        iy = age_group % 3
        dfp = dfout[dfout['age_group'] == age_group]
        ax[ix, iy].bar(dfp['protection_rounded'], dfp['proportion'], width=0.1, align='edge')
        ax[ix, iy].set_title('Age group '+str(age_group))
    plt.tight_layout()
    plt.savefig("plots/"+scenario_name +
                "-immunity_distribution_by_age_group.jpg", dpi=300, transparent=False)
    if plots:
        plt.show()

    # Write to csv
    dfout = dfout[['age_group', 'protection_rounded', 'proportion']]
    # Add missing combinations for csv file consistency (will be 0)
    if len(dfout) < 9 * 10: # Not all combinations of age group and protection percentile present
        dfout_missing = {
            'age_group': [],
            'protection_rounded': [],
            'proportion': []
        }
        for age_group in range(9):
            for protection_group in range(10):
                protection = 0.1 * protection_group
                if np.min(np.abs(dfout[dfout['age_group'] == age_group]['protection_rounded'] - protection)) > 1e-3: # Float comparison between each category...
                    dfout_missing['age_group'].append(age_group)
                    dfout_missing['protection_rounded'].append(protection)
                    dfout_missing['proportion'].append(0)
        dfout_missing = pd.DataFrame(dfout_missing)
        dfout = pd.concat([dfout, dfout_missing]).sort_values(['age_group', 'protection_rounded'])
    dfout.to_csv(
        "scenarios/"+scenario_name + "-proportion_in_protection_compartments_by_age_group.csv", index = False)
    # Also estimate time of transition between each compartment...based on plot of mean VE?

    # %%
    if plots:
        # Understanding transition rate from one compartment to next
        # Plot waning curve together with percentiles
        f, ax = plt.subplots(1)
        ax.plot(tt, A*np.exp(-w*tt))
        for yval in 0.1*np.linspace(0,9,10):
            ax.plot(tt, np.ones(len(tt))*yval, color="grey")
        plt.show()


    # %%
    if plots:
        # Plot avg VE for one age group together with percentiles, and try to fit an exponential function 
        # by eye to see if individual-level decay rate matches pop avg. It does, as it should if things have been done correctly.
        age_group = 6
        dfp = df_avg_over_time_by_age_group[df_avg_over_time_by_age_group['age_group'] == age_group]
        f, ax = plt.subplots(1)
        ax.plot(dfp['time'], dfp['protection']['mean'], label='age_group {:d}'.format(age_group))
        for yval in 0.1*np.linspace(0,9,10):
            ax.plot(dfp['time'], np.ones(len(dfp['time']))*yval, color="grey")
        # Try to superimpose an exponential curve to see if it matches
        tp = dfp['time']
        ax.plot(tp, np.where(tp > 300, 0.35*np.exp(-w*(tp-300)), 0), color='red')
        plt.show()
    # %%
    # Calculating time diff between VE value q1 and q2 based on formula
    # delta_t = (ln(q1)-ln(q2))/w
    dict_waning_time = {
        "compartment_upper": [],
        "compartment_lower": [],
        "delta_t": []
    }
    if prints:
        print("Given decay rate w =", w, "days^-1")
    for q1 in np.arange(1, 0.1, -0.1):
        q2 = q1 - 0.1
        delta_t = (np.log(q1 - 0.05) - np.log(q2 - 0.05)) / w # Using mid-point value in each compartment
        dict_waning_time['compartment_upper'].append(q1)
        dict_waning_time['compartment_lower'].append(q2)
        dict_waning_time['delta_t'].append(delta_t)
        if prints:
            print("Mean waning time out of VE compartment ({:.1f}, {:.1f}) = {:f}".format(q2, q1, delta_t))
    # Write to df
    df_waning_time = pd.DataFrame(dict_waning_time)
    df_waning_time.to_csv("scenarios/"+scenario_name + "-waning_times.csv")

    # %%


if __name__ == "__main__":

    prob_infected_by_age_group = {  # key is 1/10 of lower bound of interval, e.g. 0 means 0-9 years
        0: 0.55,  # 0-9   # Assuming same as 16-29 which is the lowest age category in Symptometer
        1: 0.55,  # 10-19 # Assuming same as 16-29 which is the lowest age category in Symptometer
        2: 0.55,  # 20-29 # Assuming same as 16-29 which is the lowest age category in Symptometer
        3: 0.725,  # 30-39
        4: 0.7,   # 40-49
        5: 0.625,  # 50-59
        6: 0.55,  # 60-69
        7: 0.46,  # 70-79 # Same value for all 70+
        8: 0.46,  # 80+ # Same value for all 70+
        # 9: 0.46   # 90+   # Same value for all 70+
    }
    # Renormalise overall probability to allow fraction frac of population to be infected in total
    frac_pop_infected = 0.8

    t_out = 304  # October 1st counting from 2022-12-01
    simple_immunity_model(
        scenario_name = "test_scenario",  # Name will be prepended to output files
        age_file="age_FHI_2022.txt",  # Population file
        A=0.9,  # Max value of protection (VE or infection-induced)
        w=0.005,  # Decay rate of protection
        prob_infected_by_age_group = None,
        frac_pop_infected = 0.8,
        t_omicron_wave1_mean = 90,  # days -- 1st of March 2022 counting from 1st december 2021
        sd_omicron_wave1 = 30,  # days
        t_omicron_wave2_mean = 212,  # days -- 1st of July 2022 counting from 1st december 2021
        sd_omicron_wave2 = 30,  # days
        # Day to evaluate (Default October 1st counting from 2022-12-01)
        t_out = 304,
        prints = True,  # Whether to print info as you go
        plots = True  # Whether to plot
    )