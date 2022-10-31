import pandas as pd
import numpy as np
import re


filenames = [
    "scenarios/sc1_medium-w0p005_A0p9_fracinfected0p9_wsecondwave0p5-proportion_in_protection_compartments_by_age_group.csv",
    "scenarios/sc2_low-w0p0075_A0p9_fracinfected0p8_wsecondwave0p35-proportion_in_protection_compartments_by_age_group.csv",
    "scenarios/sc3_high-w0p003_A0p9_fracinfected0p9_wsecondwave0p6-proportion_in_protection_compartments_by_age_group.csv"
]

for filename in filenames:
    df = pd.read_csv(filename)
    # Add proportion-weighted protection value for summing below:
    df['proportion_weighted_protection'] = df['protection_rounded'] * df['proportion']

    # Collapse to single, average value each age by summing proportion-weighted protection:
    dfc = df.groupby("age_group", as_index=False).agg({"proportion_weighted_protection": 'sum'})

    # Fix names and add necessary columns
    dfc = dfc.rename(columns={"proportion_weighted_protection": "protection_rounded"})
    dfc['protection_rounded'] = np.round(dfc['protection_rounded'], decimals=1) # TODO should this rounding be done on 0.05-shifted values to preserve midpoint of bin value?
    dfc['proportion'] = 1

    # Add missing combinations for csv file consistency (will be 0)
    if len(dfc) < 9 * 10:  # Not all combinations of age group and protection percentile present
        dfc_missing = {
            'age_group': [],
            'protection_rounded': [],
            'proportion': []
            }
        for age_group in range(9):
            for protection_group in range(10):
                protection = 0.1 * protection_group
                # Float comparison between each category...
                if np.min(np.abs(dfc[dfc['age_group'] == age_group]['protection_rounded'] - protection)) > 1e-3:
                    dfc_missing['age_group'].append(age_group)
                    dfc_missing['protection_rounded'].append(protection)
                    dfc_missing['proportion'].append(0)
        dfc_missing = pd.DataFrame(dfc_missing)
        dfc = pd.concat([dfc, dfc_missing]).sort_values(
            ['age_group', 'protection_rounded'])

    dfc.to_csv(re.sub(".csv", "", filename)+"-homogenous_immunity.csv", index = False)
