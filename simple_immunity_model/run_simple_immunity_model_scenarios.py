#%%
from simple_immunity_model import simple_immunity_model



# Set up scenarios
scenarios = {
    "sc3_high-w0p003_A0p9_fracinfected0p9_wsecondwave0p6": {
        "frac_pop_infected": 0.95,
        "A": 0.9,
        "w": 0.003,
        "weight_wave2": 0.6
    },
    "sc1_medium-w0p005_A0p9_fracinfected0p9_wsecondwave0p5": {
        "frac_pop_infected": 0.9,
        "A": 0.9,
        "w": 0.005,
        "weight_wave2": 0.5
    },
    "sc2_low-w0p0075_A0p9_fracinfected0p8_wsecondwave0p35": {
        "frac_pop_infected": 0.8,
        "A": 0.9,
        "w": 0.0075,
        "weight_wave2": 0.35
    }
}

for scname, scdict in scenarios.items():
    simple_immunity_model(
        scenario_name=scname,
        A = scdict["A"],
        w = scdict["w"],
        frac_pop_infected=scdict["frac_pop_infected"],
        weight_wave2 = scdict["weight_wave2"],
        plots = False,
        prints = True
    )