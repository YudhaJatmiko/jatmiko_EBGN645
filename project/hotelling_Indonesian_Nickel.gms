* HOTELLING MODEL FOR INDONESIAN NICKEL EXTRACTION
* Author: Yudha Jatmiko
* Course: EBGN 645 - Computational Economics (Fall2025)

* Sets
Sets
    t           'time periods (years)' /1*30/
    countries   'nickel producing countries' /
                    us, australia, brazil, canada, china,
                    indonesia, newcaledonia, philippines, russia, other
                /
    scenario    'price scenarios' /
                    oversupply, baseline, supply_constraint, ev_boom
                /
;

Alias(t, t_prev);

* Production and reserves data
Table prod_reserve(countries,*)
                 prod2022    prod2023    reserves
us                 17500       17000      340000
australia         155000      160000    24000000
brazil             88500       89000    16000000
canada            143000      180000     2200000
china             114000      110000     4200000
indonesia        1580000     1800000    55000000
newcaledonia      200000      230000     7100000
philippines       345000      400000     4800000
russia            222000      200000     8300000
other             404000      380000     9100000
;

* Price forecast data from CSV
$ondelim
Table price_forecast(t,scenario)
$include nickel_prices_gams.csv
$offdelim

* Cost function parameters (calibrated to Indonesia)
Parameters
    beta1       'base extraction cost (USD per tonne)' /13500/
    b_tech      'technology improvement factor' /0.92/
    gamma       'stock depletion cost penalty (USD per tonne)' /6750/
    disc_rate   'annual discount rate' /0.11/
    p0          'initial nickel price (USD per tonne)' /15003.11/
;

* Calculated parameters
Parameters
    r0              'total global initial reserves'
    q0              'global monthly extraction rate'
    r0_INA          'Indonesian initial reserves'
    q0_INA          'Indonesian annual production'
    disc(t)         'discount factor'
    tech_factor(t)  'technology cost multiplier'
    price_path(scenario,t)      'price evolution by scenario'
    hotelling_price(t)          'Hotelling price path for global model'
;

* Initialize parameters
r0 = sum(countries, prod_reserve(countries,'reserves'));
q0 = sum(countries, prod_reserve(countries,'prod2023')) / 12;
r0_INA = prod_reserve('indonesia','reserves');
q0_INA = prod_reserve('indonesia','prod2023');
disc(t) = 1 / power(1 + disc_rate, ord(t) - 1);
tech_factor(t) = power(b_tech, ord(t) - 1);
hotelling_price(t) = p0 * power(1 + disc_rate, ord(t) - 1);

* Load price paths from CSV
price_path(scenario,t) = price_forecast(t,scenario);

* Variables
Positive Variables
    q(t)                        'global extraction'
    reserve(t)                  'global remaining reserves'
    global_unit_cost(t)         'global marginal extraction cost'
    q_INA(scenario,t)           'Indonesian extraction by scenario'
    reserve_INA(scenario,t)     'Indonesian reserves by scenario'
    unit_cost(scenario,t)       'marginal extraction cost'
    cost_variable(scenario,t)   'variable cost component'
    cost_depletion(scenario,t)  'stock depletion cost'
    scarcity_rent(scenario,t)   'scarcity rent'
;

Free Variables
    NPV                 'global net present value'
    NPV_INA(scenario)   'Indonesian NPV by scenario'
    NPV_INA_total       'total NPV across all scenarios'
;

* Equations - Global model
Equations
    obj_global              'global NPV maximization'
    global_cost_calc(t)     'global cost calculation'
    res_init_global         'global initial reserve'
    res_bal_global(t)       'global reserve balance'
;

global_cost_calc(t).. global_unit_cost(t) =e= beta1 * tech_factor(t);
obj_global.. NPV =e= sum(t, disc(t) * (hotelling_price(t) - global_unit_cost(t)) * q(t));
res_init_global.. reserve('1') =e= r0 - q('1');
res_bal_global(t)$(ord(t) > 1).. reserve(t) =e= reserve(t-1) - q(t);

* Equations - Indonesia model
Equations
    cost_var_calc(scenario,t)   'variable cost with technology improvement'
    cost_depl_calc(scenario,t)  'stock depletion cost penalty'
    cost_total_calc(scenario,t) 'total marginal cost'
    rent_calc(scenario,t)       'scarcity rent calculation'
    obj_INA(scenario)           'Indonesian NPV by scenario'
    obj_INA_total               'total NPV across all scenarios'
    res_init_INA(scenario)      'Indonesian initial reserves'
    res_bal_INA(scenario,t)     'Indonesian reserve balance'
;

cost_var_calc(scenario,t).. cost_variable(scenario,t) =e= beta1 * tech_factor(t);
cost_depl_calc(scenario,t).. cost_depletion(scenario,t) =e= gamma * (1 - reserve_INA(scenario,t) / r0_INA);
cost_total_calc(scenario,t).. unit_cost(scenario,t) =e= cost_variable(scenario,t) + cost_depletion(scenario,t);
rent_calc(scenario,t).. scarcity_rent(scenario,t) =e= price_path(scenario,t) - unit_cost(scenario,t);
obj_INA(scenario).. NPV_INA(scenario) =e= sum(t, disc(t) * scarcity_rent(scenario,t) * q_INA(scenario,t));
obj_INA_total.. NPV_INA_total =e= sum(scenario, NPV_INA(scenario));
res_init_INA(scenario).. reserve_INA(scenario,'1') =e= r0_INA - q_INA(scenario,'1');
res_bal_INA(scenario,t)$(ord(t) > 1).. reserve_INA(scenario,t) =e= reserve_INA(scenario,t-1) - q_INA(scenario,t);

* Bounds and initial values
q.lo(t) = 0;
q.up(t) = q0 * 2;
q.l(t) = q0;
reserve.lo(t) = 0;
reserve.l(t) = r0 - q0 * ord(t);
global_unit_cost.lo(t) = 0;
global_unit_cost.l(t) = beta1;

q_INA.lo(scenario,t) = q0_INA * 0.8;
q_INA.up(scenario,t) = q0_INA * 2;
q_INA.l(scenario,t) = q0_INA;
reserve_INA.lo(scenario,t) = 0;
reserve_INA.l(scenario,t) = r0_INA - q0_INA * ord(t);
unit_cost.lo(scenario,t) = 1000;
unit_cost.up(scenario,t) = 30000;
unit_cost.l(scenario,t) = 13500;
scarcity_rent.lo(scenario,t) = 0;
scarcity_rent.l(scenario,t) = p0 - 13500;

* Model definition
Model hotelling_global /
    obj_global, global_cost_calc, res_init_global, res_bal_global
/;

Model hotelling_indonesia /
    cost_var_calc, cost_depl_calc, cost_total_calc, rent_calc,
    obj_INA, obj_INA_total, res_init_INA, res_bal_INA
/;

* Solver options
option nlp = conopt;
option decimals = 2;
option limrow = 0;
option limcol = 0;

* Solve global model
Solve hotelling_global using nlp maximizing NPV;

* Solve Indonesia model for all scenarios
Solve hotelling_indonesia using nlp maximizing NPV_INA_total;

* Post-processing parameters
Parameters
    cum_extraction(t)           'cumulative global extraction'
    global_depletion            'percentage of global reserves depleted'
    avg_global_extraction       'average global extraction'
    total_extraction_INA(scenario)      'total extraction over 30 years'
    avg_annual_extraction(scenario)     'average annual extraction'
    depletion_rate(scenario)            'percentage of reserves extracted'
    depletion_year(scenario)            'year when reserves < 1%'
    early_extraction_INA(scenario)      'extraction years 1-10'
    mid_extraction_INA(scenario)        'extraction years 11-20'
    late_extraction_INA(scenario)       'extraction years 21-30'
    timing_ratio(scenario)              'early vs late extraction ratio'
    initial_price(scenario)             'price in year 1'
    final_price(scenario)               'price in year 30'
    avg_price(scenario)                 'average price over 30 years'
    price_appreciation(scenario)        'total price change'
    avg_unit_cost(scenario)             'average unit cost'
    min_unit_cost(scenario)             'minimum cost achieved'
    max_unit_cost(scenario)             'maximum cost reached'
    cost_at_year1(scenario)             'cost in year 1'
    cost_at_year30(scenario)            'cost in year 30'
    cost_change(scenario)               'percentage cost change'
    avg_variable_cost(scenario)         'average variable cost'
    avg_depletion_cost(scenario)        'average depletion cost'
    variable_cost_share(scenario)       'variable cost as % of total'
    depletion_cost_share(scenario)      'depletion cost as % of total'
    avg_scarcity_rent(scenario)         'average scarcity rent'
    rent_at_year1(scenario)             'rent in year 1'
    rent_at_year30(scenario)            'rent in year 30'
    rent_growth_rate(scenario)          'average annual rent growth'
    implied_discount_rate(scenario,t)   'implied r from rent growth'
    avg_implied_rate(scenario)          'average implied discount rate'
    hotelling_test(scenario)            'deviation from theoretical 5%'
    current_extraction_rate             'current Indonesian rate' /1800000/
    optimal_extraction_avg(scenario)    'optimal average extraction'
    extraction_deviation(scenario)      'deviation from current'
    npv_values(scenario)                'NPV values for comparison'
    npv_range                           'difference between best and worst NPV'
;

* Calculate metrics - Global model
cum_extraction(t) = sum(t_prev$(ord(t_prev) <= ord(t)), q.l(t_prev));
global_depletion = 100 * sum(t, q.l(t)) / r0;
avg_global_extraction = sum(t, q.l(t)) / 30;

* Calculate metrics - Indonesia extraction
total_extraction_INA(scenario) = sum(t, q_INA.l(scenario,t));
avg_annual_extraction(scenario) = total_extraction_INA(scenario) / 30;
depletion_rate(scenario) = 100 * total_extraction_INA(scenario) / r0_INA;

loop((scenario,t),
    if(reserve_INA.l(scenario,t) < 0.01 * r0_INA,
        depletion_year(scenario) = ord(t);
    );
);

* Calculate metrics - Extraction timing
early_extraction_INA(scenario) = sum(t$(ord(t) <= 10), q_INA.l(scenario,t));
mid_extraction_INA(scenario) = sum(t$(ord(t) >= 11 and ord(t) <= 20), q_INA.l(scenario,t));
late_extraction_INA(scenario) = sum(t$(ord(t) >= 21), q_INA.l(scenario,t));
timing_ratio(scenario) = early_extraction_INA(scenario) / (late_extraction_INA(scenario) + 0.001);

* Calculate metrics - Price analysis
initial_price(scenario) = price_path(scenario,'1');
final_price(scenario) = price_path(scenario,'30');
avg_price(scenario) = sum(t, price_path(scenario,t)) / 30;
price_appreciation(scenario) = 100 * (final_price(scenario) - initial_price(scenario)) / initial_price(scenario);

* Calculate metrics - Cost analysis
avg_unit_cost(scenario) = sum(t, unit_cost.l(scenario,t)) / 30;
min_unit_cost(scenario) = smin(t, unit_cost.l(scenario,t));
max_unit_cost(scenario) = smax(t, unit_cost.l(scenario,t));
cost_at_year1(scenario) = unit_cost.l(scenario,'1');
cost_at_year30(scenario) = unit_cost.l(scenario,'30');
cost_change(scenario) = 100 * (cost_at_year30(scenario) - cost_at_year1(scenario)) / cost_at_year1(scenario);

avg_variable_cost(scenario) = sum(t, cost_variable.l(scenario,t)) / 30;
avg_depletion_cost(scenario) = sum(t, cost_depletion.l(scenario,t)) / 30;
variable_cost_share(scenario) = 100 * avg_variable_cost(scenario) / avg_unit_cost(scenario);
depletion_cost_share(scenario) = 100 * avg_depletion_cost(scenario) / avg_unit_cost(scenario);

* Calculate metrics - Scarcity rent and Hotelling validation
avg_scarcity_rent(scenario) = sum(t, scarcity_rent.l(scenario,t)) / 30;
rent_at_year1(scenario) = scarcity_rent.l(scenario,'1');
rent_at_year30(scenario) = scarcity_rent.l(scenario,'30');
rent_growth_rate(scenario)$(rent_at_year1(scenario) > 0.01 and rent_at_year30(scenario) > 0) =
    100 * (rpower(rent_at_year30(scenario) / rent_at_year1(scenario), 1/29) - 1);
rent_growth_rate(scenario)$(rent_at_year1(scenario) <= 0.01 or rent_at_year30(scenario) <= 0) = 0;

loop((scenario,t)$(ord(t) < 30),
    if(scarcity_rent.l(scenario,t) > 100,
        implied_discount_rate(scenario,t) = scarcity_rent.l(scenario,t+1) / scarcity_rent.l(scenario,t) - 1;
    else
        implied_discount_rate(scenario,t) = 0;
    );
);

avg_implied_rate(scenario) =
    sum(t$(ord(t) < 30 and implied_discount_rate(scenario,t) > -0.20
          and implied_discount_rate(scenario,t) < 0.30),
        implied_discount_rate(scenario,t)) /
    sum(t$(ord(t) < 30 and implied_discount_rate(scenario,t) > -0.20
          and implied_discount_rate(scenario,t) < 0.30), 1);

hotelling_test(scenario) = avg_implied_rate(scenario) - disc_rate;

* Calculate metrics - Policy comparison
optimal_extraction_avg(scenario) = avg_annual_extraction(scenario);
extraction_deviation(scenario) = 100 * (optimal_extraction_avg(scenario) - current_extraction_rate) / current_extraction_rate;

npv_values(scenario) = NPV_INA.l(scenario);
npv_range = smax(scenario, npv_values(scenario)) - smin(scenario, npv_values(scenario));

* Display results
Display
    hotelling_global.modelstat, hotelling_global.solvestat,
    NPV.l, avg_global_extraction, global_depletion,
    NPV_INA.l, npv_range,
    total_extraction_INA, avg_annual_extraction, depletion_rate,
    early_extraction_INA, mid_extraction_INA, late_extraction_INA, timing_ratio,
    initial_price, final_price, avg_price, price_appreciation,
    cost_at_year1, cost_at_year30, avg_unit_cost, cost_change,
    avg_variable_cost, avg_depletion_cost, variable_cost_share, depletion_cost_share,
    avg_scarcity_rent, rent_at_year1, rent_at_year30, rent_growth_rate,
    avg_implied_rate, hotelling_test,
    current_extraction_rate, optimal_extraction_avg, extraction_deviation
;

* Export results
execute_unload "indonesia_hotelling_complete.gdx";
execute_unload "results_for_excel.gdx"
    q_INA, reserve_INA, unit_cost, scarcity_rent,
    NPV_INA, price_path, implied_discount_rate;
execute_unload "summary_stats.gdx"
    NPV_INA, avg_annual_extraction, timing_ratio,
    avg_implied_rate, hotelling_test, extraction_deviation;
