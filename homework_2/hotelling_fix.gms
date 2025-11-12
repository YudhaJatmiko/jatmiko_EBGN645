*Hotelling Model for Nickel Extraction

Sets
t 'time period' /1*30/
countries 'nickel production per country' /us, australia, brazil, canada, china, indonesia, newcaledonia, philippines, russia, other/
scenario 'price scenarios' /oversupply, balanced, supply_constraint, ev_boom/
;

*production and reserve data
table prod_reserve(countries,*)
                 prod2022  prod2023  reserves
us                 17500     17000    340000
australia         155000    160000  24000000
brazil             88500     89000  16000000
canada            143000    180000   2200000
china             114000    110000   4200000
indonesia        1580000   1800000  55000000
newcaledonia      200000    230000   7100000
philippines       345000    400000   4800000
russia            222000    200000   8300000
other             404000    380000   9100000
;

*parameter declaration
parameter
p0 'initial nickel price (USD per tonne)'
r0 'total initial reserves (tonnes)'
q0 'initial monthly extraction rate (tonnes)'
c 'extraction cost (USD per tonne)'
disc_rate 'discount rate'
disc(t) 'discount factor by period'
;

*assumptions
*cost assumption: 1000 USD per tonne
*discount rate: 5% per period
*initial price: latest nickel price from June 2025 data
p0 = 15003.11;
r0 = sum(countries, prod_reserve(countries,'reserves'));
q0 = sum(countries, prod_reserve(countries,'prod2023')) / 12;
c = 1000;
disc_rate = 0.05;

*discount factor calculation
disc(t) = 1 / power(1 + disc_rate, ord(t) - 1);
display p0, r0, q0;

*new parameter (indonesia and scenarios)
parameter
r0_INA 'Indonesian reserves only'
q0_INA 'Indonesian annual production'
growth_rate(scenario) 'annual price growth rates'
prob_scenario(scenario) 'scenario probabilities'
price_path(scenario,t) 'price evolution by scenario'
;

*Indonesia specific parameters
r0_INA = prod_reserve('indonesia','reserves');
q0_INA = prod_reserve('indonesia','prod2023');


*price scenario (based on IEA mineral outlook), I dont find nickel price forecast
growth_rate('oversupply') = -0.03;
growth_rate('balanced') = 0.02;
growth_rate('supply_constraint') = 0.04;
growth_rate('ev_boom') = 0.06;

*probability of scenario
prob_scenario('oversupply') = 0.25;
prob_scenario('balanced') = 0.40;
prob_scenario('supply_constraint') = 0.20;
prob_scenario('ev_boom') = 0.15;

*price path based on the scenario
price_path(scenario,t) = p0 * power(1 + growth_rate(scenario), ord(t) - 1);


*variables
positive variables
q(t) 'quantity of extraction in period t (tonnes)'
reserve(t) 'remaining reserves by the end of period t (tonnes)'
q_INA(scenario,t) 'Indonesian extraction by scenario'
reserve_INA(scenario,t) 'Indonesian reserves by scenario'
;

free variables
NPV 'net present value (USD)'
NPV_INA(scenario) 'Indonesian NPV by scenario'
expected_NPV_INA 'expected Indonesian NPV'
;


*equations
equations
obj 'objective function - maximize NPV'
res_init 'initial reserve constraint'
res_bal(t) 'reserve balance equation'
*Indonesia Equations
obj_INA(scenario) 'Indonesian NPV by scenario'
exp_obj_INA 'expected Indonesian NPV'
res_init_INA(scenario) 'Indonesian initial reserves'
res_bal_INA(scenario,t) 'Indonesian reserve balance'
;


*NPV maximization objective function
obj.. NPV =e= sum(t, disc(t) * (p0 * q(t) - c * q(t)));

*initial reserve constraint
res_init.. reserve('1') =e= r0 - q('1');

*reserve balance equation
res_bal(t)$(ord(t) > 1).. reserve(t) =e= reserve(t-1) - q(t);


*Indonesia scenario
obj_INA(scenario).. NPV_INA(scenario) =e= sum(t, disc(t) * (price_path(scenario,t) - c) * q_INA(scenario,t));
exp_obj_INA.. expected_NPV_INA =e= sum(scenario, prob_scenario(scenario) * NPV_INA(scenario));
res_init_INA(scenario).. reserve_INA(scenario,'1') =e= r0_INA - q_INA(scenario,'1');
res_bal_INA(scenario,t)$(ord(t) > 1).. reserve_INA(scenario,t) =e= reserve_INA(scenario,t-1) - q_INA(scenario,t);


*bounds and initial values
q.lo(t) = 0;
q.up(t) = q0 * 2;
q.l(t) = q0;
reserve.lo(t) = 0;
reserve.l(t) = r0 - q0 * ord(t);

*Indonesian bounds
q_INA.lo(scenario,t) = 0;
q_INA.up(scenario,t) = q0_INA * 2;
q_INA.l(scenario,t) = q0_INA;
reserve_INA.lo(scenario,t) = 0;
reserve_INA.l(scenario,t) = r0_INA - q0_INA * ord(t);

*solve model
model hotelling /all/;
solve hotelling using nlp maximizing NPV;

*solve indonesia model
model hotelling_INA /obj_INA, exp_obj_INA, res_init_INA, res_bal_INA/;
solve hotelling_INA using nlp maximizing expected_NPV_INA;

*post-processing calculations
alias(t, t_prev);

parameter
cum_extraction(t) 'cumulative extraction by period t (tonnes)'
depletion 'percentage of reserves depleted (%)'
;

*Indonesia analysis
parameter
early_extraction_INA(scenario) 'Indonesian extraction years 1-10'
late_extraction_INA(scenario) 'Indonesian extraction years 11-20'
timing_ratio(scenario) 'early vs late extraction ratio'
final_price(scenario) 'price in year 30'
;


cum_extraction(t) = sum(t_prev$(ord(t_prev) <= ord(t)), q.l(t_prev));
depletion = 100 * sum(t, q.l(t)) / r0;
*new timing
early_extraction_INA(scenario) = sum(t$(ord(t) <= 10), q_INA.l(scenario,t));
late_extraction_INA(scenario) = sum(t$(ord(t) >= 11), q_INA.l(scenario,t));
timing_ratio(scenario) = early_extraction_INA(scenario) / (late_extraction_INA(scenario) + 0.001);
final_price(scenario) = price_path(scenario,'30');


*display results
display
hotelling.modelstat, hotelling.solvestat,
NPV.l, q.l, reserve.l, cum_extraction, depletion,
expected_NPV_INA.l,
NPV_INA.l,
timing_ratio,
final_price
;


*export results
execute_unload "result.gdx";
execute_unload "indonesia_results.gdx" q_INA, NPV_INA, timing_ratio;

*result
*VARIABLE NPV.L                 =  1.34334E+11  net present value (USD)
*PARAMETER depletion            =       13.607  percentage of reserves depleted (%)
*VAR q  quantity of extraction in period t (tonnes) 594333.3333
*Time value money 3402.0041 in period 30th
