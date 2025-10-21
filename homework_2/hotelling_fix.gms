*Hotelling Model for Nickel Extraction

Sets
t 'time period' /1*30/
countries 'nickel production per country' /us, australia, brazil, canada, china, indonesia, newcaledonia, philippines, russia, other/
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

*variables
positive variable
q(t) 'quantity of extraction in period t (tonnes)'
reserve(t) 'remaining reserves by the end of period t (tonnes)';

free variable
NPV 'net present value (USD)';

*equations
equations
obj 'objective function - maximize NPV'
res_init 'initial reserve constraint'
res_bal(t) 'reserve balance equation';

*NPV maximization objective function
obj.. NPV =e= sum(t, disc(t) * (p0 * q(t) - c * q(t)));

*initial reserve constraint
res_init.. reserve('1') =e= r0 - q('1');

*reserve balance equation
res_bal(t)$(ord(t) > 1).. reserve(t) =e= reserve(t-1) - q(t);

*bounds and initial values
q.lo(t) = 0;
q.up(t) = q0 * 2;
q.l(t) = q0;
reserve.lo(t) = 0;
reserve.l(t) = r0 - q0 * ord(t);

*solve model
model hotelling /all/;
solve hotelling using nlp maximizing NPV;

*post-processing calculations
alias(t, t_prev);

parameter
cum_extraction(t) 'cumulative extraction by period t (tonnes)'
depletion 'percentage of reserves depleted (%)'
;

cum_extraction(t) = sum(t_prev$(ord(t_prev) <= ord(t)), q.l(t_prev));
depletion = 100 * sum(t, q.l(t)) / r0;

*display results
display
hotelling.modelstat, hotelling.solvestat,
NPV.l, q.l, reserve.l, cum_extraction, depletion
;

*export results
execute_unload "result.gdx";


*result
*VARIABLE NPV.L                 =  1.34334E+11  net present value (USD)
*PARAMETER depletion            =       13.607  percentage of reserves depleted (%)
*VAR q  quantity of extraction in period t (tonnes) 594333.3333
*Time value money 3402.0041 in period 30th
