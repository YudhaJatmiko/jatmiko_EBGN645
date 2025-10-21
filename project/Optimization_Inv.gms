set i /ABX "Barrick Gold", FCX "Freeport McMoRan"/ ;

parameter
c(i) "score investment attractiveness"
/
ABX 0.60,
FCX 0.70
/ ;

parameter r(i) "expected annual return"
/
ABX 0.10,
FCX 0.20
/ ;

scalar budget "Total budget to be invested" /500000/;

positive variable X(i) "weight investment ratio"; 
variable portofolioreturn; 
variable objectivemaxvalue;

equation
eq_budget "Investment sum 1"
eq_portofolioreturn "expected return"
eq_objective "maximize score weighted to return"
;

eq_budget.. sum(i, X(i)) =e= 1 ;

eq_portofolioreturn.. =e= sum(i, X(i) * r(i)) ;

eq_objective.. objectivemaxvalue =e= sum(i, c(i) * X(i) * r(i)) ;

model simple_mining /all/ ;

solve simple_mining using lp maximizing objectivemaxvalue ;

execute_unload 'miningmaxreturn.gdx' ; 

$exit