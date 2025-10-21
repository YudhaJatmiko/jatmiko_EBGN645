
set i /hamburgers,hotdogs,frenchfries/ ; 

parameter 
r(i) "--$/item-- revenue per unit sold"
/
hamburgers 1.5,
hotdogs 0.75,
frenchfries 0.25
/ ;

parameter c(i) "--$/item-- cost per unit sold"
/
hamburgers 0.75,
hotdogs 0.1,
frenchfries 0.05
/ ;

parameter h(i) "--hours/item-- hours required per item sold"
/
hamburgers 0.5,
hotdogs 0.2,
frenchfries 0.1
/ ;

scalar hbar "--hours-- total hours in a week" /40/; 

positive variable X(i) "--units-- production of units"; 
variable profit ; 

equation
eq_objfn "target of our optimization", 
eq_hourlimit "cant work more than 40 hours per week"
;

eq_objfn.. profit =e= sum(i,(r(i)-c(i)) * X(i)) ; 
* =e= "equal to"
* =g= "greater than"
* =l= "less than"

eq_hourlimit.. hbar =g= sum(i,h(i) * X(i)) ; 

scalar sw_combo /%combo%/ ; 

equation eq_combo;
eq_combo$sw_combo.. X("frenchfries") =g= X("hotdogs") ;

model ripley /all/ ; 

solve ripley using lp maximizing profit ;

execute_unload 'ripleysdata.gdx' ; 

$exit


parameter rep ; 
rep(i) = X.l(i) ; 

