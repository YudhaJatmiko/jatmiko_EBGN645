
set i /roll,croissant,bread/ ; 

parameter 
r(i) "--$/item-- revenue per unit sold"
/
roll 2.25,
croissant 5.5,
bread 10
/ ;

parameter c(i) "--$/item-- cost per unit sold"
/
roll 1.5,
croissant 2,
bread 5
/ ;

parameter h(i) "--hours/item-- hours required per item sold"
/
roll 1.5,
croissant 2.25,
bread 5
/ ;

scalar hbar "--hours-- total hours in a week" /40/; 

positive variable X(i) "--units-- production of units"; 
variable profit ; 

equation
eq_objfn "target of our optimization", 
eq_hourlimit "cant work more than 40 hours per week"
;

eq_objfn.. profit =e= sum(i,(r(i)-c(i)) * X(i)) ;

eq_hourlimit.. sum(i, h(i) * X(i)) =l= hbar ; 
* =e= "equal to"
* =g= "greater than"
* =l= "less than"

scalar sw_combo /0/ ; 

equation eq_combo;
eq_combo$sw_combo.. X("croissant") =e= X("roll") ;

model benny /all/ ; 

sw_combo = 0 ;
solve benny using lp maximizing profit ;
display "COMBO 0 RESULTS:", X.l, profit.l ;

sw_combo = 1 ;
solve benny using lp maximizing profit ;
display "COMBO 1 RESULTS:", X.l, profit.l ;

execute_unload 'benny.gdx' ; 

parameter rep(i) ;
rep(i) = X.l(i) ;

$exit 

