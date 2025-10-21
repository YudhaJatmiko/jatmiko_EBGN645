set m /X1, X2/ ;
set b /yellow, blue, green, orange, purple/ ;

parameter r(b) "--$/beans-- revenue per beans"
/
yellow 1,
blue 1.05,
green 1.07,
orange 0.95,
purple 0.9
/ ;

scalar hours_per_week "--hours/week-- available hours per week" /40/ ;
scalar beans_per_hour "--beans/hour-- production rate" /100/ ;
scalar capacity "--beans/week-- capacity per machine" ;
capacity = hours_per_week * beans_per_hour ;

scalar sw_equal /0/ ;
scalar sw_restrict /0/ ;

set m_b(m,b) "valid machine-bean combinations" ;
m_b(m,b) = yes ;

positive variable Q(m,b) "--beans/week-- production quantity" ;
variable profit ;

equation
eq_objfn "objective function",
eq_capacity(m) "machine capacity constraint"
;

eq_objfn.. profit =e= sum((m,b)$m_b(m,b), r(b) * Q(m,b)) ;

eq_capacity(m).. sum(b$m_b(m,b), Q(m,b)) =l= capacity ;

alias(b, bb) ;

equation eq_upper(b,bb), eq_lower(b,bb) ;

eq_upper(b,bb)$(ord(b) lt ord(bb) and sw_equal)..
sum(m$m_b(m,b), Q(m,b)) =l= 1.06 * sum(m$m_b(m,bb), Q(m,bb)) ;

eq_lower(b,bb)$(ord(b) lt ord(bb) and sw_equal)..
sum(m$m_b(m,b), Q(m,b)) =g= 0.94 * sum(m$m_b(m,bb), Q(m,bb)) ;

model base /eq_objfn, eq_capacity/ ;
model equal /eq_objfn, eq_capacity, eq_upper, eq_lower/ ;

solve base using lp maximizing profit ;
display Q.l, profit.l ;

sw_equal = 1 ;
solve equal using lp maximizing profit ;
display Q.l, profit.l ;

sw_restrict = 1 ;
m_b(m,b) = no ;
m_b("X1","yellow") = yes ;
m_b("X1","blue") = yes ;
m_b("X1","green") = yes ;
m_b("X2","yellow") = yes ;
m_b("X2","orange") = yes ;
m_b("X2","purple") = yes ;

model restrict /eq_objfn, eq_capacity, eq_upper, eq_lower/ ;
solve restrict using lp maximizing profit ;
display Q.l, profit.l ;
