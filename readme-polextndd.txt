Replication code for "Extending Policies to Control Groups or Past Periods in Nonlinear Difference in Differences" by Eugene Lee, Myoung-jae Lee, and Goeun Lee (Applied Economics, 2025, Advance online publication)

The code file called "polextndd-04aug2025.do" (for Stata) can replicate "Probit Estimates" and "Actual and Counterfactual Treatment Effects". It includes a program called "polextndd" to present the latter results. Syntax and example of the program are provided below.

syntax  : polextndd depvar s q d covariates [if] [in]
example :
. global y "y"
. global s "s"
. global q "q"
. global d "d"
. global w "age agesq medicost female white msa diabetes copd chf educoll employed income"

. polextndd ${y} ${s} ${q} ${d} ${w}

Actual, Treat Q = 1 at t = 1 (tv):  0.0141  ( 7.39 )
    Treat everybody at t = 0 (tv):  0.1611  ( 8.91 )
    Treat everybody at t = 1 (tv):  0.0312  ( 7.81 )
        Treat Q = 1 at t = 0 (tv):  0.0402  ( 8.21 )
        Treat Q = 0 at t = 0 (tv):  0.1208  ( 9.09 )
        Treat Q = 0 at t = 1 (tv):  0.0171  ( 7.85 )
