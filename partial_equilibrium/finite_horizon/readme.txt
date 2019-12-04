please run main_partial_eq_finite_horizon from line 1 to line 138

and then change line 31 to
[vIncomeShocks, mTransition]=rouwenhorstFunction(ddelta,ssigmaErrorHigh,nGridShocks); % TRY DIFFERENT VALUES of ssigmaError, low and high

and then change line 138 to
save('finite_horizon_high_shock')

then run the whole file


