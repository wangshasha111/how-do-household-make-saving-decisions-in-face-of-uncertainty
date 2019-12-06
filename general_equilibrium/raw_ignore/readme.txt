please run main_partial_eq_infinite_horizon from line 1 to line 147

and then change line 31 to
[vIncomeShocks, mTransition]=rouwenhorstFunction(ddelta,ssigmaErrorHigh,nGridShocks); % TRY DIFFERENT VALUES of ssigmaError, low and high

and then change line 147 to
save('infinite_horizon_high_shock')

then run the whole file

Question: I don't understand why the lagrangian multiplier is not always positive
