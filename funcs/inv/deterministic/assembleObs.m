%%% =======================================================================
%%% = assembleObs.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Read the observations and assemble them in a vector.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): in -- Structure with the observations.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- Vector with the observations.
%%% =======================================================================

function [ out ] = assembleObs( in )
global use_strat_N2O % fit for N2O strat obs?
out = [in.nh_ch4;in.sh_ch4;in.nh_ch4c13;in.sh_ch4c13;in.nh_mcf;in.sh_mcf;in.nh_n2o;in.sh_n2o;in.nh_c2h6;in.sh_c2h6;in.nh_co;in.sh_co];

if use_strat_N2O
  out = [out; in.nh_n2o_strat; in.sh_n2o_strat];
end

end


%%% =======================================================================
%%% = END
%%% =======================================================================
