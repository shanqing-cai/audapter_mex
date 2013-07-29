function ComputeF0_helper
%COMPUTEF0_HELPER  - wrapper for mexified YIN functions
%
%	usage:  varargout = ComputeF0_helper(selector, varargin)
%
%    where SELECTOR is one of
%    	CUMNORM_INPLACE
%    	DFTOPERIOD
%    	DFTOPERIOD2
%    	INTERP_INPLACE
%    	MININRANGE
%    	MINPARABOLIC
%    	RDIFF_INPLACE
%    	RSMOOTH
%    	RSUM_INPLACE
%    	
% These originally standalone functions supporting Alain de Cheveigne's
% YIN approach to pitch estimation are gathered together here to simplify
% compilation and porting.  Note that many of these modify their input
% arguments in place, which Mathworks deprecates...
    
% mkt 01/08
