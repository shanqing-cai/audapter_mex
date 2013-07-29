#include "pvocWarpAtom.h"

pvocWarpAtom::pvocWarpAtom(){
	tBegin = 0;
	rate1 = 0.5;	
	dur1 = 0.5; 
	durHold = 1;
	rate2 = 2;
		
	dur2 = (1 - rate1) / (rate2 - 1) * dur1;
}