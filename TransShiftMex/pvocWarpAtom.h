class pvocWarpAtom{
public:
	pvocWarpAtom();

	mytype tBegin;
	mytype rate1;	// Should aways be in the interval of (0, 1)
	mytype dur1;
	mytype durHold;
	mytype rate2;	// Should aways be > 1
	mytype dur2;

	
};