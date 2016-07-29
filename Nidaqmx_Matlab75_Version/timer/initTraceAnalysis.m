function initTraceAnalysis
	global state gh
	gh.traceAnalyzer=guihandles(traceAnalyzer);
	openini('traceAnalysis.ini');
	
	state.analysis.setup={1 0 1 1 0 100 [1 0 100] {}};
	state.analysis.usedPC={};
	state.analysis.usedIC={};
	
	