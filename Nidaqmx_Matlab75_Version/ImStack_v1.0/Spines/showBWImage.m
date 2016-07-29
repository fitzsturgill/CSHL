function showBWImage
global gh state
try
	[XStart, XEnd, YStart, YEnd] = getCurrentAxisLimits(gh.spineGUI.initialaxis);
	if state.imageProc.spine.topImage
		if state.imageProc.spine.maxFlag == 0
			b=state.imageProc.spine.initialImage;
		else
			b=state.imageProc.spine.maxProjection;
		end		
	end
catch
	[XStart, XEnd, YStart, YEnd] = getCurrentAxisLimits(gh.spineGUI.initialaxis2);
	if state.imageProc.spine.bottomImage
		if state.imageProc.spine.maxFlag == 0
			b=state.imageProc.spine.initialImage2;
		else
			b=state.imageProc.spine.maxProjection2;
		end		
	end
	
end
	
image2bw(b(YStart:YEnd,XStart:XEnd,:),state.imageProc.spine.spineThreshold);