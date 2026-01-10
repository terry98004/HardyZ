// -------------------------------------------------------------------
// Program last modified November 30, 2025. 
// Copyright (c) 2024-2026 Terrence P. Murphy
// MIT License -- see hardyZ.c for details.
// -------------------------------------------------------------------

struct HZ {
	char	tBuf[HGT_MAX_CMDLINE_STRLEN + 2];		// holds entered '-t' value
	char	incrBuf[HGT_MAX_CMDLINE_STRLEN + 2];	// holds entered '-i' value
	int 	Count;    	 							// number of t values to check
	bool	Verbose;								// verbose report? true or false
	bool	ShowSeconds;							// report compute second taken T/F
	int		DebugFlags;								// used for debugging
	int		OutputDPt;								// digits after '.' in reported t value 
	int		OutputDPz;								// digits after '.' in reported HardyZ value 
	int		DefaultBits;							// Bits for MPFR floating point 
	int		Threads;								// Number of threads to use 
}; 

struct HZ_RPT {
	bool	Verbose;								// verbose report? true or false
	int		OutputDPt;								// digits after '.' in reported t value 
	int		OutputDPz;								// digits after '.' in reported HardyZ value 							
}; 

int ComputeAllHardyZ(struct HZ hz);
int HardyZCallback(mpfr_t t, mpfr_t HardyZ, int i, int CallerID);
