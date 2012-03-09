/*functions for converting from my format to Vienna program.*/

void Initialize_Params();
double HP_Energy(int i, int j, short *S0, char* sequence);
double IL_Energy(int i, int j, int ip, int jp, short *S0);
double AU_Penalty(int i, int j, short *S0);
double MLbasepairAndAUpenalty(int i, int j, short *S0);
