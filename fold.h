/* function from fold.c */
extern float  fold(const char *sequence, char *structure);
/* calculate mfe-structure of sequence */
extern float  energy_of_struct(const char *string, const char *structure);
/* calculate energy of string on structure */
extern void   free_arrays(void);           /* free arrays for mfe folding */
extern void   initialize_fold(int length); /* allocate arrays for folding */
extern void   update_fold_params(void);    /* recalculate parameters */
extern char  *backtrack_fold_from_pair(char *sequence, int i, int j);
extern float energy_of_circ_struct(const char *string, const char *structure);
int HairpinE(int size, int type, int si1, int sj1, const char *string);
int LoopEnergy(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1);

