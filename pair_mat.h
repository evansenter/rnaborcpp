#include <ctype.h>
#define NBASES 8
/*@notnull@*/

#define MAXALPHA 20       /* maximal length of alphabet */

extern short alias[MAXALPHA+1];
extern int pair[MAXALPHA+1][MAXALPHA+1];
/* rtype[pair[i][j]]:=pair[j][i] */
extern int rtype[8];

/* for backward compatibility */
#define ENCODE(c) encode_char(c)

int encode_char(char c);

/*@+boolint +charint@*/
/*@null@*/
extern char *nonstandards;
extern void   nrerror(const char message[]);
void make_pair_matrix(void);
