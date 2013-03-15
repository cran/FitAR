void asBinary(int *n, int *m, int *k)
/* Input parameters:
* n, non-negative integer for which base 2 representation is needed
* m, integer vector of length k>0
* Output parameters
* m contains the base 2 representation as integer vector
*/
{
int i, j0, j=*n, kk=*k;
for (i=0; i<kk; i++ ) {
	j0 = j;
	j = j/2;
	m[i] = j0 -j*2;
}
return;
}
