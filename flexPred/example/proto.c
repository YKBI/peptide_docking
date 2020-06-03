#include         <stdio.h>
#include         <math.h>
#include         <string.h>
#include         <stdlib.h>

/*----------------------------------------------------------*/
int     main(argc, argv)
int argc;
char *argv[];
{
	FILE    *infile,*out1,*out2,*out3,*fopen();
        char line[12000],fname[120],fname1[120],fname2[120];
        int  pubid, i,k,l,j,n,m,count;

        if(argc==2) {
                sscanf(*++argv, "%s", fname); 
	}
        else {
                printf(" nothing read from you \n");
                return 0 ;
        }
	infile = fopen(fname,"r");
        if (!infile) {
                printf(" Failed to open : %s\n",fname);
                exit(0);
        }
/* 
ATOM   3236  O   VAL B  10      24.285  22.172   4.053  1.00 18.19           O
0123456789012345678901234567890123456789012345678901234567890123456789
*/
	n = 0; int tmp;
	char aa[1200], bb[1200], cc[1200];
       	while (1) {
                if (!fgets(line, 12000,infile)) break;
		line[61]='0';
		line[62]='0';
		line[64]='0';
                line[65]='0';
		printf("%s",line);

		n++;
	}
	
/*
	infile = fopen(fname,"r"); 
        n = 0;
        while(1){
                if(EOF==fscanf(infile,"%*d%*d%*d%d",&nl)) break;
                for(i=0;i<nl;i++) {
                        fscanf(infile,"%d",&k); 
                }
                n++; 
        }
        fclose(infile);
*/
	return 0;
}
