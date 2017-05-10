#include <stdio.h> 
#include <math.h>
#include "mt64.h"

#define MAXDEG  6   /* Maximum degree of network, for memory alloc */
#define MAXNODES 100000 /* Maximum number of nodes, for memory alloc */
#define MYWORDSIZE 63 /*size of word for storage purposes */
#define MAXWORDS (MAXNODES/MYWORDSIZE)+1  /* Maximum number of words, for mem alloc */
#define SEED 5113533914238270004ULL /* Random seed value */
/* Note -- very important to include ULL suffix for 64bit unsigned long longs. omission denotes a 32bit int */
#define HITCOUNTPTS 4 /* number of attractors whose fraction of phase space should be printed */
#define NR_END 1 /* Some weird thing from Numerical Recipes */
#define FREE_ARG char* /* ditto */

/* This structure contains all the data for a single node */
typedef struct {
  int state;  /* Binary -- is this node on or off */
  int inputid[MAXDEG];  /* id numbers of the input nodes , each from 0 to N-1 */
  int inputstate;  /* aggregate state of the input nodes, from 0 to 2^K-1 */
  unsigned long long strategy;  /* id number of the strategy function that this node uses, from 0 to 2^(2^k)-1 */
  } Node;

/* This structure contains all the data for a single attractor */
typedef struct {
  int length; /* length of attractor */
  unsigned long long maxstate[MAXWORDS];  /* The state in the loop with the largest numerical value, starting with [1] */
  unsigned long long frozen[MAXWORDS];  /* Contains a 1 for each bit frozen in this attractor, 0 otherwise */
  int Nfrozen; /* Number of frozen nodes in this attractor */
  int hits; /* number of times this attractor has been found */
  } Attractor;

/* These are the pointers to the output files */
/* Only global variables in the entire program */
FILE *outfile,*detailfile;

/*----------------------------------Main Loop--------------------------*/
void main(int argc,char *argv[])  /* command line input unused */
{ 

void initsystem(Node agent[],int *,int *,int *,float *,float,float,double,float *);
void sharestateinfo();
void updatenodes();
void printstate();
void storestateinfo(unsigned long long **,Node agent[],int ,int ,int , int);
int lookforattractor(unsigned long long **,int , int, int);
int countfrozen(Node agent[],int ,int ,int );
void newinitcond(Node agent[],int );
int classifyattractor(Attractor attractorlib[],Node agent[],int *, int, int, int, int, int,int);
void readparamfile(int *, int *, int *, int *, int *, int *,float *,float *,int *,int *,int *);
unsigned long long **llmatrix(int, int, int, int);
Node *Nvector(int,int);
Attractor *Avector(int,int);
void nrerror();
int frozenheld(Attractor attractorlib[],int, int,int);
void systemstats(int, Attractor attractorlib[], int , int , int , int , int , int ,float, float,int);

Node *agent;  /* holds current state of system */
Attractor *attractorlib; /* stores all of the attractors found */
unsigned long long **state;  /* saves state at given interval in a set of words */
int Nagents, deg; /* stores N and (nominal) K */
int time,i,j,l,n; /* index variables */
int length;  /* Attractor length */
int frozennodes; /* total number of nodes in the attractor's frozen core */
int randomizations; /* total number of random initial configurations to try */
int systems; /* total number of initial system realizations, i.e. strats & inputs, "the wiring", to try */
int printflag;  /* Indicates whether initialization and state data will printf */
int diagflag; /* Indicates whether to print diagnostics to screen (i.e. every strategy & state in long form) & recount the number of frozen elements */
int Nattractors; /* Current number of attractors in library */
int attractorid; /* Id number of located attractor */
int numheldfrozen; /* number of nodes which are always frozen for a given system */
int words,overflow;  /* number of words needed for storage of a state and number of significant bits in last word */
float meanhom,meank; /* mean internal homogeneity of strategies, mean effective connectivity */
int maxsteps,maxsaves,saveint; /* Maximum number of iterations to run each initialization, maximum number of saves, and the interval */
int attrcapacity; /* maximum number of attractors, read from input file */
float lowp,highp;  /* even samples over this interval of desired meanhom, read from input file */
int chaoscount; /* keeps track of number of initial configurations which do not find an attractor */
double randiniter; /* drawn once for each system realization, determines the target mean homogeneity */

/* Here we read from a parameters file */

readparamfile(&deg,&Nagents,&randomizations,&diagflag,&printflag,&systems,&lowp,&highp,&maxsteps,&saveint,&attrcapacity);

/* Determine the number of state words, maximum number of state saves */
words=(int)((Nagents-1)/MYWORDSIZE)+1;
if(!(Nagents%MYWORDSIZE)) overflow=MYWORDSIZE;
else overflow=Nagents%MYWORDSIZE;
maxsaves=(int)(maxsteps/saveint)+1;

/* Allocate memory for the current state vector (somewhat deceptively named agent), the attractor library, and the saved states */
agent=Nvector(0,Nagents);
attractorlib=Avector(0,attrcapacity);
state=llmatrix(0,maxsaves,0,words);

init_genrand64(SEED);  /* Initialize the Mersenne twister */

/* Open output files. Detail file only prints if printflag is on -- will contain lines approx = realizations*configurations */
if((outfile = fopen("RBNout","w"))==NULL)
    nrerror("Error opening RBNout.");
fprintf(outfile,"Nagents=%d deg=%d randomizations=%d maxsteps=%d\n",Nagents,deg,randomizations,maxsteps);
fprintf(outfile,"meank meanhom Nattractors numalwaysfrozen avgattrlength stdevlength fracchaos fracpop-HITCOUNTPTS\n");
if(printflag) if((detailfile = fopen("RBNdetail","w"))==NULL)
    nrerror("Error opening RBNdetail.");
if(diagflag) printf("Random seed = %llu\n",SEED);

for(j=0;j<systems;j++){  /* This loop runs over system realizations */

  randiniter=genrand64_real2();  /* this random number draw will control the p of this system realization */
  initsystem(agent,&Nagents,&deg,&diagflag,&meanhom,lowp,highp,randiniter,&meank);  /* this subroutine establishes "the wiring" */
  if(printflag) fprintf(detailfile,"Initialized random system with N=%d meank %f and meanhom %f\n",Nagents,meank,meanhom);

  Nattractors=0;  /* start these tallies at zero */
  chaoscount=0;

  for(i=0;i<randomizations;i++){  /* This loop runs over random initial configurations */
	newinitcond(agent,Nagents);  /* This subroutine does not touch the wiring but chooses a new initial configuration */
	if(diagflag) printf("Random initial condition %d\n",i);
	for(time=0;time<maxsteps;time++){
 		if(diagflag) printstate(agent,Nagents,time);
		storestateinfo(state,agent,Nagents,time,saveint,words);
		if(time) length=lookforattractor(state,time,words,saveint);  /* returns zero if no attractor is found */
		else length=0;
		if(length){  /* If an attractor has been found, see if it's in the attractor library, if not add it */
		   attractorid=classifyattractor(attractorlib,agent,&Nattractors,Nagents,length,deg,words,overflow,attrcapacity);
      		   if(printflag) fprintf(detailfile,"Attr %d iterations %d length %d Nfrozen %d\n",attractorid,time,length,attractorlib[attractorid].Nfrozen);
		   break;
		}
		sharestateinfo(agent,Nagents,deg);  /* If no attractor is found, distribute current state information */
		updatenodes(agent,Nagents);  /* And update */
	}

	if(length&&diagflag){  /* This is a redundant calculation, and is only used if the diagnostic flag is on */
	   frozennodes=countfrozen(agent,Nagents,length,deg);
	   printf("Frozen core of %d nodes out of %d\n",frozennodes,Nagents);
 	}

        if(!length){  /* If no attractor is found, tick up the chaoscount */
	   if(printflag) fprintf(detailfile,"None %d iterations\n",time);
	   chaoscount+=1;
	}
  }

	/* This routine works out the relevant statistics for the system realization and prints to the outfile */
	systemstats(Nattractors,attractorlib,words,overflow,randomizations,chaoscount,Nagents,deg,meank,meanhom,printflag);
	fflush(outfile);
	fflush(detailfile);
}

	fclose(outfile);
	fclose(detailfile);
}

/*---------------------------nrerror-----------------------------------------*/
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
printf("Numerical Recipes run-time error...\n");
printf("%s\n",error_text);
printf("...now exiting to system...\n");
exit(1);
}

/*---------------------------------------------------------------------------*/
 
/*------------------------------------MT---------------------------------------*/
#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */


/* The array for the state vector */
static unsigned long long mt[NN];
/* mti==NN+1 means mt[NN] is not initialized */
static int mti=NN+1;

/* initializes mt[NN] with a seed */
void init_genrand64(unsigned long long seed)
{
    mt[0] = seed;
    for (mti=1; mti<NN; mti++)
        mt[mti] =  (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
}

/* generates a random number on [0, 2^64-1]-interval */
unsigned long long genrand64_int64(void)
{
    int i;
    unsigned long long x;
    static unsigned long long mag01[2]={0ULL, MATRIX_A};

    if (mti >= NN) { /* generate NN words at one time */

        /* if init_genrand64() has not been called, */
        /* a default initial seed is used     */
        if (mti == NN+1)
            init_genrand64(5489ULL);

        for (i=0;i<NN-MM;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        for (;i<NN-1;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        x = (mt[NN-1]&UM)|(mt[0]&LM);
        mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];

        mti = 0;
    }

    x = mt[mti++];

    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);

    return x;
}

/* generates a random number on [0,1)-real-interval */
double genrand64_real2(void)
{
    return (genrand64_int64() >> 11) * (1.0/9007199254740992.0);
}

/*------------------------------------------------------------------------*/

/*---------------------------initsystem-----------------------------------*/
/* passed the set of all nodes, initializes it with input ids, */
/* strategy, and state. */
void initsystem(Node agent[],int *Nagents,int *Kdegree,int *diagflag,float *meanhom,float lowp,float highp,double randiniter,float *meank)
{	
   int determinek(int,unsigned long long);

   int i,j,l;
   int N,K;
   unsigned long long potstrat;
   int potinput;
   int digits,necshift; /* how many binary digits is the strategy */
   float chooser;
   int digsone;
   float oneprob;
   int digszero;
   float homcomponent;
   int effk;
   int tryagain;

   N=*Nagents;
   K=*Kdegree;
   *meanhom=0.0;
   *meank=0.0;
   digits=1 << K;
   necshift=64-digits;

 /* first decide what the nodes inputs will be */

for(i=0;i<N;i++){
   for(j=0;j<K;j++){
	do{
	  tryagain=0;
	  potinput=(int)(N*genrand64_real2());
	  if(!((potinput>=0)&&(potinput<N))) tryagain=1;  /* check to be sure it's on [0,N) */
	  for(l=0;l<j;l++) if(potinput==agent[i].inputid[l]) tryagain=1; /* check to be sure we don't duplicate inputs -creates O(1/N) effects cf. Socolar */
        }while(tryagain); 
	agent[i].inputid[j]=potinput;
	if(*diagflag) printf("agent %d, input %d = %d ; ",i,j,potinput);
   }

 /* then decide what the strategy will be*/

   homcomponent=0.0;  /* this will count the p as we go */
   chooser=randiniter*(highp-lowp)+lowp;  /* the "chooser" is the p we are shooting for in this realization */
   if(genrand64_real2()<0.5) chooser-=(2.0*(chooser-0.5));  /* half the time we want a majority of zeroes */

   agent[i].strategy=0ULL;  /* Start with all zeros, then flip on accordingly */
   for(j=0;j<digits;j++){
        agent[i].strategy=agent[i].strategy << 1; /* Scoot all bits over to the left one */
   	if(genrand64_real2() < chooser){
        	agent[i].strategy+=1ULL;
        	homcomponent+=1.0;
        }
   }

   homcomponent/=(float)(digits);
   if(homcomponent<0.5) homcomponent+=2.0*(0.5-homcomponent);  /* If zeroes predominate, put in range [0.5,1.0] */
   *meanhom+=homcomponent;
   effk=determinek(K,agent[i].strategy);  /* This will determine the effective K of this node based on the strategy */
   *meank+=(float)effk;
   if(*diagflag) printf("strategy %d = %llu , effk = %d\n",i,agent[i].strategy,effk);

}

   *meanhom/=(float)N;
   *meank/=(float)N;
   if((*meanhom>1.0)|(*meanhom<0.5)) nrerror("meanhom calc messed up");

}

/*----------------------------------------------------------------------*/

/*-------------------------determinek-----------------------------------*/
/* This routine looks at the symmetries of the 2^(2^K) bit strategy to determine */
/* the effective K of the strategy */
int determinek(int maxk, unsigned long long strategy)
{
   int effectivek;
   int j,l,m,n;
   int slices,comparisons,digsperslice;
   int issame;
   unsigned long long usestrat,compare[1<<MAXDEG];

effectivek=maxk;  /* Start out assuming the effective K is the nominal K */
for(j=1;j<=maxk;j++){
	slices=1<<j;
	comparisons=1<<(j-1);
	digsperslice=1<<(maxk-j);
	issame=1;
	usestrat=strategy;
	for(l=0;l<slices;l++){
		compare[l]=0ULL;
		for(n=0;n<digsperslice;n++){
			compare[l]=compare[l] << 1;
			if(usestrat%2) compare[l]+=1ULL;
                	usestrat=usestrat >> 1;
		}
	}
	for(m=0;m<comparisons;m++){
		issame=(issame&&(compare[2*m]==compare[2*m+1]));
		if(!issame) break;
	}
	if(issame) effectivek--;  /* If all of the required comparisons reveal symmetries, decrement k */
}

return effectivek;

}
/*----------------------------------------------------------------------*/

/*-------------------------sharestateinfo-------------------------------*/
/* passed the set of all nodes, updates correct inputstate */
void sharestateinfo(Node agent[], int Nagents, int Kdegree)
{
   int i,j;
   int incomingnode;

for(i=0;i<Nagents;i++){
	agent[i].inputstate=0;
   for(j=0;j<Kdegree;j++){
	agent[i].inputstate=agent[i].inputstate << 1;
	incomingnode=agent[i].inputid[j];
        agent[i].inputstate+=agent[incomingnode].state;
        }
   }

}

/*----------------------------------------------------------------------*/

/*---------------------------updatenodes---------------------------------*/
/* passed the set of all nodes, updates state of node */
void updatenodes(Node agent[], int Nagents)
{
   int i;
   unsigned long long tag;
   unsigned long long result;

for(i=0;i<Nagents;i++){
   tag=1ULL << agent[i].inputstate;
   result=tag & agent[i].strategy;
   agent[i].state=(!(!result));
   }

}

/*----------------------------------------------------------------------*/

/*---------------------------printstate---------------------------------*/
/* prints all state data to screen */
void printstate(Node agent[], int Nagents, int time)
{
   int i;

   for(i=0;i<Nagents;i++)
   printf("%d",agent[i].state);
   printf(" at time=%d\n",time);
   
}

/*-----------------------------------------------------------------------*/

/*-------------------------storestateinfo-----------------------------*/
void storestateinfo(unsigned long long **state,Node agent[],int Nagents, int time, int saveint, int words)
{
   int i,j,word;
   int saveid;

   saveid=(int)(time/saveint);

   for(j=0;j<words;j++) state[saveid][j]=0ULL;  /* Initialize everything to zero */

   for(i=0;i<Nagents;i++){
        word=(int)(i/MYWORDSIZE); /* determine the word to which this agent belongs */
        state[saveid][word]=state[saveid][word] << 1; /* Scoot all the bits over to the left one */
        state[saveid][word]+=(unsigned long long)agent[i].state;  /* Then flip on the appropriate nodes */
   }

}

/*-----------------------------------------------------------------------*/

/*------------------------lookforattractor------------------------------*/
int lookforattractor(unsigned long long **state, int time,int words,int saveint)
{
   int i,j,lastsaved;
   int length=0;

   lastsaved=(int)(time/saveint);

   for(j=0;j<lastsaved;j++){
        length=0;
        for(i=0;i<words;i++){
		/* Compare the current state to all previously saved ones, working backwards */
                if(!(state[lastsaved-j-1][i]==state[lastsaved][i])) break;  /* If any word is unequal, break */
                else if(!(i==words-1)) continue;  /* If they're equal but we're not at the last word, continue */
		else length=j*saveint+time%saveint+1;  /* If they've all been equal, we've found an attractor */
        }
        if(length) break;
   }

   return length;

}

/*----------------------------------------------------------------------*/

/*--------------------------countfrozen--------------------------------*/
/*passed the current system and attractor length, counts the frozen nodes */
int countfrozen(Node agent[],int Nagents, int length, int deg)
{
   int last[Nagents];
   int frozen[Nagents];
   int i,j,result; 

   for(i=0;i<Nagents;i++){
	frozen[i]=1;
	last[i]=agent[i].state;
   }

   for(i=0;i<length;i++){
        sharestateinfo(agent,Nagents,deg);
        updatenodes(agent,Nagents);
	for (j=0;j<Nagents;j++){
	if(!(last[j]==agent[j].state)) frozen[j]=0;
	last[j]=agent[j].state;
	}
   }	

   result=0;
   for(i=0;i<Nagents;i++) if(frozen[i]) result+=1;

   return result;

}

/*---------------------------newinitcond-------------------------------*/
/*leaves the strategy the same but picks new random init conds */
void newinitcond(Node agent[], int Nagents)
{
   int i;

   for(i=0;i<Nagents;i++){
     if(genrand64_real2()<0.5)
     agent[i].state=1;
     else agent[i].state=0;
   }

}

/*----------------------------------------------------------------------*/

/*-------------------------readparamfile-------------------------------*/
void readparamfile(int *Kdegree, int *Nagents, int *randomizations, int *diagflag, int *printflag, int *systems, float *lowp,float *highp,int *maxsteps,int *saveint, int *attrcapacity)
{
    FILE *parameters;
    char trash[20];
    float lowpin,highpin;
    int Kin,Nin,Randinitsin,dflagin,pflagin,systemsin,maxstepsin,saveintin,attrcapin;

    /* open parameters file*/
    if((parameters=fopen("RBNparam","r"))==NULL)
      nrerror("Error opening parameter file");

    fscanf(parameters,"%s %d",trash,&Kin);
    if(Kin>MAXDEG) nrerror("degree exceeds maximum");
    *Kdegree=Kin;

    fscanf(parameters,"%s %d",trash,&Nin);
    if(Nin>MAXNODES) nrerror("number of nodes exceeds maximum");
    *Nagents=Nin;

    fscanf(parameters,"%s %d",trash,&Randinitsin);
    *randomizations=Randinitsin;

    fscanf(parameters,"%s %d",trash,&dflagin);
    *diagflag=dflagin;

    fscanf(parameters,"%s %d",trash,&pflagin);
    *printflag=pflagin;

    fscanf(parameters,"%s %d",trash,&systemsin);
    *systems=systemsin;

    fscanf(parameters,"%s %f",trash,&lowpin);
    *lowp=lowpin;

    fscanf(parameters,"%s %f",trash,&highpin);
    *highp=highpin;

    fscanf(parameters,"%s %d",trash,&maxstepsin);
    *maxsteps=maxstepsin;

    fscanf(parameters,"%s %d",trash,&saveintin);
    *saveint=saveintin;

    fscanf(parameters,"%s %d",trash,&attrcapin);
    /* if(Randinitsin>attrcapin) printf("Warning: attractor capacity smaller than randomizations\n"); */
    *attrcapacity=attrcapin;

    fclose(parameters);

}

/*----------------------------------------------------------------------*/

/*---------------------------classifyattractor--------------------------*/
/*passed the current system and attractor length, classifies the attractor */
int classifyattractor(Attractor attractorlib[],Node agent[],int *Nattractors, int Nagents, int length, int deg, int words, int overflow,int attrcapacity)
{
   unsigned long long **llmatrix( int, int, int, int); 
   void free_llmatrix(unsigned long long **, int, int, int, int);

   int maxstateid=0;
   unsigned long long currentfrozen[MAXWORDS];
   unsigned long long placeholder,freeze;
   unsigned long long **atstate;

   int i,j;
   int undecided;
   int notthisattractor;
   int attractorid;
   int numfrozen=0;

   atstate=llmatrix(0,length,0,words);  /* allocate memory */

   /* run thru the loop and find the maximum (defining) state and frozen core of the attractor */
   /* it's complicated only because we have multi-word states */

   for(i=0;i<length;i++){
	undecided=1;
	storestateinfo(atstate,agent,Nagents,i,1,words);
	if(i){ for(j=0;j<words;j++){
		freeze=~(atstate[i-1][j] ^ atstate[i][j]);  /* 1 for equivalence of the bit, 0 for difference */
		placeholder=currentfrozen[j];
		currentfrozen[j]=placeholder & freeze;
		if(undecided){  /* Here's where we decide whether this is the numerically maximum state */
		if(atstate[i][j]>atstate[maxstateid][j]){ 
	    	maxstateid=i;
	    	undecided=0;
	    	}
		else if(atstate[i][j]==atstate[maxstateid][j]) continue;
		else undecided=0;
	      } } }
	else for(j=0;j<words;j++) currentfrozen[j]=~0;  /* In the first go around we set all currentfrozen bits to one */
        if(!(length==i)){
		sharestateinfo(agent,Nagents,deg);
        	updatenodes(agent,Nagents);	
	}
   }

   /* does that state exist in the attractorlib? */
   if(*Nattractors){
   for(i=0;i<*Nattractors;i++){
	notthisattractor=0;
	for(j=0;j<words;j++){
		if(!(attractorlib[i].maxstate[j]==atstate[maxstateid][j])){
			notthisattractor=1;
			break;
			}
		}
	if(!notthisattractor){
	    attractorid=i;
	    break;
	    }
	else if((*Nattractors-1)==i) attractorid=i+1;
   } }

   /* if it doesn't, then add it */

   /* if(attractorid>=attrcapacity) nrerror("Attractor lib capacity exceeded"); */

   if((attractorid==*Nattractors)|(*Nattractors==0)){
	if(0==*Nattractors) attractorid=0;
	attractorlib[attractorid].length=length;
	if(words>1) {
	for(j=0;j<(words-1);j++){
	   attractorlib[attractorid].maxstate[j]=atstate[maxstateid][j];
	   attractorlib[attractorid].frozen[j]=currentfrozen[j];
	   for(i=0;i<MYWORDSIZE;i++)
		if((currentfrozen[j]>>i)%2ULL) numfrozen+=1;
	} }
	j=words-1;
	attractorlib[attractorid].maxstate[j]=atstate[maxstateid][j];
        attractorlib[attractorid].frozen[j]=currentfrozen[j];
	for(i=0;i<overflow;i++)
		if((currentfrozen[j]>>i)%2ULL) numfrozen+=1;
	attractorlib[attractorid].Nfrozen=numfrozen;
	attractorlib[attractorid].hits=0;
	*Nattractors+=1;
   }

   attractorlib[attractorid].hits+=1;  /* Increment the number of times this attractor was found by one */

   free_llmatrix(atstate,0,length,0,words);

   return attractorid;
}

/*--------------------------------------------------------------------*/

/*---------------------------llmatrix---------------------------------------*/
unsigned long long **llmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a ULL matrix with subscript range m[nrl..nrh][ncl..nch] */
{
int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
unsigned long long **m;
/* allocate pointers to rows */
m=(unsigned long long **) malloc((size_t)((nrow+NR_END)*sizeof(unsigned long long*)));
if (!m) nrerror("allocation failure 1 in matrix()");
m += NR_END;
m -= nrl;
/* allocate rows and set pointers to them */
m[nrl]=(unsigned long long *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(unsigned long long)));
if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
m[nrl] += NR_END;
m[nrl] -= ncl;
for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
/* return pointer to array of pointers to rows */
return m;
}

/*---------------------------------------------------------------------------*/

/*-----------------------------freellmatrix----------------------------------*/
void free_llmatrix(unsigned long long **m, int nrl, int nrh, int ncl, int nch)
/* free a ULL matrix allocated by llmatrix() */
{
free((FREE_ARG) (m[nrl]+ncl-NR_END));
free((FREE_ARG) (m+nrl-NR_END));
}

/*----------------------------------------------------------------------------*/

/*-----------------------------Nvector---------------------------------------*/
Node *Nvector(int nl, int nh)
/* allocate a Node vector with subscript range v[nl..nh] */
{
Node *v;
v=(Node *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(Node)));
if (!v) nrerror("allocation failure in Nvector()");
return v-nl+NR_END;
}
/*---------------------------------------------------------------------------*/

/*-----------------------------Avector---------------------------------------*/
Attractor *Avector(int nl, int nh)
/* allocate a Node vector with subscript range v[nl..nh] */
{
Attractor *v;
v=(Attractor *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(Attractor)));
if (!v) nrerror("allocation failure in Avector()");
return v-nl+NR_END;
}
/*---------------------------------------------------------------------------*/

/*------------------------------frozenheld-----------------------------------*/
int frozenheld(Attractor attractorlib[], int Nattractors, int words, int overflow)
{
	int i,n,l;
	unsigned long long heldfrozen[MAXWORDS]; /* holds common frozen component */
	int numheldfrozen=0;
	
        for(i=0;i<words;i++) {
                heldfrozen[i]=~(0ULL);
                for(n=0;n<Nattractors;n++) heldfrozen[i]=heldfrozen[i] & attractorlib[n].frozen[i];
                if(!(i==words-1)){ for(l=0;l<MYWORDSIZE;l++){
                if((heldfrozen[i]>>l)%2) numheldfrozen+=1;}
                }
                if(i==words-1){ for(l=0;l<overflow;l++){
                if((heldfrozen[i]>>l)%2) numheldfrozen+=1;}
                }
        }

   return numheldfrozen;
}

/*----------------------------------------------------------------------------*/

/*------------------------------systemstats-----------------------------------*/
void systemstats(int Nattractors, Attractor attractorlib[], int words, int overflow, int randomizations, int chaoscount, int Nagents, int deg,float meank,float meanhom,int printflag)
{

    float avglength,varlength,stdevlength;
    float frachits[HITCOUNTPTS];
    int i,l,n;
    int numheldfrozen;
    float currfracpop;
    float fracchaos;
    
    avglength=0.0;
    varlength=0.0;
    for(i=0;i<HITCOUNTPTS;i++) frachits[i]=0.0;
    if(Nattractors){
             numheldfrozen=frozenheld(attractorlib,Nattractors,words,overflow);  /* How big is the common frozen component */
             for(i=0;i<Nattractors;i++){   	/* Run through the attractor lib and find the average length and dist of phase space real estate */
                    avglength+=(float)attractorlib[i].length/(float)Nattractors;
                    currfracpop=(float)attractorlib[i].hits/(float)randomizations;
                    for(l=0;l<HITCOUNTPTS;l++){
                             if(currfracpop>=frachits[l]){
                                    for(n=(HITCOUNTPTS-1);n>l;n--) frachits[n]=frachits[n-1];
                                    frachits[l]=currfracpop;
                                    break;
                             }
                    }
             }
    }
    else numheldfrozen=0;
    if(Nattractors>1) for(i=0;i<Nattractors;i++) varlength+=((float)attractorlib[i].length-avglength)*((float)attractorlib[i].length-avglength)/(float)(Nattractors-1);
    stdevlength=sqrt(varlength);
    fracchaos=(float)chaoscount/randomizations;
    if(printflag) fprintf(detailfile,"Number of nodes which are always frozen for this system %d\n",numheldfrozen);
    fprintf(outfile,"%f %f %d %d %f %f %f",meank,meanhom,Nattractors,numheldfrozen,avglength,stdevlength,fracchaos);
    for(l=0;l<HITCOUNTPTS;l++) fprintf(outfile," %f",frachits[l]);
    fprintf(outfile,"\n");

}

/*----------------------------------------------------------------------------*/
