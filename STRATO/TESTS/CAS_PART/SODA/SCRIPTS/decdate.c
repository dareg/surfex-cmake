/***
****  @(#)decdate.c version 1.1 --- 87/01/30	
****  @(#)decdate.c version 1.2 --- modif 97/12/08 annee en 4 chiffres	
****							
****  But    : Calculer YYMMDD + decalage = YYMMDD
****  Appel  : decdate YYMMDD   +-decalage
****           defaut YYMMDD   : date du jour
****           defaut decalage : 0
****  Sortie : YYMMDD
***/
 
#include <stdio.h>
#include <time.h>
#include <ctype.h>
 
#define	BIS(y)		(!(y % 4) && (y % 100) || !(y % 4) && !(y % 1000))
#define NBJOUR(y)	(BIS(y) ? 366 : 365)
#define IN(x,a,b)	((x >= a) && (x <= b))
 
char	*ermes[]  = { "Usage: decdate YYMMDD +-decalage\n",
		      "DECDATE: mauvaise date\n",
		      "DECDATE: mauvais  decalage\n",       };
 
int	nbjour[2][12] = { {31,28,31,30,31,30,31,31,30,31,30,31},
		 	  {31,29,31,30,31,30,31,31,30,31,30,31} };
 
main(argc,argv)
int	argc;
char	*argv[];
{
	int	decal;
	char	today[10];
	char	*date,*s;
	int	year,mon,day;
	int	numday;
	int	typ;
 
/* Initialisation des defauts			*/
	get_date(today);
	date  = today;
	decal = 0;
 
/* Controle des arguments			*/
 
        while (--argc){
  	   if (*argv[argc]=='+' || *argv[argc]=='-')
		decal = atoi(argv[argc]);
	   else date = argv[argc];
	}

	s = &date[6];  day  = atoi(s);			/* Jour		*/
	*s=0;s--;s--;  mon  = atoi(s);			/* Mois		*/
	*s=0;s--;s--;s--;s--;  year = atoi(s);		/* Annee	*/
	typ = BIS(year) ? 1 : 0;
	if (!IN(mon,1,12)) erreur(1);			/* 1..12   	*/
	if (!IN(day,1,nbjour[typ][mon-1])) erreur(1);	/* 28,29,30,31  */
 
/* Calcul de numday, decalage, puis conversion numday --> yymmdd	*/
	numday = quant(year,mon,day);
	numday += decal;
	while (numday < 1){
		year--;
		numday += NBJOUR(year);
	}
	while (numday > NBJOUR(year)){
		numday -= NBJOUR(year);
		year++;
	}
	yymmdd(numday,year,&mon,&day);
	printf("%04d%02d%02d\n",year,mon,day);
}
 
/*
**  quant --- retourne le numero du jour dans l'annee
*/
quant(an,mois,jour)
int	an,mois,jour;
{
	int	typ,i,n;
 
	n = 0;
	typ = BIS(an) ? 1 : 0;
	for(i=1;i<mois;i++) n += nbjour[typ][i-1];
	n += jour;
	return(n);
}
 
/*
** yymmdd --- transforme le quantieme n de l'annee n en un jour et un mois
*/
yymmdd(n,an,mois,jour)
int	n,an;					/* n=quantieme           */
int	*mois,*jour;				/* valeur mises a jour	 */
{
	int	typ;				/* type bissextile ou non */
 
	*mois = 0;
	*jour = 0;
	typ = BIS(an) ? 1 : 0;

	while (n > 0)
		n -= nbjour[typ][(*mois)++];	/* Quel mois ?		*/
	*jour = n+nbjour[typ][(*mois) - 1];	/* Quel jour ?		*/
}
 
/*
** nombre --- retourne 1 si la chaine s est composee de chiffres
*/
nombre(s)
char *s;
{
	while( *s != 0 ) if (!isdigit(*s++)) return(0);
	return(1);
}
 
 
/*
** get_date --- retourne la date du jour sous la forme AAMMJJ
*/
get_date(p)
char	*p;
{
	long	time();
	struct tm *gmtime();
	struct tm *tm;
	long	clock;

	clock = time((long *)0);
	tm = gmtime(&clock);
	sprintf(p,"%04d%02d%02d",tm->tm_year+1900,tm->tm_mon+1,tm->tm_mday);
}



/*
** erreur --- affiche le message d'erreur n et exit
*/
erreur(n)
int	n;
{
	fputs(ermes[n],stderr);
	exit(1);
}
