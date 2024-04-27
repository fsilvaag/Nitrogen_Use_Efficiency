Proc datasets lib=work kill; quit;

proc import OUT= Nitrogene
DATAFILE = '/home/estuvar3/UEN/Row_data_NUE.xlsx' 
DBMS=xlsx REPLACE;
GETNAMES=YES;
run;



proc sort data=Nitrogene;
by planted level;
run;


ods output covparms =  covparms_UENtotal;
ods output  solutionr = blup_NUE;
ods output  solutionr = blues_NUE;


proc mixed data = Nitrogene covtest covtest lognote itdetails;
by   planted level;
class  Name Block ;
model NUE =  block/s;
random Name  /s;
run;
quit;

ods output testsfornormality=norma1_RESIDS1_NUE;
proc univariate data= RESIDS1_PS_aereo normaltest;
var StudentResid;
QQPLOT StudentResid /NORMAL(MU=EST SIGMA=EST COLOR=RED L=1);
*PPPLOT Resid StudentResid/NORMAL(MU=EST SIGMA=EST COLOR=RED L=1);
HISTOGRAM /NORMAL(COLOR=MAROON W=4) CFILL = BLUE CFRAME = LIGR;
INSET MEAN STD /CFILL=BLANK FORMAT=5.2 ;
run;