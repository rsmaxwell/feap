               
/*---------------------------------------------------------------------------*/
/* Include files                                                             */
/*---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

/*---------------------------------------------------------------------------*/
/* Macros                                                                    */
/*---------------------------------------------------------------------------*/

#define BLANK           999.0
#define LINE_LENGTH     256
#define FALSE             0
#define TRUE              1

#define Get(array,row,col,limit)        array[(row) * (limit) + col]
#define Set(array,row,col,limit,value)  array[(row) * (limit) + col] = value

/*---------------------------------------------------------------------------*/
/* Declarations                                                              */
/*---------------------------------------------------------------------------*/

typedef struct control_tag CONTROL;

struct control_tag
{
    char   s1[8];
    char   s2[8];
    double d1;
    double d2;
};

/*---------------------------------------------------------------------------*/
/* blank common                                                              */
/*---------------------------------------------------------------------------*/

int M[2000] = {0};

/*---------------------------------------------------------------------------*/
/* PSIZE common                                                              */
/*---------------------------------------------------------------------------*/

int MAX = 2000;

/*---------------------------------------------------------------------------*/
/* CDATA common                                                              */
/*---------------------------------------------------------------------------*/

char  O[LINE_LENGTH + 1]    = "\n";
char  HEAD[LINE_LENGTH + 1] = "";
int   NUMNP  = 0;
int   NUMEL  = 0;
int   NUMMAT = 0;
int   NEN    = 0;
int   NEQ    = 0;
int   IPR    = 2;

/*---------------------------------------------------------------------------*/
/* ELDATA common                                                             */
/*---------------------------------------------------------------------------*/

int DM  = 0;
int N   = 0;
int MA  = 0;
int MCT = 0;
int IEL = 0;
int NEL = 0;

/*---------------------------------------------------------------------------*/
/* ENGYS common                                                              */
/*---------------------------------------------------------------------------*/

double AENGY = 0.0;

/*---------------------------------------------------------------------------*/
/* LABEL common                                                              */
/*---------------------------------------------------------------------------*/

char *BC    = " B.C. ";
char *DI[3] = {"Displ", "Veloc", "Accel"};
char *CD    = " Coordinates";
char *TE    = " Temperature";
char *FD    = " Force/Displ";

/*---------------------------------------------------------------------------*/
/* PRLOD common                                                              */
/*---------------------------------------------------------------------------*/

double PROP = 0.0;

/*---------------------------------------------------------------------------*/
/* TDATA common                                                              */
/*---------------------------------------------------------------------------*/

double TIME = 0.0;
double DT   = 0.0;
int    C1 = 0;
int    C2 = 0;
int    C3 = 0;
int    C4 = 0;
int    C5 = 0;

int esc       = 27;
int fs        = 28;
int form_feed = 12;
int condensed = 15;
int cpi_15    = 103;
int cpi_10    = 80;
int cpi_12    = 77;
int reset     = 18;

/*---------------------------------------------------------------------------*/
/*                                                                           */
/*---------------------------------------------------------------------------*/

char trace_prefix[200]  = "";
int  trace_indent       = 0;

/*---------------------------------------------------------------------------*/
/* Prototypes                                                                */
/*---------------------------------------------------------------------------*/

void   PCONTR      (void);
void   GENVEC      (int NDM, double *X, char *title, int PRT, int *ERR);
void   PMESH       (int *IDL, int *IE, double *D, int *ID, double *X, int *IX, double *F, double *T,
                    int NDF, int NDM, int NEN1, int III);
void   SETMEM      (int J);
void   PMACR       (double *UL, double *XL, double *TL, int *LD, double *P, double *S, int *IE, double *D,
                    int *ID, double *X, int *IX, double *F, double *T, int *JDIAG, double *B, double *DR,
                    CONTROL *CT, int NDF, int NDM, int NEN1, int NST, int NEND);
void   ACTCOL      (double *A, double *B, int *JDIAG, int NEQ, int AFAC, int BACK);
void   ADDSTF      (double *A, double *B, double *C, double *S, double *P, int *JDIAG, int *LD, int NST, int NFL, int AFL, int BFL, int CFL);
double DOT         (double *A, double *B, int N);
void   ELMLIB      (double *D, double *U, double *X, int *IX, double *T, double *S, double *P, int I, int J, int K, int ISW);
void   EUPDAT      (double *DR, double *U, double *V, double *A, double *XM, double DT, int NEQ);
void   NORM        (double *X, double *Y, int N);
int    PCOMP       (double A, double B);
void   PEIGS       (double *A, double *B, double *F, double *X, double *Y, double *Z, int *ID, int *IX,
                    int *JDIAG, int NDF, int NDM, int NEN1, int DFL);
void   PFORM       (double *UL, double *XL, double *TL, int *LD, double *P, double *S, int *IE,
                    double *D, int *ID, double *X, int *IX, double *F, double *T, int *JDIAG,
                    double *B, double *A, double *C, int NDF,
                    int NDM, int NEN1, int NST, int ISW, double *U, double *UD, int AFL, int BFL, int CFL, int DFL);
void   PLOAD       (int *ID, double *F, double *B, int NN, double P);
void   PROFIL      (int *JDIAG, int *ID, int *IX, int NDF, int NEN1, int *NAD);
void   PROMUL      (double *A, double *B, double *C, int *JDIAG, int NEQ);
double PROPLD      (double T, int J);
void   PRTDIS      (int *ID, double *X, double *B, double *F, int NDM, int NDF);
void   PRTREA      (double *R, int NDF);
void   PSETM       (int *NA, int *NE, int NJ, int *AFL);
void   PZERO       (double *V, int NN);
void   UACTCL      (double *A, double *C, double *B, int *JDIAG, int NEQ, int AFAC, int BACK);
void   ELMT01      (double *D, double *UL, double *XL, int *IX, double *TL, double *S, double *P, int NDF, int NDM, int NST, int ISW);
void   PGAUSS      (int l, int *LINT, double *R, double *Z, double *W);
void   PSTRES      (double *SIG, double *P1, double *P2, double *P3);
void   SHAPE       (double SS, double TT, double *X, double *SHP, double *XSJ, int NDM, int NEL, int *IX, int FLG);
void   SHAP2       (double S, double T, double *SHP, int *IX, int NEL);
void   ELMT02      (double *D, double *UL, double *XL, int *IX, double *TL, double *S, double *P, int NDF, int NDM, int NST, int ISW);
void   ELMT03      (double *D, double *UL, double *XL, int *IX, double *TL, double *S, double *P, int NDF, int NDM, int NST, int ISW);
int    read_int    (char *line, int lo, int hi);
double read_double (char *line, int lo, int hi);
char  *read_string (char *line, int lo, int hi, char *buffer);
char  *strip       (char *line);
int    ISIGN       (int a, int b);
double AMOD        (double a, double b);
void   dspMatrix   (int n, int *JDIAG, double *A, double *B, char *format, ...);
void   print2Darray(char *name,  void *pointer, char *format, int row, int col);
void   print1Darray(char *name,  void *pointer, char *format, int size);
void   print2Dvalue(int integer, void *pointer, char *format, int i, int j, int col);
void   print1Dvalue(int integer, void *pointer, char *format, int i);
void   traceIN     (char *format, ...);
void   traceOUT    (char *format, ...);

/******************************************************************************
*******************************************************************************
*******************************************************************************
*                                                                             *
* CONTROL AND DATA INPUT MODULES                                              *
*                                                                             *
*******************************************************************************
*******************************************************************************
******************************************************************************/
                                                                            
                                                                            
/*****************************************************************************/
/* Finite Element Analysis Program (FEAP) for solution of general            */
/* problem classes using the finite element method. Problem size             */
/* is controlled by the dimension of blank common and value of MAX           */
/* as set in main program.  ALL arrays must reside in central memory         */
/*                                                                           */
/* Programmed in Fortran IV by Prof. R.L. Taylor, Department of Civil        */
/* Engineering, University of California, Berkeley, California 94720, U.S.A  */
/*                                                                           */
/* Converted to C by Richard Maxwell, IBM Hursley Park UK                    */
/*****************************************************************************/

int main(void)
{
    CONTROL *CT = NULL;
    double *B = NULL, *D = NULL, *UL = NULL, *TL = NULL, *DR = NULL, *P = NULL, *S = NULL, *XL = NULL, *X = NULL, *F = NULL, *T = NULL;
    int     *ID = NULL, *IDL = NULL, *IE = NULL, *IX = NULL, NE = 0;
    int     NDM = 0, NDF = 0, NEN1 = 0, NST = 0;
    int    *JDIAG = NULL, more = TRUE;

    printf("%c",   reset);
    printf("%c%c", esc, cpi_10);
    printf("\n");

    traceIN("main");

    while (more)
    {
        /*-------------------------------------------------------------------*/
        /* Read a line and compare first 4 chars with macro list             */
        /*-------------------------------------------------------------------*/

        char line[LINE_LENGTH + 1]= "";
        char *TITL = gets(line);

        if (TITL == NULL)
        {
            printf("     **FATAL ERROR 50** Unexpected end-of-file\n");
            more = FALSE;
        }
        else
        {
            TITL = strip(TITL);
           
            if (strncmp(TITL, "FEAP", 4) == 0)
            {
                int N0 = 0, N1 = 0, N2 = 0, N3 = 0, N4 = 0, N5 = 0, N6 = 0, N7 = 0, N8 = 0, N9 = 0, N10 = 0, N11 = 0, N12 = 0, N13 = 0, N14 = 0;
                int NAD = 0, III = 0;
                char *data = NULL;
           
                /*-----------------------------------------------------------*/
                /* Read and print control information                        */
                /*-----------------------------------------------------------*/
           
                strcpy(HEAD, TITL);
                data = gets(line);

                if (data == NULL)
                {
                    printf("     **FATAL ERROR 51** Unexpected end-of-file\n");
                    more = FALSE;
                }
                else
                {
                    NUMNP  = read_int(line,  0,  5);
                    NUMEL  = read_int(line,  5, 10);
                    NUMMAT = read_int(line, 10, 15);
                    NDM    = read_int(line, 15, 20);
                    NDF    = read_int(line, 20, 25);
                    NEN    = read_int(line, 25, 30);
                    NAD    = read_int(line, 30, 35);
                 
                    printf("%s\n\n", HEAD);
                    printf("     Number of nodal points        = %d\n", NUMNP);
                    printf("     Number of elements            = %d\n", NUMEL);
                    printf("     Number of material sets       = %d\n", NUMMAT);
                    printf("     Dimension of coordinate space = %d\n", NDM);
                    printf("     Degree of freedom/node        = %d\n", NDF);
                    printf("     Nodes per element (maximum)   = %d\n", NEN);
                    printf("     Extra D.O.F. to element       = %d\n", NAD);
                 
                    /*-------------------------------------------------------*/
                    /* Set pointers for allocation of data arrays            */
                    /*-------------------------------------------------------*/
                 
                    NEN1 = NEN + 1;
                    NST = NEN*NDF + NAD;         UL    = (double*) (void*) &M[0];
                    N0  = 0  + NST*2*IPR;        XL    = (double*) (void*) &M[N0];
                    N1  = N0 + NEN*NDM*IPR;      TL    = (double*) (void*) &M[N1];
                    N2  = N1 + NEN*IPR;          IDL   = (int*)    (void*) &M[N2];
                    N3  = N2 + NST;              P     = (double*) (void*) &M[N3];
                    N4  = N3 + NST*IPR;          S     = (double*) (void*) &M[N4];
                    N5  = N4 + NST*NST*IPR;      IE    = (int*)    (void*) &M[N5];
                    N6  = N5 + NUMMAT;           D     = (double*) (void*) &M[N6];
                    N7  = N6 + 10*NUMMAT*IPR;    ID    = (int*)    (void*) &M[N7];
                    N8  = N7 + NDF*NUMNP;        X     = (double*) (void*) &M[N8];
                    N9  = N8 + NDM*NUMNP*IPR;    IX    = (int*)    (void*) &M[N9];
                    N10 = N9 + NEN1*NUMEL;       F     = (double*) (void*) &M[N10];
                    N11 = N10 + NDF*NUMNP*IPR;   T     = (double*) (void*) &M[N11];
                    N12 = N11 + NUMNP*IPR;       JDIAG = (int*)    (void*) &M[N12];
                    N13 = N12 + NDF*NUMNP;                                         
                 
                    /*-------------------------------------------------------*/
                    /* Check that sufficient memory exists                   */
                    /*-------------------------------------------------------*/
                 
                    SETMEM(N13);
                    PZERO(UL, N12);
                 
                    /*-------------------------------------------------------*/
                    /* Call mesh input subroutine to read and print all mesh */
                    /* data                                                  */
                    /*-------------------------------------------------------*/
                 
                    III = 0;
                    PMESH(IDL, IE, D, ID, X, IX, F, T, NDF, NDM, NEN1, III);
                 
                    /*-------------------------------------------------------*/
                    /* Establish profile of resulting equations for          */
                    /* stiffness, mass, etc                                  */
                    /*-------------------------------------------------------*/
                 
                    PROFIL(JDIAG, ID, IX, NDF, NEN1, &NAD);
                 
                    /*-------------------------------------------------------*/
                    /* Set pointers for solution arrays & check for          */
                    /* sufficient memory                                     */
                    /*-------------------------------------------------------*/
                 
                    N13 = N12 + NEQ;                     B     = (double*)  (void*) &M[N13];
                    N14 = N13 + NEQ*IPR;                 DR    = (double*)  (void*) &M[N14];
                    NE  = N14 + NUMNP*NDF*IPR;           CT    = (CONTROL*) (void*) &M[NE];
                 
                    SETMEM(NE);
                    PZERO(B, NEQ);
                }
            }
           
            else if (strncmp(TITL, "MACR", 4) == 0)
            {
                /*-----------------------------------------------------------*/
                /* Call macro solution module for establishing solution      */
                /* algorithm                                                 */
                /*-----------------------------------------------------------*/
           
                PMACR(UL, XL, TL, IDL, P, S, IE, D, ID, X, IX, F, T, JDIAG, B,
                      DR, CT, NDF, NDM, NEN1, NST, NE);
            }
           
            else if (strncmp(TITL, "STOP", 4) == 0)
            {
                more = FALSE;
            }
           
            else if (*TITL != 0)
            {
                printf("     **FATAL ERROR 52** Unexpected keyword '%s'\n", TITL);
                more = FALSE;
            }
        }
    }

    traceOUT("main");
    return 0;
}

/*****************************************************************************/
/*                                                                           */
/* Generate real data arrays by linear interpolation                         */
/*                                                                           */
/*****************************************************************************/

void GENVEC(int NDM, double *X, char *title, int PRT, int *ERR)
{
    int  more_input = TRUE, more = TRUE;
    double XL[7] = {0.0};
    int NN = 0;
    int NG = 0;

    traceIN("GENVEC");

    *ERR = FALSE;

    more_input = TRUE;
    while (more_input)
    {
        char line[LINE_LENGTH + 1] = "", *data = NULL;
        int  L = 0, LG = 0, LI = 0, i = 0, l = 0, n = 0;

        L  = NN;                         l = L - 1;
        LG = NG;

        data = gets(line);
        if (data == NULL)
        {
            printf("     **FATAL ERROR 55** Unexpected end-of-file\n");
            more_input = FALSE;
        }
        else
        {
            NN = read_int(line,  0,  5);     n = NN - 1;
            NG = read_int(line,  5, 10);
            for (i = 0; i < 7; i++)
                XL[i] = read_double(line,  (10+10*i), (20+10*i));
           
            if ((NN <= 0) || (NN > NUMNP))
                more_input = FALSE;
            else
            {
                for (i = 0; i < NDM; i++)
                    Set(X, n, i, NDM, XL[i]);
               
                if (LG != 0)
                {
                    LG = ISIGN(LG,NN-L);
                    LI =(abs(NN-L+LG)-1) / abs(LG);
               
                    for (i = 0; i < NDM; i++)
                        XL[i] = ( Get(X, n, i, NDM) - Get(X, l, i, NDM) ) / LI;
               
                    more = TRUE;
                    while (more)
                    {
                        L = L + LG;          l = L - 1;
                       
                        if ((NN-L)*LG > 0)
                        {
                            if ((L <= 0) || (L > NUMNP))
                            {
                                printf("     **FATAL ERROR 02** Attempt to generate node %d in %s\n", L, title);
                                *ERR  = TRUE;
                                more = FALSE;
                            }
                            else
                            {
                                for (i = 0; i < NDM; i++)
                                    Set(X, l, i, NDM, Get(X, l-LG, i, NDM) + XL[i]);
                            }
                        }
                        else
                            more = FALSE;
                    }
                }
            }
        }
    }

    if (PRT && !(*ERR))
    {
        int I = 0, J = 0, j = 0, L = 0, l = 0;

        for (I = 1; I <= NUMNP; I += 50)
        {
            printf(O);
            printf("%s\n\n     NODAL%s\n\n      NODE", HEAD, title);
            for (L = 1; L <= NDM; L++) printf("%7d%6.6s", L, title);
            printf("\n");
       
            NN = min(NUMNP, I + 49);
            
            for (J = I, j = J - 1; J <= NN; J++, j = J - 1)
            {
                if (Get(X, j, 0, NDM) == BLANK)
                    printf("%d Has not been input or generated\n", J);
                else
                {
                    printf("%10d", J);
                    for (l = 0; l < NDM; l++) printf("%13.4lf", Get(X, j, l, NDM));
                    printf("\n");
                }
            }
        }
    }

    traceOUT("GENVEC");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Data input routine for mesh description                                   */
/*                                                                           */
/*****************************************************************************/

void PMESH(int *IDL, int *IE, double *D, int *ID, double *X, int *IX, double *F, double *T,
           int NDF, int NDM, int NEN1, int III)
{
    int    PRT = TRUE, ERR = FALSE, loop = TRUE;
    double XL[3] = {0.0}, FL[6] = {0.0};
    int    I = 0, i = 0, J = 0, K = 0, k = 0, L = 0, l = 0, NN = 0, n = 0;
    char   buffer[LINE_LENGTH + 1] = "", line[LINE_LENGTH + 1] = "", *XHED = NULL;

    traceIN("PMESH");

    /*-----------------------------------------------------------------------*/
    /* Initialise arrays                                                     */
    /*-----------------------------------------------------------------------*/

    if (III >= 0)
    {
        for (n = 0; n < NUMNP; n++)
        {
            T[n] = 0.0;
            for (i = 0; i < NDF; i++)
            {
                Set(ID, n, i, NDF, 0);
                Set(F,  n, i, NDF, 0.0);
            }
            for (i = 0; i < NDM; i++)
                Set(X,  n, i, NDM, BLANK);
        }
    }

    loop = TRUE;
    while (loop)
    {
        char *CC = gets(line);

        if (CC == NULL)
        {
            printf("     **FATAL ERROR 53** Unexpected end-of-file\n");
            loop = FALSE;
        }
        else
        {
            CC = strip(CC);

            /*---------------------------------------------------------------*/
            /* Nodal coordinate data input                                   */
            /*---------------------------------------------------------------*/
           
            if (strncmp(CC, "COOR", 4) == 0)
            {
                GENVEC(NDM, X, CD, PRT, &ERR);
            }
           
            /*---------------------------------------------------------------*/
            /* Element data input                                            */
            /*---------------------------------------------------------------*/
           
            else if (strncmp(CC, "ELEM", 4) == 0)
            {
                int LK = 0, LX = 0, NX = 0;
           
                L = 0;      l = L - 1;
                for (I = 1, i = I - 1; (I <= NUMEL) && loop; I += 50, i = I - 1)
                {
                    if (PRT)
                    {
                        printf(O);
                        printf("%s\n\n     Elements\n\n   Element  Material", HEAD);
                        for (K = 1; K <= NEN; K++) printf("%3d Node", K);
                        printf("\n");
                    }
           
                    J = min(NUMEL, I + 49);
           
                    for (NN = I, n = NN - 1; (NN <= J) && loop; NN++, n = NN - 1)
                    {
                        if ((L-NN) < 0)
                        {
                            char *data = gets(line);

                            if (data == NULL)
                            {
                                printf("     **FATAL ERROR 56** Unexpected end-of-file\n");
                                loop = FALSE;
                            }
                            else
                            {
                                L  = read_int(line,  0,  5);         l = L - 1;
                                LK = read_int(line,  5, 10);
                                for (k = 0; k < NEN; k++) IDL[k] = read_int(line, (10+5*k), (15+5*k));
                                LX = read_int(line, (10+5*NEN), (15+5*NEN));
                               
                                if (L  == 0) { L  = NUMEL + 1;       l = L - 1; }
                                if (LX == 0)   LX = 1;
                                if ((L-NN) < 0)
                                {
                                    printf("     **ERROR 03** Element %d appears after element %d\n", L, NN);
                                    ERR = TRUE;
                                }
                            }
                        }
                       
                        if (loop && ((L-NN) == 0))
                        {
                            NX = LX;
                           
                            for (k = 0; k < NEN; k++)
                                Set(IX, l, k, NEN1, IDL[k]);
           
                            Set(IX, l, NEN1-1, NEN1, LK);
                        }
           
                        else if (loop && ((L-NN) > 0))
                        {
                            Set(IX, n, NEN1-1, NEN1, Get(IX, n-1, NEN1-1, NEN1));
           
                            for (k = 0; k < NEN; k++)
                            {
                                if (Get(IX, n-1, k, NEN1) == 0)
                                    Set(IX, n, k, NEN1, 0);
                                else
                                    Set(IX, n, k, NEN1, Get(IX, n-1, k, NEN1) + NX);
                            }
                        }
           
                        if (loop && PRT && !ERR)
                        {
                            printf("%10d%10d", NN, Get(IX, n, NEN1-1, NEN1));
                            for (k = 0; k < NEN; k++) printf("%8d", Get(IX, n, k, NEN1));
                            printf("\n");
                        }
                    }
                }
            }
           
            /*---------------------------------------------------------------*/
            /* Material data input                                           */
            /*---------------------------------------------------------------*/
           
            else if (strncmp(CC, "MATE", 4) == 0)
            {
                printf(O);
                printf("%s\n\n     Material Properties\n", HEAD);
           
                for (n = 0; n < NUMMAT; n++)
                {
                    int ma = 0, NST = 0;
                    double DUM = 0.0, S = 0.0, P = 0.0;
           
                    gets(line);
                    MA      = read_int   (line,  0,  5);   ma  = MA  - 1;
                    IEL     = read_int   (line,  5, 10);
                    XHED    = read_string(line, 10, 80, buffer);
           
                    printf("\n     Material set %d for element type %d - %s\n\n", MA, IEL, strip(XHED));
           
                    IE[ma] = IEL;
                    ELMLIB(&D[ma], &DUM, X, IX, T, &S, &P, NDF, NDM, NST, 1);
                }
            }
           
            /*---------------------------------------------------------------*/
            /* Read in the restraint conditions for each node                */
            /*---------------------------------------------------------------*/
           
            else if (strncmp(CC, "BOUN", 4) == 0)
            {
                int more = TRUE, more_input = TRUE, found = TRUE, NG = 0, LG = 0;
           
                if (PRT)
                {
                    printf(O);
                    printf("%s\n\n     NODAL B.C.\n\n      NODE", HEAD);
                    for (I = 1; I <= NDF; I++) printf("%7d%6.6s", I, BC);
                    printf("\n");
                }
           
                III = 1;
                NN   = 0;
                NG  = 0;
           
                more_input = TRUE;
                while (more_input)
                {
                    L  = NN;          l = L - 1;
                    LG = NG;
           
                    gets(line);
                    NN  = read_int(line,  0,  5);   n = NN - 1;
                    NG  = read_int(line,  5, 10);
                    for (i = 0; i < 6; i++)
                        IDL[i] = read_int(line, 10 + 5*i, 15 + 5*i);
           
                    if ((NN <= 0) || (NN > NUMNP))
                        more_input = FALSE;
                    else
                    {
                        for (i = 0; i < NDF; i++)
                        {
                            if ((L != 0) && (IDL[i] == 0) && (Get(ID, l, i, NDF) < 0))
                                Set(ID, n, i, NDF, -1);
                            else
                                Set(ID, n, i, NDF, IDL[i]);
                        }
                       
                        LG = ISIGN(LG, NN-L);
                       
                        more = TRUE;
                        while (more)
                        {
                            L = L + LG;          l = L - 1; 
           
                            if ((NN-L)*LG <= 0)
                                more = FALSE;
                            else
                            {
                                for (i = 0; i < NDF; i++)
                                    if (Get(ID, l-LG, i, NDF) < 0)    Set(ID, l, i, NDF, -1);
                            }
                        }
                    }
                }
                   
                for (NN = 1, n = NN - 1; NN <= NUMNP; NN++, n = NN - 1)
                {
                    found = FALSE;
                    for (i = 0; (i < NDF) && !found; i++)
                        if (Get(ID, n, i, NDF) != 0) found = TRUE;
                    
                    if (PRT && found)
                    {
                        printf("%10d", NN);
                        for (i = 0; i < NDF; i++)
                            printf("%13d", Get(ID, n, i, NDF));
                        printf("\n");
                    }
                }
            }
           
            /*---------------------------------------------------------------*/
            /* Force/displacement data input                                 */
            /*---------------------------------------------------------------*/
           
            else if (strncmp(CC, "FORC", 4) == 0)
            {
                GENVEC(NDF, F, FD, PRT, &ERR);
            }
           
            /*---------------------------------------------------------------*/
            /* Temperature data input                                        */
            /*---------------------------------------------------------------*/
           
            else if (strncmp(CC, "TEMP", 4) == 0)
            {
                GENVEC(1, T, TE, PRT, &ERR);
            }
           
            /*---------------------------------------------------------------*/
            /* "END "                                                        */
            /*---------------------------------------------------------------*/
           
            else if (strncmp(CC, "END", 3) == 0)
            {
                if (ERR) exit(1);
                loop = FALSE;
            }
           
            /*---------------------------------------------------------------*/
            /* "PRIN",                                                       */
            /*---------------------------------------------------------------*/
           
            else if (strncmp(CC, "PRIN", 4) == 0)
            {
                PRT = TRUE;
            }
           
            /*---------------------------------------------------------------*/
            /* "NOPR",                                                       */
            /*---------------------------------------------------------------*/
           
            else if (strncmp(CC, "NOPR", 4) == 0)
            {
                PRT = FALSE;
            }
           
            /*---------------------------------------------------------------*/
            /* "PAGE",                                                       */
            /*---------------------------------------------------------------*/
           
            else if (strncmp(CC, "PAGE", 4) == 0)
            {
                strcpy(O, strip(gets(line)));
            }
           
            else if (*CC != 0)
            {
                printf("     **FATAL ERROR 54** Unexpected keyword '%s'\n", CC);
                loop = FALSE;
            }
        }
    }

    traceOUT("PMESH");
}

/*****************************************************************************/
/*                                                                           */
/* Monitor available memory in blank common                                  */
/*                                                                           */
/*****************************************************************************/

void SETMEM(int J)
{
    int K = J;

    if (K > MAX)
    {
        printf("**ERROR 01** Insufficient storage in blank common\n");
        printf("                 required  = %d\n", K);
        printf("                 available = %d\n", MAX);
        exit(1);
    }

    return;
}

/******************************************************************************
*******************************************************************************
*******************************************************************************
*                                                                             *
* SOLUTION AND OUTPUT MODULES                                                 *
*                                                                             *
*******************************************************************************
*******************************************************************************
******************************************************************************/
 



/*****************************************************************************/
/*                                                                           */
/* Macro instruction subprogram                                              */
/*                                                                           */
/* Controls problem solution and output algorithms by order of specifying    */
/* macro commands in array wd.                                               */
/*                                                                           */
/*****************************************************************************/

void PMACR(double *UL, double *XL, double *TL, int *LD, double *P, double *S, int *IE, double *D,
           int *ID, double *X, int *IX, double *F, double *T, int *JDIAG, double *B, double *DR,
           CONTROL *CT, int NDF, int NDM, int NEN1, int NST, int NEND)
{
    int    AFR = 0, BFR = 0, CFR = 0, AFL = 0, BFL = 0, CFL = 0, DFL = 0, EFL = 0, FFL = 0, GFL = 0;
    char   line[LINE_LENGTH + 1] = "";
    double RNMAX = 0.0, TOL = 0.0, UN = 0.0;
    int    LL = 0, ll = 0, NNEQ = 0, NPLD = 0, LMAX = 0, more  = 0;
    int    LX = 0, L = 0, l = 0, I = 0, K = 0, J = 0, i = 0, LV = 0, lv = 0, lx = 0;
    int    NA = 0, NC = 1, NM = 0, NN = 0, NV = 0;
    int    NQ = 0, NR = 0, NE = 0;
    int    LVE[9] = {0}, LVS[9] = {0};

    traceIN("PMACR");

    /*-----------------------------------------------------------------------*/
    /* Set initial values of parameters                                      */
    /*-----------------------------------------------------------------------*/

    DT    = 0.0;
    PROP  = 1.0;
    RNMAX = 0.0;
    TIME  = 0.0;
    TOL   = 1.E-9;
    UN    = 0.0;
    AFL   = TRUE;
    AFR   = FALSE;
    BFL   = TRUE;
    BFR   = FALSE;
    CFL   = TRUE;
    CFR   = FALSE;
    DFL   = TRUE;
    EFL   = TRUE;
    FFL   = FALSE;
    GFL   = TRUE;
    NE    = NEND;
    NNEQ  = NDF * NUMNP;
    NPLD  = 0;

    /*-----------------------------------------------------------------------*/
    /* Read macro cards                                                      */
    /*-----------------------------------------------------------------------*/

    printf(O);
    printf("%s\n\n", HEAD);
    printf("     Macro instructions\n\n");
    printf("     Macro statement     Variable 1     Variable 2\n");

    LL = 1;                  ll = LL - 1;
    LMAX = 16;
    SETMEM(NE + LMAX*4*IPR);

    strcpy(CT[ll].s1, "LOOP");
    CT[ll].d1 = 1.0;

    more = TRUE;
    while (more)
    {
        char *data = NULL;
        LL = LL + 1;         ll = LL - 1;

        if (LL >= LMAX)
        {
            LMAX = LMAX + 16;
            SETMEM(NE + LMAX*4*IPR);
        }

        memset(line , 0, sizeof(line));
        data = gets(line);

        if (data == NULL)
        {
            printf("     **FATAL ERROR 57** Unexpected end-of-file\n");
            more = FALSE;
        }
        else
        {
            memset(&CT[ll], 0, sizeof(CONTROL));
            strncpy(CT[ll].s1, &line[0], 4);
            strncpy(CT[ll].s2, &line[5], 4);
            CT[ll].d1 = read_int(line, 10, 25);
            CT[ll].d2 = read_int(line, 25, 40);
           
            printf("          %-4.4s %-4.4s %15.5lf %15.5lf\n", CT[ll].s1, CT[ll].s2, CT[ll].d1, CT[ll].d2);
           
            if (strncmp(CT[ll].s1, "END", 3) == 0)   more = FALSE;
        }
    }

    strcpy(CT[ll].s1, "NEXT");

    /*-----------------------------------------------------------------------*/
    /* Set loop markers                                                      */
    /*-----------------------------------------------------------------------*/

    NE = NE + LMAX * 4 * IPR;
    LX = LL - 1;
    for (L = 1, l = L - 1; L <= LX; L++, l = L - 1)
    {
        if (strncmp(CT[l].s1, "LOOP", 4) == 0)
        {
            J = 1;
            K = L + 1;
            for (I = K, i = I - 1; (I <= LL) && (J != 0); I++, i = I - 1)
            {
                if      (strncmp(CT[i].s1, "LOOP", 4) == 0) J = J + 1;
                else if (strncmp(CT[i].s1, "NEXT", 4) == 0) J = J - 1;
               
                if (J > 9)
                {
                    printf("     **FATAL ERROR 11** Loops nested deeper than 8\n");
                    return;
                }
            }

            if (J != 0)
            {
                printf("     **FATAL ERROR 10** Unbalanced loop/next macros\n");
                return;
            }

            CT[i].d2 = (double) L;
            CT[l].d2 = (double) I;
        }
    }

    J = 0;

    for (l = 0; l < LL; l++)
    {
        if      (strncmp(CT[l].s1, "LOOP", 4) == 0) J = J + 1;
        else if (strncmp(CT[l].s1, "NEXT", 4) == 0) J = J - 1;
    }

    if (J != 0)
    {
        printf("     **FATAL ERROR 10** Unbalanced loop/next macros\n");
        return;
    }

    /*-----------------------------------------------------------------------*/
    /* Execute macro instruction program                                     */
    /*-----------------------------------------------------------------------*/

    LV = 0;
    L  = 1;                            l = L - 1;

    more = TRUE;
    while (more)
    {
        I = L - 1;                     i = I - 1;
      
        if ((L != 1) && (L != LL))
        {
            printf("  **Macro instruction %4d executed**  ", I);
            printf("  %-4.4s %-4.4s %15.5lf %15.5lf\n", CT[l].s1, CT[l].s2, CT[l].d1, CT[l].d2);
        }

        /*-------------------------------------------------------------------*/
        /* Set solution tolerance                                            */
        /*-------------------------------------------------------------------*/
      
        if (strncmp(CT[l].s1, "TOL", 3) == 0)
        {
            TOL = CT[l].d1;
        }
      
        /*-------------------------------------------------------------------*/
        /* Set time increment                                                */
        /*-------------------------------------------------------------------*/

        else if (strncmp(CT[l].s1, "DT", 2) == 0)
        {
            DT = CT[l].d1;
        }
      
        /*-------------------------------------------------------------------*/
        /* Print stress values                                               */
        /*-------------------------------------------------------------------*/
      
        else if (strncmp(CT[l].s1, "STRE", 4) == 0)
        {
            LX = LVE[lv];              lx = LX - 1;

            if (AMOD(CT[lx].d1, max(CT[l].d1, 1.0)) == 0.0)
            {
                PFORM(UL, XL, TL, LD, P, S, IE, D, ID, X, IX, F, T, JDIAG, DR, DR, DR,
                      NDF, NDM, NEN1, NST, 4, B, (double*) (void*) &M[NV], FALSE, FALSE, FALSE, FALSE);
            }
        }
      
        /*-------------------------------------------------------------------*/
        /* Print displacement                                                */
        /*-------------------------------------------------------------------*/
      
        else if (strncmp(CT[l].s1, "DISP", 4) == 0)
        {
            LX = LVE[lv];              lx = LX - 1;

            if (AMOD(CT[lx].d1, max(CT[l].d1, 1.0)) == 0.0)
            {
                printf(O);
                printf("%s\n\n", HEAD);
                printf("          Time %13.5lf\n\n", TIME);
                printf("     Proportional Load %13.5lf\n", PROP);

                PRTDIS(ID, X, B, F, NDM, NDF);
            }
        }
      
        /*-------------------------------------------------------------------*/
        /* Form tangent stiffness                                            */
        /*-------------------------------------------------------------------*/
      
        else if (strncmp(CT[l].s1, "TANG", 4) == 0)
        {
            if (J == 5) CFR = FALSE;
            if (GFL)    PSETM(&NA, &NE, JDIAG[NEQ-1]*IPR, &GFL);
            if (NPLD > 0) PROP = PROPLD(TIME, 0);
            PZERO((double*) (void*) &M[NA], JDIAG[NEQ-1]);
            PFORM(UL, XL, TL, LD, P, S, IE, D, ID, X, IX, F, T, JDIAG, DR, (double*) (void*) &M[NA], (double*) (void*) &M[NC],
                  NDF, NDM, NEN1, NST, 3, B, (double*) (void*) &M[NV], TRUE, FALSE, CFR, FALSE);
            AFR = TRUE;
        }
      
        else if (strncmp(CT[l].s1, "UTAN", 4) == 0)
        {
            if (CFL)    PSETM(&NC, &NE, JDIAG[NEQ-1]*IPR, &CFL);
            PZERO((double*) (void*) &M[NC], JDIAG[NEQ-1]);
            CFR = TRUE;
            if (J == 5) CFR = FALSE;
            if (GFL)    PSETM(&NA, &NE, JDIAG[NEQ-1]*IPR, &GFL);
            if (NPLD > 0) PROP = PROPLD(TIME, 0);
            PZERO((double*) (void*) &M[NA], JDIAG[NEQ-1]);
            PFORM(UL, XL, TL, LD, P, S, IE, D, ID, X, IX, F, T, JDIAG, DR, (double*) (void*) &M[NA], (double*) (void*) &M[NC],
                  NDF, NDM, NEN1, NST, 3, B, (double*) (void*) &M[NV], TRUE, FALSE, CFR, FALSE);
            AFR = TRUE;
        }
      
        /*-------------------------------------------------------------------*/
        /* Form out of balance force for time step/iteration                 */
        /*-------------------------------------------------------------------*/
      
        else if (strncmp(CT[l].s1, "FORM", 4) == 0)
        {
            double RN = 0.0;
            int n = 0;

            if (NPLD > 0) PROP = PROPLD(TIME, 0);
            PLOAD(ID, F, DR, NNEQ, PROP);
            PFORM(UL, XL, TL, LD, P, S, IE, D, ID, X, IX, F, T, JDIAG, DR, DR, DR,
                  NDF, NDM, NEN1, NST, 6, B, (double*) (void*) &M[NV], FALSE, TRUE, FALSE, FALSE);
            BFR = TRUE;

            RN = 0.0;
            for (n = 0; n < NEQ; n++)
                RN = RN + DR[n] * DR[n];
            RN = sqrt(RN);
               
            RNMAX = max(RNMAX, RN);

            printf("     Force Convergence Test\n");
            printf("          RNMAX = %15.5lf     RN    = %15.5lf     TOL   = %15.5lf\n", RNMAX,RN,TOL);

            if (RN < RNMAX*TOL)
            {
                int L0 = 0, l0 = 0;

                LX = LVE[lv];          lx = LX - 1;
                L0 = LVS[lv];          l0 = L0 - 1;
                CT[lx].d1 = CT[l0].d1;
                L = LX - 1;            l  = L  - 1;
            }
        }
      
        /*-------------------------------------------------------------------*/
        /* Set loop start indicators                                         */
        /*-------------------------------------------------------------------*/
      
        else if (strncmp(CT[l].s1, "LOOP", 4) == 0)
        {
            LV = LV + 1;                         lv = LV - 1;
            LX = (int) (0.5 + CT[l].d2);         lx = LX - 1;
            LVS[lv] = L;
            LVE[lv] = LX;
            CT[lx].d1 = 1.0;
        }
      
        /*-------------------------------------------------------------------*/
        /* Loop terminator control                                           */
        /*-------------------------------------------------------------------*/
      
        else if (strncmp(CT[l].s1, "NEXT", 4) == 0)
        {
            int n = 0;

            NN = (int) CT[l].d2;                           n  = NN  - 1;
            CT[l].d1 = CT[l].d1 + 1.0;
            if (CT[l].d1 >  CT[n].d1) LV = LV - 1;         lv = LV - 1;
            if (CT[l].d1 <= CT[n].d1) L  = NN;             l  = L  - 1;
        }
      
        /*-------------------------------------------------------------------*/
        /* Input proportional load table                                     */
        /*-------------------------------------------------------------------*/
      
        else if (strncmp(CT[l].s1, "PROP", 4) == 0)
        {
            NPLD = (int) CT[l].d1;
            PROP = PROPLD(0.0, NPLD);
        }
      
        /*-------------------------------------------------------------------*/
        /* Read commands                                                     */
        /*-------------------------------------------------------------------*/
      
        else if (strncmp(CT[l].s1, "DATA", 4) == 0)
        {
            CONTROL CTL = {{0}};

            gets(line);
            memset(&CTL, 0, sizeof(CONTROL));
            strncpy(CTL.s1, &line[0], 4);
            strncpy(CTL.s2, &line[5], 4);
            CTL.d1 = read_int(line, 10, 25);
            CTL.d2 = read_int(line, 25, 40);

            if (strcmp(CT[l].s2, CTL.s1) != 0)
            {
                printf("**FATAL ERROR 12** Macro label mismatch on a read command\n");
                return;
            }

            if (strncmp(CTL.s1, "TOL", 3) == 0) TOL = CTL.d1;
            if (strncmp(CTL.s1, "DT",  2) == 0) DT  = CTL.d1;
        }
      
        /*-------------------------------------------------------------------*/
        /* Increment time                                                    */
        /*-------------------------------------------------------------------*/
      
        else if (strncmp(CT[l].s1, "TIME", 4) == 0)
        {
            TIME  = TIME + DT;
            RNMAX = 0.0;
            UN    = 0.0;
        }
      
        /*-------------------------------------------------------------------*/
        /* Compute convergence test                                          */
        /*-------------------------------------------------------------------*/
      
        else if (strncmp(CT[l].s1, "CONV", 4) == 0)
        {
            double RN = 0.0, CN = 0.0;
            int n = 0, L0 = 0, l0 = 0;

            for (n = 0; n < NEQ; n++)
            {
                UN = UN + B[n]  * B[n];
                RN = RN + DR[n] * DR[n];
            }

            UN = max(UN, RN);
            CN = sqrt(UN);
            RN = sqrt(RN);

            printf("     Displacement Convergence Test\n");
            printf("          UNMAX = %15.5lf     UN    = %15.5lf      TOL   = %15.5lf\n", CN, RN, TOL);

            LX = LVE[lv];                                  lx = LX - 1;
            L0 = LVS[lv];                                  l0 = L0 - 1;
            if (RN < CN*TOL) CT[lx].d1 = CT[l0].d1;
        }

        /*-------------------------------------------------------------------*/
        /* Solve the equations                                               */
        /*-------------------------------------------------------------------*/
      
        else if (strncmp(CT[l].s1, "SOLV", 4) == 0)
        {
            if (CFR)
                UACTCL((double*) (void*) &M[NA], (double*) (void*) &M[NC], DR, JDIAG, NEQ, AFR, BFR);
            else
                ACTCOL((double*) (void*) &M[NA], (double*) DR, JDIAG, NEQ, AFR, BFR);

            AFR = FALSE;

            if (BFR)
            {
                int n = 0;

                BFR = FALSE;

                for (n = 0; n < NEQ; n++)
                    B[n] = B[n] + DR[n];
            }
        }

        /*-------------------------------------------------------------------*/
        /* Form a lumped mass approximation                                  */
        /*-------------------------------------------------------------------*/
      
        else if (strncmp(CT[l].s1, "LMAS", 4) == 0)
        {

            AFL = FALSE;
            BFL = TRUE;
            if (EFL) PSETM(&NN, &NE, NEQ*IPR, &EFL);
            PZERO((double*) (void*) &M[NN], NEQ);
            PFORM(UL, XL, TL, LD, P, S, IE, D, ID, X, IX, F, T, JDIAG, (double*) (void*) &M[NN], (double*) (void*) &M[NM],
                  (double*) (void*) &M[NM], NDF, NDM, NEN1, NST, 5, B, (double*) (void*) &M[NV], AFL, BFL, FALSE, FALSE);
        }

        /*-------------------------------------------------------------------*/
        /* Form a consistent mass approximation                              */
        /*-------------------------------------------------------------------*/
      
        else if (strncmp(CT[l].s1, "CMAS", 4) == 0)
        {
            AFL = TRUE;
            BFL = FALSE;
            if (DFL)     PSETM(&NM, &NE, JDIAG[NEQ-1]*IPR, &DFL);
            PZERO((double*) (void*) &M[NM], JDIAG[NEQ-1]);
            PFORM(UL, XL, TL, LD, P, S, IE, D, ID, X, IX, F, T, JDIAG, (double*) (void*) &M[NN], (double*) (void*) &M[NM],
                  (double*) (void*) &M[NM], NDF, NDM, NEN1, NST, 5, B, (double*) (void*) &M[NV], AFL, BFL, FALSE, FALSE);
        }

        else if (strncmp(CT[l].s1, "MESH", 4) == 0)
        {
            I = -1;
            PMESH(LD, IE, D, ID, X, IX, F, T, NDF, NDM, NEN1, I);

            if (I > 0)
            {
                printf("     **FATAL ERROR 14** Attempt to change boundary restraint codfs during macro execution\n");
                return;
            }
        }

        else if (strncmp(CT[l].s1, "EIGE", 4) == 0)
        {
            J = NM;
            if (DFL) J = NN;

            PEIGS((double*) (void*) &M[NA], (double*) (void*) &M[J], F, X, B, DR, ID, IX, JDIAG, NDF, NDM, NEN1, DFL);
        }

        else if (strncmp(CT[l].s1, "EXCD", 4) == 0)
        {
            if (!FFL)
            {
                /*-----------------------------------------------------------*/
                /* Macro *EXCD* explicit integration of equations of motion  */
                /*-----------------------------------------------------------*/
               
                NQ = NE;
                NV = NQ + NEQ*IPR;
                NR = NV + NDF*NUMNP*IPR;
                NE = NR + NEQ*IPR;
                SETMEM(NE+1);
                PZERO((double*) (void*) &M[NQ], NE-NQ);
                FFL = TRUE;
            }
            else if (!BFR || EFL)
            {
                printf("     **FATAL ERROR 13** Macro exed must be preceded by lmas and form\n");
                return;
            }
            else
                EUPDAT(DR, B, (double*) (void*) &M[NQ], (double*) (void*) &M[NR], (double*) (void*) &M[NN], DT, NEQ);
        }

        /*-------------------------------------------------------------------*/
        /* Compute reactions and print                                       */
        /*-------------------------------------------------------------------*/
      
        else if (strncmp(CT[l].s1, "REAC", 4) == 0)
        {
            PZERO(DR, NNEQ);
            PFORM(UL, XL, TL, LD, P, S, IE, D, ID, X, IX, F, T, JDIAG, DR, DR, DR,
                  NDF, NDM, NEN1, NST, 6, B, (double*) (void*) &M[NV], FALSE, TRUE, FALSE, TRUE);
            PRTREA(DR, NDF);
        }

        else if (strncmp(CT[l].s1, "CHEC", 4) == 0)
        {
            PFORM(UL, XL, TL, LD, P, S, IE, D, ID, X, IX, F, T, JDIAG, DR, DR, DR,
               NDF, NDM, NEN1, NST, 2, B, F, FALSE, FALSE, FALSE, FALSE);
        }
      
        /*-------------------------------------------------------------------*/
        /* Output the JDIAG array                                            */
        /*-------------------------------------------------------------------*/
      
        else if (strncmp(CT[l].s1, "OUTJ", 4) == 0)
        {
            int j = 0;

            printf(O);
            printf("%s\n\n     The JDIAG array\n\n", HEAD);

            for (j = 0; j < NEQ; j++)
            {
                printf("    %4d: %7d\n", j+1, JDIAG[j]);
            }
        }
      
        /*-------------------------------------------------------------------*/
        /* Output the assembled matrix                                       */
        /*-------------------------------------------------------------------*/
      
        else if (strncmp(CT[l].s1, "OUTP", 4) == 0)
        {
            printf("%c",     form_feed);
            printf("%s\n", HEAD);

            dspMatrix(NEQ, JDIAG, (double*) (void*) &M[NA], (double*) DR, "Macro command 'OUTP'");
        }
      
        L = L + 1;                  l = L - 1;
        if (L > LL) more = FALSE;
    }

    traceOUT("PMACR");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Active column profile symmetric equation solver                           */
/*                                                                           */
/*****************************************************************************/

void ACTCOL(double *A, double *B, int *JDIAG, int NEQ, int AFAC, int BACK)
{
    double D = 0.0;
    int jh = 0, ie = 0, ir = 0, ih = 0;
    int j = 0, k = 0, is = 0, id = 0, i = 0, jd = 0, jr = 0;

    traceIN("ACTCOL");

    /*-----------------------------------------------------------------------*/
    /* Factor A to UT*D*U, reduce B                                          */
    /*-----------------------------------------------------------------------*/

    printf("-----------------------------------------------------------------\n");

    AENGY = 0.0;
    jr = -1;
    for (j = 0; j < NEQ; j++)
    {
        printf("j = %d\n\n", j);

        jd = JDIAG[j] - 1;
        jh = jd - jr;
        is = j - jh + 2;

        if ((jh-2 > 0) && AFAC)
        {
            ie = j;
            k = jr + 2;
            id = JDIAG[is - 1] - 1;
       
            /*---------------------------------------------------------------*/
            /* Reduce all equations except diagonal                          */
            /*---------------------------------------------------------------*/
         
            for (i = is; i < ie; i++)
            {
                ir = id;
                id = JDIAG[i] - 1;
                ih = min(id-ir-1, i-is+1);
         
                if (ih > 0)
                {
                    double dot = DOT(&A[k-ih], &A[id-ih], ih);
                    double x = A[k] - dot;

                    printf("A[%d] = A[%d] - DOT(&A[%d], &A[%d], %d)\n", k, k, k-ih, id-ih, ih);
                    printf("     = %7.3lf - %7.3lf = %7.3lf\n\n", A[k], dot, x);

                    A[k] = x;
                }
         
                k = k + 1;
            }
        }

        /*-------------------------------------------------------------------*/
        /* Reduce diagonal term                                              */
        /*-------------------------------------------------------------------*/

        if ((jh-2 >= 0) && AFAC)
        {
            ir = jr + 1;
            ie = jd;
            k  = j - jd - 1;
            for (i = ir; i < ie; i++)
            {
                id = JDIAG[k+i+1] - 1;

                if (A[id] != 0.0)
                {
                    double x = A[i] / A[id];
                    D = A[i];

                    printf("A[%d] = A[%d] / A[%d]\n", i, i, id);
                    printf("     = %7.3lf / %7.3lf = %7.3lf\n\n", A[i], A[id], x);

                    A[i] = x;
                    x = A[jd] - D*A[i];

                    printf("A[%d] = A[%d] - D*A[%d]\n", jd, jd, i);
                    printf("     = %7.3lf - (%7.3lf * %7.3lf) = %7.3lf\n\n", A[jd], D, A[i], x);

                    A[jd] = x;
                }
            }
        }

        /*-------------------------------------------------------------------*/
        /* Reduce RHS                                                        */
        /*-------------------------------------------------------------------*/

        if (BACK)
        {
            double dot = DOT(&A[jr+1], &B[is-1], jh-1);
            double x   = B[j] - dot;

            printf("--------------------\n");

            printf("B[%d] = B[%d] - DOT(&A[%d], &B[%d], %d)\n", j, j, jr+1, is-1, jh-1);
            printf("     = %7.3lf - %7.3lf = %7.3lf\n\n", B[j], dot, x);

            B[j] = x;
        }

        jr = jd;
        printf("--------------------\n");
        dspMatrix(NEQ, JDIAG, A, B, "Subroutine 'ACTCOL'");
        printf("-----------------------------------------------------------------\n");
    }

    if ( ! BACK) return;

    printf("Divide by diagonal pivots\n\n");

    /*-----------------------------------------------------------------------*/
    /* Divide by diagonal pivots                                             */
    /*-----------------------------------------------------------------------*/

    for (i = 0; i < NEQ; i++)
    {
        id = JDIAG[i] - 1;

        if (A[id] != 0.0)
        {
            double x = B[i] / A[id];

            printf("B[%d] = B[%d] / A[%d]\n", i, i, id);
            printf("     = %7.3lf / %7.3lf = %7.3lf\n\n", B[i], A[id], x);

            B[i] = x;
        }

        AENGY = AENGY + B[i] * B[i] * A[id];
    }

    /*-----------------------------------------------------------------------*/
    /* Backsubstitute                                                        */
    /*-----------------------------------------------------------------------*/

    j = NEQ - 1;  
    jd = JDIAG[j] - 1;

    printf("Backsubstitute\n\n");

    while (j >= 0)
    {
        D = B[j];
        j = j - 1;

        if (j >= 0)
        {
            jr = JDIAG[j] - 1;
           
            if (jd-jr > 1)
            {
                is = j - jd + jr + 2;
                k = jr - is;
                for (i = is; i <= j; i++)
                {
                    double x = B[i] - A[i+k+1] * D;

                    printf("B[%d] = B[%d] - A[%d] * D\n", i, i, i+k+1);
                    printf("     = %7.3lf - (%7.3lf * %7.3lf)  = %7.3lf\n\n", B[i], A[i+k+1], D , x);

                    B[i] = x;
                }
            }
           
            jd = jr;
        }
    }

    traceOUT("ACTCOL");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Assemble global arrays                                                    */
/*                                                                           */
/*****************************************************************************/

void ADDSTF(double *A, double *B, double *C, double *S, double *P,
            int *JDIAG, int *LD, int NST, int NFL, int AFL, int BFL, int CFL)
{
    int i = NFL, j = 0, k = 0;

    traceIN("ADDSTF: NST = %d, NFL = %d, AFL = %d, BFL = %d, CFL = %d", NST, NFL, AFL, BFL, CFL);

    for (j = 0; j < NST; j++)
    {
        k = LD[j] - 1;

        if (k >= 0)
        {
            if (BFL)     B[k] = B[k] + P[j];

            if (AFL || CFL)
            {
                int l = JDIAG[k] - k - 1;

                for (i = 0; i < NST; i++)
                {
                    int m = LD[i] - 1;

                    if ((m <= k) && (m >= 0))
                    {
                        int n = l + m;
                        int height = k - m + 1;
                        int max_height = JDIAG[k] - l;

                        if (height > max_height)
                            printf("----> ERROR: ");
                        else
                            printf("----> OK:    ");

                        printf("j = %2d, k = %2d, i = %2d, m = %2d, n = %2d", j, k, i, m, n);
                        printf(", height = %2d",                              k - m + 1);
                        printf(", max_height = JDIAG[%2d] - l",               k);
                        printf(" = %3d - %3d",                                JDIAG[k], l);
                        printf(" = %3d",                                      JDIAG[k] - l);
                        printf(": S[%d][%d] = %12.6lf\n",                     j, i, Get(S, j, i, NST));

                        if (AFL) A[n] = A[n] + Get(S, j, i, NST);
                        if (CFL) C[n] = C[n] + Get(S, i, j, NST);
                    }
                }
            }
        }
    }

    dspMatrix(NEQ, JDIAG, A, B, "ADDSTF");

    traceOUT("ADDSTF");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Vector dot product                                                        */
/*                                                                           */
/*****************************************************************************/

double DOT(double *A, double *B, int NN)
{
    int i = 0;
    double dot = 0.0;

    for (i = 0; i < NN; i++)
      dot += A[i] * B[i];

    return dot;
}

/*****************************************************************************/
/*                                                                           */
/* Element library                                                           */
/*                                                                           */
/*****************************************************************************/

void ELMLIB(double *D, double *U, double *X, int *IX, double *T, double *S, double *P, int I, int J, int K, int ISW)
{
    traceIN("ELMLIB: ISW = %d, IEL = %d", ISW, IEL);

    if (ISW >= 3)
    {
        int l = 0;

        for (l = 0; l < K; l++)
        {
            int m = 0;

            P[l] = 0.0;
           
            for (m = 0; m < K; m++)
                Set(S, m, l, K, 0.0);
        }
    }

    switch (IEL)
    {
    case 1:    ELMT01(D,U,X,IX,T,S,P,I,J,K,ISW);    break;
    case 2:    ELMT02(D,U,X,IX,T,S,P,I,J,K,ISW);    break;
    case 3:    ELMT03(D,U,X,IX,T,S,P,I,J,K,ISW);    break;
    default:  
        printf("    **FATAL ERROR 04** Element class number %d  input\n", IEL);
        exit(0);
    }

    traceOUT("ELMLIB");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Update solution using explicit central differences                        */
/*                                                                           */
/*****************************************************************************/

void EUPDAT(double *DR, double *U, double *V, double *A, double *XM, double DT, int NEQ)
{
    double DTHP = 0.0;
    double DTH  = DT / 2.0;
    double DTAV = DTH + DTHP;
    int n = 0;

    traceIN("EUPDAT");
    DTHP = DTH;

    for (n = 0; n < NEQ; n++)
    {
        A[n] = DR[n] / XM[n];
        V[n] = V[n] + DTAV * A[n];
        U[n] = U[n] + DT   * V[n];
    }

    traceOUT("EUPDAT");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Normalize vector Y to unit vector X                                       */
/*                                                                           */
/*****************************************************************************/

void NORM(double *X, double *Y, int NN)
{
    double SCALE = sqrt(DOT(Y, Y, NN));
    int i = 0;

    traceIN("NORM");

    for (i = 0; i < NN; i++)
        X[i] = Y[i] / SCALE;

    traceOUT("NORM");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* It may be necessary to replace the following alphanumeric                 */
/* comparison statement if computer produces an overflow                     */
/*                                                                           */
/*****************************************************************************/

int PCOMP(double A, double B)
{
    return (A == B);
}

/*****************************************************************************/
/*                                                                           */
/* Compute dominant eiganvalue by inverse iteration                          */
/*                                                                           */
/*****************************************************************************/

void PEIGS(double *A, double *B, double *F, double *X, double *Y, double *Z, int *ID, int *IX,
           int *JDIAG, int NDF, int NDM, int NEN1, int DFL)
{
    int    ITS  = 100;
    double TOL  = 1.E-9;
    double EIGP = 0.0;
    double EIG  = 0.0;
    int    converged = FALSE;
    int    i = NEN1;
    int    j = *IX;

    traceIN("PEIGS");

    /*-----------------------------------------------------------------------*/
    /* Get start vector from diagonal of mass matrix                         */
    /*-----------------------------------------------------------------------*/

    for (i = 0; i < NEQ; i++)
    {
        if (DFL)
            j = 0;
        else
            j = JDIAG[i] - 1;

        Y[i] = B[j];
    }

    EIGP = 0.0;
    ACTCOL(A, Z, JDIAG, NEQ, TRUE, FALSE);

    for (i = 0; (i < ITS) && !converged; i++)
    {
        PZERO(Z, NEQ);
        PROMUL(B, Y, Z, JDIAG, NEQ);

        /*-------------------------------------------------------------------*/
        /* Rayleigh quotient                                                 */
        /*-------------------------------------------------------------------*/
    
        EIG = AENGY / DOT(Y, Z, NEQ);

        if (fabs(EIG-EIGP) < TOL * fabs(EIG))
            converged = TRUE;
        else
        {
            NORM(Y, Z, NEQ);
            EIGP = EIG;
    
            /*---------------------------------------------------------------*/
            /* Inverse iteration                                             */
            /*---------------------------------------------------------------*/

            ACTCOL(A, Y, JDIAG, NEQ, FALSE, TRUE);
        }
    }

    if (converged)
    {
        printf(O);
        printf("%s\n\n", HEAD);
        printf("     Eigenvalue =  %13.4lf\n", EIG);
        printf("     Iterations =  %d\n\n", i);
        NORM(Z, Y, NEQ);
        PRTDIS(ID, X, Z, F, NDM, NDF);
    }
    else
        printf("     **FATAL ERROR 09** No convergence in eigenvalues, ITS =  %d\n", ITS);

    traceOUT("PEIGS");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Compute element arrays and assemble global arrays                         */
/*                                                                           */
/*****************************************************************************/

void PFORM(double *UL, double *XL, double *TL, int *LD, double *P, double *S, int *IE,
           double *D, int *ID, double *X, int *IX, double *F, double *T, int *JDIAG,
           double *B, double *A, double *C, int NDF,
           int NDM, int NEN1, int NST, int ISW, double *U, double *UD,
           int AFL, int BFL, int CFL, int DFL)
{
    /*-----------------------------------------------------------------------*/
    /* Loop on elements                                                      */
    /*-----------------------------------------------------------------------*/

    int n = 0, i = 0, j = 0;

    traceIN("PFORM");

    for (n = 0; n < NUMEL; n++)
    {
        /*-------------------------------------------------------------------*/
        /* Set up local arrays                                               */
        /*-------------------------------------------------------------------*/
       
        for (i = 0; i < NEN; i++)
        {
            int II = Get(IX, n, i, NEN1), ii = II - 1;

            if (II == 0)
            {
                TL[i] = 0.0;
              
                for (j = 0; j < NDM; j++)
                    Set(XL, i, j, NDM, 0.0);
              
                for (j = 0; j < NDF; j++)
                {
                    Set(UL, i,       j, NDF, 0.0);
                    Set(UL, i + NEN, j, NDF, 0.0);
                    Set(LD, i,       j, NDF, 0  );
                }
            }
            else
            {
                int IID = II * NDF - NDF;
                NEL = i + 1;
                TL[i] = T[ii];
               
                for (j = 0; j < NDM; j++)
                    Set(XL, i, j, NDM, Get(X, ii, j, NDM));
               
                for (j = 0; j < NDF; j++)
                {
                    int K = labs(Get(ID, ii, j, NDF));
                    Set(UL, i,       j, NDF, Get(F,  ii, j, NDF) * PROP);
                    Set(UL, i + NEN, j, NDF, Get(UD, ii, j, NDF));
                    if (K > 0) Set(UL, i, j, NDF, U[K]);
                    if (DFL) K = IID + j + 1;
                    Set(LD, i, j, NDF, K);
                }
            }
        }

        /*-------------------------------------------------------------------*/
        /* Dump the local arrays for element "1"                             */
        /*-------------------------------------------------------------------*/
       
        if (n == 0)
        {
            printf("ISW = %d, n = %d\n", ISW, n);
            print2Darray("D",  D,  "%12.2lf", 10,    NUMMAT);
            print2Darray("UL", UL, "%12.2lf", NEN,   NDF);
            print2Darray("XL", XL, "%12.2lf", NEN,   NDM);
            print2Darray("IX", IX, "%2d",     NUMEL, NEN1);
            print1Darray("TL", TL, "%lf",     NEN);
            print2Darray("LD", LD, "%2d",     NEN,   NDF);

            for (i = 0; i < NEN; i++)
            {
                int ii = Get(IX, n, i, NEN1);

                printf("         i = %d, ii = %d", i, ii);
                if (ii != 0) printf(", T[%d] = %10.4lf\n", ii, T[ii]);
            }
        }

        /*-------------------------------------------------------------------*/
        /* Form element array                                                */
        /*-------------------------------------------------------------------*/
       
        MA = Get(IX, n, NEN1-1, NEN1);
        if (IE[MA-1] != IEL) MCT = 0;
        IEL = IE[MA-1];
        ELMLIB(&Get(D, MA-1, 0, 10), UL, XL, &Get(IX, n, 0, NEN1), TL, S, P, NDF, NDM, NST, ISW);
       
        /*-------------------------------------------------------------------*/
        /* Add to total array                                                */
        /*-------------------------------------------------------------------*/
       
        if (AFL || BFL || CFL)   ADDSTF(A, B, C, S, P, JDIAG, LD, NST, NEL*NDF, AFL, BFL, CFL);

        if (n == 0)   dspMatrix(NEQ, JDIAG, A, B, "PFORM");
    }

    traceOUT("PFORM");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Form load vector in compact form                                          */
/*                                                                           */
/*****************************************************************************/

void PLOAD(int *ID, double *F, double *B, int NN, double P)
{
    int n = 0;

    traceIN("PLOAD");
    for (n = 0; n < NN; n++)
    {
        int j = ID[n] - 1;

        if (j >= 0)    B[j] = F[n] * P;
    }

    traceOUT("PLOAD");

    return;
}

/*****************************************************************************/
/*                                                                           */
/* Compute profile of global arrays                                          */
/*                                                                           */
/*****************************************************************************/

void PROFIL(int *JDIAG, int *ID, int *IX, int NDF, int NEN1, int *NAD)
{
    int i = 0, j = 0, k = 0, l = 0, n = 0, ii = 0, jj = 0;

    traceIN("PROFIL");

    /*-----------------------------------------------------------------------*/
    /* Set up the equation numbers                                           */
    /*-----------------------------------------------------------------------*/

    NEQ = 0;

    for (n = 0; n < NUMNP; n++)
    {
        for (i = 0; i < NDF; i++)
        {
            j = Get(ID, n, i, NDF);
             
            if (j == 0)
            {
                NEQ = NEQ + 1;
                Set(ID, n, i, NDF, NEQ);
                JDIAG[NEQ-1] = 0;
            }
            else
                Set(ID, n, i, NDF, 0);
        }
    }
                          
    /*-----------------------------------------------------------------------*/
    /* Compute column heights                                                */
    /*-----------------------------------------------------------------------*/

    for (n = 0; n < NUMEL; n++)
    {
        for (i = 0; i < NEN; i++)
        {
            int II = Get(IX, n, i, NEN1);                         ii = II - 1;
       
            if (II != 0)
            {
                for (k = 0; k < NDF; k++)
                {
                    int KK = Get(ID, ii, k, NDF);
                
                    if (KK != 0)
                    {
                        for (j = i; j < NEN; j++)
                        {
                            int JJ = Get(IX, n, j, NEN1);         jj = JJ - 1;
                        
                            if (JJ != 0)
                            {
                                for (l = 0; l < NDF; l++)
                                {
                                    int LL = Get(ID, jj, l, NDF);
                                
                                    if (LL != 0)
                                    {
                                        int MM = max(KK, LL),     mm = MM - 1;
                                        JDIAG[mm] = max(JDIAG[mm], abs(KK-LL));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /*-----------------------------------------------------------------------*/
    /* Compute diagonal pointers for profile                                 */
    /*-----------------------------------------------------------------------*/

    *NAD = 1;
    JDIAG[0] = 1;

    if (NEQ > 1)
    {
        for (n = 1; n < NEQ; n++)
            JDIAG[n] = JDIAG[n] + JDIAG[n-1] + 1;

        *NAD = JDIAG[NEQ -1];
    }

    traceOUT("PROFIL");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Routine to form C = C + A*B where A is a symmetric square matrix          */
/* stored in profile form, B,C are vectors, and JDIAG locates the            */
/* diagonals in A.                                                           */
/*                                                                           */
/*****************************************************************************/

void PROMUL(double *A, double *B, double *C, int *JDIAG, int NEQ)
{
    int j = 0, js = 0;

    traceIN("PROMUL");

    for (j = 0; j < NEQ; j++)
    {
        int jd = JDIAG[j] - 1;

        if (js <= jd)
        {
            double BJ = B[j];
            double AB = A[jd] * BJ;

            if (js != jd)
            {
                int jb = j - jd;
                int je = jd + 2;
                int jj = 0;
               
                for (jj = js; jj < je; jj++)
                {
                    AB       += A[jj] * B[jj+jb];
                    C[jj+jb] += A[jj] * BJ;
                }
               
                C[j] = C[j] + AB;
            }

            js = jd + 1;
        }
    }

    traceOUT("PROMUL");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Proprotional load table (one load card only)                              */
/*                                                                           */
/*****************************************************************************/

double PROPLD(double T, int J)
{
    static int    K     = 0;
    static int    L     = 0;
    static double value = 0.0;
    static double TMIN  = 0.0;
    static double TMAX  = 0.0;
    static double A[5]  = {0.0};

    traceIN("PROPLD");

    if (J <= 0)
    {
        /*-------------------------------------------------------------------*/
        /* Compute value at time T                                           */
        /*-------------------------------------------------------------------*/

        if ((T >= TMIN) && (T <= TMAX))
        {                                          
            L = max(L, 1);
            value = A[0] + A[1] * T + pow(A[2]*(sin(A[3]*T+A[4])), (double) L);
        }
    }
    else
    {
        /*-------------------------------------------------------------------*/
        /* Input table of proportional loads                                 */
        /*-------------------------------------------------------------------*/

        int  I = 1;
        char line[LINE_LENGTH + 1] = "";

        gets(line);
        K    = read_int   (line,   0,  5);
        L    = read_int   (line,   5, 10);
        TMIN = read_double(line,  10, 20);
        TMAX = read_double(line,  20, 30);
        A[0] = read_double(line,  30, 40);
        A[1] = read_double(line,  40, 50);
        A[2] = read_double(line,  50, 60);
        A[3] = read_double(line,  60, 70);
        A[4] = read_double(line,  70, 80);
                                  
        printf("     Proportional load table\n\n");
        printf("  Number    Type    Exp.   Minimum time    Maximum time      A0             A1             A2             A3             A4\n");
        printf("%8ld %8ld %8ld %15.5lf %15.5lf %15.5lf %15.5lf %15.5lf %15.5lf %15.5lf\n", I, K, L, TMIN, TMAX, A[0], A[1], A[2], A[3], A[4]);
    }

    traceOUT("PROPLD");
    return value;
}

/*****************************************************************************/
/*                                                                           */
/* Output nodal values                                                       */
/*                                                                           */
/*****************************************************************************/

void PRTDIS(int *ID, double *X, double *B, double *F, int NDM, int NDF)
{
    int i = 0, ii = 0, j = 0, jj = 0, k = 0, n = 0;
    double UL[6] = {0};

    traceIN("PRTDIS");

    for (ii = 0; ii < NUMNP; ii += 50)
    {
        printf(O);
        printf("%s\n\n", HEAD);
        printf("     Nodal Displacements     Time %13.5lf\n\n         Node", TIME);

        for (j = 0; j < NDM; j++)    printf("%7d %s",  j+1, CD);
        for (j = 0; j < NDF; j++)    printf("%14d %s", j+1, DI[0]);
        printf("\n");

        jj = min(NUMNP, ii+49);

        for (n = ii; n < jj; n++)
        {                                    
            if ( ! PCOMP(Get(X, n, 0, NDM), BLANK))
            {
                for (i = 0; i < NDF; i++)
                {
                    UL[i] = Get(F, n , i, NDF) * PROP;
                    k = abs(Get(ID, n, i, NDF));
                    if (k > 0) UL[i] = B[k];
                }
               
                printf("%13d", n+1);
                for (i = 0; i < NDM; i++)    printf("%20.4lf", Get(X, n, i, NDM));
                for (i = 0; i < NDF; i++)    printf("%20.4lf", UL[i]);
                printf("\n");
            }
        }
    }

    traceOUT("PRTDIS");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Print Nodal Reactions                                                     */
/*                                                                           */
/*****************************************************************************/

void PRTREA(double *R, int NDF)
{
    int k = 0, n = 0, i = 0, j = 0;
    double RSUM[6] = {0.0};
    double ASUM[6] = {0.0};

    traceIN("PRTREA");

    for (k = 0; k < NDF; k++)
    {
        RSUM[k] = 0.0;
        ASUM[k] = 0.0;
    }

    for (n = 0; n < NUMNP; n += 50)
    {
        j = min(NUMNP, n + 50);

        printf(O);
        printf("%s\n\n     Nodal Reactions\n\n      Node", HEAD);

        for (k = 0; k < NDF; k++)
            printf("%9d DOF", k+1);

        for (i = n; i < j; i++)
        {
            for (k = 0; k < NDF; k++)
            {
                Set(R, i, k, NDF, -Get(R, i, k, NDF));
                RSUM[k] = RSUM[k] + Get(R, i, k, NDF);
                ASUM[k] = ASUM[k] + fabs(Get(R, i, k, NDF));
            }

            printf("%10d", i+1);

            for (k = 0; k < NDF; k++)
                printf("%13.4lf", Get(R, i, k, NDF));
        }
    }

    /*-----------------------------------------------------------------------*/
    /* Print Statics Check                                                   */
    /*-----------------------------------------------------------------------*/

    for (k = 0; k < NDF; n++)
        printf("       Sum %13.4\n", RSUM[k]);

    for (k = 0; k < NDF; n++)
        printf("       Abs sum %13.4\n", ASUM[k]);


    traceOUT("PRTREA");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Set Pointer For Arrays                                                    */
/*                                                                           */
/*****************************************************************************/

void PSETM(int *NA, int *NE, int NJ, int *AFL)
{
    traceIN("PSETM");

    *NA  = *NE;
    *NE  = *NE + NJ;
    *AFL = FALSE;
    SETMEM(*NE);

    traceOUT("PSETM");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Zero Real Array                                                           */
/*                                                                           */
/*****************************************************************************/

void PZERO(double *V, int NN)
{
    int n = 0;

    traceIN("PZERO");

    for (n = 0; n < NN; n++)
        V[n] = 0.0;

    traceOUT("PZERO");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Unsymmetric, Active Column Profile Equation Solver                        */
/*                                                                           */
/*****************************************************************************/

void UACTCL(double *A, double *C, double *B, int *JDIAG, int NEQ, int AFAC, int BACK)
{
    int i = 0, j = 0, k = 0;
    int id = 0, ie = 0, ih = 0, ir = 0, is = 0, jd = 0, jh = 0, jr = 0;

    traceIN("UACTCOL");

    /*-----------------------------------------------------------------------*/
    /* Factor A to UT*D*U, reduce B to Y                                     */
    /*-----------------------------------------------------------------------*/

    jr = 0;

    for (j = 0; j < NEQ; j++)
    {
        jd = JDIAG[j];
        jh = jd - jr;
       
        if (jh > 1)
        {
            is = j + 1 - jh;
            ie = j - 1;
           
            if (AFAC)
            {
                k = jr + 1;
                id = 0;
              
                /*-----------------------------------------------------------*/
                /* Reduce all equations except diagonal                      */
                /*-----------------------------------------------------------*/
              
                for (i = is; i <= ie; i++)
                {
                    ir = id;
                    id = JDIAG[i];
                    ih = min(id - ir - 1, i - is);
              
                    if (ih != 0)
                    {
                        A[k] = A[k] - DOT(&A[k-ih], &C[id-ih], ih);
                        C[k] = C[k] - DOT(&C[k-ih], &A[id-ih], ih);
              
                        if (A[id] != 0.0) C[k] = C[k] / A[id];
                    }
              
                    k = k + 1;
                }
              
                /*-----------------------------------------------------------*/
                /* Reduce diagonal term                                      */
                /*-----------------------------------------------------------*/
              
                  A[jd] = A[jd] - DOT(&A[jr+1], &C[jr+1], jh-1);
            }
           
            /*---------------------------------------------------------------*/
            /* Forward reduce the R.H.S                                      */
            /*---------------------------------------------------------------*/
           
            if (BACK) B[j] = B[j] - DOT(&C[jr+1], &B[is], jh-1);
        }
       
        jr = jd;
    }

    if (BACK)
    {
        double D = 0.0;

        /*-------------------------------------------------------------------*/
        /* Back substitution                                                 */
        /*-------------------------------------------------------------------*/
       
        j  = NEQ - 1;
        jd = JDIAG[j];
       
        while (j >= 0)
        {
            if (A[jd] != 0.0)   B[j] = B[j] / A[jd];
       
            D = B[j];
            j = j - 1;
           
            if (j >= 0)
            {
                jr = JDIAG[j];
               
                if (jd - jr > 1)
                {
                    is = j - jd + jr + 2;
                    k = jr - is + 1;
               
                    for (i = is; i <= j; i++)
                        B[i] = B[i] - A[i+k] * D;
               
                    jd = jr;
                }
            }
        }
    }

    traceOUT("UACTCOL");
}

/******************************************************************************
*******************************************************************************
*******************************************************************************
*                                                                             *
C ELEMENT MODULES
*                                                                             *
*******************************************************************************
*******************************************************************************
******************************************************************************/
 



/*****************************************************************************/
/*                                                                           */
/* Plane Linear Elastic Element Routine                                      */
/*                                                                           */
/*****************************************************************************/

void ELMT01(double *D, double *UL, double *XL, int *IX, double *TL, double *S,
            double *P, int NDF, int NDM, int NST, int ISW)
{
    double SHP[9][3] = {{0.0}};

    static int    LINT = 0;
    static double SG[9] = {0.0}, TG[9] = {0.0}, WG[9] = {0.0};

    traceIN("ELMT01");

    /*-----------------------------------------------------------------------*/
    /* Go to correct array processor                                         */
    /*-----------------------------------------------------------------------*/

    switch (ISW)
    {
    /*-----------------------------------------------------------------------*/
    /* Input material properties                                             */
    /*-----------------------------------------------------------------------*/

    case 1:
        {
            int    I = 0, K = 0, L = 0;
            double E = 0.0, XNU = 0.0;
            char  *WD = NULL;

            char line[LINE_LENGTH + 1] = "";
            gets(line);
            E    = read_double(line,  0, 10);
            XNU  = read_double(line, 10, 20);
            D[3] = read_double(line, 20, 30);
            L    = read_int   (line, 30, 35);
            K    = read_int   (line, 35, 40);
            I    = read_int   (line, 35, 40);

            if (I == 0)
                { I = 2;    WD = "Plane Strain"; }
            else
                { I = 1;    WD = "Plane Stress"; }

            D[0] = E * (1.0 + (1.0 - I) * XNU) / (1.0 + XNU) / (1.0 - I * XNU);
            D[1] = XNU * D[0] / (1.0 + (1.0 - I) * XNU);
            D[2] = 0.5 * E / (1.0 + XNU);
            L = min(3, max(1, L));
            D[4] = L;
            K = min(3, max(1, K));
            D[5] = K;
            LINT = 0;

            printf("     %s Linear Elastic Element\n",     WD);
            printf("\n");
            printf("          Modulus %18.5lf\n",          E);
            printf("          Poisson Ratio     %8.5lf\n", XNU);
            printf("          Density %18.5lf\n",          D[3]);
            printf("          Gauss pts/dir %6d\n",        L);
            printf("          Stress pts    %6ld\n",       K);
            printf("\n");
            printf("          D[%d] = {%lf, %lf, %lf, %lf, %lf, %lf}\n\n",
                              6, D[0], D[1], D[2], D[3], D[4], D[5]);
            break;
        }

    /*-----------------------------------------------------------------------*/
    /* Compute Jacobiabs, etc. to check an element for errors                */
    /*-----------------------------------------------------------------------*/

    case 2:
        break;

    /*-----------------------------------------------------------------------*/
    /* Compute element (tangent) stiffness matrix                            */
    /*-----------------------------------------------------------------------*/

    case 3:
        {
            int l = 0, j = 0, k = 0, NSL = 0;
    
            /*----------------------------*/
    
            print1Darray("D",  D,  "%lf",     6);
            print2Darray("UL", UL, "%12.2lf", NEN, NDF);
            print2Darray("XL", XL, "%12.2lf", NEN, NDM);
            print1Darray("IX", IX, "%d",      NEN+1);
            print1Darray("TL", TL, "%lf",     NEN);
    
            /*----------------------------*/
    
            l = (int) D[4];
            if (l*l != LINT)     PGAUSS(l, &LINT, SG, TG, WG);
    
            /*---------------------------------------------------------------*/
            /* Fast stiffness computation,  compute integrals of shape       */
            /* functions                                                     */
            /*---------------------------------------------------------------*/
    
            for (l = 0; l < LINT; l++)
            {
                int jj = 0;
                double XSJ = *TL;                     /* dummy use of TL !!! */
    
                SHAPE(SG[l], TG[l], XL, SHP[0], &XSJ, NDM, NEL, IX, FALSE);
    
                XSJ = XSJ * WG[l];
               
                /*-----------------------------------------------------------*/
                /* Loop over rows                                            */
                /*-----------------------------------------------------------*/
               
                jj = 0;
               
                for (j = 0; j < NEL; j++)
                {
                    int k1 = 0;
                    double W11 = SHP[j][0] * XSJ;
                    double W12 = SHP[j][1] * XSJ;
                   
                    /*-------------------------------------------------------*/
                    /* Loop over columns (symmetry noted)                    */
                    /*-------------------------------------------------------*/
                   
                    k1 = jj;
                   
                    for (k = j; k < NEL; k++)
                    {
                        Set(S, k1  , jj  , NST, Get(S, k1  , jj  , NST) + W11 * SHP[k][0]);
                        Set(S, k1+1, jj  , NST, Get(S, k1+1, jj  , NST) + W11 * SHP[k][1]);
                        Set(S, k1  , jj+1, NST, Get(S, k1  , jj+1, NST) + W12 * SHP[k][0]);
                        Set(S, k1+1, jj+1, NST, Get(S, k1+1, jj+1, NST) + W12 * SHP[k][1]);
                        k1 = k1 + NDF;
                    }
                   
                    jj = jj + NDF;
                }
            }
    
            print2Darray("S", S, "%12.2lf", NST, NST);
    
            /*---------------------------------------------------------------*/
            /* Assemble the stiffness matrix from integrals and material     */
            /* props.                                                        */
            /*---------------------------------------------------------------*/
    
            NSL = NEL * NDF;
    
            printf("    NSL = %d, NEL = %d, NDF = %d, NST = %d, LINT = %d\n", NSL, NEL, NDF, NST, LINT);
    
            for (j = 0; j < NSL; j += NDF)
            {
                for (k = j; k < NSL; k += NDF)
                {
                    double W11 = Get(S, k,   j,   NST);
                    double W12 = Get(S, k+1, j,   NST);
                    double W21 = Get(S, k,   j+1, NST);
                    double W22 = Get(S, k+1, j+1, NST);
    
                    Set(S, k,   j,   NST, D[0]*W11 + D[2]*W22);
                    Set(S, k+1, j,   NST, D[1]*W12 + D[2]*W21);
                    Set(S, k,   j+1, NST, D[1]*W21 + D[2]*W12);
                    Set(S, k+1, j+1, NST, D[0]*W22 + D[2]*W11);
    
                    /*-------------------------------------------------------*/
                    /* Form lower part by symmetry                           */
                    /*-------------------------------------------------------*/
                   
                    Set(S, j,   k,   NST, Get(S, k,   j,   NST));
                    Set(S, j+1, k,   NST, Get(S, k,   j+1, NST));
                    Set(S, j,   k+1, NST, Get(S, k+1, j,   NST));
                    Set(S, j+1, k+1, NST, Get(S, k+1, j+1, NST));
                }
            }
    
            /*----------------------------*/
    
            print2Darray("S", S, "%12.2lf", NST, NST);
    
            /*----------------------------*/
    
            break;
        }

    /*-----------------------------------------------------------------------*/
    /* Compute and output element quantities (eg. stresses)                  */
    /*-----------------------------------------------------------------------*/

    case 4:
    case 6:
        {
            int l = 0, j = 0;

            l = (int) D[4];
            if (ISW == 4) l = (int) D[5];
    
            if (l*l != LINT)     PGAUSS(l, &LINT, SG, TG, WG);
    
            /*---------------------------------------------------------------*/
            /* Compute element stresses, strains, and forces                 */
            /*---------------------------------------------------------------*/
    
            for (l = 0; l < LINT; l++)
            {
                int    i = 0;
                double EPS[3] = {0.0};
                double SIG[6] = {0.0};
                double XX = 0.0, YY = 0.0;
                double XSJ = *TL;                     /* dummy use of TL !!! */
    
                /*-----------------------------------------------------------*/
                /* Compute element shape functions                           */
                /*-----------------------------------------------------------*/
               
                SHAPE(SG[l], TG[l], XL, SHP[0], &XSJ, NDM, NEL, IX, FALSE);
               
                /*-----------------------------------------------------------*/
                /* Compute strains and coordinates                           */
                /*-----------------------------------------------------------*/
               
                for (i = 0; i < 3; i++)
                    EPS[i] = 0.0;
    
                XX = 0.0;
                YY = 0.0;
    
                for (j = 0; j < NEL; j++)
                {
                    XX = XX + SHP[j][2] * Get(XL, j, 0, NDM);
                    YY = YY + SHP[j][2] * Get(XL, j, 1, NDM);
    
                    EPS[0] = EPS[0] + SHP[j][0] * Get(UL, j, 0, NDF);
                    EPS[2] = EPS[2] + SHP[j][1] * Get(UL, j, 1, NDF);
                    EPS[1] = EPS[1] + SHP[j][0] * Get(UL, j, 1, NDF) + SHP[j][1] * Get(UL, j, 0, NDF);
                }
    
                /*-----------------------------------------------------------*/
                /* Compute stresses                                          */
                /*-----------------------------------------------------------*/
               
                SIG[0] = D[0] * EPS[0] + D[1] * EPS[2];
                SIG[2] = D[1] * EPS[0] + D[0] * EPS[2];
                SIG[1] = D[2] * EPS[1];
    
                if (ISW != 6)
                {
                    PSTRES(SIG, &SIG[3], &SIG[4], &SIG[5]);
                   
                    /*-------------------------------------------------------*/
                    /* Output stresses and strains                           */
                    /*-------------------------------------------------------*/
                   
                    MCT = MCT - 2;
    
                    if (MCT <= 0)
                    {
                        printf(O);
                        printf("%s\n\n     Element stresses\n\n", HEAD);
                        printf("   Element  material      1-Coord      2-Coord    11-Stress    12-Stress    22-Stress    1-Stress    2-Stress   Angle\n");
                        printf("                                                  11-Strain    12-Strain    22-Strain\n");
    
                        MCT = 50;
                    }
    
                    printf("%10d %10d %13.4lf %13.4lf %13.4lf %13.4lf %13.4lf %13.4lf %13.4lf %8.2lf\n", N, MA, XX, YY, SIG[0], SIG[1], SIG[2], SIG[3], SIG[4], SIG[5]);
                    printf("                                              %13.4lf %13.4lf %13.4lf\n", EPS[0], EPS[1], EPS[2]);
                }
                else
                {
                    /*-------------------------------------------------------*/
                    /* Compute internal forces                               */
                    /*-------------------------------------------------------*/
                   
                    double DV = XSJ * WG[l];
                    int jj = 1;
    
                    for (j = 0; j < NEL; j++)
                    {
                        P[jj]   = P[jj]   - (SHP[j][0] * SIG[0] + SHP[j][1] * SIG[1]) * DV;
                        P[jj+1] = P[jj+1] - (SHP[j][0] * SIG[1] + SHP[j][1] * SIG[2]) * DV;
                        jj = jj + NDF;
                    }
                }
            }
            break;
        }

    /*-----------------------------------------------------------------------*/
    /* Compute constistent mass matrix                                       */
    /*-----------------------------------------------------------------------*/

    case 5:
        {
            int l = 0, j = 0, k = 0, NSL = 0;

            l = (int)  D[4];
       
            if (l*l != LINT)      PGAUSS(l, &LINT, SG, TG, WG);
       
            for (l = 0; l < LINT; l++)
            {
                double DV = 0.0;
                double XSJ = *TL;                     /* dummy use of TL !!! */
                int    jj = 0;
       
                /*-----------------------------------------------------------*/
                /* Compute shape functions                                   */
                /*-----------------------------------------------------------*/
               
                SHAPE(SG[l], TG[l], XL, SHP[0], &XSJ, NDM, NEL, IX, FALSE);
                DV = WG[l] * XSJ * D[3];
               
                /*-----------------------------------------------------------*/
                /* For each node J compute DB = RHO*SHAPE*DV                 */
                /*-----------------------------------------------------------*/
               
                jj = 0;
               
                for (j = 0; j < NEL; j++)
                {
                    double W11 = SHP[j][2] * DV;
                 
                    /*-------------------------------------------------------*/
                    /* For each node K compute mass matrix (upper triangular */
                    /* part)                                                 */
                    /*-------------------------------------------------------*/
                 
                    int k1 = jj;
                 
                    for (k = j; k < NEL; k++)
                    {
                        Set(S, k1, jj, NST, Get(S, k1, jj, NST) + SHP[j][2] * W11);
                        k1 = k1 + NDF;
                    }
                 
                    jj = jj + NDF;
                }
            }
       
            /*---------------------------------------------------------------*/
            /* Compute missing parts and lower part by symmetries            */
            /*---------------------------------------------------------------*/
       
            NSL = NEL * NDF;
       
            for (j = 0; j < NSL; j+= NDF)
            {
                for (k = j; k < NSL; k += NDF)
                {
                    Set(S, k+1, j+1, NST, Get(S, k, j, NST));
                    Set(S, j,   k,   NST, Get(S, k, j, NST));
                    Set(S, j+1, k+1, NST, Get(S, k, j, NST));
                }
            }
            break;
        }
    }

    traceOUT("ELMT01");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Gauss points and weights for two dimensions                               */
/*                                                                           */
/*****************************************************************************/

void PGAUSS(int l, int *LINT, double *R, double *Z, double *W)
{
    int LR[9] = {-1,  1,  1, -1,  0,  1,  0, -1,  0};
    int LZ[9] = {-1, -1,  1,  1, -1,  0,  1,  0,  0};
    int LW[9] = {25, 25, 25, 25, 40, 40, 40, 40, 64};
    int i = 0;
    double G = 0.0, H = 0.0;

    traceIN("PGAUSS: l = %d", l);

    *LINT = l * l;

    switch (l)
    {
    /*-----------------------------------------------------------------------*/
    /* 1X1 Integration                                                       */
    /*-----------------------------------------------------------------------*/

    case 1:
        R[0] = 0.0;
        Z[0] = 0.0;
        W[0] = 4.0;
        break;

    /*-----------------------------------------------------------------------*/
    /* 2X2 Integration                                                       */
    /*-----------------------------------------------------------------------*/

    case 2:
        G = 1.0 / sqrt(3.0);

        for (i = 0; i < 4; i++)
        {
            R[i] = G * LR[i];
            Z[i] = G * LZ[i];
            W[i] = 1.0;
        }
        break;

    /*-----------------------------------------------------------------------*/
    /* 3X3 Integration                                                       */
    /*-----------------------------------------------------------------------*/

    case 3:
        G = sqrt(0.6);
        H = 1.0 / 81.0;

        for (i = 0; i < 9; i++)
        {
            R[i] = G * LR[i];
            Z[i] = G * LZ[i];
            W[i] = H * LW[i];
        }
        break;
    }

    print1Darray("R",  R,  "%lf", *LINT);
    print1Darray("Z",  Z,  "%lf", *LINT);
    print1Darray("W",  W,  "%lf", *LINT);

    traceOUT("PGAUSS");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Compute principle stresses (2 dimensioons)                                */
/*                                                                           */
/*****************************************************************************/

void PSTRES(double *SIG, double *P1, double *P2, double *P3)
{
    /*-----------------------------------------------------------------------*/
    /* Stresses must be stored in array SIG[3] in the order                  */
    /* tau-xx, tau-XY, tau-YY                                                */
    /*-----------------------------------------------------------------------*/
 
    double XI1 = (SIG[0] + SIG[2]) / 2.0;
    double XI2 = (SIG[0] - SIG[2]) / 2.0;
    double RHO = sqrt(XI2*XI2 + SIG[1]*SIG[1]);

    traceIN("PSTRES");

    *P1 = XI1 + RHO;
    *P2 = XI1 - RHO;
    *P3 = 45.0;
                                         
    if (XI2 != 0.0)   *P3 = 22.5 * atan2(SIG[1], XI2) / atan(1.0);

    traceOUT("PSTRES");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Shape function routine for two dimensions                                 */
/*                                                                           */
/*****************************************************************************/

void SHAPE(double SS, double TT, double *X, double *SHP, double *XSJ, 
           int NDM, int NEL, int *IX, int FLG)
{
    double S[4] = {-0.5,  0.5,  0.5, -0.5};
    double T[4] = {-0.5, -0.5,  0.5,  0.5};
    double XS[2][2] = {{0.0}};
    double SX[2][2] = {{0.0}};
    int i = 0, j = 0;

    traceIN("SHAPE: SS = %lf, TT = %lf, NDM = %d, NEL = %d, FLG = %d", SS, TT, NDM, NEL, FLG);

    /*-----------------------------------------------------------------------*/
    /* Form 4-node quadrilateral shape functions                             */
    /*-----------------------------------------------------------------------*/

    for (i = 0; i < 4; i++)
    {
        Set(SHP, i, 2, 3, (0.5 + S[i] * SS) * (0.5 + T[i] * TT));
        Set(SHP, i, 0, 3, S[i] * (0.5 + T[i] * TT));
        Set(SHP, i, 1, 3, T[i] * (0.5 + S[i] * SS));
    }

    if (NEL < 4)
    {
        /*-------------------------------------------------------------------*/
        /* Form triangle by adding third and fourth together                 */
        /*-------------------------------------------------------------------*/
       
        for (i = 0; i < 3; i++)
        {
            Set(SHP, 2, i, 3, Get(SHP, 2, i, 3) + Get(SHP, 3, i, 3));
        }
    }

    /*-----------------------------------------------------------------------*/
    /* Add quadratic terms if necessary                                      */
    /*-----------------------------------------------------------------------*/

    if (NEL > 4)     SHAP2(SS, TT, SHP, IX, NEL);
    
    /*-----------------------------------------------------------------------*/
    /* Construct jacobian and its inverse                                    */
    /*-----------------------------------------------------------------------*/

    for (i = 0; i < NDM; i++)
    {
        for (j = 0; j < 2; j++)
        {
            int k = 0;

            XS[j][i] = 0.0;

            for (k = 0; k < NEL; k++)
            {
                XS[j][i] = XS[j][i] + Get(X, k, i, NDM) * Get(SHP, k, j, 3);
            }
        }
    }

    *XSJ = XS[0][0] * XS[1][1] - XS[0][1] * XS[1][0];

    print2Darray("SHP", SHP, "%12.6lf", NEL, 3);
    print2Darray("XS", XS, "%12.6lf", 2, 2);
    printf("*XSJ = %lf\n", *XSJ);

    if (!FLG)
    {
        SX[0][0] =   XS[1][1] / (*XSJ);
        SX[1][1] =   XS[0][0] / (*XSJ);
        SX[1][0] = - XS[1][0] / (*XSJ);
        SX[0][1] = - XS[0][1] / (*XSJ);
                                 
        /*-------------------------------------------------------------------*/
        /* Form global derivatives                                           */
        /*-------------------------------------------------------------------*/
       
        for (i = 0; i < NEL; i++)
        {
            double TP       = Get(SHP, i, 0, 3) * SX[0][0] + Get(SHP, i, 1, 3) * SX[0][1];
            Set(SHP, i, 1, 3, Get(SHP, i, 0, 3) * SX[1][0] + Get(SHP, i, 1, 3) * SX[1][1]);
            Set(SHP, i, 0, 3, TP);
        }

        print2Darray("SX", SX, "%12.6lf", 2, 2);
        print2Darray("SHP", SHP, "%12.6lf", NEL, 3);
    }

    traceOUT("SHAPE");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Add quadratic functions as necessary                                      */
/*                                                                           */
/*****************************************************************************/

void SHAP2(double S, double T, double *SHP, int *IX, int NEL)
{
    double S2 = (1.0 - S * S) / 2.0;
    double T2 = (1.0 - T * T) / 2.0;
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;

    traceIN("SHAP2");

    for (i = 4; i < 9; i++)
    {
        for (j = 0; j < 3; j++)
        {
            Set(SHP, i, j, 3, 0.0);
        }
    }

    /*-----------------------------------------------------------------------*/
    /* Midside nodes (serendipty)                                            */
    /*-----------------------------------------------------------------------*/

    if (IX[4] != 0)
    {
        Set(SHP, 4, 0, 3, -S * (1.0 - T));
        Set(SHP, 4, 1, 3, -S2);
        Set(SHP, 4, 2, 3, S2 * (1.0 - T));
    }

    if (NEL >= 6)
    {
        if (IX[5] != 0)
        {
            Set(SHP, 5, 0, 3, T2);
            Set(SHP, 5, 1, 3, -T * (1.0 + S));
            Set(SHP, 5, 2, 3, T2 * (1.0 + S));
        }
       
        if (NEL >= 7)
        {
            if (IX[6] != 0)
            {
                Set(SHP, 6, 0, 3, -S * (1.0 + T));
                Set(SHP, 6, 1, 3, S2);
                Set(SHP, 6, 2, 3, S2 * (1.0 + T));
            }
           
            if (NEL >= 8)
            {
                  if (IX[7] != 0)
                  {
                      Set(SHP, 7, 0, 3, -T2);
                      Set(SHP, 7, 1, 3, -T * (1.0 - S));
                      Set(SHP, 7, 2, 3, T2 * (1.0 - S));
                  }
               
                /*-----------------------------------------------------------*/
                /* Interior node (lagrangian)                                */
                /*-----------------------------------------------------------*/
               
                if ((NEL >= 9) && (IX[8] != 0))
                {
                    Set(SHP, 8, 0, 3, -4.0 * S * T2);
                    Set(SHP, 8, 1, 3, -4.0 * T * S2);
                    Set(SHP, 8, 2, 3,  4.0 * S2 * T2);
                   
                    /*-------------------------------------------------------*/
                    /* Correct edge nodes for interior node (lagrangian)     */
                    /*-------------------------------------------------------*/
                   
                    for (j = 0; j < 3; j++)
                    {
                        for (i = 0; j < 4; j++)
                        {
                            Set(SHP, i, j, 3, Get(SHP, i, j, 3) - 0.25 * Get(SHP, 8, j, 3));
                        }
                       
                        for (i = 4; i < 8; i++)
                        {
                            if (IX[i] != 0)   Set(SHP, i, j, 3, Get(SHP, i, j, 3) - 0.5 * Get(SHP, 8, j, 3));
                        }
                    }
                }
            }
        }
    }

    /*-----------------------------------------------------------------------*/
    /* Correct corner nodes for presence of midside nodes                    */
    /*-----------------------------------------------------------------------*/

    k = 8;

    for (i = 0; i < 4; i++)
    {
        l = i + 4;

        for (j = 0; j < 3; j++)
        {
            Set(SHP, i, j, 3, Get(SHP, i, j, 3) - 0.5 * (Get(SHP, k, j, 3) + Get(SHP, l, j, 3)));
        }

        k = l;
    }

    traceOUT("SHAP2");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Two dimensional heat transfer element                                     */
/*                                                                           */
/*****************************************************************************/

void ELMT02(double *D, double *UL, double *XL, int *IX, double *TL,
            double *S, double *P, int NDF, int NDM, int NST, int ISW)
{
    double XSJ = *TL;
    int i = NDF, j = 0, l = 0;
    double SG[4] = { 1.0 ,  1.0 , -1.0 , -1.0};
    double TG[4] = {-1.0 ,  1.0 ,  1.0 , -1.0};
    double G = 1.0 / sqrt(3.0);
    double SHP[9][3] = {{0.0}};

    static int    KAT = 0;

    traceIN("ELMT02");

    /*-----------------------------------------------------------------------*/
    /* Transfer to correct processor                                         */
    /*-----------------------------------------------------------------------*/

    switch (ISW)
    {

    /*-----------------------------------------------------------------------*/
    /* Input material properties                                             */
    /*-----------------------------------------------------------------------*/
 
    case 1:
        {
            char *WLAB = NULL;
            char line[LINE_LENGTH + 1] = "";

            gets(line);
            D[0] = read_double(line,  0, 10);
            D[1] = read_double(line, 10, 20);
            D[2] = read_double(line, 20, 30);
            KAT  = read_int   (line, 30, 35);
        
            printf("    Linear heat conduction element\n\n");
            printf("     Conductivity %12.5lf     Specific heat %12.5lf       Density        %12.5lf\n", D[0], D[1], D[2]);
        
            G = 1.0 / sqrt(3.0);
            D[1] = D[1] * D[2];

            if (KAT != 2)
            {
                KAT = 1;
                WLAB = "Plane";
            }
            else
            {
                KAT = 2;
                WLAB = "Axisym";
            }

            printf("          %s Analysis\n", WLAB);

            break;
        }

    /*-----------------------------------------------------------------------*/
    /* Insert check of mesh if desired                                       */
    /*-----------------------------------------------------------------------*/

    case 2:
    case 4:
        break;

    /*-----------------------------------------------------------------------*/
    /* Compute conductivity (stiffness) matrix                               */
    /*-----------------------------------------------------------------------*/

    case 3:
    case 6:
        for (l = 0; l < 4; l++)
        {
            SHAPE(SG[l] * G, TG[l] * G, XL, SHP[0], &XSJ, NDM, NEL, IX, FALSE);

            if (KAT == 2)
            {
                double RR = 0.0;

                for (i = 0; i < NEL; i++)
                    RR = RR + SHP[i][2] * Get(XL, i, 0, NDM);

                XSJ = XSJ * RR;
            }

            for (j = 0; j < NEL; j++)
            {
                double A1  = D[0] * SHP[j][0] * XSJ;
                double A2  = D[0] * SHP[j][1] * XSJ;
               
                for (i = 0; i < NEL; i++)
                {
                    Set(S, j, i, NST, Get(S, j, i, NST) + A1 * SHP[i][0] + A2 * SHP[i][1]);
                }
            }
        }

        for (i = 0; i < NEL; i++)
        {
            for (j = 0; j < NEL; j++)
            {
                P[i] = P[i] - Get(S, j, i, NST) * UL[j];
            }
        }

        break;

    /*-----------------------------------------------------------------------*/
    /* Compute heat capacity (mass) matrix                                   */
    /*-----------------------------------------------------------------------*/

    case 5:
        for (l = 0; l < 4; l++)
        {
            SHAPE(SG[l] * G, TG[l] * G, XL, SHP[0], &XSJ, NDM, NEL, IX, FALSE);

            if (KAT == 2)
            {
                double RR = 0.0;
               
                for (i = 0; i < NEL; i++)
                    RR = RR + SHP[i][2] * Get(XL, i, 0, NDM);
               
                XSJ = XSJ * RR;
            }

            for (j = 0; j < NEL; j++)
            {
                double SHJ = D[1] * SHP[j][2] * XSJ;
                P[j] = P[j] + SHJ;
               
                for (i = 0; i < NEL; i++)
                {
                    Set(S, j, i, NST, Get(S, j, i, NST) + SHJ * SHP[i][2]);
                }
            }
        }

        break;
    }

    traceOUT("ELMT02");
    return;
}

/*****************************************************************************/
/*                                                                           */
/* Two dimensional fluid element for navier-stokes equations                 */
/*                                                                           */
/*****************************************************************************/

void ELMT03(double *D, double *UL, double *XL, int *IX, double *TL,
            double *S, double *P, int NDF, int NDM, int NST, int ISW)
{
    double XSJ  = *TL, XLAM = 0.0, XMU = 0.0, XRHO = 0.0;
    double V[2] = {0.0};
    int    l = 0, i = 0, j = 0, k = 0;
    double DV[2][2] = {{0.0}};
    double SHP[9][3] = {{0.0}};

    static int    LINT = 0;
    static double SG[9] = {0.0}, TG[9] = {0.0}, WG[9] = {0.0};

    traceIN("ELMT03");

    /*-----------------------------------------------------------------------*/
    /* Transfer to correct processor                                         */
    /*-----------------------------------------------------------------------*/

    switch (ISW)
    {

    /*-----------------------------------------------------------------------*/
    /* Input/output fluid properties                                         */
    /*-----------------------------------------------------------------------*/

    case 1:
        {
            int L = 0;

            char line[LINE_LENGTH + 1] = "";
           
            gets(line);
            D[0] = read_double(line,  0, 10);
            D[1] = read_double(line, 10, 20);
            D[2] = read_double(line, 20, 30);
            L    = read_int   (line, 30, 35);
            
           
            printf("    Two dimensional fluid element\n\n");
            printf("          Viscosity  = %12.5lf\n", D[0]);
            printf("          Constraint = %12.5lf\n", D[1]);
            printf("          Density    = %12.5lf\n", D[2]);
            printf("          Gauss pt/dir = %5d\n", L);
           
            D[3] = L;
            LINT = 0;
            break;
        }

    /*-----------------------------------------------------------------------*/
    /* Insert check of mesh if desired                                       */
    /*-----------------------------------------------------------------------*/

    case 2:
        break;

    /*-----------------------------------------------------------------------*/
    /* Compute unsymmetric tangent stiffness or out of balance forces        */
    /*-----------------------------------------------------------------------*/

    case 3:
    case 6:
        l = (int) D[3];

        if (l*l != LINT)    PGAUSS(l, &LINT, SG, TG, WG);

        for (l = 0; l < LINT; l++)
        {
            SHAPE(SG[l], TG[l], XL, SHP[0], &XSJ, NDM, NEL, IX, FALSE);
            XLAM = D[1] * XSJ * WG[l];
            XMU  = D[0] * XSJ * WG[l];
            XRHO = D[2] * XSJ * WG[l];
           
            /*---------------------------------------------------------------*/
            /* Compute velocities and gradients                              */
            /*---------------------------------------------------------------*/
           
            for (i = 0; i < 2; i++)
            {
                V[i] = 0.0;
               
                for (k = 0; k < NEL; k++)
                    V[i] = V[i] + SHP[k][2] * Get(UL, k, i, NDF);
               
                for (j = 0; j < 2; j++)
                {
                    DV[j][i] = 0.0;
                   
                    for (k = 0; k < NEL; k++)
                    {
                        DV[j][i] = DV[j][i] + SHP[k][j] * Get(UL, k, i, NDF);
                    }
                }
            }
           
            if (ISW == 3)
            {
                /*-----------------------------------------------------------*/
                /* Compute tangent, loop over columns of S                   */
                /*-----------------------------------------------------------*/
               
                int k1 = 0;
                for (k = 0; k < NEL; k++)
                {
                    int jj = 0;

                    double A1 = XMU * SHP[k][0];
                    double A2 = XMU * SHP[k][1];
/*                  double A3 = XRHO * (DV[0][0] * SHP[k][2] + V[0] * SHP[k][0] + V[1] * SHP[k][1]); */
/*                  double A4 = XRHO * (DV[1][1] * SHP[k][2] + V[0] * SHP[k][0] + V[1] * SHP[k][1]); */
/*                  double A5 = XRHO * DV[1][0] * SHP[k][2];                                         */
/*                  double A6 = XRHO * DV[0][1] * SHP[k][2];                                         */
                    double B1 = XLAM * SHP[k][0];
                    double B2 = XLAM * SHP[k][1];
                   
                    /*-------------------------------------------------------*/
                    /* Loop over rows of S                                   */
                    /*-------------------------------------------------------*/
                   
                    jj = 0;
                   
                    for (j = 0; j < NEL; j++)
                    {
                        Set(S, k1,   jj,   NST, Get(S, k1,   jj,   NST) + SHP[j][0] * (A1+A1+B1) + SHP[j][1] * A2);
                        Set(S, k1+1, jj,   NST, Get(S, k1+1, jj,   NST) + SHP[j][0] * B2 + SHP[j][1] * A1);
                        Set(S, k1,   jj+1, NST, Get(S, k1,   jj+1, NST) + SHP[j][0] * A2 + SHP[j][2] * B1);
                        Set(S, k1+1, jj+1, NST, Get(S, k1+1, jj+1, NST) + SHP[j][0] * A1 + SHP[j][2] * (A2+A2+B2));
                        jj = jj + NDF;
                    }
                   
                    k1 = k1 + NDF;
                }
            }
            else
            {
                /*-----------------------------------------------------------*/
                /* Compute divergence term                                   */
                /*-----------------------------------------------------------*/
               
                double XDIV = (DV[0][0] + DV[1][1]) * XLAM;
               
                /*-----------------------------------------------------------*/
                /* Compute internal forces                                   */
                /*-----------------------------------------------------------*/
               
                for (k = 0; k < NEL; k++)
                {
                    for (j = 0; j < 2; j++)
                    {
                        double SUM = XDIV * SHP[k][j];
                       
                        for (i = 0; i < 2; i++)
                        {
                            SUM = SUM + XMU * (DV[i][j] + DV[j][i]) * SHP[k][i] + XRHO * V[i] * DV[i][j] * SHP[k][2];
                        }
                       
                        Set(P, k, j, NDF, Get(P, k, j, NDF) - SUM);
                    }
                }
            }
        }
        break;

    /*-----------------------------------------------------------------------*/
    /* Compute stresses and velocity gradients                               */
    /*-----------------------------------------------------------------------*/

    case 4:
        break;

    /*-----------------------------------------------------------------------*/
    /* Compute mass matrix                                                   */
    /*-----------------------------------------------------------------------*/

    case 5:
        break;
    }

    traceOUT("ELMT03");
}
 
/*****************************************************************************/
/*                                                                           */
/* Read an integer value from a string between margins                       */
/*                                                                           */
/*****************************************************************************/

int read_int(char *line, int lo, int hi)
{
    int value = 0, count = 0;

    if (strlen(line) <= lo)
        value = 0;
    else
    {
        char ch = line[hi];
        line[hi] = 0;
        count = sscanf(&line[lo], "%d", &value);
        line[hi] = ch;

        if (count != 1) value = 0;
    }

    return value;
}

/*****************************************************************************/
/*                                                                           */
/* Read a double value from a string between margins                         */
/*                                                                           */
/*****************************************************************************/

double read_double(char *line, int lo, int hi)
{
    double value = 0.0;
    int    count = 0;

    if (strlen(line) <= lo)
        value = BLANK;
    else
    {
        char ch = line[hi];
        line[hi] = 0;
        count = sscanf(&line[lo], "%lf", &value);
        line[hi] = ch;

        if (count != 1) value = BLANK;
    }

    return value;
}

/*****************************************************************************/
/*                                                                           */
/* Read a string from a string between margins                               */
/*                                                                           */
/*****************************************************************************/

char *read_string(char *line, int lo, int hi, char *buffer)
{
    size_t limit = (size_t) (hi - lo + 1);

    if (strlen(line) <= lo)
        *buffer = 0;
    else
        strncpy(buffer, &line[lo], limit);

    return buffer;
}

/*****************************************************************************/
/*                                                                           */
/* Left justify and strip trailing blanks from a string                      */
/*                                                                           */
/*****************************************************************************/

char *strip(char *line)
{
    char *p = line, *q = line, *end = line;

    while (*p == ' ')     p++;

    while (*p)
    {
        if ((*q++ = *p++) != ' ')     end = q;
    }

    *end = 0;

    return line;
}

/*****************************************************************************/
/*                                                                           */
/* Transfer of sign.  Return the sign of "b" times "absolute value of 'a'"   */
/*                                                                           */
/*****************************************************************************/

int ISIGN(int a, int b)
{         
    int sign = 1;

    if (b < 0)    sign = -1;
 
    return sign * abs(a);
}

/*****************************************************************************/
/*                                                                           */
/* Remaindering.  Return "a - [a/b]*b" where [x] is the integer whose        */
/* magnitude does not exceed the magnitude of "x" and whose sign is the same */
/* as "x"                                                                    */
/*                                                                           */
/*****************************************************************************/

double AMOD(double a, double b)
{
    double x     = a / b;
    double sign  = (x >= 0.0) ? 1.0 : -1.0;
    double thing = sign * floor(fabs(x));

    return (a - thing * b);
}

/*****************************************************************************/
/*                                                                           */
/* Print the matrix                                                          */
/*                                                                           */
/*****************************************************************************/

void dspMatrix(int n, int *JDIAG, double *A, double *B, char *format, ...)
{
    int i = 0;
    va_list arguments = NULL;

    printf("\n");

    va_start(arguments, format);
    vprintf(format, arguments);
    va_end(arguments);

    printf("\n");
    printf("%c%c",   esc, condensed);
    printf("%c%c",   esc, cpi_15);
    printf("\n");

    for (i = 0; i < n; i++)
    {
        int j = 0, jr = -1;

        for (j = 0; j < n; j++)
        {
            int jd = JDIAG[j] - 1;
            int jh = jd - jr;
            int k  = jh - (j - i);
            int l  = jr +  k;

            if ((j < i) || (k < 1))
                printf("       ");
            else
                printf("%7.1lf", A[l]);

            jr = jd;
        }

        printf("  : %7.3lf\n", B[i]);
    }

    printf("%c",   reset);
    printf("%c%c", esc, cpi_10);
    printf("\n");
}

/*****************************************************************************/
/*                                                                           */
/* Dump a 2 dimensional array                                                */
/*                                                                           */
/*****************************************************************************/

void print2Darray(char *name, void *pointer, char *format, int row, int col)
{
    char string[40] = "";
    int i = 0, j = 0, indent = 20, integer = 1;
    if (strchr(format, 'f')) integer = 0;

    sprintf(string, "%s[%d][%d] = {{", name, row, col);
    printf("%*.*s", indent, indent, string);
    print2Dvalue(integer, pointer, format, 0, 0, col);

    for (j = 1; j < col; j++)
    {
        printf(", ");
        print2Dvalue(integer, pointer, format, 0, j, col);
    }

    printf("}");

    for (i = 1; i < row; i++)
    {
        printf(",\n%*.*s{", indent-1, indent-1, "");
        print2Dvalue(integer, pointer, format, i, 0, col);

        for (j = 1; j < col; j++)
        {
            printf(", ");
            print2Dvalue(integer, pointer, format, i, j, col);
        }
        printf("}");
    }
    printf("}\n");
}

/*****************************************************************************/
/*                                                                           */
/* Dump a 1 dimensional array                                                */
/*                                                                           */
/*****************************************************************************/

void print1Darray(char *name, void *pointer, char *format, int size)
{
    char string[40]   = "";
    int i = 0, indent = 20, integer = 1;
    if (strchr(format, 'f')) integer = 0;

    sprintf(string, "%s[%d] = {", name, size);
    printf("%*.*s", indent, indent, string);
    print1Dvalue(integer, pointer, format, 0);

    for (i = 1; i < size; i++)
    {
        printf(", ");
        print1Dvalue(integer, pointer, format, i);
    }
    printf("}\n");
}

/*****************************************************************************/
/*                                                                           */
/* Print a value from a 2 dimensional array                                  */
/*                                                                           */
/*****************************************************************************/

void print2Dvalue(int integer, void *pointer, char *format, int i, int j, int col)
{
    if (integer)
    {
        int *array = (int*) pointer;
        int value  = Get(array, i, j, col);
        printf(format, value);
    }
    else
    {
        double *array = (double*) pointer;
        double value  = Get(array, i, j, col);
        printf(format, value);
    }
}

/*****************************************************************************/
/*                                                                           */
/* Print a value from a 1 dimensional array                                  */
/*                                                                           */
/*****************************************************************************/

void print1Dvalue(int integer, void *pointer, char *format, int i)
{
    if (integer)
    {
        int *array = (int*) pointer;
        int value  = array[i];
        printf(format, value);
    }
    else
    {
        double *array = (double*) pointer;
        double value  = array[i];
        printf(format, value);
    }
}

/*****************************************************************************/
/*                                                                           */
/* Trace into a function                                                     */
/*                                                                           */
/*****************************************************************************/

void traceIN(char *format, ...)
{
    va_list arguments = NULL;

    trace_prefix[trace_indent  ] = 'Ú';
    trace_prefix[trace_indent+1] = 0;
 
    printf("%s ", trace_prefix);
 
    va_start(arguments, format);
    vprintf(format, arguments);
    va_end(arguments);

    printf("\n");

    trace_prefix[trace_indent  ] = '³';
    trace_indent++;
}

/*****************************************************************************/
/*                                                                           */
/* Trace out of a function                                                   */
/*                                                                           */
/*****************************************************************************/

void traceOUT(char *format, ...)
{
    va_list arguments = NULL;

    if (trace_indent > 0)    trace_indent--;

    trace_prefix[trace_indent  ] = 'À';
    trace_prefix[trace_indent+1] = 0;
 
    printf("%s ", trace_prefix);
 
    va_start(arguments, format);
    vprintf(format, arguments);
    va_end(arguments);

    printf("\n");
}

