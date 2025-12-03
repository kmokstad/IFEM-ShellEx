/*!
  Simple program to generate a Nastran bulk data file for a 2D bicycle wheel.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main (int argc, const char** argv)
{
  double t = 0.1, E = 2.0e11, G = 8.0e9, nu = 0.3, rho = 5000.0;

  char useRBE3 = 'n';
  double r, theta, Ri = 1.0, Ro = 1.5;
  int i, j, k, nRef, nC = 10, nR = 2;

  const char* outputFile = NULL;

  for (i = 1; i < argc; i++)
    if      (!strcmp(argv[i],"-nC") && i+1 < argc)
      nC = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-nR") && i+1 < argc)
      nR = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-Ri") && i+1 < argc)
      Ri = atof(argv[++i]);
    else if (!strcmp(argv[i],"-Ro") && i+1 < argc)
      Ro = atof(argv[++i]);
    else if (!strcmp(argv[i],"-rbe3"))
      useRBE3 = 'y';
    else if (!outputFile)
      outputFile = argv[i];

  if (!outputFile)
  {
    fprintf(stderr,"Usage: %s [-nC <nC>] [-nR <nR>] [-Ri <Ri>] [-Ro <Ro>]"
            " [-rbe3] <Nastran-file>\n",argv[0]);
    return 1;
  }

  FILE* fd = fopen(outputFile,"w");
  if (!fd)
  {
    perror(outputFile);
    return 1;
  }

  printf("Generating bicycle wheel in domain [%g,%g]x2*pi",Ri,Ro);
  printf(" with %dx%d elements\n",nC,nR);

  fprintf(fd,"$* Bicycle wheel Ri=%g Ro=%g nC=%d nR=%d\n",Ri,Ro,nC,nR);
  fprintf(fd,"BEGIN BULK\n");
  for (j = k = 0, r = Ri; j <= nR; ++j, r += (Ro-Ri)/nR)
    for (i = 0, theta = 0.0; i < nC; ++i, theta += 2.0*M_PI/nC)
      fprintf(fd,"GRID,%d,,%g,%g\n",++k, r*cos(theta), r*sin(theta));
  fprintf(fd,"GRID,%d,,0.0,0.0\n",nRef=++k);
  for (j = k = 0; j < nR; j++)
  {
    for (i = 0; i+1 < nC; i++)
      fprintf(fd,"CQUAD4,%d,1,%d,%d,%d,%d\n",++k,
              1+i+nC*j,2+i+nC*j,2+i+nC*(1+j),1+i+nC*(1+j));
    fprintf(fd,"CQUAD4,%d,1,%d,%d,%d,%d\n",++k,
            nC*(1+j),1+nC*j,1+nC*(1+j),nC*(2+j));
  }
  fprintf(fd,"CONM2,%d,%d,,10.0\n",++k,nRef);
  if (useRBE3 == 'y')
    fprintf(fd,"RBE3,%d,,%d,123456,1.0,123",++k,nRef);
  else
    fprintf(fd,"RBE2,%d,%d,123456",++k,nRef);
  for (i = 1; i <= nC; i++) fprintf(fd,",%d",i);
  fprintf(fd,"\nPSHELL,1,1,%g\n",t);
  fprintf(fd,"MAT1,1,%g,%g,%g,%g\n",E,G,nu,rho);
  fprintf(fd,"END DATA\n");
  fclose(fd);

  return 0;
}
