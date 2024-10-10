/*!
  Simple program to generate a Nastran bulk data file for a rectangular plate.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main (int argc, const char** argv)
{
  float t = 0.1f, E = 2.0e11f, G = 8.0e9f, nu = 0.3f, rho = 5000.0f;

  float x, y, Lx = 1.0f, Ly = 1.0f;
  int i, j, k, nx = 10, ny = 10;

  const char* outputFile = NULL;

  for (i = 1; i < argc; i++)
    if      (!strcmp(argv[i],"-nx") && i+1 < argc)
      nx = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-ny") && i+1 < argc)
      ny = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-Lx") && i+1 < argc)
      Lx = atof(argv[++i]);
    else if (!strcmp(argv[i],"-Ly") && i+1 < argc)
      Ly = atof(argv[++i]);
    else if (!outputFile)
      outputFile = argv[i];

  if (!outputFile)
  {
    fprintf(stderr,"Usage: %s [-nx <nx>] [-ny <ny>] [-Lx <Lx>] [-Ly <Ly>]"
	    " <Nastran-file>\n",argv[0]);
    return 1;
  }

  FILE* fd = fopen(outputFile,"w");
  if (!fd)
  {
    perror(outputFile);
    return 1;
  }

  printf("Generating rectangular plate in domain [0,%g]x[0,%g]",Lx,Ly);
  printf(" with %dx%d elements\n",nx,ny);

  fprintf(fd,"$* Rectangular plate Lx=%g Ly=%g Nx=%d Ny=%d\n",Lx,Ly,nx,ny);
  fprintf(fd,"BEGIN BULK\n");
  for (j = k = 0, y = 0.0f; j <= ny; j++, y += Ly/(float)ny)
    for (i = 0, x = 0.0f; i <= nx; i++, x += Lx/(float)nx)
      fprintf(fd,"GRID,%d,,%g,%g\n",++k,x,y);
  for (j = k = 0; j < ny; j++)
    for (i = 0; i < nx; i++)
      fprintf(fd,"CQUAD4,%d,1,%d,%d,%d,%d\n",++k,
	      1+i+(nx+1)*j,2+i+(nx+1)*j,2+i+(nx+1)*(1+j),1+i+(nx+1)*(1+j));
  fprintf(fd,"PSHELL,1,1,%g\n",t);
  fprintf(fd,"MAT1,1,%g,%g,%g,%g\n",E,G,nu,rho);
  if (nx > 1)
  {
    fprintf(fd,"SPC1,1,123,%d,THRU,%d\n",1,nx+1);
    fprintf(fd,"SPC1,1,123,%d,THRU,%d\n",ny*(nx+1)+1,(nx+1)*(ny+1));
  }
  else
  {
    fprintf(fd,"SPC1,1,123,%d,%d\n",1,2);
    fprintf(fd,"SPC1,1,123,%d,%d\n",ny*2+1,2*(ny+1));
  }
  if (ny > 1)
  {
    fprintf(fd,"SPC1,1,123");
    for (j = 1; j < ny; j++)
      fprintf(fd,",%d",(nx+1)*j+1);
    fprintf(fd,"\nSPC1,1,123");
    for (j = 1; j < ny; j++)
      fprintf(fd,",%d",(nx+1)*(j+1));
    fprintf(fd,"\n");
  }
  fprintf(fd,"END DATA\n");

  return 0;
}
