#include <stdio.h>
#include <cstdlib>

//Convert the function to standard C to
//be called in fortran:
extern "C" void helloworld();
extern "C" void printvariable(double x);
extern "C" void printvector(double* x,int size);
extern "C" void writeheader(const char* filename,const int* ln,const char* comment1,const int* lnc1,const int* flg);
extern "C" void PrintVariable(double x);



void helloworld()
{
  printf("Hello-world from C\n");
}

void printvariable(double x)
{
  printf("Print variable:%.5f\n",x);
}

void printvector(double* x,int size)
{
  int i = 0;
    for(i=0;i<size;++i)
      {
	printf("%.4f ",x[i]);
      }
    printf("\n");
}

void PrintVariable(double x)
{
  printf("Print variable:%.5f\n",x);
}




//!> @brief Write the header file for a vtk file.
//!!
//!!
void writeheader(const char* filename,const int* ln,
                   const char* comment1,const int* lnc1,
                   const int* flg)
{
  int i;
  char *cfilename,*ccomment1,*ccomment2;
  cfilename = (char *)malloc(sizeof(char)*(*ln+1));
  ccomment1 = (char *)malloc(sizeof(char)*(*lnc1+1));
  FILE *fptr;

  /* Checking Endian-ness */
  const char *Endian[] = { "BigEndian", "LittleEndian" };
  unsigned char EndianTest[2] = {1,0};
  short tmp = *(short *)EndianTest;
  if( tmp != 1 ) tmp = 0;
  
  /* string conversion Fortran -> C */
  for (i=0; i<*ln  ; i++) cfilename[i]=filename[i];
  for (i=0; i<*lnc1; i++) ccomment1[i]=comment1[i];
  cfilename[*ln]='\0';
  ccomment1[*lnc1]='\0';
  
  if(*flg == 1){
    /* opening file */
    fptr = fopen(cfilename,"wb"); 

    /* writing header */
    fprintf(fptr, "<?xml version=\"1.0\"?>\n"
	    "<!-- %s -->\n"
	    "<VTKFile type=\"UnstructuredGrid\" "
	    "version=\"0.1\" "
	    "byte_order=\"%s\">\n",ccomment1,Endian[tmp]);
  }
  else if(*flg == 2){
    /* opening file */
    fptr = fopen(cfilename,"wb"); 

    /* writing header */
    fprintf(fptr, "<?xml version=\"1.0\"?>\n"
	    "<VTKFile type=\"Collection\" "
	    "version=\"0.1\" "
	    "byte_order=\"%s\">\n"
	    "<Collection>\n",Endian[tmp]);
  }
  else if(*flg == 3){
    /* opening file */
    fptr = fopen(cfilename,"ab"); 

    fprintf(fptr, "</Collection>\n"
	    "</VTKFile>\n");
  }
  else{
    /* opening file */
    fptr = fopen(cfilename,"ab"); 

    fprintf(fptr, "</VTKFile>\n");
  }
  
  /* closing file */
  fclose(fptr);
  free(cfilename);

  return;
};

