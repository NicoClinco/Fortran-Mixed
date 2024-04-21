#include <stdio.h>

//Convert the function to standard C to
//be called in fortran:
extern "C" void helloworld();
extern "C" void printvariable(double x);
extern "C" void printvector(double* x,int size);

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


