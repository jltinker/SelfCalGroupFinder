#include <stdlib.h>

int search(int n, float *x, float val)
{
  int first, last, middle;
  first = 1;
   last = n ;
   middle = (first+last)/2;
   while (first <= last) {
      if (x[middle] < val)
         first = middle + 1;    
      else
         last = middle - 1;
      middle = (first + last)/2;
   }
   return first;
}
