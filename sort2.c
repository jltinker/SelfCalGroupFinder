#define NRANSI
#include "nrutil.h"
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define iSWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50

void sort2(int n, float arr[], int id[])
{
	unsigned long i,ir=n,j,k,l=1;
	int jstack=0,*istack;
	float a,temp,temp_Id;
	int temp_m,itemp;

	istack=ivector(1,NSTACK);
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				temp_m=id[j];
				for (i=j-1;i>=1;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
					id[i+1]=id[i];
				}
				arr[i+1]=a;
				id[i+1]=temp_m;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1]);
			iSWAP(id[k],id[l+1]);
			if (arr[l+1] > arr[ir]) {
			  SWAP(arr[l+1],arr[ir]);
			  iSWAP(id[l+1],id[ir]);
			}
			if (arr[l] > arr[ir]) {
			  SWAP(arr[l],arr[ir]);
			  iSWAP(id[l],id[ir]);
			}
			if (arr[l+1] > arr[l]) {
			  SWAP(arr[l+1],arr[l]);
			  iSWAP(id[l+1],id[l]);
			}
			i=l+1;
			j=ir;
			a=arr[l];
			temp_m=id[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j]);
				iSWAP(id[i],id[j]);
			}
			arr[l]=arr[j];
			id[l]=id[j];
			arr[j]=a;
			id[j]=temp_m;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI
