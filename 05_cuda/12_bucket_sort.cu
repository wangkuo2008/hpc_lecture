#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void sort(int *key,int *bucket, int n, int range) {
 int i = blockIdx.x * blockDim.x + threadIdx.x;
 for(int j=0,k=0; k<=i; j++) {
   key[i]=j;
   __syncthreads();
   k+=bucket[j];
   __syncthreads();
  }
}

__global__ void bucket1(int *bucket){
 int i = blockIdx.x * blockDim.x + threadIdx.x;
 bucket[i] = 0;
}

__global__ void bucket2(int *key, int *bucket){
 int i = blockIdx.x * blockDim.x + threadIdx.x;
 atomicAdd(&bucket[key[i]],1);
}

int main() {
 int n = 50;
 int range = 5;
 int *key, *bucket;
 cudaMallocManaged(&key,n*sizeof(int));
 cudaMallocManaged(&bucket,range*sizeof(int));
 for (int i=0; i<n; i++) {
   key[i] = rand() % range;
   printf("%d ",key[i]);
  }
 printf("\n");

bucket1<<<1,range>>>(bucket);
bucket2<<<1,n>>>(key,bucket);
sort<<<1,n>>>(key, bucket, n, range);
cudaDeviceSynchronize();

 for (int i=0; i<n; i++) {
   printf("%d ",key[i]);
  }
 printf("\n");
 cudaFree(key);
 cudaFree(bucket);
}
