#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>
int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }
  for(int i=0; i<N; i++) {
    __m256 axi = _mm256_setzero_ps();
    __m256 ayi = _mm256_setzero_ps();
    __m256 xi = _mm256_set1_ps(x[i]);
    __m256 yi = _mm256_set1_ps(y[i]);
    for(int j=0; j<N; j+=8) {
      __m256 xj = _mm256_load_ps(x+j);
      xj = _mm256_sub_ps(xj,xi);
      __m256 yj = _mm256_load_ps(y+j);
      yj = _mm256_sub_ps(yj,yi);
      __m256 R2 = _mm256_setzero_ps();
      R2[i] = 1;
      R2 = _mm256_fmadd_ps(xj,xj,R2);
      R2 = _mm256_fmadd_ps(yj,yj,R2);
      __m256 mj = _mm256_load_ps(m+j);
      __m256 invR = _mm256_rsqrt_ps(R2);
      mj = _mm256_mul_ps(mj,invR);
      invR = _mm256_mul_ps(invR,invR);
      invR = _mm256_mul_ps(invR,mj);
      axi = _mm256_fmadd_ps(xj,invR,axi);
      ayi = _mm256_fmadd_ps(yj,invR,ayi);
      __m256 c = _mm256_permute2f128_ps(axi,axi,1);
      __m256 d = _mm256_permute2f128_ps(ayi,ayi,1);
      c = _mm256_add_ps(c,axi);
      c = _mm256_hadd_ps(c,c);
      c = _mm256_hadd_ps(c,c);
      d = _mm256_add_ps(d,ayi);
      d = _mm256_hadd_ps(d,d);
      d = _mm256_hadd_ps(d,d);
      _mm256_store_ps(fx,c);
      _mm256_store_ps(fy,d);
    }
    fx[i]=fx[1];
    fy[i]=fy[1];
    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
