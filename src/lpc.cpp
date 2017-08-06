class LPC {
public:
void memcof(float data[], int n, int m, float *xms, float d[])
{
  int k,j,i;
  float p=0.0,*wk1,*wk2,*wkm;

  wk1=vector(1,n);

  wk2=vector(1,n);

  wkm=vector(1,m);

  for (j=1;j<=n;j++) p += SQR(data[j]);

  *xms=p/n;

  wk1[1]=data[1];

  wk2[n-1]=data[n];

  for (j=2;j<=n-1;j++) {

    wk1[j]=data[j];

    wk2[j-1]=data[j];

  }

  for (k=1;k<=m;k++) {

    float num=0.0,denom=0.0;

    for (j=1;j<=(n-k);j++) {

      num += wk1[j]*wk2[j];

      denom += SQR(wk1[j])+SQR(wk2[j]);

    }

    d[k]=2.0*num/denom;

    *xms *= (1.0-SQR(d[k]));

    for (i=1;i<=(k-1);i++)

      d[i]=wkm[i]-d[k]*wkm[k-i];
                                                                                                        if (k == m) {

                                                                                                          free_vector(wkm,1,m);

                                                                                                          free_vector(wk2,1,n);

                                                                                                          free_vector(wk1,1,n);

                                                                                                          return;

                                                                                                          for (i=1;i<=k;i++) wkm[i]=d[i];

                                                                                                          for (j=1;j<=(n-k-1);j++) {

                                                                                                            wk1[j] -= wkm[k]*wk2[j];

                                                                                                            wk2[j]=wk2[j+1]-wkm[k]*wk1[j+1];

                                                                                                            nrerror("never get here in memcof.");

                                                                                                          }
                                                                                                          
void predic(float data[], int ndata, float d[], int m, float future[], int nfut)

{

  int k,j;

  float sum,discrp,*reg;

  reg=vector(1,m);

  for (j=1;j<=m;j++) reg[j]=data[ndata+1-j];

  for (j=1;j<=nfut;j++) {

    discrp=0.0;

                                                            sum=discrp;

 for (k=1;k<=m;k++) sum += d[k]*reg[k];

 for (k=m;k>=2;k--) reg[k]=reg[k-1]; [If you want to implement circular

                                      future[j]=reg[1]=sum;

                                      }

 free_vector(reg,1,m);

  }
