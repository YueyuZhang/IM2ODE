#include<stdio.h>
#include<string.h>
char str[200];
char tt1[50], tt2[50], tt3[50], tt4[50], tt5[50], tt6[50], tt7[50], tt8[50], tt9[50];
int n,m, i,j,k;
int tr[3][150], ty[150], ntp[150];
double etot, tmp;
double lat[3][3], pos[3][150];

int main(){
  freopen("MOVEMENT", "r", stdin);
  while(scanf("%d", &n) != EOF){
    //printf("%d\n", n);
    scanf("%s %s %s %s %s %lf %s %s %s %s %s %s %s %s", tt1, tt1, tt1, tt1, tt1, &etot, tt1, tt1, tt1, tt1, tt1, tt1, tt1, tt1);
    //printf("%.6f\n", etot);
    scanf("%s %s", tt1, tt1);
    for(i = 0; i < 3; i ++){
      for(j = 0; j < 3; j ++){
        scanf("%lf", &lat[i][j]);
      }
    }
    scanf("%s %s %s %s", tt1, tt1, tt1, tt1);
    for(i = 0; i < n; i ++){
      scanf("%d %lf %lf %lf %d %d %d", &ty[i], &pos[0][i], &pos[1][i], &pos[2][i], &tr[0][i], &tr[1][i], &tr[2][i], &tr[3][i]);
    }
    scanf("%s", tt1);
    for(i = 0; i < 8; i ++) 
      for(j = 0; j < 4; j ++)
        scanf("%lf", &tmp);
    scanf("%s", tt1);
    m = 1;
    k = ty[0];
    ntp[0] = 1;
    for(i = 1; i < n; i ++){
      if(ty[i] != ty[i - 1]){
        ty[m] = ty[i];
        m ++;
      }
      ntp[m-1] ++;
    }
  }
  freopen("CONTCAR", "w", stdout);
  puts("system");
  puts("1.000000");
  for(i = 0; i < 3; i ++){
    for(j = 0; j < 3; j ++){
      printf("     %.8f", lat[i][j]);
    }
    puts("");
  }
  for(i = 0; i < m; i ++)
    printf("   %d", ty[i]);
  puts("");
  for(i = 0; i < m; i ++)
    printf("   %d", ntp[i]);
  puts("");
  puts("Direct");
  for(i = 0; i < n; i ++){
    for(j = 0; j < 3; j ++){
      printf("   %.6f", pos[j][i]);
    }
    puts("");
  }
  freopen("OUTCAR", "w", stdout);
  printf("  free  energy   TOTEN  =       %.8f eV\n", etot);
  return 0;
}
