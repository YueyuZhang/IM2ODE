#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<algorithm>
using namespace std;
#define eps 1e-8
#define inf 0xfffffff
#define MAX 20000

struct node{
	double a, b;
	int c;
	friend bool operator<(node p, node q){
		return p.a < q.a;
	}
}ss[MAX];

int sn, fn, cnt, pn, tag;
double eng_cut, gap_cut;
int i,j,k;
int f1[MAX];
char str[100];
bool flag;
int main(){
	cnt = 0;
	pn = 1;
	freopen("input.dat", "r", stdin);
	scanf("%d", &tag);
	scanf("%lf", &eng_cut);
	if(tag == 2)
		scanf("%lf", &gap_cut);
	else
		gap_cut = eng_cut;
	scanf("%d %d", &sn, &fn);
//	printf("%d %d\n", sn, fn);
	for(i = sn; i <= fn; i ++){
		sprintf(str, "de_opt_%d", i);
		freopen(str, "r", stdin);
		while(scanf("%lf", &ss[cnt].a)!= EOF){
			if(tag == 2)
				scanf("%lf", &ss[cnt].b);
			else
				ss[cnt].b = ss[cnt].a; //cnt represents total of all the structures.
//			printf("%.6f %.6f\n", ss[cnt].a, ss[cnt].b);
			ss[cnt].c=0;
			flag = 0;
			if(ss[cnt].a < eng_cut && ss[cnt].b < gap_cut)
				flag = 1;
			cnt ++;
			if(flag) {
				ss[cnt-1].c=pn;
				sprintf(str, "POSCAR%d", pn++);//pn represents the number of wanted structures.
			}
			if(flag) freopen(str, "w", stdout);
			for(j = 0; j < 8; j ++){
				gets(str);
				if(flag && j > 0) puts(str);
			}
			int ln = strlen(str), p = 0, q = 0;
			for(j = 0; j < ln; j ++){
				if(str[j] >= '0' && str[j] <= '9'){
					q *= 10;
					q += str[j] - '0';
				}
				else{
					p += q;
					q = 0;
				}
			}
			p += q;
			for(j = 0; j < p + 1; j ++){
				gets(str);
				if(flag) puts(str);
			}
		}
	}
	sort(ss, ss+cnt);
	freopen("out.dat", "w", stdout);
	printf("Frontier   NUM_POSCAR   Energy   Fitness\n");
	memset(f1, 0, sizeof(f1));
	flag = 1;
	double gap;
	int front = 0;
	while(flag){
		front ++;
		for(i = 0; i < cnt; i ++){
			if(f1[i] == 0){
				f1[i] = 1;
				gap = ss[i].b;
				printf("%d\t%d\t%.6f\t%.6f\n", front, ss[i].c, ss[i].a, ss[i].b);
				break;
			}
		}
		for(; i < cnt; i ++){
			if(f1[i] == 0 && ss[i].b < gap){
				f1[i] = 1;
				gap = ss[i].b;
				//printf("%d %.6f %.6f\n", front, ss[i].a, ss[i].b);
				printf("%d\t%d\t%.6f\t%.6f\n", front, ss[i].c, ss[i].a, ss[i].b);
			}
		}
		flag = 0;
		for(i = 0; i < cnt; i ++)
			if(f1[i] == 0)
				flag = 1;
	}
	return 0;
}

