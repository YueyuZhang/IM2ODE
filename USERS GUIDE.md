#IM2ODE ʹ��˵��
#GUIDE FOR USERS
de.in�Ĳ�������
��δ��˵���ĳ��ȵ�λ����Angstrom��������λ��eV��

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%��Ҫ�������
%%����ṹ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SystemName				#��ϵ���ƣ������һ���ַ��������Ȳ�����20

NumberOfSpecies			#��ϵԪ�����࣬��������,��Χ[1,10]

NumberOfElements		#ÿ��Ԫ�ض�Ӧ��ԭ�Ӹ���������NumberOfSpecies������������Χ[1,100]

NameOfElements			#ÿ��Ԫ�ض�Ӧ��Ԫ�����ƣ�����NumberOfSpecies�����ַ��������Ȳ�����2

Volumn					#��ʼԭ���������Ź�һ�£������������Ժ�ÿ�����Զ�����

DistanceOfAtom			#�������NumberOfSpecies�У�ÿ����NumberOfSpecies�����������ͣ�������
##Ԫ�ؼ����С���룬ÿ��Ԫ��������Ԫ�ؼ�ľ�����һ������DIS1, DIS2��ʾ
## (���ӣ���ɫ����д��de.in���ɫ������ȥ��д)
	DistanceOfAtom=|	Ti|	O
-----------------|----|----
Ti	|DIS1=	|1.6	|1.2
O	|DIS2=	|1.2	|1.0

Population				#��Ⱥ������ÿ�������Ľṹ������������Χ[4, 500]

MaxStep					#���������������ܹ��ܼ���������

De_ratio					#ÿ���ж��ٱ����Ľṹֱ����DE�����������������������������Χ(0,1)
\##Ҫ����De_ratio*Population>4�������޷�����DE�Ĳ�����DE�㷨�������ƣ�
\##�����Ҫÿ����random������DE������ʹDe_ratio*Population=0

SelectiveDynamics 		#ֻ��һ���ֵ�ԭ�ӣ������������ļ�struct.in���������Ч��bool, default = F
\##��Tʱ��struct.in�ļ���д�������仯����ԭ������������T/F, T��ʾ��##vasp����λ���Ż���F��ʾ�̶�������������struct.in_selective_dynamics)
Symmetry				#�Ƿ���Ͽռ�Ⱥ�Ż���bool, default = T(���ڲ�������ṹʱ��Ч)
spg_front					#��ʼ�ռ�Ⱥ�ţ�default=1
spg_rear					#��ֹ�ռ�Ⱥ�ţ�default=230�����������޶������ṹѡȡ�Ŀռ�Ⱥ��Χ
1-2	��б
3-15	��б
16-74	����
75-142	�ķ�
143-167	����
168-194	����
195-230	����

Pickup					#�Ƿ��ĳһ����������, bool
\###ע����results�ļ�����һ��Ҫ��de_opt_*. ����ӵ�5����ʼ�ܣ���###Pickup_step���5, ׼����results/de_opt_5�ļ�
Pickup_step				#��������Ĵ���������

Mystruct					#�����Լ���Ϊ�ȽϺõĽṹ�ĸ���������
\##�ṹ�ļ�����mystruct1, mystruct2, ����ԭ������˳��Ҫ��de.in��һ��
\##��������гĵ׵��������д�ĵ�ԭ�ӣ���дȫ���Ż���ԭ��

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%��Ŀ��ṹ�������μ�example1��
%%�ھ���ṹ�����Ļ����ϼ������ǩ
%%Ŀǰ֧��free energy & hardness ͬʱ�Ż��� free energy & bandgap ͬʱ�Ż�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Multi-Objective			#��Ŀ�꣬bool

Hardness	#�����㷽���õ�PRB32,7988 && Mat. Sci. & Eng. A209 1996 1-4���Ƿ����##bulk modulus, bool����ͬʱ��rcut��ionicity��
Rcut						#max bond length��ȡһ����һ���ں͵ڶ����ڼ����)��������
Iconicity					#����Ԫ�ز�ͬ��ȡ����ϸ�ο�Cohen�����£�PRB41,10727��

ESflag					#�Ƿ���bandgap��bool����ҪдMatrix_elements,���ж�������ԾǨ
Es_mod					#bandgap���Ż���ʽ����������Χ[1,6]
\## 1: maximal bandgap
\## 2: minimal bandgap
\## 3: a target bandgap
\## 4: only direct gap, largest
\## 5: only direct gap, largest (recommended)
\## 6: a target direct gap
ES_Eg					#Ŀ��gap������������Es_mod=3/6ʱ��Ҫ����
ES_opt	#��es_mod=3/6ʱ�����ڴ���ֱ�Ӵ�϶�Ƿ���ԾǨ����������ҽ���vaspд##Matrix_elements��ʱ����Ч����Ҫֱ�Ӵ�϶ԾǨ��ES_opt��ֵ��ĺ�ES_Egһ�����������ES_opt������ԾǨ�Ķ�����

HSE						#��band gap�����м���HSE��bool��������run_pvasp_HSE.shͶ��HSE������
\##��һ���ܾ�̬��LDA��������ֵ��Ҫ��LDA�������ⲽ�ܵ�OUTCACR��##��OUTCAR.old��ȡ������
\## Population���LDA_population�������������²�����
HSE_population			#ÿһ���Ż����������HSE����Ľṹ��
LDA_population			#��LDA����ȫ���Ż��Ľṹ������population������ͬ
energy_cut				#�����ضϣ���ÿ��LDA�����ֻ����������energy_cut�Ľṹ����HSE����
gap_cut					#�ڶ���Ŀ�꺯���Ľض�
LDA_ES_Eg				#ͬ����ES_Eg�Ķ��壬���LDA��������gapֵ
LDA_Es_opt				#ͬ����ES_opt�Ķ��壬���LDA��������gapֵ
HSE_ES_Eg				#ͬ����ES_Eg�Ķ��壬���HSE��������gapֵ
HSE_Es_opt				#ͬ����ES_opt�Ķ��壬���HSE��������gapֵ

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%׼��ά�ṹ�������μ�example2��
%%�ھ���ṹ�����Ļ����ϼ������ϱ�ǩ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q2D					#׼��ά, bool
vacuum_layer			#��ղ㣬����������Ҫ�Ȳ��ԣ�ͨ����10?���ң���������Ľṹ��ղ���c����
Area				#׼��ά����ԭ���ڶ�άƽ����ͶӰ�������������
Layer_hight			#׼��ά���ϵ�ԭ�Ӳ��ȣ�������

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%cluster�������μ�example3��
%%�����򵥵�cluster������cluster+substrate������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ע��һ��Ҫ��fix_lat��ǩ

cluster				#����cluster��ǩ��bool
\## cluster������3��ģʽ���ֱ���ball, shell,��plate����3��ģʽ���ڲ�����ͬ���͵�cluster���������£�
model_ball			#��������cluster�ı�����(0,1)
model_shell			#������״cluster�ı�����(0,1)
model_plate			#������״cluster�ı�����(0,1)�����������������ܺ�Ϊ1���磺
\#model_ball=0.8
\#model_shell=0.1
\#model_plate=0.1

\#model_boll�в������£�
init_radius			#��ʼ��cluster�뾶�������������������ž��У�
cluster_ctr_x			# cluster����������(direct)
cluster_ctr_y
cluster_ctr_z

\#model_shell�в������£�
shell_radius_out		#��״cluster��뾶
shell_radius_in		#��״cluster�ڰ뾶
shell_ctr_x			#��״cluster����������(direct)
shell_ctr_y
shell_ctr_z

\#model_plate�в������£�
plate_radius			#��״cluster�뾶
plate_height			#��״cluster���
plate_ctr_x			#��״cluster����������(direct)
plate_ctr_y
plate_ctr_z

cluster_substrate		#�Ƿ�ӳĵױ�ǩ��bool��������һ�������ļ�struct.in�������ǳĵ׵�ԭ������
\##Ԫ��˳���SubstrateElements��Ӧ����ϵ��NumberOfElements�ǳĵ�ԭ������+clusterԭ������
SubstrateElements		#�ĵ���ÿ��Ԫ�ض�Ӧ��ԭ��������ʽ��˳���϶���NumberOfElements��д��һ��


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%�̶��������������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fix_lat			#��ǩ��bool
fix_a
fix_b
fix_c
fix_alpha
fix_beta
fix_gama			#���϶��Ǹ�������alpha, beta, gama��(0,180)�������ʾ
\##ע����Q2D����´�fix_lat��ǩ��fix_c��ʾԭ�Ӳ��ȶ����Ǿ�����

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%ȱ������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

find_defect		#��ǩ��bool
\## �ھ����в���defect����ʽ��cluster+substrateһ����֮��Ĳ�����ϸ���ü�cluster���õ��Ĳ����У�
cluster=T 		#(fixed)
model_ball=1.0 	#(fixed)
init_radius		#(�������ã�defectԭ���޶��뾶)
cluster_ctr_x 		#(�������ã�defect��������)
cluster_ctr_y		#(ͬ��)
cluster_ctr_z 		#(ͬ��)
cluster_substrate=T
SubstrateElements	#(defect�й̶���ԭ�Ӹ�������struct.in�ļ�����)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%��������(example 4)
%%��Ϊԓģ������cluster_substrate�ϸĵģ�����Ҫ���ϲ���cluster_substrate��tag��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\## ע������ԭ�Ӳ��ԭ��������struct.in�ļ���д
\## ǿ�ҽ�����SelectiveDynamics��struct.in����ԭ�Ӳ�Զ��ȫ���Ż�����F��������������T

model_surface		#��ǩ��bool
surface_height		#����ȫ���Ż���ԭ�Ӳ�ĺ�ȣ���������slab�ĺ�ȣ�������������λ��Astron
cluster=T
cluster_substrate=T
SubstrateElements	#���²��ԭ�Ӹ�������struct.in�ļ����룩

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%��������(example 5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\## ע������ղ�ģ��,����ԭ�Ӳ�����д��struct.in�ļ���

model_gb			#��ǩ��bool
bottom_height		#�ײ�ĵ׵ĺ�ȣ�����㿪ʼ���㣩
gb_height		#ȫ���Ż�����ԭ�Ӳ�ĺ�ȣ�������Զ��Ż���
top_height		#����ԭ�Ӳ��ȣ���0,0,1�㿪ʼ���㣩
transverse_a		#����ԭ��ƽ��x����λ��������ɶȣ���Χ[0,1]
transverse_b		#����ԭ��ƽ��y����λ��������ɶȣ���Χ[0,1]


de.in������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
example1 (Multi-objective)��
\# input file for de

\# search for TiO2, muti-objective, objective bandgap = 0.8
SystemName=TiO
NumberOfSpecies=2
NumberOfElements=4 8
NameOfElements=Ti O
Volumn=132
DistanceOfAtom=
DIS1=1.6 1.2
DIS2=1.2 1.0
Population=30
MaxStep=20
De_ratio=0.6
Symmetry=T

Multi-Objective=T
hardness=F
#rcut=2.5
#ionicity=0.0
ESflag=T
ES_mod=6
ES_Eg=0.8
ES_opt=1.5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Example2: (quasi-2D hexagonal && fix lattice)
\# input file for de
SystemName=Fe Se
NumberOfSpecies=2
NumberOfElements=4 4
NameOfElements=Fe Se
Volumn=160
DistanceOfAtom=
DIS1=2.2 1.5
DIS2=1.5 2.2
Population=8
MaxStep=1
De_ratio=0.6
Symmetry=T
\#spg_front=168
\#spg_rear=194

Q2D=T
vacuum_layer=10
Area=70.15

fix_lat=T
fix_a=7.8
fix_b=7.8
fix_c=0.8
fix_alpha=120
fix_beta=90
fix_gama=90

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Example3:(cluster)
SystemName=PdAu
NumberOfSpecies=2
NumberOfElements=4 1
NameOfElements=Pd Au
Volumn=3375
DistanceOfAtom=
DIS1=1.5 1.5
DIS2=1.5 1.5
Population=15
MaxStep=10
De_ratio=0.6

cluster=T
init_radius=3
cluster_ctr_x=0.5
cluster_ctr_y=0.5
cluster_ctr_z=0.5

fix_lat=T
fix_a=15.0
fix_b=15.0
fix_c=15.0
fix_alpha=90
fix_beta=90
fix_gama=90

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Example4:(surface)
SystemName=TiO
NumberOfSpecies=2
NumberOfElements=16 32
NameOfElements=Ti O
Volumn=954.87
DistanceOfAtom=
DIS1=2.2 1.8
DIS2=1.8 1.2
Population=30
MaxStep=20
De_ratio=0.6
SelectiveDynamics=T
Pickup=F

cluster=T
cluster_substrate=T
SubstrateElements=12 24

model_surface=T
surface_height=3

fix_lat=T
fix_a=10.2913
fix_b=3.8015
fix_c=24.4223
fix_alpha=90
fix_beta=90
fix_gama=90


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Example4:(interface)
SystemName=STOSi
NumberOfSpecies=5
NumberOfElements=28 8 10 28 8
NameOfElements=Si Sr Ti O H
Volumn=1983.885
DistanceOfAtom=
DIS1=1.0 1.5 1.5 1.0 1.0
DIS2=1.5 2.0 2.0 1.5 1.0
DIS3=1.0 2.0 2.0 1.5 1.0
DIS4=1.0 1.5 1.5 1.5 1.0
DIS5=1.0 1.0 1.0 1.0 1.0
Population=30
MaxStep=30
De_ratio=0.6
SelectiveDynamics=T
cluster_substrate=T
SubstrateElements=28 8 8 24 8

\#for grain boundary generation
model_gb=T
\#the following three are cartisian coordinations
bottom_height=14.964
gb_height=4.084
top_height=14.102
\#transverse freedom for the top layer (direct)
transverse_a=1.0
transverse_b=1.0

fix_lat=T
fix_a=7.736
fix_b=7.736
fix_c=35.15              
fix_alpha=90
fix_beta=90
fix_gama=90



���룬���м����ڴ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
���ڳ����������
���ڷ����ĳ��������DE-package_����.tar.gz�����µ��ȶ��汾��DE-package-20140627.tar.gz
ʹ��tar -xzvf DE-package_����.tar.gz�����
װ��gfortran�Ļ������ܱ��룬makeһ�¾�����

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
��������
��de.pbs��de.xͶ��������ȥ�ܣ�Ȼ��de.x�����run_pvasp.sh������ű��ǰ����е�POSCAR*��vasp.pbsͶ�����ϲ������û�����ꡣ
���ȣ�de.x
һ����ڲ�ͬ�Ļ�����һ���⼸���ű������ˡ�������д�¼����ű���ô���еģ�

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
�����ļ�
���к�ÿһ���Ľ������results�ļ�����
de_ini_*��ÿһ����ʼ���Ľṹ
de_opt_*���Ż��õĽṹ��ÿ���ṹ�ļ�ǰ��������(eV/atom)�������muti-objective��ǰ����һ��fitness��ֵ

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
������������tools�ļ����µ�outPOSCAR.cpp��input.dat����results�ļ�����
input.dat���һ�е��Ǹ����Ǹ���ǩ����һ���de��1, mode��2
�ڱ�ǩ=1������£��ڶ�����һ��ʵ ������ʾenergy/atom�����Ǹ����Ľṹȫ��ȡ����
�ڱ�ǩ=2������£��ڶ���������ʵ������һ��ͬ�ϣ��ڶ�����mode����һ��Ŀ�꺯�������ֵ
�����е�������a,b��ʾ�����a-b�����еĽ����һ����a=1, b��de.in�����һ�¡�
���к�out.dat�������нṹ������multi-objective����ǰ�frontier��������
���POSCAR*�Ƿ���input.dat���趨�����Ľṹ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


