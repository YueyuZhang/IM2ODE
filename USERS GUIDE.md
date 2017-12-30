#IM2ODE 使用说明
#GUIDE FOR USERS
de.in的参数设置
（未加说明的长度单位都是Angstrom，能量单位是eV）

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%主要输入参数
%%晶体结构搜索
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SystemName				#体系名称，后面加一个字符串，长度不超过20

NumberOfSpecies			#体系元素种类，输入整数,范围[1,10]

NumberOfElements		#每种元素对应的原子个数，输入NumberOfSpecies个，整数，范围[1,100]

NameOfElements			#每种元素对应的元素名称，输入NumberOfSpecies个，字符串，长度不超过2

Volumn					#初始原胞体积（大概估一下），浮点数，以后每代会自动更新

DistanceOfAtom			#下面跟着NumberOfSpecies行，每行有NumberOfSpecies个数（浮点型），输入
##元素间的最小距离，每种元素与其他元素间的距离用一个向量DIS1, DIS2表示
## (例子，黑色部分写在de.in里，灰色部分隐去不写)
	DistanceOfAtom=|	Ti|	O
-----------------|----|----
Ti	|DIS1=	|1.6	|1.2
O	|DIS2=	|1.2	|1.0

Population				#种群数，即每代产生的结构数，整数，范围[4, 500]

MaxStep					#最大代数，即程序总共跑几代，整数

De_ratio					#每代有多少比例的结构直接由DE产生，其余随机产生，浮点数，范围(0,1)
\##要满足De_ratio*Population>4，否则无法进行DE的操作（DE算法本身限制）
\##如果想要每代都random，不用DE，可以使De_ratio*Population=0

SelectiveDynamics 		#只动一部分的原子（仅在有输入文件struct.in的情况下有效，bool, default = F
\##在T时，struct.in文件的写法有所变化，在原子坐标后面加上T/F, T表示由##vasp进行位置优化，F表示固定（见程序包里的struct.in_selective_dynamics)
Symmetry				#是否加上空间群优化，bool, default = T(仅在产生块体结构时有效)
spg_front					#起始空间群号，default=1
spg_rear					#终止空间群号，default=230，这两个数限定产生结构选取的空间群范围
1-2	三斜
3-15	单斜
16-74	正交
75-142	四方
143-167	三角
168-194	六角
195-230	立方

Pickup					#是否从某一代继续计算, bool
\###注：在results文件夹里一定要有de_opt_*. 比如从第5代开始跑，把###Pickup_step设成5, 准备好results/de_opt_5文件
Pickup_step				#继续计算的代数，整数

Mystruct					#输入自己认为比较好的结构的个数，整数
\##结构文件放在mystruct1, mystruct2, 其中原子种类顺序、要和de.in的一致
\##如果遇到有衬底的情况，先写衬底原子，再写全局优化的原子

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%多目标结构搜索（参见example1）
%%在晶体结构搜索的基础上加下面标签
%%目前支持free energy & hardness 同时优化和 free energy & bandgap 同时优化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Multi-Objective			#多目标，bool

Hardness	#（计算方法用的PRB32,7988 && Mat. Sci. & Eng. A209 1996 1-4）是否最大化##bulk modulus, bool，需同时加rcut和ionicity）
Rcut						#max bond length，取一个第一近邻和第二近邻间的数)，浮点数
Iconicity					#根据元素不同来取，详细参考Cohen的文章（PRB41,10727）

ESflag					#是否考虑bandgap，bool，需要写Matrix_elements,来判断哪里有跃迁
Es_mod					#bandgap的优化方式，整数，范围[1,6]
\## 1: maximal bandgap
\## 2: minimal bandgap
\## 3: a target bandgap
\## 4: only direct gap, largest
\## 5: only direct gap, largest (recommended)
\## 6: a target direct gap
ES_Eg					#目标gap，浮点数，当Es_mod=3/6时需要设置
ES_opt	#在es_mod=3/6时，用于处理直接带隙是否有跃迁的情况，当且仅当vasp写##Matrix_elements的时候有效，需要直接带隙跃迁把ES_opt的值设的和ES_Eg一样，其余的在ES_opt以下有跃迁的都可以

HSE						#在band gap运算中加入HSE，bool，需增加run_pvasp_HSE.sh投递HSE的任务
\##第一步跑静态的LDA，算能量值。要在LDA跑完后把这步跑的OUTCACR拷##成OUTCAR.old读取能量。
\## Population设成LDA_population，另外增加以下参数：
HSE_population			#每一代优化完后用来做HSE计算的结构数
LDA_population			#用LDA来做全局优化的结构数，与population含义相同
energy_cut				#能量截断，在每次LDA计算后只有能量低于energy_cut的结构进行HSE计算
gap_cut					#第二个目标函数的截断
LDA_ES_Eg				#同上面ES_Eg的定义，针对LDA计算所得gap值
LDA_Es_opt				#同上面ES_opt的定义，针对LDA计算所得gap值
HSE_ES_Eg				#同上面ES_Eg的定义，针对HSE计算所得gap值
HSE_Es_opt				#同上面ES_opt的定义，针对HSE计算所得gap值

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%准二维结构搜索（参见example2）
%%在晶体结构搜索的基础上加上以上标签
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q2D					#准二维, bool
vacuum_layer			#真空层，浮点数，需要先测试（通常在10?左右，这里产生的结构真空层沿c方向）
Area				#准二维材料原胞在二维平面内投影的面积，浮点数
Layer_hight			#准二维材料的原子层厚度，浮点数

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%cluster搜索（参见example3）
%%包括简单的cluster搜索和cluster+substrate的搜索
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
注：一定要打开fix_lat标签

cluster				#搜索cluster标签，bool
\## cluster下面有3种模式，分别是ball, shell,和plate。这3种模式用于产生不同构型的cluster。定义如下：
model_ball			#产生球形cluster的比例，(0,1)
model_shell			#产生笼状cluster的比例，(0,1)
model_plate			#产生盘状cluster的比例，(0,1)，上述三个加起来总和为1，如：
\#model_ball=0.8
\#model_shell=0.1
\#model_plate=0.1

\#model_boll有参数如下：
init_radius			#初始的cluster半径，浮点数（这里给个大概就行）
cluster_ctr_x			# cluster的中心坐标(direct)
cluster_ctr_y
cluster_ctr_z

\#model_shell有参数如下：
shell_radius_out		#笼状cluster外半径
shell_radius_in		#笼状cluster内半径
shell_ctr_x			#笼状cluster的中心坐标(direct)
shell_ctr_y
shell_ctr_z

\#model_plate有参数如下：
plate_radius			#盘状cluster半径
plate_height			#盘状cluster厚度
plate_ctr_x			#盘状cluster的中心坐标(direct)
plate_ctr_y
plate_ctr_z

cluster_substrate		#是否加衬底标签，bool，需增加一个输入文件struct.in，里面是衬底的原子坐标
\##元素顺序和SubstrateElements对应，体系中NumberOfElements是衬底原子总数+cluster原子总数
SubstrateElements		#衬底中每种元素对应的原子数，形式和顺序上都和NumberOfElements的写法一致


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%固定晶格参数的搜索
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fix_lat			#标签，bool
fix_a
fix_b
fix_c
fix_alpha
fix_beta
fix_gama			#以上都是浮点数，alpha, beta, gama用(0,180)间的数表示
\##注：在Q2D情况下打开fix_lat标签，fix_c表示原子层厚度而不是晶格常数

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%缺陷搜索%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

find_defect		#标签，bool
\## 在晶体中产生defect的形式与cluster+substrate一样，之后的参数详细设置见cluster，用到的参数有：
cluster=T 		#(fixed)
model_ball=1.0 	#(fixed)
init_radius		#(自行设置，defect原子限定半径)
cluster_ctr_x 		#(自行设置，defect中心坐标)
cluster_ctr_y		#(同上)
cluster_ctr_z 		#(同上)
cluster_substrate=T
SubstrateElements	#(defect中固定的原子个数，由struct.in文件读入)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%表面搜索(example 4)
%%因为該模块是在cluster_substrate上改的，所以要加上部分cluster_substrate的tag：
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\## 注：下面原子层的原子坐标在struct.in文件里写
\## 强烈建议用SelectiveDynamics，struct.in里面原子层远离全局优化层设F，靠近的两层用T

model_surface		#标签，bool
surface_height		#表面全局优化洒原子层的厚度（不是整个slab的厚度），浮点数，单位：Astron
cluster=T
cluster_substrate=T
SubstrateElements	#（下层的原子个数，由struct.in文件读入）

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%界面搜索(example 5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\## 注：用真空层模型,上下原子层坐标写在struct.in文件里

model_gb			#标签，bool
bottom_height		#底层衬底的厚度（从零点开始计算）
gb_height		#全局优化的洒原子层的厚度（程序会自动优化）
top_height		#顶层原子层厚度（从0,0,1点开始计算）
transverse_a		#顶层原子平移x方向单位晶格的自由度，范围[0,1]
transverse_b		#顶层原子平移y方向单位晶格的自由度，范围[0,1]


de.in的样例
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
example1 (Multi-objective)：
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



编译，运行及后期处理
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
关于程序包及编译
现在发布的程序包都是DE-package_日期.tar.gz，最新的稳定版本是DE-package-20140627.tar.gz
使用tar -xzvf DE-package_日期.tar.gz对这个
装有gfortran的机器都能编译，make一下就行了

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
运行流程
用de.pbs把de.x投到机器上去跑，然后de.x会调用run_pvasp.sh，这个脚本是把所有的POSCAR*用vasp.pbs投机器上并检测有没有跑完。
首先，de.x
一般对于不同的机器改一下这几个脚本就行了。（具体写下几个脚本怎么运行的）

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
生成文件
运行后每一步的结果都在results文件夹里
de_ini_*是每一代初始化的结构
de_opt_*是优化好的结构，每个结构文件前面是能量(eV/atom)，如果用muti-objective，前面多加一个fitness的值

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
主程序跑完后把tools文件夹下的outPOSCAR.cpp和input.dat拷到results文件夹下
input.dat里第一行的那个数是个标签，跑一般的de用1, mode用2
在标签=1的情况下，第二行设一个实 数，表示energy/atom低于那个数的结构全部取出来
在标签=2的情况下，第二行设两个实数，第一个同上，第二个是mode的另一个目标函数的最大值
第三行的两个数a,b表示处理第a-b代运行的结果。一般设a=1, b和de.in里面的一致。
运行后out.dat里是所有结构能量（multi-objective情况是按frontier）的排序
输出POSCAR*是符合input.dat中设定条件的结构
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


