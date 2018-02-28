#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <string.h>
#define PI 301415926535
#define RAND ((float)rand()/RAND_MAX)

//items数据类型
//#define BirchType int
#define BirchType double

using namespace std;
//命令行解析函数   参数-r 数据文件名称  -n  items项数   -a 属性数   -d 半径阈值
int parseOptions(int argc, char** argv,char filename[256],long  &item_num,int &attribute_num,int &fraction,long &radius)
{
    int i = 0;
    if(argc == 1)
    {
        return -1;
    }
    for (i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-r") == 0)
        {
          i++;
          memcpy(filename,argv[i],strlen(argv[i])+1);
          cout << "filename_len:" << strlen(argv[i]) << endl;
          cout << "filename1:" << filename <<endl;
        }
        else if (strcmp(argv[i], "-n") == 0)
        {
         item_num = (long)atol(argv[++i]);
        }
        else if (strcmp(argv[i], "-a") == 0)
        {
            attribute_num = (int)atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-c") == 0)
        {
            fraction = (int)atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-d") == 0)
        {
         radius = (long)atol(argv[++i]);
        }
        else
        {
             cout << argv[i] <<endl;
             return -1;
        }
    }
    return 0;
}

//运算符重载  + - * k* 用于vector运算
vector<BirchType> operator+(vector<BirchType>aa, vector<BirchType> bb){
	for (int i = 0; i < aa.size(); i++)
		aa[i] += bb[i];
	return aa;
}


vector<BirchType> operator*(vector<BirchType>aa, vector<BirchType>bb){
	for (int i = 0; i < aa.size(); i++)
		aa[i] *= bb[i];
	return aa;
}

vector<BirchType> operator-(vector<BirchType>aa, vector<BirchType>bb){
	for (int i = 0; i < aa.size(); i++)
		aa[i] -= bb[i];
	return aa;
}

vector<BirchType> operator*(vector<BirchType>aa, double k){

	for (int i = 0; i < aa.size(); i++)
		aa[i] = double(aa[i])* k;
	return aa;
}
vector<BirchType> operator*(int k, vector<BirchType>aa){

	for (int i = 0; i < aa.size(); i++)
		aa[i] *= k;
	return aa;
}

class birch//CF树类
{
public:
	struct Attribute//表示一条item数据，属性维度dim，值为data[dim]
	{
		unsigned int dim;
		vector<BirchType>data;
		Attribute(unsigned int d) :dim(d)
		{
			data.resize(dim);
		}
	};
	struct CF//CF节点，表示一个簇，feature <N,LS,SS>
	{
		unsigned int N;
		vector<BirchType> LS;
		vector<BirchType> SS;
		vector<BirchType> data;
		CF(unsigned int N,
			vector<BirchType> LS,
			vector<BirchType>SS) :N(N), LS(LS), SS(SS),data(LS){}
        CF(unsigned int N,
			vector<BirchType> LS,vector<BirchType>SS,vector<BirchType> data)
			:N(N), LS(LS), SS(SS),data(data){}
		CF(unsigned int dim){
			N = 0;
			LS.resize(dim);
			SS.resize(dim);
		};
		CF(){};
	};

	struct Leaf;
	struct MinCluster
	{
		CF cf;
		Leaf*parent;
		MinCluster()
		{
			parent = NULL;
		}
		MinCluster(CF cf)
		{
			parent = NULL;
			this->cf = cf;
		}
	};

	struct Leaf//节点结构  包括叶子节点和非叶子节点
	/*非叶子节点使用 parent child CF
	叶子几点使用pre next parent cluster[] CF*/
	{
		Leaf*pre, *next;//用于创建叶子节点链表
		Leaf*parent;
		vector<Leaf*>*child;
		vector<MinCluster>*cluster;
		CF cf;
		Leaf()
		{
			parent = pre = next = NULL;
			child = NULL;
			cluster = NULL;
		}
	};
	void init_generate_data(int num, int dim, vector<int>&span)
	{
		this->dim = dim;
		for (int i = 0; i < num; i++)
		{
			Attribute att(dim);
			for (int j = 0; j < dim; j++)
				att.data[j] = span[j] * double(rand()) / double(RAND_MAX + 1.0);
			dataset.push_back(att);
		}
	}
 
	//一种数据集生成方法  生成center_x,y为圆心 radius为半径的数据
	void generate_data(int dim,float center_x,float center_y,float radius)
	{
	        this->dim = dim;
			Attribute att(dim);
            float u = sqrt(RAND) * radius;
            float v = RAND *2 * PI;

            att.data[0] = center_x + u * cos(v);
            att.data[1] = center_y + u * sin(v);

            dataset.push_back(att);
	}

	void generate_data_from_file(int dim,vector<BirchType>&value)//用items生成dataset[]
	{
        this->dim = dim;
	    Attribute att(dim);
        for (int j = 0; j < dim; j++)
            att.data[j] = value[j];
        dataset.push_back(att);
	}

	vector<Attribute>dataset;

public:
	birch(unsigned int b, unsigned int l, unsigned int t)//构造函数，给出参数B L T
	/*B 为非叶子节点最多孩子个数  L为叶子几点最多CF簇个数  T为CF簇半径阈值*/
		:B(b), L(l), T(t){
		root = NULL;
		time_t tt;
		srand(time(&tt));
	}
	~birch();
	void insert(Attribute att);


private:

	unsigned int B; //为非叶子节点最多孩子个数
	unsigned int L;//L为叶子几点最多CF簇个数
	//unsigned int T;//CF簇半径阈值
	double T;
	Leaf*root;
	Leaf*head;//叶子节点链表头  便于直接遍历叶子节点
	int dim;


private:
	double cal_inter_cluster_dis(CF &cf1, CF &cf2);
	vector<BirchType>updateSS(vector<BirchType>&LS, vector<BirchType>&SS)
	{
		for (int i = 0; i < LS.size(); i++)
			SS[i] += pow(LS[i], 2.0);
		return SS;
	}
	CF updateCF(CF &c1, CF &c2)//更新CF 合并c1 c2得到新CF 替换cf1
	{
	    for(int i = 0 ; i < c2.data.size() ; i++)
            c1.data.push_back(c2.data[i]);
		return CF(c1.N + c2.N, c1.LS + c2.LS, c1.SS + c2.SS,c1.data);
	}
	void updateCF(Leaf*leaf)
	{
		CF cf(dim);
		if (leaf->cluster != NULL)
		{

			for (int i = 0; i < leaf->cluster->size(); i++)
			{
				cf.N = cf.N + (*leaf->cluster)[i].cf.N;
				cf.LS = cf.LS + (*leaf->cluster)[i].cf.LS;
				cf.SS = cf.SS + (*leaf->cluster)[i].cf.SS;
			}
		}
		else if (leaf->child != NULL)
		{
			for (int i = 0; i < leaf->child->size(); i++)
			{
				cf.N = cf.N + (*leaf->child)[i]->cf.N;
				cf.LS = cf.LS + (*leaf->child)[i]->cf.LS;
				cf.SS = cf.SS + (*leaf->child)[i]->cf.SS;
			}
		}
		leaf->cf = cf;
	}

	MinCluster create_mincluster(Attribute att)
	{
		vector<BirchType>aa;
		aa.resize(att.dim);
		return MinCluster(CF(1, att.data, updateSS(att.data, aa)));
	}

	void insert(Leaf*close, bool &split, MinCluster &clu);
public:
	void dfs(FILE * out,Leaf* p,Leaf* fa);
    int num = 1;
    Leaf* getroot(){
        return root;
    }
    int sum = 0;
    int res = 0;

};

birch::~birch()//析构，回收child和叶子节点链表
{
	Leaf*plist = head;
	while (plist != NULL)
	{
		delete plist->cluster;
		plist = plist->next;
	}
	vector<Leaf*>aa, bb;
	aa.push_back(root);
	while (!aa.empty())
	{
		Leaf*pleaf = aa.back();
		aa.pop_back();
		bb.push_back(pleaf);
		if (pleaf->child != NULL)
			aa.insert(aa.end(), pleaf->child->begin(), pleaf->child->end());
	}
	for (int i = 0; i < bb.size(); i++)
	{
		if (bb[i]->child != NULL)
			delete bb[i]->child;
		delete bb[i];
	}
}

void birch::insert(Attribute att)//建树 逐条插入item
{
	if (root == NULL)//第一个节点
	{
		root = new Leaf;
		root->cluster = new vector < MinCluster > ;
		(*root->cluster).push_back(create_mincluster(att));//用这条数据创建一条Mincluster
		root->cf = CF((*root->cluster)[0].cf);
		head = root;
		head->pre = NULL;
		head->next = NULL;
		return;
	}
	MinCluster clu = create_mincluster(att);//否则 用这条数据创建一条新的Mincluster
	Leaf*leaf = root;

	vector<int>path;

	while (leaf->cluster == NULL)//从root开始  遍历非叶子节点
	{
		int k = -1;
		double mindis = 10000000000000;
		double dd;
		for (int i = 0; i < (*leaf->child).size(); i++)//所有child数组节点
		{
			double dis = cal_inter_cluster_dis(clu.cf, (*leaf->child)[i]->cf);//找与新clu的cf最近的leaf->child[i]
			if (dis < mindis)
			{
				mindis = dis;
				k = i;//记录leaf->child[i]下标到path[]
			}
			dd = dis;
		}

		path.push_back(k);
		leaf = (*leaf->child)[k];//继续循环，找这个中间节点的child数组节点中与clu距离最近的leaf->child[i]，下标放入path[]
		//到叶子节点停止
	}


	int k = -1;
	double mindis = 100000;
	for (int i = 0; i < (*leaf->cluster).size(); i++)//最终到path的叶子节点，对于节点的cluster数组
	{
		double dis = cal_inter_cluster_dis(clu.cf, (*leaf->cluster)[i].cf);//找与当前clu的cf最近的leaf->cluster[i]，记录i
		if (dis < mindis)
		{
			mindis = dis;
			k = i;
		}
	}
	double ttt = cal_inter_cluster_dis(clu.cf, (*leaf->cluster)[k].cf);//当前clu的cf与最近的leaf->cluster[i].cf的距离
	if (ttt < T)
	{
		(*leaf->cluster)[k].cf = updateCF((*leaf->cluster)[k].cf, clu.cf);//小于T，更新leaf->cluster[i].cf（放入新clu的cf）
		//absorbnum++;
	}
	else//大于T，新clu的cf放入这个leaf的cluster数组
	{
		(*leaf->cluster).push_back(clu);
	}
	//更新加入clu的这个path的所有节点的cf值（把clu.cf加进去）
	Leaf*lea = root;//路径开始是root
	(*lea).cf = updateCF((*lea).cf, clu.cf);
	for (int i = 0; i < path.size(); i++)
	{
		(*lea->child)[path[i]]->cf = updateCF((*lea->child)[path[i]]->cf, clu.cf);//path中root下面那层的中间节点
		lea = (*lea->child)[path[i]];//继续走到child层
	}

	if ((*leaf->cluster).size() > L)//叶子节点cluster数组长度超L，split
	{
		double maxdis = 0;
		int th1 = -1;
		int th2 = -1;
		double**dismatrix = new double*[(*leaf->cluster).size()];
		for (int i = 0; i < (*leaf->cluster).size(); i++)
			dismatrix[i] = new double[(*leaf->cluster).size()];
		//在当前leaf的cluster数组中找到距离最远的两个簇
		for (int i = 0; i < (*leaf->cluster).size() - 1; i++)
			for (int h = i + 1; h < (*leaf->cluster).size(); h++)
			{
				double dis = cal_inter_cluster_dis((*leaf->cluster)[i].cf, (*leaf->cluster)[h].cf);
				dismatrix[i][h] = dis;
				dismatrix[h][i] = dis;
				if (dis > maxdis)
				{
					maxdis = dis;
					th1 = i; th2 = h;
				}
			}
		Leaf*new_leaf = new Leaf;
		new_leaf->cluster = new vector < MinCluster > ;
		new_leaf->cluster->push_back((*leaf->cluster)[th2]);//把2簇放到new_leaf
		int len = (*leaf->cluster).size();
		(*leaf->cluster)[th2].parent = new_leaf;//th2创建的new_leaf做child

		//根据各簇与两个新簇的距离将原节点分配到两个新簇中
		for (int i = 0; i < len; i++)
		{
			if (i == th1 || i == th2)
				continue;
			if (dismatrix[i][th2] < dismatrix[i][th1])
			{
				(*leaf->cluster)[i].parent = new_leaf;
				new_leaf->cluster->push_back((*leaf->cluster)[i]);

			}
		}
		for (int i = 0; i < (*leaf->cluster).size(); i++)
			delete[] dismatrix[i];
		delete[]dismatrix;

		vector < MinCluster >::iterator it, it1;
		it = (*leaf->cluster).begin();
		while (it != (*leaf->cluster).end())
		{
			if (it->parent == new_leaf)
				it = (*leaf->cluster).erase(it);
			else
			{
				it++;
			}
		}
		//更新leaf和new_leaf的cf值
		updateCF(leaf);
		updateCF(new_leaf);
		//将new_leaf加入到叶子节点链表中
		Leaf*next = leaf->next;
		leaf->next = new_leaf;
		new_leaf->pre = leaf;
		new_leaf->next = next;
		if (next)
			next->pre = new_leaf;
		if (leaf->parent != NULL)
		{
			leaf->parent->child->push_back(new_leaf);
			new_leaf->parent = leaf->parent;
		}
		else//root是叶子节点，创建新root
		{
			Leaf*new_root = new Leaf;
			new_root->child = new vector < Leaf* > ;
			new_root->child->push_back(leaf);
			new_root->child->push_back(new_leaf);
			leaf->parent = new_root;
			new_leaf->parent = new_root;
			updateCF(new_root);
			root = new_root;
			return;
		}
	}
	//叶子节点超L问题解决，接着判断child数组是否超B
	Leaf*cur = leaf->parent;//找叶子节点的parent，加入clu后是否超过B
	while (cur != NULL&&cur->child->size() > B)//中间节点child数组超过B，split，循环找parent
	{
		double maxdis = 0;
		int th1 = -1;
		int th2 = -1;
		double**dismatrix = new double*[cur->child->size()];
		for (int i = 0; i < cur->child->size(); i++)
			dismatrix[i] = new double[cur->child->size()];
		//找到距离最远的两个child
		for (int i = 0; i < cur->child->size() - 1; i++)
			for (int h = i + 1; h < cur->child->size(); h++)
			{
				double dis = cal_inter_cluster_dis((*cur->child)[i]->cf, (*cur->child)[h]->cf);
				dismatrix[i][h] = dis;
				dismatrix[h][i] = dis;
				if (dis > maxdis)
				{
					maxdis = dis;
					th1 = i; th2 = h;
				}
			}

		Leaf*new_leaf1 = new Leaf;
		new_leaf1->child = new vector < Leaf* > ;
		(*cur->child)[th2]->parent = new_leaf1;
		(*new_leaf1->child).push_back((*cur->child)[th2]);
		int len = (*cur->child).size();

		//将其他叶子节点分到th1 th2做child
		for (int i = 0; i < len; i++)
		{
			if (i == th1 || i == th2)
				continue;
			if (dismatrix[i][th2] < dismatrix[i][th1])
			{
				(*cur->child)[i]->parent = new_leaf1;
				new_leaf1->child->push_back((*cur->child)[i]);

			}
		}
		for (int i = 0; i < (*cur->child).size(); i++)
			delete[] dismatrix[i];
		delete[]dismatrix;

		vector < Leaf* >::iterator it;
		it = (*cur->child).begin();
		while (it != (*cur->child).end())
		{
			if ((*it)->parent == new_leaf1)
				it = (*cur->child).erase(it);
			else
				it++;
		}
		//更新cur和new_leaf1的cf值
		updateCF(cur);
		updateCF(new_leaf1);

		//root是叶子节点，创建新root
		if (cur->parent == NULL)
		{
			Leaf*new_root = new Leaf;
			new_root->child = new vector < Leaf* > ;
			new_root->child->push_back(cur);
			new_root->child->push_back(new_leaf1);
			cur->parent = new_root;
			new_leaf1->parent = new_root;
			updateCF(new_root);
			root = new_root;
			return;
		}

		//将new_leaf1加入cur的父亲节点的child
		cur->parent->child->push_back(new_leaf1);
		new_leaf1->parent = cur->parent;
		cur = cur->parent;
	}

}

void birch::dfs(FILE * out,Leaf * p,Leaf * fa){//建树完成  遍历打印结果
    Leaf* it ;
    //int old_clu_num;
    if( p->child != NULL){  //&& ((*p->child)[0])->child != NULL ){
        for(int i = 0 ; i < (*p->child).size() ; i++)
        {
            if(((*p->child)[i])->child != NULL)
                dfs(out,(*p->child)[i],p);
            else{
                for(int j = 0 ; j < (*(*p->child)[i]->cluster).size() ; j++ ){
                    cout<<" 第"<<num<<"类"<<"包含个数为:"<<(*(*p->child)[i]->cluster)[j].cf.N<<endl;
                    fprintf(out,"第%d类包含个数为:%d\n",num,(*(*p->child)[i]->cluster)[j].cf.N);
                    num++;
                    for(int k = 0 ; k < (*(*p->child)[i]->cluster)[j].cf.data.size() ; k+=this->dim){
                        /*BirchType x = (*(*p->child)[i]->cluster)[j].cf.data[k];
                        BirchType y = (*(*p->child)[i]->cluster)[j].cf.data[k+1] ;
                        old_clu_num = 10 * (int)((x + 2.5)/5 - 1) + (int)((y + 2.5)/5);
                        cout <<old_clu_num<<":"<< (*(*p->child)[i]->cluster)[j].cf.data[k] << " "
                            << (*(*p->child)[i]->cluster)[j].cf.data[k+1] <<endl;

                        fprintf(out,"old_clu_num:%d  ,%f  %f\n",old_clu_num,(*(*p->child)[i]->cluster)[j].cf.data[k],(*(*p->child)[i]->cluster)[j].cf.data[k+1]);
                        */
                        for(int l = 0;l < this->dim;l++)
                          {
                              cout << (*(*p->child)[i]->cluster)[j].cf.data[l] << " ";
                              fprintf(out,"%f  ",(*(*p->child)[i]->cluster)[j].cf.data[l]);
                          }
                        cout << endl;
                        fprintf(out,"\n");
                    }
                }
            }
        }
    }
}


double birch::cal_inter_cluster_dis(CF &cf1, CF &cf2)//计算2个cf间距离（各维分别求差取平方和）
{
	double dis = 0;
	double temp;
	for (int i = 0; i < dim; i++)
	{
		double t1 = double(cf1.LS[i]) / double(cf1.N);
		double t2 = double(cf2.LS[i]) / double(cf2.N);
		temp = t1 - t2;
		dis += temp*temp;
	}

	return sqrt(dis);
}

int main(int argc,char **argv)
//-r 1.txt -n 10 -a 2 -d 5
//int main(void)
{
    //程序开始计时
    clock_t start, finish;
	double duration;
	//start = clock();

	//解析命令行参数   数据集文件名   数据条数   半径阈值
    char filename[256];
	long item_num = 10;
	int attribute_num = 2;
    int fraction = 0;
	long radius = 10;


    if(-1 == parseOptions(argc, argv, filename, item_num,attribute_num,fraction, radius))
    {
        cout << "Get parameter failed!" << endl;
        exit(-1);
    }

    cout << "filename:" << filename << endl;
    cout << "item_num:" << item_num << endl;
    cout << "attribute_num:" << attribute_num << endl;
    //cout << "fraction:" << fraction << endl;
    cout <<"radius:" << radius << endl;


    FILE * out;
    out = fopen("out.txt","w");
    if(!out)
    {
        printf("file open error!\n");
        return -1;
    }

    //读文件获取数据集  读到bir.dataset[]
    birch bir(5, 6, radius);//B默认5 L默认4

    //const char * file = "1.txt";
    //memcpy(filename,file,sizeof(filename));
	FILE * data;
	data = fopen(filename,"r");
	if(!data)
	{
		cout << "file open error" << endl;
		return -1;
	}

	vector<BirchType>value;
	char  temp [1024];
	size_t length = 1024;
	int len;
	while(fgets(temp,length,data)){//按行读取数据文件  空格分割属性
       len = 0;
       value.clear();
       //cout << "temp:" << temp <<endl;
	   char * token = strtok(temp," ");
       while(token!= NULL)
       {
          value.push_back(atof(token));
          //cout << "value:" << value[len] << " ";
          token = strtok(NULL," ");
          len++;
       }
       bir.generate_data_from_file(attribute_num,value);
       //cout << endl;
	}

	//打印当前birch的dataset[]
	for (int i = 0; i < item_num; i++)
    {
      cout <<i<<":";
      for (int j = 0; j < attribute_num; j++)
      {
          cout << bir.dataset[i].data[j] << "   ";
      }
      cout << endl;
    }

    /////////////////////////////////////////////////////////////////////////////////////
    //开始用bir.dataset[]建树
    start = clock();//从开始建树计时 不计读文件时间
	for (int i = 0; i < item_num; i++)
		bir.insert(bir.dataset[i]);

	bir.dfs(out,bir.getroot(),NULL);
	cout <<"一共分了" << bir.num-1 << "个类" << endl;
	fprintf(out,"一共分了%d个类\n",bir.num-1);

	//程序结束计时
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout<< "running " << duration << " seconds" << endl;
	fprintf(out,"running%lf seconds\n",duration);

	system("pause");
	return 0;
}
