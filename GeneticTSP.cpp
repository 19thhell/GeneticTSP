#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <memory.h>
using namespace std;

//Particular Const
#define NUM_OF_CITY 48

//Basic Const
#define POP_SIZE 50
#define MAX_GENERATION 100000
#define NUM_OF_GENE NUM_OF_CITY
#define P_CROSSOVER 0.7
#define P_MUTATION 0.01
#define P_SUBSTITUTE 0.02
#define FACTOR 100000000
#define MAX_STABLE MAX_GENERATION * 2 / 5

//Particular Variable
double dist[NUM_OF_CITY + 1][NUM_OF_CITY + 1] = {0};
int order_table[NUM_OF_GENE] = {0};
int max_duration;

//Basic Variable
int generation;
int best_record;
int best_mem;
int worst_mem;

//Mutation adjust
double adjust[MAX_GENERATION + 1];

class Gene
{
public:
	int order;
};

class Genetype
{
public:
	Genetype():fitness(0),rfitness(0),cfitness(0){}
    double max_load;
	//Basic Property
	Gene gene[NUM_OF_GENE];
	double fitness;
	double rfitness;
	double cfitness;
	bool chosen;
    void display(){
        for (int i = 0;i < NUM_OF_GENE;i++)
            cout << this->gene[i].order << " ";
        cout << endl;
    }
};
Genetype population[POP_SIZE + 1];
Genetype newpopulation[POP_SIZE + 1];
Genetype history_best;

void initialize(int,int);
void evaluate();
void keep_best();
void elitist();
void select();
void crossover();
void crossover(int,int);
void substitute();
void substitute(int,int);
int compete(int,int);
void mutate();
void rollback();
void swirl(int,Genetype &);
void report(int,string);
void trace(int,string);

void initialize()
{
    ifstream fin("data.txt");
    int order;
    double x[NUM_OF_GENE],y[NUM_OF_GENE];
    for (int i = 1;i <= NUM_OF_CITY;i++)
        fin >> order >> x[i] >> y[i];
    fin.close();
    for (int i = 1;i <= NUM_OF_CITY;i++)
        for (int j = 1;j <= NUM_OF_CITY;j++)
            dist[i][j] = dist[j][i] = (double)sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));
    for (int i = 0;i < NUM_OF_GENE;i++)
        order_table[i] = NUM_OF_GENE - i;
	//初始化种群
	for (int i = 0;i < POP_SIZE;i++)
	{
        population[i].gene[0].order = rand() % NUM_OF_CITY + 1;
        bool visit[NUM_OF_CITY + 1];
        for (int j = 1;j <= NUM_OF_CITY;j++)
            visit[j] = false;
        visit[population[i].gene[0].order] = true;
        for (int j = 1;j < NUM_OF_GENE;j++){
            double minimum = 999999;
            int index = 1;
            for (int k = 1;k <= NUM_OF_CITY;k++){
                if (k == population[i].gene[j - 1].order)
                    continue;
                if (dist[population[i].gene[j - 1].order][k] < minimum && !visit[k]){
                    index = k;
                    minimum = dist[population[i].gene[j - 1].order][k];
                }
            }
            population[i].gene[j].order = index;
            visit[index] = true;
        }
        for (int j = 0;j < NUM_OF_GENE;j++)
          for (int k = j + 1;k < NUM_OF_GENE;k++)
            if (population[i].gene[k].order > population[i].gene[j].order)
              population[i].gene[k].order--;
	}
	population[POP_SIZE] = population[0];
}

//结合Grefenstetee解码进行适应性评价
void evaluate()
{
    double sum,t;
    for (int i = 0;i < POP_SIZE;i++){
        sum = 0;
        int origin[NUM_OF_GENE];
        for (int j = 0;j < NUM_OF_GENE;j++)
            origin[j] = population[i].gene[j].order;
        for (int j = NUM_OF_GENE - 1;j >= 0;j--){
            for (int k = j - 1;k >= 0;k--){
                if (origin[k] <= origin[j])
                  origin[j]++;
            }
        }
        for (int j = 0;j < NUM_OF_GENE;j++)
            sum += dist[origin[j]][origin[(j + 1) % NUM_OF_GENE]];
        population[i].fitness = FACTOR / sum;
    }
}

//2_opt局部搜索
void opt_2(){
    double sum;
    bool flag;
    for (int i = 0;i < POP_SIZE;i++){
        flag = false;
        int origin[NUM_OF_GENE];
        for (int j = 0;j < NUM_OF_GENE;j++)
            origin[j] = population[i].gene[j].order;
        for (int j = NUM_OF_GENE - 1;j >= 0;j--){
            for (int k = j - 1;k >= 0;k--){
                if (origin[k] <= origin[j])
                    origin[j]++;
            }
        }
        int a,b,tmp;
        sum = 0;
        a = rand() % NUM_OF_GENE;
        b = rand() % NUM_OF_GENE;
        while (a + 1 == b || b + 1 == a)
            b = rand() % NUM_OF_GENE;
        tmp = origin[(a + 1) % NUM_OF_GENE];
        origin[(a + 1) % NUM_OF_GENE] = origin[b];
        origin[b] = tmp;
        for (int j = 0;j < NUM_OF_CITY;j++)
            sum += dist[origin[j]][origin[(j + 1) % NUM_OF_CITY]];
        if (FACTOR / sum > population[i].fitness)
            flag = true;
        if (flag) {
            for (int j = 0;j < NUM_OF_GENE;j++)
              for (int k = j + 1;k < NUM_OF_GENE;k++)
                if (origin[k] > origin[j])
                    origin[k]--;
            for (int j = 0;j < NUM_OF_GENE;j++)
                population[i].gene[j].order = origin[j];
        }
    }
}

//保留当前代最优个体
void keep_best()
{
	int mem,i;
	best_record = 0;
	for (mem = 0;mem < POP_SIZE;mem++)
	{
		if (population[mem].fitness > population[POP_SIZE].fitness)
		{
			best_record = mem;
			population[POP_SIZE].fitness = population[mem].fitness;
		}
	}
    population[POP_SIZE] = population[best_record];
}

//精炼
void elitist()
{
	int i;
	double best,worst;
	best = population[0].fitness;
	worst = population[0].fitness;
	for (i = 0;i < POP_SIZE - 1;i++)
	{
		if (population[i].fitness > population[i + 1].fitness)
		{
			if (population[i].fitness >= best)
			{
				best = population[i].fitness;
				best_mem = i;
			}
			if (population[i + 1].fitness <= worst)
			{
				worst = population[i + 1].fitness;
				worst_mem = i + 1;
			}
		}
		else
		{
			if (population[i].fitness <= worst)
			{
				worst = population[i].fitness;
				worst_mem = i;
			}
			if (population[i + 1].fitness >= best)
			{
				best = population[i + 1].fitness;
				best_mem = i + 1;
			}
		}
	}
	if (best > population[POP_SIZE].fitness)
        population[POP_SIZE]= population[best_mem];
	else
		population[worst_mem] = population[POP_SIZE];
}

//轮盘选择
void select()
{
	int first = 0,mem,i,j,one;
	double sum = 0,p,x;
	for (mem = 0;mem < POP_SIZE;mem++)
		sum += population[mem].fitness;
	for (mem = 0;mem < POP_SIZE;mem++)
	{
		population[mem].rfitness = population[mem].fitness / sum;
		population[mem].chosen = false;
	}
	population[0].cfitness = population[0].rfitness;
	for (mem = 1;mem < POP_SIZE;mem++)
		population[mem].cfitness = population[mem - 1].cfitness + population[mem].rfitness;
	for (i = 0;i < POP_SIZE;i++)
	{
		p = rand() % 1000 / 1000.0;
		if (p < population[0].cfitness)
			newpopulation[i] = population[0];
		else
			for (j = 0;j < POP_SIZE;j++)
				if (p >= population[j].cfitness && p < population[j + 1].cfitness)
					newpopulation[i] = population[j + 1];
	}
}

//两点交叉
void crossover()
{
	int first = 0,mem,one;
	double x;
	for (mem = 0;mem < POP_SIZE;mem++)
	{
		x = rand() % 1000 / 1000.0;
		if (x < P_CROSSOVER)
		{
			first++;
			if (first % 2 == 0)
				crossover(one,mem);
			else one = mem;
		}
	}
}

void crossover(int one,int two)
{
	int i,a,b;
	a = rand() % (NUM_OF_GENE - 1) + 1;
	b = rand() % (NUM_OF_GENE - a) + a;
	for (i = a;i < b;i++)
		swap(population[one].gene[i],population[two].gene[i]);
}

//修补
void substitute()
{
	int first = 0,mem,one;
	double x;
	for (mem = 0;mem < POP_SIZE;mem++)
	{
		x = rand() % 1000 / 1000.0;
		if (x < P_CROSSOVER)
		{
			first++;
			if (first % 2 == 0)
				substitute(one,mem);
			else one = mem;
		}
	}
}

void substitute(int one,int two)
{
	int i;
    double x;
	for (i = 0;i < NUM_OF_GENE;i++)
	{
        x = rand() % 1000 / 1000.0;
        if (x < P_SUBSTITUTE)
            population[one].gene[i] = population[two].gene[i];
	}
}

//竞争选择
int compete(int one,int two)
{
	double p;
	if (population[one].fitness < population[two].fitness)
		return one;
	else
	{
		p = rand() % 1000 / 1000.0;
		if (p > 1 / (1 + exp((population[one].fitness + population[two].fitness) / generation)))
			return two;
		else return one;
	}
}

//单点变异
void mutate()
{
	int i,j;
	double x;
	for (i = 0;i < POP_SIZE;i++)
		for (j = 0;j < NUM_OF_GENE - 1;j++)
		{
			x = (rand() % 1000) / 1000.0;
			if (x < P_MUTATION * adjust[generation])
                population[i].gene[j].order = rand() % order_table[j] + 1;
		}
}

//数据记录
void report(int t)
{
	ofstream fout("output_no_opt" + to_string(t) + ".txt",ios::app);
	fout << fixed << setprecision(4) << FACTOR / history_best.fitness << "\n";
	fout.close();
}

int main()
{
    cout << "Genetic TSP demo. \nData set: 48 city\nGeneration: 100000\n";
	int count,j;
	clock_t start,end;
	string progress_bar,prefix;
	double complete_percent;
	for (int i = 1;i <= MAX_GENERATION;i++)
		adjust[i] = 0.3 + exp(float(-i));
	srand((unsigned int)time(0));
	start = clock();
    count = 0;
    complete_percent = 0;
    for (int t = 1;t <= 5;t++){
    ofstream fout("route_no_opt" + to_string(t) + ".txt");
    generation = 0;
    initialize();
    history_best = population[POP_SIZE];
    ifstream fsin("answer.txt");
    int a[NUM_OF_CITY];
    for (int i = 0;i < NUM_OF_CITY;i++)
      fsin >> a[i];
    fsin.close();
    double sum = 0;
    for (int j = 0;j < NUM_OF_GENE;j++)
        sum += dist[a[j]][a[(j + 1) % NUM_OF_CITY]];
    cout << "Best solution: " << sum << endl;
    
    evaluate();
    keep_best();
    while (generation++ < MAX_GENERATION)
    {
        elitist();
        if (population[POP_SIZE].fitness > history_best.fitness)
            history_best = population[POP_SIZE];
        select();
        report(t);
        crossover();
        mutate();
        evaluate();
//        opt_2();
    }
    for (int j = NUM_OF_GENE - 1;j >= 0;j--){
        for (int k = j - 1;k >= 0;k--){
            if (history_best.gene[k].order <= history_best.gene[j].order)
              history_best.gene[j].order++;
        }
    }
    for (int i = 0;i < NUM_OF_GENE;i++){
        fout << history_best.gene[i].order << endl;
    }
    fout.close();
    }

    end = clock();
    cout << "Current solution: " << FACTOR / history_best.fitness << endl;
    printf("\nCompleted.\nTime used: %.2f min\n",(double)(end - start) / (CLOCKS_PER_SEC * 60));
    
    return 0;
}
