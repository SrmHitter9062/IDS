
//-----------------------------------------project_train_main.cpp--------------------------------------


#include <iostream>
#include<vector>
#include<set>
#include<map>
#include<queue>
#include<stack>
#include<string>
#include<algorithm>
#include<functional>
#include<iomanip>
#include<cstdio>
#include<cmath>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<climits>
#include <fstream>
#include <utility>
#include <time.h>
#define R_M  0x7FFFU
#define CROSSOVER_RATE            0.9
#define MUTATION_RATE             0.01
#define TOTAL_OBJECT              24
#define CHROMO_LENGTH             42
#define GENE_LENGTH               4
#define MAX_ALLOWABLE_GENERATIONS 1
#define DISCRITE_COF              4
#define DECISION_ABILITY_COFF     0.9
#define RULE_THERSOLD             4
#define RANDOM_NUM      ((float)rand()/(R_M+1))

using namespace std;

void Discritization(float &input_value, float min, float max){
        float width;
        width = (max+1-min)/DISCRITE_COF;

        int l;
        for (l= 1; l <= DISCRITE_COF ; l++) {
            if (input_value >= min + (l-1)*width && input_value < min + l*width) {
                input_value = l;
                break;
            }
        }
}

void Discritization_algorithm(vector < vector <float> > & input ){

    for(int i=0;i<TOTAL_OBJECT;i++) {
		for(int j=0;j<CHROMO_LENGTH-1 ;j++) {
			//input_value = input[i][j];

			switch(j) {
				case 0:
                    Discritization(input[i][j],0,58329);
                    break;
                case 1:
                    Discritization(input[i][j],1,3);
                    break;
                case 2:
                    Discritization(input[i][j],1,70);
                    break;
                case 3:
                    Discritization(input[i][j],1,10);
                    break;
                case 4:
                    Discritization(input[i][j],0,90000000);
                    break;
                case 5:
                    Discritization(input[i][j],0,90000000);
                    break;
                case 6:
                    Discritization(input[i][j],0,1);
                    break;
                case 7:
                    Discritization(input[i][j],0,3);
                    break;
                case 8:
                    Discritization(input[i][j],0,14);
                    break;
                case 9:
                    Discritization(input[i][j],0,101);
                    break;
                case 10:
                    Discritization(input[i][j],0,5);
                    break;
                case 11:
                    Discritization(input[i][j],0,1);
                    break;
                case 12:
                    Discritization(input[i][j],0,7479);
                    break;
                case 13:
                    Discritization(input[i][j],0,1);
                    break;
                case 14:
                    Discritization(input[i][j],0,1);
                    break;
                case 15:
                    Discritization(input[i][j],0,7468);
                    break;
                case 16:
                    Discritization(input[i][j],0,100);
                    break;
                case 17:
                    Discritization(input[i][j],0,5);
                    break;
                case 18:
                    Discritization(input[i][j],0,9);
                    break;
                case 19:
                    Discritization(input[i][j],0,1);
                    break;
                case 20:
                    Discritization(input[i][j],0,1);
                    break;
                case 21:
                    Discritization(input[i][j],0,1);
                    break;
                case 22:
                    Discritization(input[i][j],0,511);
                    break;
                case 23:
                    Discritization(input[i][j],0,511);
                    break;
                case 24:
                    Discritization(input[i][j],0,1);
                    break;
                case 25:
                    Discritization(input[i][j],0,1);
                    break;
                case 26:
                    Discritization(input[i][j],0,1);
                    break;
                case 27:
                    Discritization(input[i][j],0,1);
                    break;
                case 28:
                    Discritization(input[i][j],0,1);
                    break;
                case 29:
                    Discritization(input[i][j],0,1);
                    break;
                case 30:
                    Discritization(input[i][j],0,1);
                    break;
                case 31:
                    Discritization(input[i][j],0,255);
                    break;
                case 32:
                    Discritization(input[i][j],0,255);
                    break;
                case 33:
                    Discritization(input[i][j],0,1);
                    break;
                case 34:
                    Discritization(input[i][j],0,1);
                    break;
                case 35:
                    Discritization(input[i][j],0,1);
                    break;
                case 36:
                    Discritization(input[i][j],0,1);
                    break;
                case 37:
                    Discritization(input[i][j],0,1);
                    break;
                case 38:
                    Discritization(input[i][j],0,1);
                    break;
                case 39:
                    Discritization(input[i][j],0,1);
                    break;
                case 40:
                    Discritization(input[i][j],0,1);
                    break;
			}
		}
    }
}

/* build the discernibility Matrix */

void Discernibility_Matrix(vector <vector <float> > & input,vector <vector <vector <int> > > & Discern_Matrix){
    for (int i = 0; i < TOTAL_OBJECT ; i++ ) {
        for (int j = 0 ; j < i; j++){
            if (input[i][CHROMO_LENGTH-1] != input[j][CHROMO_LENGTH-1]){            
                for (int k =0 ;k < CHROMO_LENGTH-1; k++ ) {
                    if (input[i][k] != input[j][k])                         
                    Discern_Matrix[i][j].push_back(k+1);
                }
            }
        }
    }
}

/* generate_Multi_Set */

void Generate_Multi_Set(vector <vector <int> > & Discernibility_Function,vector <vector <int> > & Multi_Set){
    for (int i = 0; i < Discernibility_Function.size(); i++) {
        for (int j = 0; j < Discernibility_Function[i].size();j++){
            Multi_Set[i][Discernibility_Function[i][j]-1] = 1;      
        }
    }
}


bool cmp1(pair< int , int > &a,pair< int , int > &b)
{
    return a.second<=b.second;
}




double min(double a,double b){
    return (a>b?b:a);
}

//---------------------------------Fitness-----------------------------------------


double Fitness(vector <vector <int> > & Multi_Set, vector <int> &v){
    int size_b = 0;
    int size_Multi_Set = Multi_Set.size();
    int reduct_Check_Count=0;
    double fitness_Value;

	for (int i = 0; i < v.size(); i++ ){
        if (v[i] == 1) size_b+=1;
    }

    for (int i = 0; i < Multi_Set.size(); i++ ){
        int cnt = 0;

        for (int j = 0; j < Multi_Set[i].size(); j++ ) {
            if(Multi_Set[i][j] == v[j]){cnt =1;break;}
        }
        if (cnt==1) reduct_Check_Count++;
    }
    fitness_Value = (double)(CHROMO_LENGTH - 1-size_b)/(CHROMO_LENGTH - 1) + min(DECISION_ABILITY_COFF,(double)reduct_Check_Count/size_Multi_Set);
    cout << CHROMO_LENGTH - 1 << " " << size_b << " "<< reduct_Check_Count << " fitness " << size_Multi_Set << " " << fitness_Value << endl;
    //cout << CHROMO_LENGTH - 1 << " " << size_b << " "<< reduct_Check_Count << " fitness " << size_Multi_Set << " " << fitness_Value << endl;
    return fitness_Value;
}

//---------------------------------Print_data-----------------------------------------


void Print_data(vector <vector <int> > & Multi_Set){
    ofstream outfile1;
	outfile1.open("output_multiset_before_ga");
    for (int i=0;i<Multi_Set.size();i++) {
        for (int j = 0; j < Multi_Set[i].size(); j++ ){
        outfile1<< Multi_Set[i][j] << " ";
        cout << Multi_Set[i][j] << " ";
        }
        outfile1 << endl;
	cout << endl;
    }
}

bool cmp(const pair<vector <int>, double > a, const pair<vector <int> ,double > b)
{
    return a.second<=b.second;
}

//---------------------------------Crossover-----------------------------------------

void Crossover(vector < int > & offspring1, vector < int > &offspring2){
    if (RANDOM_NUM < CROSSOVER_RATE) {

    int crossover = (int) (RANDOM_NUM * (CHROMO_LENGTH-1));

    vector <int> t1(offspring1.size());
    vector <int> t2(offspring2.size());

    for (int i = 0; i < CHROMO_LENGTH-1; i++ ){
        if (i<= crossover){
            t1[i] = offspring1[i];
            t2[i] = offspring2[i];
        } else {
            t1[i] = offspring2[i];
            t2[i] = offspring1[i];
        }
    }
    offspring1 = t1; offspring2 = t2;
  }
}

//---------------------------------Mutate-----------------------------------------

void Mutate(vector < int > & offspring){

    for (int i=0; i<offspring.size(); i++) {
        if (RANDOM_NUM < MUTATION_RATE) {
            if (offspring[i] == 1)
                offspring[i] = 0;
            else
                offspring[i] = 1;
        }
    }
}

bool cmp3(const pair< int, double > &a, const pair< int , double > &b)
{
    return a.second<=b.second;
}

//---------------------------------Decision_Rule-----------------------------------------

void Decision_Rule( vector <vector <float> > & input,vector <pair <int , double > > &best_attr, vector <int> &final_reduct){
    ofstream outfile4,outfile5;
	outfile4.open("rules_class1");
	outfile5.open("rules_class1_binary");
	ofstream outfile9;
	outfile9.open("reduct_attribute_index_class1");

    for (int i =1; i <=RULE_THERSOLD; i++ ) {
        final_reduct.push_back(best_attr[best_attr.size()-i].first);
    }
    sort(final_reduct.begin(),final_reduct.end());
    for (int i =0; i <final_reduct.size(); i++ ) {
        outfile9 << final_reduct[i] << " ";
    }
    outfile9 << endl;

    for (int i = 0; i < input.size(); i++ ) {
        outfile4 << "IF " ;
            for (int j = 1; j <= RULE_THERSOLD ;j++ ) {
                outfile5 << input[i][final_reduct[j-1]] << " ";
                if (j == RULE_THERSOLD)
                    outfile4 <<"ATTR " << final_reduct[j-1]+1<< " = "<<input[i][final_reduct[j-1]]<<" " ;
                else
                    outfile4 <<"ATTR " << final_reduct[j-1] +1<< " = "<<input[i][final_reduct[j-1]] << " AND " ;
            }
            outfile5 << input[i][CHROMO_LENGTH-1] << endl;
            outfile4 << " THEN " << input[i][CHROMO_LENGTH-1] << endl;
    }

}

//---------------------------------Decision_Rule_Generation-----------------------------------------


void Decision_Rule_Generation( vector <vector <float> > & input,vector <int> &reduct , vector <int> &final_reduct){
    int reduct_attr_cnt=0;
    vector <int> index;
    for (int i=0; i < reduct.size(); i++ )
        if(reduct[i] == 1){
            reduct_attr_cnt++;
            index.push_back(i);
        }

    vector < pair <int , double > > best_attr;
    vector < vector < int > > store(reduct_attr_cnt,vector <int> (2*DISCRITE_COF));
    vector < vector < int > > store_cnt(reduct_attr_cnt,vector <int> (DISCRITE_COF));

    for (int i=0; i < input.size(); i++ ) {
        for (int j =0 ; j < index.size(); j++) {
            for (int l = 1; l <= DISCRITE_COF ; l++ ){
                if (input[i][index[j]] == l  ) {

                    store_cnt[j][l-1] +=1;
                    if ( input[i][CHROMO_LENGTH-1] == 1 ){
                        store[j][l-1]+=1;
                    } else {
                        store[j][l - 1 + DISCRITE_COF] +=1;
                    }

                }
            }
        }
    }

    double sum,total_instance;
    for (int i = 0; i < reduct_attr_cnt; i++ ) {
        sum = 0.0;
        total_instance = 0;
        for (int j = 0; j < DISCRITE_COF; j++ ) {

            if (store[i][j] > 0){
                total_instance++;
                sum += ((double)store[i][j]/store_cnt[i][j])*100;
            }
            if (store[i][j+DISCRITE_COF] > 0){
                total_instance++;
                sum +=  ((double)store[i][j+DISCRITE_COF]/store_cnt[i][j])*100;
            }
        }
       // cout << "sum = " << sum << "instance = " << total_instance << endl;
        best_attr.push_back(make_pair(index[i],sum/total_instance));
    }
    sort(best_attr.begin(),best_attr.end(),cmp3);
    cout << "\n\nbest attribute with strength values\n\n";
     for (int i = 0; i < reduct_attr_cnt; i++ )
        cout << "Attributet No. " << best_attr[i].first<<"  strength value  "<< best_attr[i].second << endl;

    Decision_Rule(input,best_attr,final_reduct);


}

//---------------------------------Genetic_Algorithm-----------------------------------------


void Genetic_Algorithm(vector <vector <int> > & Multi_Set,vector <vector <float> > & input , vector <int> &final_reduct){
    double ft;
    bool flag = true;
    int GenerationsRequiredToFindASolution = 0; 
	

	int rb = 0;
    while(flag){
	cout << " ** loop "<< ++rb << endl;
        vector < pair <double, vector <int> > > vec;
        vector <vector <int> > Multi_Set_copy(Multi_Set.size(),vector<int> (Multi_Set[0].size()));
        vector <int> offspring1;
        vector <int>  offspring2;


	    //cout <<"checking1  at" << rb << endl;
        for (int i = 0 ; i < Multi_Set.size(); i++ ){
            ft = Fitness(Multi_Set,Multi_Set[i]);
            //vec.push_back(make_pair(Multi_Set[i],ft));
            vec.push_back(make_pair(ft,Multi_Set[i]));
        }
	    cout <<"checking1  at" << rb << endl;
        sort(vec.begin(),vec.end());
	   
	  cout <<"checking1 aqq at" << rb << endl;

        if (GenerationsRequiredToFindASolution == MAX_ALLOWABLE_GENERATIONS){
            flag = false;
	    cout <<"checking2  at" << rb << endl;
            Decision_Rule_Generation(input,vec[0].second,final_reduct);
	    cout <<"checking3  at" << rb << endl;
            break;
        }else GenerationsRequiredToFindASolution++;

cout << "selected off1 - " <<vec[vec.size()-1].first <<"  selected off2 - " << vec[vec.size()-2].first<< endl;
        offspring1 = vec[vec.size()-1].second;
        offspring2 = vec[vec.size()-2].second;

        Crossover(offspring1, offspring2);

        Mutate(offspring1);
        Mutate(offspring2);

       
        Multi_Set_copy.clear();
        Multi_Set_copy.push_back(offspring1);
        Multi_Set_copy.push_back(offspring2);

        for (int i =2; i < vec.size(); i++ ){
            Multi_Set_copy.push_back(vec[i].second);
        }

        Multi_Set.clear();
        for (int i = 0; i < Multi_Set_copy.size(); i++ ){
            Multi_Set.push_back(Multi_Set_copy[i]);

        }


    }
}

//---------------------------------MAIN-----------------------------------------


int main()
{
	ifstream infile;
	infile.open("training_data_for_class1_dos_raw.txt");
	ofstream outfile,outfile1,outfile3;
	outfile.open("output_multiset_after_ga");
	outfile1.open("output_multiset_before_ga");
	outfile3.open("output_input_after_discretization");

/*traning data reading */
    vector < vector <float> > input (TOTAL_OBJECT,vector<float> (CHROMO_LENGTH));
    for (int i=0;i<TOTAL_OBJECT;i++) {
        for (int j=0;j<CHROMO_LENGTH;j++) {
            infile >> input[i][j];
        }
   }
/*traning data printing */
 
   cout << "input data\n";
    for (int i=0;i<TOTAL_OBJECT;i++) {
        for (int j=0;j<CHROMO_LENGTH;j++) {
            cout << input[i][j] << " ";
        }
        cout << endl;
   }

    vector <int> final_reduct;


    Discritization_algorithm(input); 

/*output after Discrization stored into outfile3 */

    for (int i=0;i<TOTAL_OBJECT;i++){
        for (int j=0;j<CHROMO_LENGTH;j++)
        outfile3 << input[i][j] << " ";
        outfile3 << endl;
    }

    vector < vector < vector <int> > > Discern_Matrix(TOTAL_OBJECT,vector < vector <int> >(TOTAL_OBJECT,vector<int> () ));
    vector <vector <int > > Discernibility_Function;

    Discernibility_Matrix(input,Discern_Matrix);

//Function to calculate Discern Matrix
//****wrote**********************************************//
    cout << " Discernibility_Matrix\n";
    for(int i=0;i<TOTAL_OBJECT;i++) {
            for (int j=0;j<TOTAL_OBJECT;j++) {
		if(Discern_Matrix[i][j].size() > 0) {
			cout << "block["<<i <<"]["<<j<< "] => "; 
               	 	for(int k = 0;k < Discern_Matrix[i][j].size();k++){
                       	 	cout << Discern_Matrix[i][j][k] << " ";
                	}
                	cout << endl;
		}
            }
        }

/* push to vector */
    for (int i=0;i<TOTAL_OBJECT;i++) {
        for (int j=0;j<TOTAL_OBJECT;j++) {
            if(Discern_Matrix[i][j].size()){
                    Discernibility_Function.push_back(Discern_Matrix[i][j]);
            }
        }
    }

    cout << "Discernibility_Function is \n\n";

    for (int i=0;i<Discernibility_Function.size();i++) {
        for (int j=0;j<Discernibility_Function[i].size();j++) {
                cout << Discernibility_Function[i][j] << " ";
        }
        cout << endl;
    }
    cout <<  "\n\n";

/****************************************************************************/
    int size = Discernibility_Function.size();

    vector< vector <int> > Multi_Set(size,vector<int>(CHROMO_LENGTH-1,0));

    Generate_Multi_Set(Discernibility_Function,Multi_Set);

    Print_data(Multi_Set);

    Genetic_Algorithm(Multi_Set,input,final_reduct);
    
   // Rule_Testing(width_attr,final_reduct);

    for (int i=0;i<Multi_Set.size();i++) {
        for (int j = 0; j < Multi_Set[i].size(); j++ ){
        outfile<< Multi_Set[i][j] << " ";
        }
        outfile << endl;
    }
   //cout <<Discernibility_Function.size()<<endl;

	return 0;
}
