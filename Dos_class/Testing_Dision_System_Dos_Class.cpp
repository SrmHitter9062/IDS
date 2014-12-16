//-----------------------------------------project_test_main.cpp--------------------------------------
//
//  Code is use to taking object which we want to classify and after that dscretization happen and
//  using the generated rule we decide which class our new object is related

saaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa




dsdsdsdsds

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

// hello this is  dfdfdfIDS project //

#define CROSSOVER_RATE            0.9               //Crossover rate for crossover method
#define MUTATION_RATE             0.01              //mutation rate for mutation method
#define TOTAL_OBJECT              24                // no of object inlude in training purpose
#define CHROMO_LENGTH             42                //Number of attribute for a particular chromosome
#define MAX_ALLOWABLE_GENERATIONS 1                //maximum allowed no of repetation of GA
#define DISCRITE_COF              4                 //Cofficient for discretization
#define DECISION_ABILITY_COFF     0.9               //Decision cofficient
#define RULE_THERSOLD             4               //Number of attribute choosen
#define TOTAL_TEST_OBJECT		  136                 //Number f testing object

using namespace std;


//---------------------------------Discritization_algorithm-----------------------------------------
//
//  This function returns a real value attribute value as in eqaul width interwal .
//
//-----------------------------------------------------------------------------------------

void Discritization(float &input_value, float min, float max){
        float width;
        width = (max+1-min)/DISCRITE_COF;               //Calculate width for attribute j

        int l;
        for (l= 1; l <= DISCRITE_COF ; l++) {
            if (input_value >= min + (l-1)*width && input_value < min + l*width) {
                input_value = l;
                break;
            }
        }

}
void Discritization_algorithm(vector < vector <float> > & input ){
    int i,j;
    for(i=0;i<TOTAL_TEST_OBJECT;i++) {
		for(j=0;j<CHROMO_LENGTH-1 ;j++) {
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









int main(){
	 ifstream infile1,infile2,infile3,infile4;      //Declare objects for open file in reading mode
	infile1.open("DoS_Class_data250_raw.txt");              //This file contans attribute data in real value formate
	infile2.open("rules_class1_binary");            //Generated Rule by training dataset using hybrid GA
	//infile3.open("discretization_usefull_data");    //Usefull data for Dscretization of new objects
	infile4.open("reduct_attribute_index_class1");         //Optimize reduct attribute index

    ofstream outfile4,outfile5;                              //Declare objects for open file in writting mode
	outfile4.open("testing_result_class1");                //Store final classified result for testing object
outfile5.open("check1");

    vector < vector <float> > test_input (TOTAL_TEST_OBJECT,vector<float> (CHROMO_LENGTH)); //Store test input as marix formate
    vector < vector <int> > rule_matrix (TOTAL_OBJECT,vector<int> (RULE_THERSOLD + 1));     //Store generated rule in matrix formate
    vector <int> final_reduct;                                                             //Store attribute index
    vector < pair <float , float > > width_data;
    for (int i = 0; i < CHROMO_LENGTH-1; i++ ) {
        float max,min;
        infile3 >> min >> max;
        width_data.push_back(make_pair(min,max));
    }

	for (int i=0;i<TOTAL_TEST_OBJECT;i++) {
        for (int j=0;j<CHROMO_LENGTH;j++) {
            infile1 >> test_input[i][j];
		}
    }
    Discritization_algorithm(test_input);
for (int i=0;i<TOTAL_TEST_OBJECT;i++) {
        for (int j=0;j<CHROMO_LENGTH-1;j++) {
            outfile5 <<  test_input[i][j]<< " ";
        }outfile5 << endl;
}
    int x;
    for (int i = 0; i < RULE_THERSOLD; i++ ){
        infile4 >> x;
        final_reduct.push_back(x);                  //Puting index of attribute in vector
    }

    for (int i=0;i<TOTAL_OBJECT;i++) {
        for (int j=0;j<RULE_THERSOLD + 1;j++) {
            infile2 >> rule_matrix[i][j];           //Puting rules iin matrix formate
        }
    }
     for (int i=0;i<TOTAL_TEST_OBJECT;i++) {
        int flag;
        for (int j=0;j<TOTAL_OBJECT;j++) {
            flag = 0;
            for (int k = 0; k < RULE_THERSOLD; k++  ){
        cout << test_input[i][final_reduct[k]] <<" ";
                 if (test_input[i][final_reduct[k]] == rule_matrix[j][k]){  //Testing attribute value of new objects with generated rules to classify object
                    flag = 1;
                }else {
                    flag = 0;
                    break;
                }

            }
            cout <<endl;
            if (flag == 1){                          // Condition to find class for object
                outfile4 << "class of object "<< i<<" is calculated by our model is -- " << rule_matrix[j][RULE_THERSOLD]<< " and actual by dataset is " << test_input[i][CHROMO_LENGTH-1] <<endl;
                break;
            }
        }
        if (flag == 0)                              //Condition to not find class for object
            outfile4 << "No result for this object --" << endl;

  }

}


