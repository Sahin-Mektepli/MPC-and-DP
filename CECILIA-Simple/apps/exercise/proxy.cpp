//
// Created by Seyma Selcan on 16.10.23.
//
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <tuple>
#include <iomanip>
#include <bitset>
#include "../../core/core.h"
#include <fstream>//ben ekledim bunları!!!
#include <sstream>//
#include <string>//

#include <cerrno> // for errno
#include <cstring> // for strerror






using namespace std;

void Regression(Party *const proxy, uint64_t **X, uint64_t **X_test,uint64_t *y, uint64_t *y_test, int rows, int columns,int rows_test, int iterations){
    double learning_rate = 0.01;//0.01 idi
    uint32_t param[2] = {(uint32_t)rows, (uint32_t)columns};
    uint32_t param_test[2] = {(uint32_t)rows_test, (uint32_t)columns};

    //bu yukarıdakiler ve epsilon değeri ve learning rate elde olduktan sonra çağırmala başlayacak.

    //ne kadar epsilon istediğimize göre alttakilerden birini yoruma alıyoruz

    vector<double> epsilons = {1.0}; 
    //vector<double> epsilons = {0.01, 0.05, 0.1, 0.2,0.5,0.7,1.0};

    int run_time = 1;
    for(auto epsilon: epsilons){
        double total_time_taken = 0; //to follow how long Gradient_Descent takes in total.


        vector<double> accuracies = {}; //python kodunu takip ederek ilerliyorum.
        for(int i=0; i<run_time; ++i){ //makalede 100 olarak verilen "run time"

            //başlangıç thetalarının tayini de bir mesele tabii.
            uint64_t *thetas = new uint64_t[columns]; 

            double* zeros_for_initial_theta = new double[columns];
            for(int j=0; j<columns; ++j){
                zeros_for_initial_theta[j] = .0; //önce bir sıfırlar array'i
            }
            thetas = proxy->CreateShare(zeros_for_initial_theta, columns);
            //Proxylere doğrudan 0 vermek aklıma yatmıyor ve içime sinmiyor
            //Bu "doğru" yaklaşım olmayabilir...

            std::cout << "gradient'ı çağırıyoruz" << endl; //some debugging nonsense


            auto start = chrono::high_resolution_clock::now();

            uint64_t *new_thetas = GradientDescent(proxy, X, y, thetas, rows, columns, learning_rate, epsilon, iterations); //X_train ve y_train ile

            auto end = chrono::high_resolution_clock::now();
            double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
            time_taken *= 1e-9;
            total_time_taken += time_taken; //so I sum up the durations of GradientDescent

            
            std::cout << "thetas:\n"; //some debugging nonsense

            uint64_t *new_theta_double = Reconstruct(proxy, new_thetas, columns);

            proxy->SendBytes(lgPredict, param_test, 2);            
            uint64_t *prediction = Predict(proxy, X_test, new_thetas, rows_test, columns); //bu tabii ki X_test ile çalışmalı

            proxy->SendBytes(equals, param_test, 1);            
            uint64_t *equals_shares = Equals(proxy, prediction, y_test, rows_test);
            uint64_t *equals = Reconstruct(proxy, equals_shares, rows_test);
            double acc = 0;
            for(int j=0; j<rows_test; ++j){
                if(equals[j]){
                    ++acc;
                }
            }
            acc /= rows_test;
            accuracies.push_back(acc);
        }
        double acc_sum = 0;
        for(auto x: accuracies){
            acc_sum += x;
        }
        double acc = acc_sum/run_time; 
        std::cout << "for epsilon: "<< epsilon <<" the accuracy is: " << acc << "\n";
        std::cout << "total time taken during descents for this e in seconds: " << total_time_taken << "\n";
    }
}

int main(int argc, char* argv[]) {

   

    uint8_t role = atoi(argv[1]);
    uint16_t cport = atoi(argv[2]);
    string caddress(argv[3]);
    uint16_t hport = atoi(argv[4]);
    string haddress(argv[5]);



    Party *proxy;
    if (role == 0)
        proxy = new Party(proxy1, hport, haddress, cport, caddress);
    else
        proxy = new Party(proxy2, hport, haddress, cport, caddress);
 


    uint32_t params[1] = {(uint32_t)4};//neden 32 bit?? doesn't really matter now does it?
    



 
  std::ifstream infile; 
  //infile.open("/Users/dulumrae/Downloads/geçiciCecilia/MPC-and-DP/CECILIA-Simple/apps/exercise/datasets/random_dataset_10000.csv");  
  infile.open("/Users/dulumrae/Downloads/datasets/diabet.txt");  
    std::cout << "diabettt\n"; 
    if (!infile.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        
        std::cerr << "Error code: " << errno << " - " << std::strerror(errno) << std::endl; 

        return 1;
    }
    
    //random_dataset'ler için bunları kullanmıyoruz.
    //const int numRows = 30162; // number of rows
    //const int numCols = 15;  // number of columns (target column dahil) 16 idi

    const int numRows = 101766; // number of rows
    const int numCols = 50;
    
    //const int numRows = 10000; // satır sayısı.
    //const int numCols = 250;   // feautre sayısı 144

    double** X = new double*[numRows]; // Matrix X
    for (int i = 0; i < numRows; ++i) {
        X[i] = new double[numCols - 1]; 
    }
    double* y = new double[numRows]; // vector y
    std::string line;
    int row = 0; //works like a for loop

    while (std::getline(infile, line) && row < numRows) {
        std::stringstream ss(line);
        std::string token;
        int col = 0;

        
        while (std::getline(ss, token, ',') && col < numCols - 1) { // , ile ayrılıyorlar
            if(!token.compare("")){token = "0";}

            X[row][col] = std::stod(token); // Convert string to double

            col++;
        }

        // target için bu kısım (sondaki sütun)
        if (col == numCols - 1 && std::getline(ss, token, ',')) {
            if(!token.compare("")){token = "0";}

            y[row] = std::stod(token);
        }

        row++;
    }

    infile.close();

  

//test ve train olarak bölüyoruz buradan sonra
    int num_samples = numRows;
    int num_features = (numCols-1);

    const int train_size = static_cast<int>(num_samples * 0.7);
    const int test_size = num_samples - train_size;

    // Create arrays for X_train, X_test, y_train, y_test
    double** X_train = new double*[train_size];
    double** X_test = new double*[test_size];
    double* y_train = new double[train_size];
    double* y_test = new double[test_size];

    for (int i = 0; i < train_size; ++i) {
        X_train[i] = new double[num_features];
    }
    for (int i = 0; i < test_size; ++i) {
        X_test[i] = new double[num_features];
    }

    // Fill X_train and y_train arrays
    for (int i = 0; i < train_size; ++i) {
        for (int j = 0; j < num_features; ++j) {
            X_train[i][j] = X[i][j];
        }
        y_train[i] = y[i];
    }

    // Fill X_test and y_test arrays
    for (int i = train_size; i < num_samples; ++i) {
        for (int j = 0; j < num_features; ++j) {
            X_test[i - train_size][j] = X[i][j];
        }
        y_test[i - train_size] = y[i];
        //std::cout << y[i] << " ";
    }
   
uint64_t** X_te=proxy->CreateShare(X_test,test_size,numCols-1); 
uint64_t* y_te=proxy ->CreateShare(y_test,test_size);
uint64_t** X_tr=proxy->CreateShare(X_train,train_size,numCols-1);
uint64_t* y_tr=proxy ->CreateShare(y_train,train_size);




//auto start = chrono::high_resolution_clock::now();
Regression(proxy, X_tr, X_te,y_tr, y_te, train_size, numCols-1, test_size, 1); //Bu sondaki iterations = 10 olamlı ben test için 1 yaptım
//auto end = chrono::high_resolution_clock::now();
//double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
//time_taken *= 1e-9;
  
     // Clean up dynamically allocated memory
    for (int i = 0; i < train_size; ++i) {
        delete[] X_train[i];
    }
    for (int i = 0; i < test_size; ++i) {
        delete[] X_test[i];
    }
    // Clean up
    for (int i = 0; i < numRows; ++i) {
        delete[] X[i];
    }
    delete[] X;
    delete[] y; 
    delete[] X_train;
    delete[] X_test;
    delete[] y_train;
    delete[] y_test;

    
   proxy->SendBytes(coreEnd);


    return 0;
}

