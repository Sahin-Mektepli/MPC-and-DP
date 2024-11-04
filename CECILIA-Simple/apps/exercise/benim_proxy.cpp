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
#include <fstream>
#include <sstream>

#include "../../booleancore/core.h" //bu bool için lazım oldu

constexpr uint64_t sz = 4;

using namespace std;


void Regression(Party *const proxy, uint64_t **X, uint64_t **X_test,uint64_t *y, uint64_t *y_test, int rows, int columns,int rows_test, int iterations){
    double learning_rate = 0.1;
    //int rows = 3; // r
    //int columns = 2; // n
    uint32_t param[2] = {(uint32_t)rows, (uint32_t)columns};
    uint32_t param_test[2] = {(uint32_t)rows_test, (uint32_t)columns};

    //bu yukarıdakiler ve epsilon değeri ve learning rate elde olduktan sonra çağırmala başlayacak.
    
    vector<double> epsilons = {0.01, 0.05, 0.1, 0.2,0.5,0.7,1.0}; //epsilonların kullanmaya alışkın olduğumuz bir konfigürasyonu
    epsilons = {0, 0.5, 1};
    for(auto epsilon: epsilons){
        vector<double> accuracies = {}; //python kodunu takip ederek ilerliyorum.
        for(int i=0; i<1; ++i){ //makalede 100 olarak verilen "run time"
            //model.fit() FONKSİYONUNU AYRI YAZMAK NE KADAR GEREKLİ BİLMİYORUM, ŞİMDİLİK BU ALT KISMA YAZMAKLA YETİNECEĞİM.

            //başlangıç thetalarının tayini de bir mesele tabii.
            uint64_t *thetas = new uint64_t[columns]; 
            //thetaları rastgele oluşturmak için şu da mümkün:
            // random_device rd{};
            // mt19937 gen{rd()};
            // normal_distribution d{.0,0.1}; //mean, std_dev
            // double theta1 = d(gen);
            // double theta2 = d(gen);  BU RASTGELE ATAMALAR ASIL İMPLEMANTASYONDA TUTULAN TARİK OLACAKTIR, ŞU ANSA TUTARLILIK İÇİN TERK EDİLDİLER.
            for(int j=0; j<columns; ++j){
                thetas[j] = proxy->CreateShare(.0); //bu yavaş la!!!!!!!!!!!!!!!!!!!!!!!
            }

            proxy->SendBytes(lgGradient, param, 2); //NORMALİ BU, TÜRKMEN İÇİN ÇIKARDIM
            //proxy->SendBytes(lgGradientTurkmen, param, 2);

            uint64_t *new_thetas = GradientDescent(proxy, X, y, thetas, rows, columns, learning_rate, epsilon, iterations); //X_train ve y_train ile NORMALİ BU TÜRKMEN İÇİN ÇIKARDIM.
            //uint64_t *new_thetas = GradientDescent_Turkmen(proxy, X, y, thetas, rows, columns, iterations); //Türkmen'den kopya.
            //thetas =  GradientDescent(proxy, X, y, thetas, rows, columns, 0.5, epsilon, iterations); //X_train ve y_train ile

            cout << "************************************" << endl;

            uint64_t *new_theta_double = Reconstruct(proxy, new_thetas, columns);
            for(int j=0; j<columns; ++j){
                cout << ConvertToDouble(new_theta_double[j]) << endl;
            }


            //model.fit() tamam, yeni theta değerleriyle önsöylemeye başlıyoruz.
            //tek satır olarak predict() çalışıyor 

            proxy->SendBytes(lgPredict, param_test, 2);
            
            uint64_t *prediction = Predict(proxy, X_test, new_thetas, rows_test, columns); //bu tabii ki X_test ile çalışmalı


            //accuracy'yi hesap için elde olduğu varsayılan bir y_test array'i kullanılacak.

            proxy->SendBytes(equals, param, 1);            
            uint64_t *equals_shares = Equals(proxy, prediction, y_test, rows);
            uint64_t *equals = Reconstruct(proxy, equals_shares, rows);
            double acc = 0;
            for(int j=0; j<rows_test; ++j){
                if(equals[j]){
                    ++acc;
                }
            }
            acc /= rows_test;
            accuracies.push_back(acc);
            //çirkin görünüyor, bunu toparlamak lazım.

            // cout << ConvertToDouble(Reconstruct(proxy, prediction[0])) << endl;
            // cout << ConvertToDouble(Reconstruct(proxy, prediction[1])) << endl;
            // cout << ConvertToDouble(Reconstruct(proxy, prediction[2])) << endl;
            // cout << ConvertToDouble(Reconstruct(proxy, prediction[3])) << endl;
            // cout << ConvertToDouble(Reconstruct(proxy, prediction[4])) << endl;
            // cout << ConvertToDouble(Reconstruct(proxy, prediction[5])) << endl;
            // cout << ConvertToDouble(Reconstruct(proxy, prediction[6])) << endl;
            // cout << ConvertToDouble(Reconstruct(proxy, prediction[7])) << endl;
            // cout << ConvertToDouble(Reconstruct(proxy, prediction[8])) << endl; //bu çıktıların anlamı azalmakta
            // cout << ConvertToDouble(Reconstruct(proxy, prediction[9])) << endl; //bu çıktıların anlamı azalmakta
            // cout << ConvertToDouble(Reconstruct(proxy, prediction[10])) << endl; //bu çıktıların anlamı azalmakta
           
        }
        double acc_sum = 0;
        for(auto x: accuracies){
            acc_sum += x;
        }
        double acc = acc_sum/10; 
        cout << "for epsilon: "<< epsilon <<" the accuracy is: " << acc << endl;

    }
}


int main(int argc, char* argv[]) {

    uint8_t role = atoi(argv[1]);
    uint16_t cport = atoi(argv[2]);
    string caddress(argv[3]);
    uint16_t hport = atoi(argv[4]);
    string haddress(argv[5]);


    Party *proxy; //direkt kopyaladım demodan.
    if (role == 0)
        proxy = new Party(proxy1, hport, haddress, cport, caddress);
    else
        proxy = new Party(proxy2, hport, haddress, cport, caddress);

    

    std::ifstream infile; 
    //infile.open("/Users/dulumrae/Desktop/preprocessed_data.txt");  //sayın kullanıcı, bunu lütfen değiştirin.
    infile.open("/Users/dulumrae/Desktop/pre.txt");  //sayın kullanıcı, bunu lütfen değiştirin.


    if (!infile.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        
        std::cerr << "Error code: " << errno << " - " << std::strerror(errno) << std::endl; 

        return 1;
    }

    const int numRows = 30162; // number of rows
    const int numCols = 15;  // number of columns (target column dahil)

    double** X = new double*[numRows]; // Matrix X
    for (int i = 0; i < numRows; ++i) {
        X[i] = new double[numCols - 1]; 
    }
    double* y = new double[numRows]; // vector y
    std::string line;
    int row = 0;

    while (std::getline(infile, line) && row < numRows) {
        std::stringstream ss(line);
        std::string token;
        int col = 0;

        
        while (std::getline(ss, token, ',') && col < numCols - 1) { // , ile ayrılıyorlar
            X[row][col] = std::stod(token); // Convert string to double
            col++;
        }

        // target için bu kısım (sondaki sütun)
        if (col == numCols - 1 && std::getline(ss, token, ',')) {
            y[row] = std::stod(token);
        }

        row++;
    }

    infile.close();

  


//test ve train olarak bölüyorum buradan sonra
    int num_samples = numRows;
    int num_features = (numCols-1);

    const int train_size = static_cast<int>(num_samples * 0.7);//0.7 olmalı aslında ama 0.7 yapınca acc=0.3 geldi
    const int test_size = num_samples - train_size;
    cout<<"test size:"<<test_size<<endl;
    cout<<"train size:"<<train_size<<endl;

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
    }
//burası çıktı alıp test etmek için
    /*// Output the arrays for verification
    std::cout << "X_train:" << std::endl;
    for (int i = 0; i < train_size; ++i) {
        for (int j = 0; j < num_features; ++j) {
            std::cout << X_train[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "y_train:" << std::endl;
    for (int i = 0; i < train_size; ++i) {
        std::cout << y_train[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "X_test:" << std::endl;
    for (int i = 0; i < test_size; ++i) {
        for (int j = 0; j < num_features; ++j) {
            std::cout << X_test[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "y_test:" << std::endl;
    for (int i = 0; i < test_size; ++i) {
        std::cout << y_test[i] << " ";
    }
    std::cout << std::endl;

   */
    //!!!!numCols=14 bunu senin değiştirmen gerekebilir sanki sen farklı yazmıştın
uint64_t** X_te=proxy->CreateShare(X_test,test_size,numCols-1);
uint64_t* y_te=proxy ->CreateShare(y_test,test_size);
uint64_t** X_tr=proxy->CreateShare(X_train,train_size,numCols-1);
uint64_t* y_tr=proxy ->CreateShare(y_train,train_size);


Regression(proxy, X_tr, X_te,y_tr, y_te, train_size,numCols-1, test_size, 3); //!
  
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


    ///////////////////////////
/*
int rows = 3;
int columns = 3;

double thetas[columns];
thetas[0] = 0;
thetas[1] = 0;
thetas[2] = 0;
uint64_t *theta_shares = proxy->CreateShare(thetas,columns);

double ys[rows];
ys[0] = 13;
ys[1] = 24;
ys[2] = 36;
uint64_t *y_shares = proxy->CreateShare(ys,columns);

double **Xs = new double*[rows]; //r rows
Xs[0] = new double[columns]();
Xs[1] = new double[columns]();
Xs[2] = new double[columns]();


Xs[0][0] = 12; 
Xs[0][1] = 10; 
Xs[0][2] = 1;

Xs[1][0] = 16; 
Xs[1][1] = 15; 
Xs[1][2] = 17;

Xs[2][0] = 12; 
Xs[2][1] = 18; 
Xs[2][2] = 16;


uint64_t **X_shares = proxy->CreateShare(Xs, rows, columns);

uint32_t param[2] = {(uint32_t)rows, (uint32_t)columns};


proxy->SendBytes(lgGradientTurkmen, param, 2);
auto result = GradientDescent_Turkmen(proxy, X_shares, y_shares, theta_shares, rows,columns,2);

cout << "**************************" << endl;
for(int i=0; i<columns; ++i){
    cout << ConvertToDouble(Reconstruct(proxy, result[i])) << endl;;
}  

*/
    
   proxy->SendBytes(coreEnd);





    return 0;
}
