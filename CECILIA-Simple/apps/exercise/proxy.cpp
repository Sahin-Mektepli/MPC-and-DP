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


constexpr uint64_t sz = 4;




double generateLaplaceNoise_old(Party *const proxy,double scale) {
    //scale double mı olmalı yoksa uint64_t mi olmalı emin olamadım şu an kullanacağımız yere göre değişebilir
    //orijinal halini test ettim 0 ile bir arasında sayılar döndürüyor ama şu anki hali öyle değil
    
    
    
    //NORMALDE ALTTAKİ SATIRLARI KULLANIYORDUK RANDOM SAYI ÜRETMEK İÇİN
    std::random_device rd; //bu seed oluşturmak için
    std::mt19937 gen(rd()); //bu o seed'i kullanarak random bir sayı olışturuyor
    

    std::uniform_real_distribution<> dis(0.0, 1.0);
    
   
    double u = dis(gen) - 0.5;//bu da o random sayıyı kullanıyor
    //auto n= proxy->GenerateCommonRandom();//onlarınkini kullandım burada random sayı oluşturmak için
    //ama ufak bir sorun var bu double değil de uint64 döndürüyor biz bunu nerede çağıracaksak ona göre değiştirmeliyiz
 
    //double u = dis(n) - 0.5;
    // Apply the inverse CDF of the Laplace distribution
    double noise = scale * ((u >= 0) ? -std::log(1 - 2*u) : std::log(1 + 2*u));
    
    return noise;
}


using namespace std;

void Regression(Party *const proxy, uint64_t **X, uint64_t **X_test,uint64_t *y, uint64_t *y_test, int rows, int columns,int rows_test, int iterations){
    double learning_rate = 0.01;//0.01 idi
    //int rows = 3; // r
    //int columns = 2; // n
    uint32_t param[2] = {(uint32_t)rows, (uint32_t)columns};
    uint32_t param_test[2] = {(uint32_t)rows_test, (uint32_t)columns};


    

    

    //bu yukarıdakiler ve epsilon değeri ve learning rate elde olduktan sonra çağırmala başlayacak.
    
    vector<double> epsilons = {1.0}; //epsilonların kullanmaya alışkın olduğumuz bir konfigürasyonu
    //vector<double> epsilons = {0.01, 0.05, 0.1, 0.2,0.5,0.7,1.0};
    int run_time = 1;
    for(auto epsilon: epsilons){
        double total_time_taken = 0; //to follow how long Gradient_Descent takes in total.


        vector<double> accuracies = {}; //python kodunu takip ederek ilerliyorum.
        for(int i=0; i<run_time; ++i){ //makalede 100 olarak verilen "run time"
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
                thetas[j] = proxy->CreateShare(.0); //bu yavaş la!!!!!!!!!!!!!!!!!!!!!!! BUNU SIFIR'A DÖNDÜR, ÇEVİR, DEBUG
            }

            std::cout << "gradient'ı çağırıyoruz" << endl;


            auto start = chrono::high_resolution_clock::now();

            uint64_t *new_thetas = GradientDescent(proxy, X, y, thetas, rows, columns, learning_rate, epsilon, iterations); //X_train ve y_train ile
            //uint64_t *new_thetas = thetas;
           
           //uint64_t *new_thetas = GradientDescent_Turkmen(proxy, X, y, thetas, rows, columns, iterations); //X_train ve y_train ile

            auto end = chrono::high_resolution_clock::now();
            double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
            time_taken *= 1e-9;
            total_time_taken += time_taken; //so I sum up the durations of GradientDescent

            
            std::cout << "thetas:\n";

            uint64_t *new_theta_double = Reconstruct(proxy, new_thetas, columns);
            for(int j=0; j<columns; ++j){
                cout << ConvertToDouble(new_theta_double[j]) << "\n"; //bu ne la
            } 




            proxy->SendBytes(lgPredict, param_test, 2);            
            uint64_t *prediction = Predict(proxy, X_test, new_thetas, rows_test, columns); //bu tabii ki X_test ile çalışmalı


            //accuracy'yi hesap için elde olduğu varsayılan bir y_test array'i kullanılacak.

            //proxy->SendBytes(equals, param, 1);            
            //uint64_t *equals_shares = Equals(proxy, prediction, y_test, rows); BUNLARIN ROW-ROW_TRAIN SORUNU VAR
            //uint64_t *equals = Reconstruct(proxy, equals_shares, rows);

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
            //çirkin görünüyor, bunu toparlamak lazım.

            /* cout << ConvertToDouble(Reconstruct(proxy, prediction[0])) << endl;
            cout << ConvertToDouble(Reconstruct(proxy, prediction[1])) << endl;
            cout << ConvertToDouble(Reconstruct(proxy, prediction[2])) << endl;
            cout << ConvertToDouble(Reconstruct(proxy, prediction[3])) << endl;
            cout << ConvertToDouble(Reconstruct(proxy, prediction[4])) << endl;
            cout << ConvertToDouble(Reconstruct(proxy, prediction[5])) << endl;
            cout << ConvertToDouble(Reconstruct(proxy, prediction[6])) << endl;
            cout << ConvertToDouble(Reconstruct(proxy, prediction[7])) << endl;
            cout << ConvertToDouble(Reconstruct(proxy, prediction[8])) << endl; //bu çıktıların anlamı azalmakta
            cout << ConvertToDouble(Reconstruct(proxy, prediction[9])) << endl; //bu çıktıların anlamı azalmakta
            cout << ConvertToDouble(Reconstruct(proxy, prediction[10])) << endl; //bu çıktıların anlamı azalmakta
            cout << "the accuracy is: " << acc << endl; */
        }
        double acc_sum = 0;
        for(auto x: accuracies){
            acc_sum += x;
        }
        double acc = acc_sum/run_time; 
        std::cout << "for epsilon: "<< epsilon <<" the accuracy is: " << acc << endl;
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
 
    


    uint32_t params[1] = {(uint32_t)4};//neden 32 bit??


/*     uint64_t z_shared = proxy->CreateShare(1.75);
    proxy->SendBytes(lgSingleSigmoid);
    uint64_t s=Sigmoid(proxy,z_shared);
    double sr=ConvertToDouble(Reconstruct(proxy,s));
    cout<<"Sigmoid:"<<sr<<endl;

    
    params[0]=4;
    double v1[4] = {1.75, 2.0, 1.0, 0.75};
    uint64_t* v1_secret_shared = proxy->CreateShare(v1, 4);
    proxy->SendBytes(lgSigmoid,params,1);
    uint64_t* s1=Sigmoid(proxy,v1_secret_shared,4);
    double* rv=ConvertToDouble(Reconstruct(proxy, s1, 4), 4);
    for(int i=0;i<4;i++){
        cout<<rv[i]<<endl;
    } */
    
    



 
  std::ifstream infile; 
  infile.open("/Users/dulumrae/Downloads/random_dataset_30000.csv");  

    if (!infile.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        
        std::cerr << "Error code: " << errno << " - " << std::strerror(errno) << std::endl; 

        return 1;
    }
    
    //const int numRows = 30162; // number of rows
    //const int numCols = 15;  // number of columns (target column dahil) 16 idi

    //const int numRows = 101766; // number of rows
    //const int numCols = 50;
    
    const int numRows = 50000; // bu mütegayyir olacak.
    const int numCols = 140;

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

  


//test ve train olarak bölüyorum buradan sonra
    int num_samples = numRows;
    int num_features = (numCols-1);

    const int train_size = static_cast<int>(num_samples * 0.7);//0.7 olmalı aslında ama 0.7 yapınca acc=0.3 geldi
    const int test_size = num_samples - train_size;
    //cout<<"test size:"<<test_size<<endl;
    //cout<<"train size:"<<train_size<<endl;

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

uint64_t** X_te=proxy->CreateShare(X_test,test_size,numCols-1); ///bunu test için yaptım şu anda train test diye ayırmıyor
uint64_t* y_te=proxy ->CreateShare(y_test,test_size);
uint64_t** X_tr=proxy->CreateShare(X_train,train_size,numCols-1);
uint64_t* y_tr=proxy ->CreateShare(y_train,train_size);

// uint64_t** X_te=proxy->CreateShare(X,numRows,numCols-1);
// uint64_t* y_te=proxy ->CreateShare(y,numRows);
// uint64_t** X_tr=proxy->CreateShare(X,numRows,numCols-1);
// uint64_t* y_tr=proxy ->CreateShare(y,numRows);



auto start = chrono::high_resolution_clock::now();
Regression(proxy, X_tr, X_te,y_tr, y_te, train_size, numCols-1, test_size, 1); //!10 olamlı ben test için 1 yaptım
auto end = chrono::high_resolution_clock::now();
double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
time_taken *= 1e-9;
//std::cout << "time taken during regression in seconds: " << time_taken << "\n";
  
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


    
   proxy->SendBytes(coreEnd);





    return 0;
}

