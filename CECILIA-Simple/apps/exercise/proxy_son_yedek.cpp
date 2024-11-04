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





void Regression(Party *const proxy, uint64_t **X, uint64_t *y, uint64_t *y_test, int rows, int columns, int iterations){
    
    double learning_rate = 0.01; //bunu elle yazalım, nolcak yani.

    uint32_t param[2] = {(uint32_t)rows, (uint32_t)columns};
    
    vector<double> epsilons = {1}; //epsilonların kullanmaya alışkın olduğumuz bir konfigürasyonu
    for(auto epsilon: epsilons){
        vector<double> accuracies = {}; //python kodunu takip ederek ilerliyorum.
        for(int i=0; i<1; ++i){ //makalede 100 olarak verilen "run time"

            //başlangıç thetalarının tayini de bir mesele tabii.
            uint64_t *thetas = new uint64_t[columns]; 
            //thetaları rastgele oluşturmak için şu da mümkün:
            // random_device rd{};
            // mt19937 gen{rd()};
            // normal_distribution d{.0,0.1}; //mean, std_dev
            // double theta1 = d(gen);
            // double theta2 = d(gen);  BU RASTGELE ATAMALAR ASIL İMPLEMANTASYONDA TUTULAN TARİK OLACAKTIR, ŞU ANSA TUTARLILIK İÇİN TERK EDİLDİLER.
            for(int j=0; j<columns; ++j){
                thetas[j] = proxy->CreateShare(.0); //bu daha mı yavaş acaba bir double[colums] = [0,0,...]'den.
            }

            proxy->SendBytes(lgGradient, param, 2);

            uint64_t *new_thetas = GradientDescent(proxy, X, y, thetas, rows, columns, learning_rate, epsilon, iterations); //X_train ve y_train ile

            cout << "NEW_THETAS************************************" << endl;
            uint64_t *new_theta_double = Reconstruct(proxy, new_thetas, columns);
            for(int j=0; j<columns; ++j){
                cout << ConvertToDouble(new_theta_double[j]) << endl;
            }
            cout << "NEW_THETAS************************************" << endl;


            //model.fit() tamam, yeni theta değerleriyle önsöylemeye başlıyoruz.
            cout << "gradient'tan çıktık" << endl;
            //tek satır olarak predict() çalışıyor 
            proxy->SendBytes(lgPredict, param, 2);
            
            uint64_t *prediction = Predict(proxy, X, new_thetas, rows, columns); //bu tabii ki X_test ile çalışmalı


            //accuracy'yi hesap için elde olduğu varsayılan bir y_test array'i kullanılacak.

            proxy->SendBytes(equals, param, 1);            
            uint64_t *equals_shares = Equals(proxy, prediction, y_test, rows);
            uint64_t *equals = Reconstruct(proxy, equals_shares, rows);
            double acc = 0;
            for(int j=0; j<rows; ++j){
                if(equals[j]){
                    ++acc;
                }
            }
            acc /= rows;
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
            cout << "the accuracy is: " << acc << endl;
        }

    }
}


int main(int argc, char* argv[]) {

    uint8_t role = atoi(argv[1]);
    uint16_t cport = atoi(argv[2]);
    string caddress(argv[3]);
    uint16_t hport = atoi(argv[4]);
    string haddress(argv[5]);

    const int numRows = 30162; // number of rows
    const int numCols = 15;  // number of columns (target column dahil)
    double** X_read = new double*[numRows]; // Matrix X
    for (int i = 0; i < numRows; ++i)   {X_read[i] = new double[numCols - 1]; }
    double* y_read = new double[numRows]; // vector y


    std::ifstream infile; 
    infile.open("/Users/dulumrae/Desktop/preprocessed_data.txt");  //sayın kullanıcı, bunu lütfen değiştirin.

    if (!infile.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        std::cerr << "Error code: " << errno << " - " << std::strerror(errno) << std::endl; 
        return 1;
    }
    
    std::string line;
    int row = 0;
    while (std::getline(infile, line) && row < numRows) {
        std::stringstream ss(line);
        std::string token;
        int col = 0;

        while (std::getline(ss, token, ',') && col < numCols - 1) { // , ile ayrılıyorlar
            X_read[row][col] = std::stod(token); // Convert string to double
            col++;
        }
        // target için bu kısım (sondaki sütun)
        if (col == numCols - 1 && std::getline(ss, token, ',')) {
            y_read[row] = std::stod(token);
        }
        row++;
    }
    infile.close();

 


    Party *proxy; //direkt kopyaladım demodan.
    if (role == 0)
        proxy = new Party(proxy1, hport, haddress, cport, caddress);
    else
        proxy = new Party(proxy2, hport, haddress, cport, caddress);

    
/* **************************************************************
    BU AŞAĞIDAKİ KISIM KÜÇÜK MATRİSLERLE DENEMEK İÇİN HAZIRLANDI
 */

    int row_number = numRows;
    int column_number = numCols-1; //14 olacak yav bu.
    uint64_t **X = new uint64_t*[row_number]; //r rows
    X[0] = new uint64_t[2]();
    X[1] = new uint64_t[2]();
    X[2] = new uint64_t[2]();

    X[0][0] = proxy->CreateShare(1.0);
    X[0][1] = proxy->CreateShare(2.0);
    X[1][0] = proxy->CreateShare(3.0);
    X[1][1] = proxy->CreateShare(4.0);
    X[2][0] = proxy->CreateShare(2.0);
    X[2][1] = proxy->CreateShare(1.0);

    uint64_t *y = new uint64_t[row_number]; 
    y[0] = proxy->CreateShare(.0);
    y[1] = proxy->CreateShare(1.0);
    y[2] = proxy->CreateShare(0.0);

    uint64_t *y_test = new uint64_t[row_number];  //test için oluşturdum.
    y_test[0] = proxy->CreateShare(1.0);
    y_test[1] = proxy->CreateShare(0.0);
    y_test[2] = proxy->CreateShare(0.0);


    uint64_t **X_read_shares = new uint64_t*[row_number];
    for(int i=0; i<row_number; i++){
        X_read_shares[i] = proxy->CreateShare(X_read[i], column_number); 
    }
    uint64_t *y_read_shares = proxy->CreateShare(y_read, row_number); 
    double ones[row_number];
    for(int i=0; i<row_number; i++){
        ones[i] = 1;
    }
    uint64_t *y_test_shares = proxy->CreateShare(ones, row_number); //tamamen uydurdum :p ASLA BÖYLE OLMAMALI ASLINDA!!


    Regression(proxy, X_read_shares, y_read_shares, y_read_shares, row_number, column_number, 1); //

    
    
    int rr = 3;
    int nn = 2;
    //Regression(proxy, X, y, y_test, rr, nn, 10); // DENEME İÇİN KÜÇÜK MATRİSLERLE ÇALIŞAN REGRESYON.


    return 0;
}