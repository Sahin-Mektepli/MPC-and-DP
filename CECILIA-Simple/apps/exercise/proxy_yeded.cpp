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
    double learning_rate = 0.1;
    //int rows = 3; // r
    //int columns = 2; // n
    uint32_t param[2] = {(uint32_t)rows, (uint32_t)columns};

    


    //evvela X'in yani ki train matrisimizin hazırlanması gerekiyor
    // !! Kodun bu hâlinde talim ve test setleri ayrıştırılmış değil !!
    // uint64_t **X = new uint64_t*[rows]; //r rows
    // X[0] = new uint64_t[2]();
    // X[1] = new uint64_t[2]();
    // X[2] = new uint64_t[2]();

    // X[0][0] = proxy->CreateShare(1.0);
    // X[0][1] = proxy->CreateShare(2.0);
    // X[1][0] = proxy->CreateShare(3.0);
    // X[1][1] = proxy->CreateShare(4.0);
    // X[2][0] = proxy->CreateShare(2.0);
    // X[2][1] = proxy->CreateShare(1.0);
    //X tamam ve lalettayin değerler atandı.

    //şimdi y, yani ki hedef değerler, hazırlanıyor.
    // uint64_t *y = new uint64_t[rows]; 
    // y[0] = proxy->CreateShare(.0);
    // y[1] = proxy->CreateShare(1.0);
    // y[2] = proxy->CreateShare(0.0);
    //y'ler tamam

    // uint64_t *y_test = new uint64_t[rows];  //test için oluşturdum.
    // y_test[0] = proxy->CreateShare(1.0);
    // y_test[1] = proxy->CreateShare(0.0);
    // y_test[2] = proxy->CreateShare(0.0);

    //rows = 30162;
    //columns = 14; //çünkü 14 tane kolon var...
    
    //thetas[2] = proxy->CreateShare(1.0);

    //bu yukarıdakiler ve epsilon değeri ve learning rate elde olduktan sonra çağırmala başlayacak.
    
    vector<double> epsilons = {1}; //epsilonların kullanmaya alışkın olduğumuz bir konfigürasyonu
    for(auto epsilon: epsilons){
        vector<double> accuracies = {}; //python kodunu takip ederek ilerliyorum.
        for(int i=0; i<10; ++i){ //makalede 100 olarak verilen "run time"
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
            //thetas[0] = proxy->CreateShare(.0);
            //thetas[1] = proxy->CreateShare(1.0);
            cout << "gradient'ı çağırıyoruz" << endl;
            proxy->SendBytes(lgGradient, param, 2);


            // uint64_t *theta_double = Reconstruct(proxy, thetas, columns);
            // for(int j=0; j<columns; ++j){
            //     cout << ConvertToDouble(theta_double[j]) << endl;
            // }

            uint64_t *new_thetas = GradientDescent(proxy, X, y, thetas, rows, columns, learning_rate, epsilon, iterations); //X_train ve y_train ile
            //thetas =  GradientDescent(proxy, X, y, thetas, rows, columns, 0.5, epsilon, iterations); //X_train ve y_train ile

            cout << "************************************" << endl;

            uint64_t *new_theta_double = Reconstruct(proxy, new_thetas, columns);
            for(int j=0; j<columns; ++j){
                cout << ConvertToDouble(new_theta_double[j]) << endl;
            }


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

            cout << ConvertToDouble(Reconstruct(proxy, prediction[0])) << endl;
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
            cout << "the accuracy is: " << acc << endl;
        }

    }
}


int main(int argc, char* argv[]) {
    //İPEK BURDAN İTİBAREN SANA LAZIM.
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
    infile.open("/Users/dulumrae/Desktop/preprocessed_data.txt");   //burayı senin değiştirmen lazım!!!!!!!!!

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

    //işleri büyüttük
    cout << "başlıyoruz" << endl;


    


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
        X_read_shares[i] = proxy->CreateShare(X_read[i], column_number); //BUNLAR TEHLİKELİ
    }
    uint64_t *y_read_shares = proxy->CreateShare(y_read, row_number); //TEHLİKEEEE
    double ones[row_number];
    for(int i=0; i<row_number; i++){
        ones[i] = 1;
    }
    //uint64_t *y_test_shares = proxy->CreateShare(ones, row_number); //tamamen uydurdum :p ASLA BÖYLE OLMAMALI ASLINDA!!

    


    //row_number = 30000;
    // uint32_t param_rows[1] = {(uint32_t)row_number};
    // proxy->SendBytes(lgSigmoid, param_rows, 1);
    // cout << "lan nerede" << endl;
    // uint64_t *sigsss = Sigmoid(proxy, y_test_shares, row_number);
    // cout <<"çıktık buradan" << endl;
    // for(int i=0; i<row_number; ++i){
    //     cout << ConvertToDouble(Reconstruct(proxy, sigsss[i])) << endl;
    // }
    
    
    cout << "regresyona geldik" << endl;

    int ulan = 10;
    uint32_t parapara[1] = {uint32_t(ulan)};

    double *dividenda = new double[ulan];
    double *dividendori = new double[ulan];

    for(int i=0; i<ulan; ++i){
        dividenda[i] = i; //0'dan 10'a sayılar
        dividendori[i] = 2;
    }
    uint64_t *dividenda_shares = proxy->CreateShare(dividenda, ulan);
    uint64_t *dividendori_shares = proxy->CreateShare(dividendori, ulan);

    uint32_t div[1] = {uint32_t(ulan)};
    //proxy->SendBytes(coreVectorisedMultiply, div, 1);
    //uint64_t *division = Multiply(proxy, dividenda_shares, dividendori_shares, ulan);

    //proxy->SendBytes(lgSigmoid, div, 1);
    //uint64_t *divienda_sigmoid_shares = Sigmoid(proxy, dividenda_shares, ulan);

    // for(int i=0; i<ulan; ++i){
    //     cout << ConvertToDouble(Reconstruct(proxy, divienda_sigmoid_shares[i])) << endl;
    // }



    /*
    double scale = 1.0;
    uint32_t noisior[2] = {uint32_t(scale), uint32_t(ulan)};
    proxy->SendBytes(lgAddNoiseNew, noisior, 2);
    auto stuff = add_noise_new(proxy, dividenda_shares, 1.0, ulan); //scale = 1 dedim çünkü meh
    for(int i=0; i<ulan; ++i){
        cout << ConvertToDouble(Reconstruct(proxy, stuff[i])) << endl;
    }

    cout << "bunu göreyim" << endl;
    */



    //Regression(proxy, X_read_shares, y_read_shares, y_test_shares, row_number, column_number, 10); //
    //Regression(proxy, X_read_shares, y_read_shares, y_read_shares, row_number, column_number, 3); //

    int rr = 3;
    int nn = 2;
    Regression(proxy, X, y, y_test, rr, nn, 1); //


    //İPEK ZATEN BURADA RETURN EDİYORUZ, AŞAĞIDAKİ ZIKKIMATI OKUMA BİLE.
    return 0;





















    //BAMBAŞKA İŞLER PEŞİNDEYİM:
    int const three = 5; //five is the new three
    uint32_t para[1] = {three};
    double pos[three] = {-3,0,0.5,0.75,4};
    
    double pos2[three] = {3,1,43};

    uint64_t *pos_shares = proxy->CreateShare(pos, three);
    uint64_t *pos2_shares = proxy->CreateShare(pos2, three);

    proxy->SendBytes(lgSigmoid, para, 1);
    uint64_t *is_poses = Sigmoid(proxy, pos_shares, three);

    for(int i=0; i<three; i++){
        cout<<" sig of " << pos[i] <<" =  " << ConvertToDouble(Reconstruct(proxy, is_poses[i])) << endl;
    }

    return 1;
    /*SIGMOID FONKSİYONUMU DENİYORUM:
    */
    //double z = 0.5;
    //uint64_t z_share =  proxy->CreateShare(z); //double olması lazımmış.

    



    // double xs[3] = {1,2,3};
    // double ys[3] = {10,10,10};
    // uint64_t *x_shares = proxy->CreateShare(xs, 3);
    // uint64_t *y_shares = proxy->CreateShare(ys, 3);
    // proxy->SendBytes(coreVectorisedCompare, param, 1); //buna para ve 1 eklemek lazım sanırım ama anlamıyorum şu an.
    // uint64_t *results = Compare(proxy, x_shares, y_shares, size);
    // for(int i=0; i<size; i++){
    //     cout << ConvertToDouble(Reconstruct(proxy, results[i], size)) << endl;
    // }


    uint32_t const this_size = 5;
    //double say[this_size] = {2,6,69,14,99};
    //uint64_t *say_share = proxy->CreateShare(say, this_size); //şu an lazım değil aslında bu.

    double unimportant = 10;
    uint64_t unimp_share = proxy->CreateShare(unimportant);

    uint32_t parametre[1] = {this_size};
    proxy->SendBytes(lgfonk, parametre, 1);
    uint64_t *results = fonk(proxy, unimp_share, this_size);


    for(int i=0; i<this_size; i++){
        cout << ConvertToDouble(Reconstruct(proxy, results[i])) << endl;
    }


    
   
    

    


    

    for(int i=0; i<this_size; i++){
        cout << ConvertToDouble(Reconstruct(proxy, results[i])) << endl;
    }






    
    double vector1[3] = {-2,0.5,4};

    uint64_t* v1_share = proxy->CreateShare(vector1, 3);


    uint64_t* sigmoid_term1_shares = Sigmoid(proxy, v1_share, 3);

    for(int i=0; i<3; i++){
        cout << ConvertToDouble(Reconstruct(proxy,sigmoid_term1_shares[i])) << endl;
    }

    cout <<"çalışıyor lan sanki: " << endl;
    double z = 0.25;
    uint64_t z_share = proxy->CreateShare(z);
    //uint64_t sig_result = Sigmoid(proxy, z_share);
    //cout << "sigmoid(0.5) = " << ConvertToDouble(Reconstruct(proxy, sig_result)) << endl;

    

    int const boyut = 3;
    uint32_t param[1] = {boyut};
    double somes[boyut] = {-10,1,20};
    uint64_t *some_shares = proxy->CreateShare(somes, boyut);
    proxy->SendBytes(lgSigmoid, param, 1); //param muhabbeti gerekmiyro ama olsun
    uint64_t *sigs = Sigmoid(proxy, some_shares, boyut);
    for(int i=0; i<boyut; i++){
        cout << somes[i] << " > 0: " << ConvertToDouble(Reconstruct(proxy,sigs[i])) << endl;
        cout << somes[i] << " > 0: " << (Reconstruct(proxy,sigs[i])) << endl;
    }



    //CompareAll()

    /*
    BUNLARI SİGMOİD FONKSİYONUNU VEKTÖRİZE ETME SÜRECİNDE YORUMA ALIYORUM.
    vector<double> zs = {-2,0,0.5,1,2};
    vector<uint64_t> zs_shares = {};
    // for(auto z: zs){
    //     uint64_t share = proxy->CreateShare(z);
    //     zs_shares.push_back(share);
    // }
    zs_shares = proxy->CreateShare(zs, 5);
    vector<uint64_t> sigs = Sigmoid(proxy,zs_shares);
    for(int i=0; i<sigs.size(); i++){
        cout << "for z = " << zs[i] << " sigmoid(z) is " << ConvertToDouble(Reconstruct(proxy, sigs[i])) << endl;
    }*/



    
    //uint64_t sigmoid_term3 = 
    /*
    cout <<"boolean işini anlayana kadar paydos" <<endl;
    uint8_t dene1 = 1;
    uint8_t dene0 = 0;
    uint8_t d1_sec[] = {(dene1)}; //bool'lar 8'e sığıyor sanki.  []
    uint8_t d0_sec[] = {(dene0)};

    proxy->SendBytes(boolAnd);
    uint8_t* d0_and_d1_share = And(proxy,d0_sec,d1_sec,1); //anam ağladı tipleri doğru yazana kadar
    cout << "gizli hâli şöyle: " << d0_and_d1_share[0] << endl;

    bool d0_and_d1 = ConvertToDouble(Reconstruct(proxy,d0_and_d1_share[0]));
    cout << "1 && 0 işleminin sonucu: " << d0_and_d1 << " geldi valla ve 0 olmalıydı" << endl;
    */
    
    return 0;

    /* TODO: Create the random number x:
     */
    int x = 10;
    /*
     * Secret sharing is a fundamental component of multi-party computation where multiple parties work collaboratively
     * to compute a desired function using their secret shares without revealing their secrets to each other and the
     * result of the computation. In literature, there are different secret sharing techniques. In CECILIA2, we use
     * 2-out-of-2 additive secret sharing. In this secret sharing scheme, we generate two values such that their
     * summation over a ring N gives the secret value. What we meant by this is this: Assume that the ring is 16.
     * We sum two values and take their modulus 16. We use unsigned 64-bit integers to hold the shares. To create the
     * shares, we can call CreateShare function of Party instances.
     */
    //TODO: Create secret shares of our sample data:
    uint64_t x_share =  proxy->CreateShare((double)x); //double olması lazımmış.


     //TODO: Now calculate the function 5x+18


    uint64_t multiplied = 5 * x_share; //skalerle çarpma işi lokal yapılabilir.
    uint64_t result = multiplied;
    if(proxy->GetPRole() == proxy1){
        result = multiplied + ConvertToUint64(18); //kendi elemanıysa ekleme yapacak, o da lokal
        }

   // cout << "x'in paylaşılmış hâli: " << x_share << endl;
   // cout << "reconsturcted x: " << ConvertToDouble(Reconstruct(proxy, x_share)) << endl;
   // cout << "Reconstructed result (5x+18): " << ConvertToDouble(Reconstruct(proxy, result)) << endl;
   // cout << "yukarıdaki 68 olmalı..." << endl;


    //TODO: Let's reconstruct the resulting share
    double result_reconstructed = ConvertToDouble(Reconstruct(proxy, result));
    cout << "Reconstructed result (5x+18): " << ConvertToDouble(Reconstruct(proxy, result)) << endl;

    /**
     * TODO: Now let us now lets calculate the square of our variable x
     * Note that this is not a local operation we also need to
     * inform the helper, the third computing party of CECILIA2, about the operation that we would like to perform.
     * */

    proxy->SendBytes(coreMultiply); //çarpım yapacağımızın haberini veriyoruz.
    uint64_t square = Multiply(proxy, x_share, x_share); //simple yet effective (?)
    double square_reconstructed = ConvertToDouble(Reconstruct(proxy, square));
    cout << "the square of " << x << " is " << x*x << endl;
    cout << "this code calculates it as " << square_reconstructed << endl;


    /*
     * Great, we can do operations on a single variable.
     * But usually we have more than one data point.
     * TODO: Create 2 random arrays of length 4 each
     * TODO: And create their shares
    */
    double arr1[4] = {4,2,7,12};
    double arr2[4] = {1,5,2,3};

    uint64_t* arr1_share = proxy->CreateShare(arr1,4); //iki nokta: 1)CreateShare size'ı da alıyor, 2)uint64_t ARRAYİ oluşturuyoruz; pointer bu
    uint64_t* arr2_share = proxy->CreateShare(arr2,4);

    /*
     * TODO: Calculate the element-wise multiplication of these two arrays using Multiply function
     * Try the vectorized version, do not multiply them one by one
     * */
    uint32_t params[1] = {4}; //bu kısım sıkıntılı biraz...
    proxy->SendBytes(coreVectorisedMultiply,params,1);
    uint64_t* result_of_multiplication = Multiply(proxy, arr1_share, arr2_share, 4); //4 tane eleman olduğunu falan da söylemek lazım.

    double* mult_result_reconstructed = ConvertToDouble(Reconstruct(proxy, result_of_multiplication,4) ,4);

    for(int i=0; i<4; i++){
        cout << arr1[i] << " * " << arr2[i] << " = " << arr1[i]*arr2[i] << " =?= " << mult_result_reconstructed[i] << endl;
    }
    cout << endl << endl;

    /*
     * If we can multiply these two arrays element-wise, we can also calculate the dot product of them!
     * TODO: Calculate the dot product of these two arrays.
     * Reconstruct and check if it is correct
     */

    double vec1[3] = {4,6,3};
    double vec2[3] = {1,2,8};

    uint64_t* vec1_share = proxy->CreateShare(vec1, 3);
    uint64_t* vec2_share = proxy->CreateShare(vec2, 3);

    //tell the helper that we are multiplying vectors:
    uint32_t params_dot_product[1] = {3}; //bunu hâlâ pek anlamıyorum...
    proxy->SendBytes(coreVectorisedMultiply,params_dot_product,1);

    uint64_t* product_vector_shared = Multiply(proxy, vec1_share, vec2_share, 3);

    uint64_t sum = 0;

    for(int i=0; i<3; i++){
        sum = Add(proxy, product_vector_shared[i], sum); //normal toplama da olabilirdi.
    }
    double dot_product_result = ConvertToDouble(Reconstruct(proxy, sum));

    double* product_vector = ConvertToDouble(Reconstruct(proxy,product_vector_shared,3),3);

    double true_result = 0;
    for(int i=0; i<3; i++){
        true_result += product_vector[i]; 
    }

    cout << "the dot product of these two vectors should be: "<< true_result << endl;
    cout << "shared shit found: " << dot_product_result << endl;

    //TODO: Create another array of length $5$. This time add some negative integers into the list. Use the MostSignificantBit tp detect if those numbers are negative or not.

    /*
     * Once the proxies are done with the computation, they need to let the Helper know that the process is over, and
     * it can also terminate. For this, the proxies need to send coreEnd signal to the Helper.
     * Remove the comment from the following line
     */
    proxy->SendBytes(coreEnd);

    return 0;
}