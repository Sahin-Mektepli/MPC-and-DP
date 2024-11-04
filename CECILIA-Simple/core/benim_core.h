//
// Created by Mete Akgun on 28.12.21.
//

#ifndef CORE_H
#define CORE_H

#include "Party.h"
#include "../utils/test_functions.h"
#include <thread>
#include <mutex>
#include <bitset>

#include <vector> //bana lazım oldu

double recn_time = 0;
double mul_triple_gen = 0;
double mul_ef_calc = 0;
double mul_z_calc = 0;
double rec_calc = 0;
double rec_transfer = 0;
double rec_read = 0;

/**
 * Perform the truncation operation which we use to keep the number of fractional bit consistent after Multiply operation
 * @param proxy
 * @param z: value we want to Truncate
 * @return truncated z is returned
 */
uint64_t Truncate(Party *const proxy, uint64_t z, int shift = FRACTIONAL_BITS) {
    switch (proxy->GetPRole()) {
        case proxy1:
            z = ArithmeticShift(z, shift);
            break;
        case proxy2:
            z = -1 * ArithmeticShift(-1 * z, shift);
            break;
        case helper:
            break;
    }
    return z;
}

uint64_t *Reconstruct(Party *const proxy, const uint64_t *const a, uint32_t sz, uint64_t mask= RING_SIZE) {
    if (proxy->GetPRole() != helper) {
        uint64_t *b = new uint64_t[sz];
        if (proxy->GetPRole() == proxy1 ) {
            unsigned char *ptr = proxy->GetBuffer1();
            for (int i = 0; i < sz; i++) {
                AddValueToCharArray(a[i], &ptr);
            }
            thread thr1 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer1(), sz * 8);
            thread thr2 = thread(Receive, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 8);
            thr1.join();
            thr2.join();

            ptr = proxy->GetBuffer2();
            for (int i = 0; i < sz; i++) {
                b[i] = ConvertToLong(&ptr);
            }

        } else if (proxy->GetPRole() == proxy2) {
            unsigned char *ptr = proxy->GetBuffer1();
            for (int i = 0; i < sz; i++) {
                AddValueToCharArray(a[i], &ptr);
            }
            thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 8);
            thread thr2 = thread(Receive, proxy->GetSocketP1(), proxy->GetBuffer2(), sz * 8);
            thr1.join();
            thr2.join();
            ptr = proxy->GetBuffer2();
            for (int i = 0; i < sz; i++) {
                b[i] = ConvertToLong(&ptr);
            }
        }
        for (int i = 0; i < sz; i++) {
            b[i] = (a[i] + b[i]) & mask;
        }
        return b;
    } else {
        return nullptr;
    }
}

uint64_t Reconstruct(Party *const proxy, uint64_t a, uint64_t mask= RING_SIZE) {
    uint64_t *result_array = Reconstruct(proxy, &a, 1, mask);
    uint64_t result = result_array[0];
    delete[] result_array;
    return result;
}

uint64_t *ReconstructNarrow(Party *const proxy, const uint64_t *const a, uint32_t sz, uint64_t ringbits) {
    auto mask = (1<< ringbits)-1;
    auto bsz = (uint32_t)ceil(ringbits/8.0);
    uint64_t *b = new uint64_t[sz];
    if ( proxy->GetPRole() == proxy1 ) {
        unsigned char *ptr = proxy->GetBuffer1();
        //auto start = chrono::high_resolution_clock::now();
        WriteToBuffer(a, ptr, sz, bsz);
        thread thr1 = thread(Send,proxy->GetSocketP2(), proxy->GetBuffer1(), sz*bsz);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), sz*bsz);
        thr1.join();
        thr2.join();
       // auto end = chrono::high_resolution_clock::now();
       // rec_transfer +=
       //         chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
       // start = chrono::high_resolution_clock::now();
        ptr = proxy->GetBuffer2();
        ReadBuffer(b, ptr, sz, bsz);
       // end = chrono::high_resolution_clock::now();
       // rec_read +=
       //         chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
       // start = chrono::high_resolution_clock::now();
#pragma omp parallel for num_threads(4)
        for (int i = 0; i < sz; i++) {
            b[i] = (a[i] + b[i]) & mask;
        }
        //end = chrono::high_resolution_clock::now();
        //rec_calc +=
        //        chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
    } else if ( proxy->GetPRole() == proxy2) {
        unsigned char *ptr = proxy->GetBuffer1();
        WriteToBuffer(a, ptr, sz, bsz);
        //auto start = chrono::high_resolution_clock::now();
        thread thr1 = thread(Send,proxy->GetSocketP1(), proxy->GetBuffer1(), sz*bsz);
        thread thr2 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer2(), sz*bsz);
        thr2.join();
        thr1.join();
        //auto end = chrono::high_resolution_clock::now();
        //rec_transfer +=
        //               chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
        //start = chrono::high_resolution_clock::now();
        ptr = proxy->GetBuffer2();
        ReadBuffer(b, ptr, sz, bsz);
        //end = chrono::high_resolution_clock::now();
        //rec_read +=
        //          chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
        //start = chrono::high_resolution_clock::now();
#pragma omp parallel for num_threads(4)
        for (int i = 0; i < sz; i++) {
            b[i] = (a[i] + b[i]) & mask;
        }
        //end = chrono::high_resolution_clock::now();
	//rec_calc +=
        //          chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
    }
    return b;
}

/**Reconstruct a secret shared 2D array.*/
uint64_t** Reconstruct(Party *const proxy, const uint64_t *const *const a, uint32_t n_row, uint32_t n_col) {
    if (proxy->GetPRole() != helper) {
        uint64_t **b = new uint64_t*[n_row];
        if (proxy->GetPRole() == proxy1) {
            unsigned char *ptr = proxy->GetBuffer1();
            for (int i = 0; i < n_row; i++) {
                b[i] = new uint64_t[n_col];
                for( int j = 0; j < n_col; j++) {
                    AddValueToCharArray(a[i][j], &ptr);
                }
            }
            thread thr1 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer1(), n_row * n_col * 8);
            thread thr2 = thread(Receive, proxy->GetSocketP2(), proxy->GetBuffer2(), n_row * n_col * 8);
            thr1.join();
            thr2.join();
            ptr = proxy->GetBuffer2();
            for (int i = 0; i < n_row; i++) {
                for(int j = 0; j < n_col; j++) {
                    b[i][j] = ConvertToLong(&ptr);
                }
            }
        } else if (proxy->GetPRole() == proxy2) {
            unsigned char *ptr = proxy->GetBuffer1();
            for (int i = 0; i < n_row; i++) {
                for( int j = 0; j < n_col; j++) {
                    AddValueToCharArray(a[i][j], &ptr);
                }
            }
            thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), n_row * n_col * 8);
            thread thr2 = thread(Receive, proxy->GetSocketP1(), proxy->GetBuffer2(), n_row * n_col * 8);
            thr1.join();
            thr2.join();
            ptr = proxy->GetBuffer2();
            for (int i = 0; i < n_row; i++) {
                b[i] = new uint64_t[n_col];
                for( int j = 0; j < n_col; j++) {
                    b[i][j] = ConvertToLong(&ptr);
                }
            }
        }
        for (int i = 0; i < n_row; i++) {
            for( int j = 0; j < n_col; j++) {
                b[i][j] = (a[i][j] + b[i][j]);
            }
        }
        return b;
    } else {
        return nullptr;
    }
}

/**Reconstruct a secret shared 3D array.*/
uint64_t*** Reconstruct(Party *const proxy, const uint64_t *const *const *const a, uint32_t n_matrices, uint32_t n_row, uint32_t n_col) {
    if(proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t ***b = new uint64_t**[n_matrices];
        unsigned char *ptr = proxy->GetBuffer1();
        // add the values of a to the buffer and initialize b
        for(int n = 0; n < n_matrices; n++) {
            b[n] = new uint64_t * [n_row];
            for (int i = 0; i < n_row; i++) {
                b[n][i] = new uint64_t[n_col];
                for( int j = 0; j < n_col; j++) {
                    AddValueToCharArray(a[n][i][j], &ptr);
                }
            }
        }
        // determine the
        int *partner_socket;
        if (proxy->GetPRole() == proxy1) {
            partner_socket = proxy->GetSocketP2();
        }
        else {
            partner_socket = proxy->GetSocketP1();
        }
        thread thr1 = thread(Send, partner_socket, proxy->GetBuffer1(), n_matrices * n_row * n_col * 8);
        thread thr2 = thread(Receive, partner_socket, proxy->GetBuffer2(), n_matrices * n_row * n_col * 8);
        thr1.join();
        thr2.join();
        ptr = proxy->GetBuffer2();

        for(int n = 0; n < n_matrices; n++) {
            for (int i = 0; i < n_row; i++) {
                for (int j = 0; j < n_col; j++) {
                    b[n][i][j] = ConvertToLong(&ptr);
                }
            }
        }

        for(int n = 0; n < n_matrices; n++) {
            for (int i = 0; i < n_row; i++) {
                for (int j = 0; j < n_col; j++) {
                    b[n][i][j] += a[n][i][j];
                }
            }
        }

        return b;
    }
    return nullptr;
}


uint64_t Add(Party *const proxy, uint64_t a, uint64_t b) {
    return a + b;
}

/**
 * Adds values of a and b at equal position.
 * @param proxy
 * @param a
 * @param b
 * @param size length of vectors a and b
 * @return vector of length size containing the sum of according values in a and b.
 */
uint64_t* Add(Party *const proxy, const uint64_t *const a, const uint64_t *const b, uint32_t size) {
    uint64_t* sum = new uint64_t[size];
    for(int i = 0; i<size; i++){
        sum[i] = a[i] + b[i];
    }
    return sum;
}

/** For smaller ring size
 * @param mt1 3-by-size array whose rows will be a_i, b_i and c_i, respectively
 * @param mt2 3-by-size array whose rows will be a_i, b_i and c_i, respectively
 * @param size the number of multiplication triples that will be generated
 */
void GenerateMultiplicationTriple(Party *const proxy, uint64_t *const c1, uint32_t size, uint64_t mask) {

    for (int i = 0; i < size; i++) {
        uint64_t a0 = proxy->GenerateCommonRandom()&mask;
        uint64_t a1 = proxy->GenerateCommonRandom2() & mask;
        uint64_t b0 = proxy->GenerateCommonRandom()&mask;
        uint64_t b1 = proxy->GenerateCommonRandom2() & mask;
        uint64_t c0=  proxy->GenerateCommonRandom()&mask;
        c1[i] = (((a0+a1)*(b0+b1)) - c0)&mask;
        }
}

/**
 * Adds values of all vectors in a at equal position in a row to calculate their sum (sum over column where each row is one vector).
 * @param proxy
 * @param a matrix containing several vectors of length size. Values of all vectors at same position shall be summed up.
 * @param n_vectors number of vectors in a
 * @param size length of each vector in a
 * @return vector of length size
 */
uint64_t* Add(Party *const proxy, const uint64_t *const *const a, int n_vectors, int size) {
    uint64_t* res = new uint64_t [size];
    for(int i = 0; i<size; i++){
        res[i] = 0;
        for(int v = 0; v<n_vectors; v++){
            res[i] += a[v][i];
        }
    }
    return res;
}


/** Strictly a function for the helper to generate Beaver's multiplication triples.
 * It is used to compute proxy2's shares of c, which then need to be sent to it.
 * @param c1 array in which proxy2's shares of C will be stored
 * @param size the number of multiplication triples that will be generated
 */
 void GenerateMultiplicationTriple(Party *const proxy, uint64_t *const c1, uint32_t size) {
     for (int i = 0; i < size; i++) {
        uint64_t a0 = proxy->GenerateCommonRandom();
        uint64_t a1 = proxy->GenerateCommonRandom2();
        uint64_t b0 = proxy->GenerateCommonRandom();
        uint64_t b1 = proxy->GenerateCommonRandom2();
        uint64_t c0=  proxy->GenerateCommonRandom();
        c1[i] = ((a0+a1)*(b0+b1)) - c0;
    }
}


/**
 * Adds values of all matrices in a at equal position to calculate their sum (sum over all matrices in a).
 * @param proxy
 * @param a 3-dmatrix containing several 2-d matrices in dimension rows x cols. Values of all matrices at same position shall be summed up.
 * @param n_matrices number of matrices in a
 * @param rows height of each matrix in a
 * @param cols width of each matrix in a
 * @return 2-d matrix of shape rows x cols with the summed up values.
 */
uint64_t** Add(Party *const proxy, const uint64_t *const *const *const a, int n_matrices, int rows, int cols) {
    uint64_t** res = new uint64_t *[rows];
    for(int r = 0; r<rows; r++){
        res[r] = new uint64_t [cols];
        for(int c = 0; c<cols; c++){
            res[r][c] = 0;
            for (int m = 0; m < n_matrices; ++m) {
                res[r][c] += a[m][r][c];
            }
        }
    }
    return res;
}

uint64_t* Multiplex(Party *const proxy, const uint64_t *const x, const uint64_t *const y, const uint64_t *const b, uint32_t sz, int shift = FRACTIONAL_BITS) {
    if (proxy->GetPRole() == proxy1){
        unsigned char *ptr = proxy->GetBuffer1();
        uint64_t *res = new uint64_t[sz];
        uint64_t *m1 = new uint64_t[sz];
        for (uint32_t i = 0; i < sz; i++) {
            uint64_t r1= proxy->GenerateCommonRandom(), r2= proxy->GenerateCommonRandom(), r3= proxy->GenerateCommonRandom(), r4= proxy->GenerateCommonRandom();

            m1[i] = (b[i] * (x[i] - y[i])) - (r2*b[i]) - (r3*(x[i] - y[i])) - (r3*r4);
            uint64_t m2 = b[i] + r1;
            uint64_t m3 = x[i] - y[i] + r4;

            AddValueToCharArray(m2, &ptr);
            AddValueToCharArray(m3, &ptr);
        }
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 16);
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 8);
        ptr = proxy->GetBuffer1();
        for (uint32_t i = 0; i < sz; i++) {
            res[i] = m1[i] + ConvertToLong(&ptr);
            res[i] = res[i] >> shift;
            res[i] = x[i] - res[i];
        }
        delete [] m1;
        return res;

    }else if (proxy->GetPRole() == proxy2){
        unsigned char *ptr = proxy->GetBuffer1();
        uint64_t *res = new uint64_t[sz];
        uint64_t *m1 = new uint64_t[sz];
        for (uint32_t i = 0; i < sz; i++) {
            uint64_t r1= proxy->GenerateCommonRandom(), r2= proxy->GenerateCommonRandom(), r3= proxy->GenerateCommonRandom(), r4= proxy->GenerateCommonRandom();

            m1[i] = (b[i] * (x[i] - y[i])) - (r1*(x[i] - y[i])) - (r1*r2) - (r4*b[i]);
            uint64_t m2 = x[i] - y[i] + r2;
            uint64_t m3 = b[i] + r3;

            AddValueToCharArray(m2, &ptr);
            AddValueToCharArray(m3, &ptr);

        }
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 16);
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 8);
        ptr = proxy->GetBuffer1();
        for (uint32_t i = 0; i < sz; i++) {
            res[i] = m1[i] + ConvertToLong(&ptr);
            res[i] = -1 * ((-1 * res[i]) >> shift);
            res[i] = x[i] - res[i];
        }
        delete [] m1;
        return res;

    }else if (proxy->GetPRole() == helper){
        thread thr1 = thread(Receive, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 16);
        thread thr2 = thread(Receive, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 16);
        thr1.join();
        thr2.join();
        unsigned char *ptr = proxy->GetBuffer1();
        unsigned char *ptr_out = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();
        unsigned char *ptr_out2 = proxy->GetBuffer2();
        for (uint32_t i = 0; i < sz; i++) {
            uint64_t m2 = ConvertToLong(&ptr);
            uint64_t m3 = ConvertToLong(&ptr);

            uint64_t m5 = ConvertToLong(&ptr2);
            uint64_t m6 = ConvertToLong(&ptr2);

            uint64_t m = (m2 * m5) + (m3 * m6);
            m2 = proxy->GenerateRandom();
            m3 = m-m2;
            AddValueToCharArray(m2, &ptr_out);
            AddValueToCharArray(m3, &ptr_out2);
        }
        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 8);
        thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 8);
        thr1.join();
        thr2.join();
        return NULL;
    }
    return NULL;
}

uint64_t Multiplex(Party *const proxy, uint64_t x, uint64_t y, uint64_t b, int shift = FRACTIONAL_BITS) {
    if (proxy->GetPRole() == helper) {
        Multiplex(proxy, nullptr, nullptr, nullptr, 1, shift);
        return 0;
    } else {
        uint64_t *result_vector = Multiplex(proxy, &x, &y, &b, 1, shift);
        uint64_t result = result_vector[0];
        delete[] result_vector;
        return result;
    }
}

/**
 * Private Compare Bool: Check b>a
 *
 * @param a reconstructed value
 * @param b boolean share
 * @param L1
 * @return
 */
uint8_t *PrivateCompareBool(Party *const proxy, const uint64_t *const a, const uint8_t *const b, uint32_t sz, int L1) {
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        for (int j = 0; j < sz; j++) {
            int jk = j * L1;
            uint8_t w_sum = 0;
            for (int i = L1 - 1; i >= 0; i--) {
                uint8_t a_bit = Bit(a[j], i);
                int k = jk + i;
                uint8_t w = Mod((b[k] + proxy->GetPRole() * a_bit - 2 * a_bit * b[k]) % LP, LP);
                proxy->GetBuffer1()[k] =
                        (Mod((proxy->GetPRole() * a_bit - b[k] + proxy->GetPRole() + w_sum), LP) * (proxy->GenerateCommonRandom() % (LP - 1) + 1)) %
                        LP;
                w_sum = (w_sum + w) % LP;
            }
            for (int i = 0; i < L1; i++) {
                int ind1 = (proxy->GenerateCommonRandom() % L1) + jk;
                int ind2 = (proxy->GenerateCommonRandom() % L1) + jk;
                uint8_t tmp = proxy->GetBuffer1()[ind1];
                proxy->GetBuffer1()[ind1] = proxy->GetBuffer1()[ind2];
                proxy->GetBuffer1()[ind2] = tmp;
            }
        }
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * L1);
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz);
        uint8_t *r = new uint8_t[sz];
        for (int i = 0; i < sz; i++)
            r[i] = proxy->GetBuffer1()[i];
        return r;

    } else if (proxy->GetPRole() == helper) {
        Receive(proxy->GetSocketP1(), proxy->GetBuffer1(), sz * L1);
        Receive(proxy->GetSocketP2(), proxy->GetBuffer2(), sz * L1);
        unsigned char *ptr_out = proxy->GetBuffer1();
        unsigned char *ptr_out2 = proxy->GetBuffer2();
        for (int j = 0; j < sz; j++) {
            int jk = j * L1;
            uint8_t res = 0;
            for (int i = 0; i < L1; i++) {
                proxy->GetBuffer1()[jk + i] = (proxy->GetBuffer1()[jk + i] + proxy->GetBuffer2()[jk + i]) % LP;
                if (((int) proxy->GetBuffer1()[jk + i]) == 0) {
                    res = 1;
                    break;
                }
            }
            uint8_t res1 = proxy->GenerateRandom() % 2;
            uint8_t res2 = res ^res1;

            AddValueToCharArray(res1, &ptr_out);
            AddValueToCharArray(res2, &ptr_out2);
        }
        Send(proxy->GetSocketP1(), proxy->GetBuffer1(), sz);
        Send(proxy->GetSocketP2(), proxy->GetBuffer2(), sz);
        return nullptr;
    }
    return nullptr;
}

/** Check @p b>@p a
 *
 * @param a reconstructed value
 * @param b boolean share
 * @param L1
 * @return
 */
uint8_t PrivateCompareBool(Party *const proxy, uint64_t a, const uint8_t *const b, int L1) {
    if (proxy->GetPRole() == helper) {
        PrivateCompareBool(proxy, nullptr, nullptr, 1, L1);
        return 0;
    } else {
        auto result_vector = PrivateCompareBool(proxy, &a, b, 1, L1);
        uint64_t result = result_vector[0];
        delete[] result_vector;
        return result;
    }
}

/** Multiple modular conversions.
 *
 * @param x an array of values in the ring 2^63
 * @param sz the length of @p x
 * @return an array of values in the ring 2^64
 */
uint64_t *ModularConversion(Party *const proxy, const uint64_t *const x, uint32_t sz) {
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *z_1 = new uint64_t[sz];
        uint64_t *ya = new uint64_t[sz];
        uint8_t *yb = new uint8_t[sz * (L_BIT - 1)];
        uint8_t *w = new uint8_t[sz];
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * (8 + L_BIT));
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < sz; i++) {
            ya[i] = ConvertToLong(&ptr);
            ConvertToArray(&ptr, &yb[i * (L_BIT - 1)], L_BIT - 1);
            w[i] = (*ptr);
            ptr++;
            z_1[i] = (x[i] + ya[i]) & N1_MASK;
        }
        uint64_t *z = Reconstruct(proxy, z_1, sz, N1_MASK);
        uint8_t *wc = PrivateCompareBool(proxy, z, yb, sz, L_BIT - 1);

        for (int i = 0; i < sz; i++) {
            w[i] = w[i] ^ wc[i];
            if (proxy->GetPRole() == proxy1 && z_1[i] > z[i])
                z_1[i] = z_1[i] + N1;
            z_1[i] = (z_1[i] - (ya[i] + w[i] * N1));
        }
        delete[] ya;
        delete[] yb;
        delete[] w;
        return z_1;
    }
    else if (proxy->GetPRole() == helper) {
        unsigned char *ptr_out = proxy->GetBuffer1();
        unsigned char *ptr_out2 = proxy->GetBuffer2();
        for (int i = 0; i < sz; i++) {
            uint64_t y = proxy->GenerateRandom() & N1_MASK;
            uint64_t ya_1 = proxy->GenerateRandom() & N1_MASK;
            uint64_t ya_2 = (y - ya_1) & N1_MASK;
            AddValueToCharArray(ya_1, &ptr_out);
            AddValueToCharArray(ya_2, &ptr_out2);
            for (int j = 0; j < L_BIT - 1; j++) {
                uint8_t k = (y >> j) & 0x1;
                uint8_t yb_1 = proxy->GenerateRandom() % LP;
                uint8_t yb_2 = Mod(k - yb_1, LP);
                AddValueToCharArray(yb_1, &ptr_out);
                AddValueToCharArray(yb_2, &ptr_out2);
            }
            uint8_t w = 0;
            if (ya_1 > y)
                w = 1;
            uint8_t w_1 = proxy->GenerateRandom() % 2;
            uint8_t w_2 = w ^w_1;
            AddValueToCharArray(w_1, &ptr_out);
            AddValueToCharArray(w_2, &ptr_out2);
        }

        thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * (8 + L_BIT));
        thread thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * (8 + L_BIT));
        thr1.join();
        thr2.join();
        // proxy1 and proxy2 will call MPrivateCompareBool
        PrivateCompareBool(proxy, 0, 0, sz, L_BIT - 1);
        return NULL;
    }
    return NULL;
}

/** Modular conversion.
 *
 * @param x a value in the ring 2^63
 * @return
 */
uint64_t ModularConversion(Party *const proxy, uint64_t x) {
    if (proxy->GetPRole() == helper) {
        ModularConversion(proxy, nullptr, 1);
        return 0;
    } else {
        uint64_t *result_vector = ModularConversion(proxy, &x, 1);
        uint64_t result = result_vector[0];
        delete[] result_vector;
        return result;
    }
}

void MostSignificantBitSubroutine(
    Party *const proxy,
    const uint64_t *const x,
    const uint64_t *const z,
    uint64_t *const z_1,
    const uint8_t *const yb,
    const uint64_t *const ya,
    uint8_t f,
    uint8_t rnd,
    int start_index,
    int end_index
){
    int L1 = L_BIT - 1;
    unsigned char *ptr_out = proxy->GetBuffer1();
    ptr_out += (start_index * (L1+16));
    int buffer_index = (start_index * (L1+16));
    int y_index = (start_index * L1);
    for (int i = start_index; i < end_index; i++) {
        uint8_t w_sum = 0;
        for (int t = L1 - 1; t >= 0; t--) {
            uint8_t a_bit = Bit(z[i], t);
            int bi = buffer_index + t;
            int yi = y_index + t;
            uint8_t w = Mod((yb[yi] + proxy->GetPRole() * a_bit - 2 * a_bit * yb[yi]) % LP, LP);
            proxy->GetBuffer1()[bi] = (Mod((proxy->GetPRole() * a_bit - yb[yi] + proxy->GetPRole() + w_sum), LP) * ((rnd % (LP - 1)) + 1)) % LP;
            rnd += 7;
            w_sum = (w_sum + w) % LP;
        }
        buffer_index += L1;
        y_index += L1;
        ptr_out += L1;


        uint8_t isWrap = 0;
        if (z[i]<z_1[i])
            isWrap = 1;
        z_1[i] =  z_1[i] + proxy->GetPRole() * isWrap * N1;
        AddValueToCharArray(proxy->GetPRole() * f * N1 - x[i] + z_1[i] - ya[i], &ptr_out);
        AddValueToCharArray(proxy->GetPRole() * (1 - f) * N1 - x[i] + z_1[i] - ya[i], &ptr_out);
        buffer_index +=16;
    }
}


/**
 * @brief Computes the most significant bit.
 *
 * MSB has 4 communication round. ModularConversion and PC are hardcoded in MostSignificantBit to reduce the number of
 * communication rounds of MostSignificantBit calls.
 * @param proxy p_proxy:
 * @param x p_x:
 * @param sz p_sz: The number of elements in x.
 * @param format p_format: Whether to convert the result to the regular representation. Defaults to true.
 * @return uint64_t* x < 0
 */
uint64_t *MostSignificantBit(Party *const proxy, const uint64_t *const x, uint32_t sz, size_t shift=FRACTIONAL_BITS, bool format = true) {
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint8_t f = proxy->GenerateCommonRandomByte() & 0x1;
        uint64_t *z_1 = new uint64_t[sz];
        uint64_t *ya = new uint64_t[sz];
        uint8_t *yb = new uint8_t[sz * (L_BIT - 1)];

        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * (8 + L_BIT - 1));

        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < sz; i++) {
            uint64_t dk = x[i] & N1_MASK;
            ya[i] = ConvertToLong(&ptr);
            ConvertToArray(&ptr, &yb[i * (L_BIT - 1)], L_BIT - 1);
            z_1[i] = (dk + ya[i]) & N1_MASK;
        }

        uint64_t *z = Reconstruct(proxy, z_1, sz, N1_MASK);
        int block_size = (int)ceil(sz * 1.0 / SOCKET_NUMBER);
        if (block_size == 0)
            block_size = sz;

        thread thr[SOCKET_NUMBER];
        int start_index = 0;
        int end_index = block_size;
        int thr_num = 0;
        for (int i = 0; i < SOCKET_NUMBER; i++) {
            uint8_t rnd = proxy->GenerateCommonRandomByte();
            thr[i] = thread(MostSignificantBitSubroutine, proxy, x, z, z_1, yb, ya, f, rnd, start_index, end_index);
            thr_num +=1;
            start_index += block_size;
            end_index += block_size;
            if (start_index >= sz)
                break;
            if (end_index > sz)
                end_index = sz;
        }
        for (int i = 0; i < thr_num; i++) {
            thr[i].join();
        }

        delete [] yb;
        delete [] ya;
        delete [] z_1;

        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * (16 + L_BIT - 1));
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer2(), sz * 16);

        ptr = proxy->GetBuffer2();
        uint64_t *m = new uint64_t[sz];
        uint64_t val[2];
        for (int i = 0; i < sz; i++) {
            val[0] = ConvertToLong(&ptr);
            val[1] = ConvertToLong(&ptr);
            m[i] = val[f];
        }
        return m;

    }
    else if (proxy->GetPRole() == helper) {
        unsigned char *ptr_out = proxy->GetBuffer1();
        unsigned char *ptr_out2 = proxy->GetBuffer2();
        uint8_t *w = new uint8_t [sz];
        for (int i = 0; i < sz; i++) {
            uint64_t y = proxy->GenerateRandom() & N1_MASK;
            uint64_t ya_1 = proxy->GenerateRandom() & N1_MASK;
            uint64_t ya_2 = (y - ya_1) & N1_MASK;
            AddValueToCharArray(ya_1, &ptr_out);
            AddValueToCharArray(ya_2, &ptr_out2);
            for (int j = 0; j < L_BIT - 1; j++) {
                uint8_t k = (y >> j) & 0x1;
                uint8_t yb_1 = proxy->GenerateRandomByte() % 0x3f;
                uint8_t yb_2 = LP - yb_1 + k;
                AddValueToCharArray(yb_1, &ptr_out);
                AddValueToCharArray(yb_2, &ptr_out2);
            }
            w[i] = 0;
            if (y<ya_1)
                w[i] = 1;
        }
        thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * (8 + L_BIT - 1));
        thread thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * (8 + L_BIT - 1));
        thr1.join();
        thr2.join();

        thr1 = thread(Receive, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * (16 + L_BIT - 1));
        thr2 = thread(Receive, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * (16 + L_BIT - 1));
        thr1.join();
        thr2.join();


        unsigned char *ptr = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();
        ptr_out = proxy->GetBuffer1();
        ptr_out2 = proxy->GetBuffer2();

        int L1 = L_BIT-1;
        int jk = 0;
        for (int j = 0; j < sz; j++) {
            uint8_t res = 0;
            for (int i = 0; i < L1; i++) {
                uint8_t tmp = (proxy->GetBuffer1()[jk + i] + proxy->GetBuffer2()[jk + i]) % LP;
                if (((int) tmp) == 0) {
                    res = 1;
                    break;
                }
            }
            jk += L1;
            ptr += L1;
            ptr2 += L1;


            uint64_t val1 = (ConvertToLong(&ptr) + ConvertToLong(&ptr2)-(w[j]^res)*N1)/N1;
            uint64_t val2 = (ConvertToLong(&ptr) + ConvertToLong(&ptr2)-(w[j]^res)*N1)/N1;
            jk += 16;
            if(format) {
                val1 = ConvertToUint64((double) val1, shift);
                val2 = ConvertToUint64((double) val2, shift);
            }
            uint64_t vs_1 = proxy->GenerateRandom();
            uint64_t vs_2 = (val1 - vs_1);
            AddValueToCharArray(vs_1, &ptr_out);
            AddValueToCharArray(vs_2, &ptr_out2);
            vs_1 = proxy->GenerateRandom();
            vs_2 = (val2 - vs_1);
            AddValueToCharArray(vs_1, &ptr_out);
            AddValueToCharArray(vs_2, &ptr_out2);
        }

        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 16);
        thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 16);
        thr1.join();
        thr2.join();

        delete [] w;

        return 0;
    }
    return 0;
}


/** Most significant bit: Returns the first (=left-most) bit of @p x.
 *
 * @param x
 * @return The first bit of @p x
 */
uint64_t MostSignificantBit(Party *const proxy, uint64_t x, int shift=FRACTIONAL_BITS) {
    if (proxy->GetPRole() == helper) {
        MostSignificantBit(proxy, nullptr, 1, shift);
        return 0;
    } else {
        auto result_array = MostSignificantBit(proxy, &x, 1, shift);
        uint64_t result = result_array[0];
        delete[] result_array;
        return result;
    }
}

/**
 * @brief Compares all values of x to all values of y.
 *
 * In the resulting array, the first y_size results are all comparisons between the first x value and y and so on.
 * @param proxy p_proxy
 * @param x p_x: an array of shared values
 * @param y p_y:an array of shared values
 * @param x_size p_x_size: the amount of values in x
 * @param y_size p_y_size: the amount of values in y
 * @param shift p_shift: How many bits are reserved for the decimal places. Defaults to FRAC.
 * @return uint64_t* all comparison results between x and y
 */
uint64_t *CompareAll(
    Party *const proxy,
    const uint64_t *const x,
    const uint64_t *const y,
    int x_size,
    int y_size,
    int shift = FRACTIONAL_BITS
) {
    int overall_size = x_size*y_size;
    if (proxy->GetPRole() != helper) {
        uint64_t* difference = new uint64_t[overall_size];
        int overall_index;
        for (int x_index = 0; x_index < x_size; x_index++) {
            for (int y_index = 0; y_index < y_size; y_index++) {
                overall_index = x_index * y_size + y_index;
                difference[overall_index] = x[x_index] - y[y_index];
            }
        }
        uint64_t *result = MostSignificantBit(proxy, difference, overall_size, shift);
        delete[] difference;
        for (int i = 0; i < overall_size; i++) {
            result[i] = (proxy->GetPRole()<<shift) - (result[i] << shift);
        }
        return result;
    } else { // HELPER
        MostSignificantBit(proxy, nullptr, overall_size, shift);
        return nullptr;
    }
}

/** Comparison between two numbers.
 *
 * @param proxy
 * @param x
 * @param y
 * @return @p x > @p y
 */
uint64_t *Compare(Party *const proxy, const uint64_t *const x, const uint64_t *const y, uint32_t sz, int shift = FRACTIONAL_BITS) {
    if ( proxy->GetPRole() != helper) {
        uint64_t* diff = new uint64_t[sz];
        for (int i = 0; i < sz; i++) {
            diff[i] = y[i] - x[i];
        }
        uint64_t* m = MostSignificantBit(proxy, diff, sz, shift);
        delete[] diff;
        return m;
    }else if ( proxy->GetPRole() == helper) {
        MostSignificantBit(proxy, nullptr, sz, shift);
        return nullptr;
    }
    return nullptr;
}

/** Comparison between two numbers.
 *
 * @param proxy
 * @param x
 * @param y
 * @return @p x > @p y
 */
uint64_t Compare(Party *const proxy, uint64_t x, uint64_t y, int shift = FRACTIONAL_BITS) {
    if ( proxy->GetPRole() != helper) {
        uint64_t *result_array = Compare(proxy, &x, &y, 1, shift);
        uint64_t result = result_array[0];
        delete[] result_array;
        return  result;
    }else { // helper
        Compare(proxy, nullptr, nullptr, 1, shift);
        return 0;
    }
}


uint64_t* Equals(Party *const proxy, const uint64_t *const x, const uint64_t *const y, uint32_t size, int shift = FRACTIONAL_BITS) {
    if (proxy->GetPRole() != helper) {
        uint64_t *xyx = new uint64_t[size*3];
        std::copy(x, x+size, xyx);
        std::copy(y, y+size, xyx+size);
        std::copy(x, x+size, xyx+2*size);
        uint64_t* greater_and_smaller = Compare(proxy, xyx, xyx+size, size*2, shift);
        delete[] xyx;
        auto m = new uint64_t[size];
        for (int i = 0; i < size; i++) {
                m[i] = ConvertToUint64(0.5, shift) -greater_and_smaller[i] -greater_and_smaller[i+size];
        }
        delete[] greater_and_smaller;
        return m;
    } else {
        Compare(proxy, nullptr, nullptr, size*2, shift);
        return nullptr;
    }
}

uint64_t Equals(Party *const proxy, uint64_t x, uint64_t y, int shift = FRACTIONAL_BITS){
    if (proxy->GetPRole() == helper) {
        Equals(proxy, nullptr, nullptr, 1, shift);
        return 0;
    } else {
        uint64_t *result_array = Equals(proxy, &x, &y, 1, shift);
        uint64_t result = result_array[0];
        delete[] result_array;
        return result;
    }
}

/** Multiplication of two arrays of numbers.
 *
 * @param a one of the vectors of shares of the multiplicands
 * @param b the other vector of shares of the multiplicands
 * @param size the size of the vectors @p a and @p b
 * @return a vector containing the share of the result of the multiplication
 */
uint64_t *Multiply(Party *const proxy, const uint64_t *const a, const uint64_t *const b, uint32_t size, int shift = FRACTIONAL_BITS) {
    if (DEBUG_FLAG >= 1)
        cout << "************************************************************\nPMNF_MUL is called" << endl;
    if (proxy->GetPRole() == helper) {
        uint64_t *c1 = new uint64_t[size];

        GenerateMultiplicationTriple(proxy, c1, size);

        unsigned char *ptr_out2 = proxy->GetBuffer2();
        for (int j = 0; j < size; j++) {
            AddValueToCharArray(c1[j], &ptr_out2);
        }

        Send( proxy->GetSocketP2(), proxy->GetBuffer2(), size * 8);

        delete[] c1;

        if (DEBUG_FLAG >= 1)
            cout << "Returning from PMNF_MUL...\n************************************************************" << endl;
        return 0;
    } else if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *mt[3];
        mt[0] = new uint64_t[size]; //a
        mt[1] = new uint64_t[size]; //b
        mt[2] = new uint64_t[size]; //c
        uint64_t *concat_e_f = new uint64_t[size * 2];
        if (proxy->GetPRole() == proxy2) {
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), size * 8);
            unsigned char *ptr = proxy->GetBuffer1();
            for (int i = 0; i < size; ++i) {
                mt[0][i] = proxy->GenerateCommonRandom2();
                mt[1][i] = proxy->GenerateCommonRandom2();
                mt[2][i] = ConvertToLong(&ptr);
                concat_e_f[i] = a[i] - mt[0][i];
                concat_e_f[i + size] = b[i] - mt[1][i];
            }
        }
        else {
            for (int i = 0; i < size; ++i) {
                mt[0][i] = proxy->GenerateCommonRandom2();
                mt[1][i] = proxy->GenerateCommonRandom2();
                mt[2][i] = proxy->GenerateCommonRandom2();

                concat_e_f[i] = a[i] - mt[0][i];
                concat_e_f[i + size] = b[i] - mt[1][i];
            }
        }

        uint64_t *e_f = Reconstruct(proxy, concat_e_f, size * 2);
        uint64_t *e = e_f;
        uint64_t *f = &e_f[size];

        uint64_t *z = new uint64_t[size];
        for (int i = 0; i < size; i++) {
            z[i] = proxy->GetPRole() * e[i] * f[i] + f[i] * mt[0][i] + e[i] * mt[1][i] + mt[2][i];
            z[i] = Truncate(proxy, z[i], shift);
        }
        delete [] e_f;
        delete [] concat_e_f;
        for (auto &i : mt) {
            delete[] i;
        }
        if(DEBUG_FLAG >= 1)
            cout << "Returning from PMNF_MUL...\n************************************************************" << endl;
        return z;
    } else {
        return nullptr;
    }
}

 /** Multiplication of two numbers.
  *
  * @param proxy
  * @param a a share of the first multiplicand
  * @param b a share of the second multiplicand
  * @return the share of the multiplication of @p a and @p b
  */
uint64_t Multiply(Party *const proxy, uint64_t a, uint64_t b, int shift = FRACTIONAL_BITS) {
    if (proxy->GetPRole() == helper) {
        Multiply(proxy, nullptr, nullptr, 1, shift);
        return 0;
    } else {
        uint64_t *result_array = Multiply(proxy, &a, &b, 1, shift);
        uint64_t result = result_array[0];
        delete[] result_array;
        return result;
    }
}


/** Multiple exponentials. Note that this function considers only the specific number of least significant bits not to
 * cause overflow. This is different for positive and negative powers.
 *
 * @param a the vector of values that will be used as the power of exp
 * @param size the length of @p a
 * @return a vector of arithmetic secret shares for each exp(@p a)
 */
uint64_t* Exp(Party *const proxy, const uint64_t *const a, uint32_t size, int shift = FRACTIONAL_BITS) {
    int p_role = proxy->GetPRole();
    int n_bits = proxy->GetNBits();
    int neg_n_bits = proxy->GetNegNBits();

    if (p_role == proxy1 || p_role == proxy2) {
        // compute the absolute of the input value
        uint64_t* msb_a = MostSignificantBit(proxy, a, size, shift);
        uint64_t* abs_a = new uint64_t[size];
        for(uint32_t i = 0; i < size; i++) {
            abs_a[i] = ((uint64_t) 0) - a[i];
        }

        // compute the possible contribution of positive and negative values
        uint64_t* pec = new uint64_t[n_bits];
        uint64_t* nec = new uint64_t[n_bits];
        if(p_role == proxy2) {
            for (int i = n_bits - 1; i >= 0; i--) {
                pec[n_bits - i - 1] = ConvertToUint64(exp(pow(2, i - shift)));
                if (i > neg_n_bits - 1) {
                    nec[n_bits - i - 1] = (((uint64_t) 1) << shift);
                } else {
                    nec[n_bits - i - 1] = ConvertToUint64(1.0 / exp(pow(2, i - shift)));
                }
            }
        }
        else {
            for (int i = n_bits - 1; i >= 0; i--) {
                pec[n_bits - i - 1] = 0;
                nec[n_bits - i - 1] = 0;
            }
        }

        // selection of the correct contribution from each bit of the input value based on the msb of the input value
        uint64_t *pos_e_contributions = new uint64_t[size * (n_bits + 1)]; // if the value is positive
        uint64_t *neg_e_contributions = new uint64_t[size * (n_bits + 1)]; // if the value is negative
        uint64_t *one_contributions = new uint64_t[size * n_bits]; // if the bit is zero regardless of the sign
        uint64_t *repeated_msb_a = new uint64_t[size * (n_bits + 1)]; // this will be used as a selection bit for the contributions of all bits

        for(uint32_t i = 0; i < size; i++) {
            pos_e_contributions[i * (n_bits + 1)] = a[i];
            neg_e_contributions[i * (n_bits + 1)] = abs_a[i];
            repeated_msb_a[i * (n_bits + 1)] = msb_a[i];
            for (int bi = 0; bi < n_bits; bi++) {
                pos_e_contributions[(i * (n_bits + 1)) + bi + 1] = pec[bi];
                neg_e_contributions[(i * (n_bits + 1)) + bi + 1] = nec[bi];
                one_contributions[(i * n_bits) + bi] = p_role * (((uint64_t) 1) << shift);
                repeated_msb_a[(i * (n_bits + 1)) + bi + 1] = msb_a[i];
            }
        }
        delete[] msb_a;
        delete[] abs_a;
        delete[] pec;
        delete[] nec;
        uint64_t *e_contributions = Multiplex(proxy, pos_e_contributions, neg_e_contributions, repeated_msb_a,
                                              size * (n_bits + 1), shift);
        delete[] pos_e_contributions;
        delete[] neg_e_contributions;
        delete[] repeated_msb_a;
        uint64_t* new_a = new uint64_t[size];
        uint64_t* selected_e_contributions = new uint64_t[size * n_bits];
        for(uint32_t i = 0; i < size; i++) {
            new_a[i] = e_contributions[i * (n_bits + 1)];
            for(uint32_t j = 0; j < n_bits; j++) {
                selected_e_contributions[(i * n_bits) + j] = e_contributions[(i * (n_bits + 1)) + j + 1];
            }
        }
        delete[] e_contributions;
        // arrange all the shifted versions of the input value for MostSignificantBit
        uint64_t *partial_a = new uint64_t[size * n_bits];
        for(uint32_t i = 0; i < size; i++) {
            for (uint32_t j = 0; j < n_bits; j++) {
                partial_a[(i * n_bits) + j] = new_a[i] << (L_BIT - n_bits + j);
            }
        }

        // get secret shared form of the bits of the values that could contribute into the result
        uint64_t *bit_shares = MostSignificantBit(proxy, partial_a, size * n_bits, shift);
        delete[] partial_a;
        // selection of the contribution of the bits of the value
        uint64_t *contributions = Multiplex(proxy, one_contributions, selected_e_contributions, bit_shares,
                                            size * n_bits, shift);
        delete[] one_contributions;
        delete[] bit_shares;
        // binary-tree-based multiplication of the contributions into the BenchmarkExp
        int cs = n_bits;
        bool flag = false;
        uint64_t* remaining = new uint64_t[size];
        uint64_t *tmp1, *tmp2;
        for (int j = 0; j < (int) ceil(log2(n_bits)); j++) {
            tmp1 = contributions;
            tmp2 = &contributions[cs / 2];

            if (cs % 2 == 1) {
                if (!flag) {
                    for(uint32_t i = 0; i < size; i++){
                        remaining[i] = contributions[(i * cs) + cs - 1];
                    }

                    tmp1 = new uint64_t[size * (cs / 2)];
                    tmp2 = new uint64_t[size * (cs / 2)];

                    for( uint32_t i = 0; i < size; i++) {
                        copy(contributions + (i * cs), contributions + (i * cs) + (cs / 2), tmp1 + (i * (cs / 2)));
                        copy(contributions + (i * cs) + (cs / 2), contributions + (i * cs) + 2 * (cs / 2), tmp2 + (i * (cs / 2)));
                    }

                    flag = true;
                } else {
                    tmp1 = new uint64_t[size * ((cs + 1) / 2)];
                    tmp2 = new uint64_t[size * ((cs + 1) / 2)];

                    for(uint32_t i = 0; i < size; i++) {
                        copy(contributions + (i * cs), contributions + (i * cs) + ((cs + 1) / 2), tmp1 + (i * ((cs + 1) / 2)));
                        copy(contributions + (i * cs) + ((cs + 1) / 2), contributions + ((i + 1) * cs), tmp2 + (i * ((cs + 1) / 2)));
                        tmp2[(i + 1) * ((cs + 1) / 2) - 1] = remaining[i];
                    }

                    cs++;
                    flag = false;
                }
            }
            else {
                tmp1 = new uint64_t[size * (cs / 2)];
                tmp2 = new uint64_t[size * (cs / 2)];

                for( uint32_t i = 0; i < size; i++) {
                    copy(contributions + (i * cs), contributions + (i * cs) + (cs / 2), tmp1 + (i * (cs / 2)));
                    copy(contributions + (i * cs) + (cs / 2), contributions + ((i + 1) * cs), tmp2 + (i * (cs / 2)));
                }
            }
            delete [] contributions;
            contributions = Multiply(proxy, tmp1, tmp2, size * (cs / 2), shift);

            delete [] tmp1;
            delete [] tmp2;

            cs /= 2;
        }

        // deleting dynamically allocated arrays
        delete [] remaining;

        return contributions;
    }
    else if ( p_role == helper) {
        MostSignificantBit(proxy, 0, size, shift);
        Multiplex(proxy, 0, 0, 0, size * (n_bits + 1), shift);
        MostSignificantBit(proxy, 0, size * n_bits, shift);
        Multiplex(proxy, 0, 0, 0, size * n_bits, shift);
        int current_size = n_bits;
        bool flag = false;
        for (int i = 0; i < (int) ceil(log2(n_bits)); i++) {
            if (current_size % 2 == 1) {
                if (!flag) {
                    flag = true;
                } else {
                    current_size++;
                    flag = false;
                }
            }
            Multiply(proxy, 0, 0, size * (current_size / 2), shift);
            current_size /= 2;
        }

        return 0;
    }
    else {
        return nullptr;
    }
}

/** Exponential. Note that this function considers only the specific number of least significant bits not to cause
 * overflow. This is different for positive and negative powers.
 *
 * @param a the value that will be used as the power of exp
 * @return Returns the arithmetic secret share of exp(@p a)
 */
uint64_t Exp(Party *const proxy, uint64_t a, int shift = FRACTIONAL_BITS) {
    if (proxy->GetPRole() == helper) {
        Exp(proxy, nullptr, 1, shift);
        return 0;
    } else {
        uint64_t* result_array = Exp(proxy, &a, 1, shift);
        uint64_t result = result_array[0];
        delete[] result_array;
        return result;
    }
}

/** PartialSum: sum the elements of each section separately.
 *
 * @param a the vector on which we perform the partial summation
 * @param size the size of @p a
 * @param d the size of the part that we will use to Partition @p a
 * @return
 */
uint64_t* PartialSum(Party *const proxy, const uint64_t *const a, uint32_t size, uint32_t d) {
    int p_role = proxy->GetPRole();
    if(p_role == proxy1 || p_role == proxy2) {
        uint64_t *ps_x = new uint64_t[size / d];
        for(uint32_t base = 0; base < size; base += d) {
            uint64_t tmp = 0;
            for(uint32_t i = 0; i < d; i++) {
                tmp += a[base + i];
            }
            ps_x[base / d] = tmp;
        }
        return ps_x;
    }
    else {
        return NULL;
    }
}

/** Computes the dot product of arithmetically shared vectors, which are formed by vectors.
 *
 * @param a vector formed by vectors of given size
 * @param b vector formed by vectors of given size
 * @param size the length of the vectors
 * @param d the size of the partial vectors forming the main vectors
 * @return Dot product of the given (@p size / @p d) vectors as a vector of (@p size / @p d)
 */
uint64_t* DotProduct(Party *const proxy, const uint64_t *const a, const uint64_t *const b, uint32_t size, uint32_t d, int shift = FRACTIONAL_BITS) {
    int p_role = proxy->GetPRole();
    if(p_role == proxy1 || p_role == proxy2) {
        // compute elementwise multiplication of vectors
        uint64_t *ew_xy = Multiply(proxy, a, b, size, shift);
        // sum the vectors in the main vector
        uint64_t *dp_shr = PartialSum(proxy, ew_xy, size, d);

        delete [] ew_xy;

        return dp_shr;
    }
    else if(p_role == helper) {
        Multiply(proxy, 0, 0, size, shift);
        return NULL;
    }
    else {
        return NULL;
    }

}

/** computes the dot product of two single arithmetically shared vectors.
 *
 * @param proxy
 * @param a vector
 * @param b vector
 * @param size the length of the vectors
 * @return
 */
uint64_t DotProduct(Party *const proxy, const uint64_t *const a, const uint64_t *const b, uint32_t size, int shift = FRACTIONAL_BITS) {
    if (proxy->GetPRole() == helper) {
        DotProduct(proxy, nullptr, nullptr, size, size, shift);
        return 0;
    } else {
        uint64_t *result_vector = DotProduct(proxy, a, b, size, size, shift);
        uint64_t result = result_vector[0];
        delete[] result_vector;
        return result;
    }
}


// TODO test
/** Get the Modular Inverse (ModularInverse) of a given number a with modulo being the specified ring size.
 * For the resulting/returned value b, it must hold ab Mod(modulo) are congruent to 1. The modulo under which a
 * multiplied with the inverse are equal to 1, will always be the ring size.
 * @param a secret share of the value for which the modular inverse shall be calculated.
 * @return the secret share of the modular inverse of a under the ring size.
 */
uint64_t ModularInverse(Party *const proxy, uint64_t a){
    cout << "searching for ModularInverse of value " << ConvertToDouble(a) << endl;
    uint64_t exchangingBit = RING_SIZE / 64;
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        for (uint16_t step = 0; step < 64; step++) {
            cout << "step " << step << endl;
            uint64_t ringProducts [exchangingBit];
            // start with 1 because 0 does not have an inverse value.
            for (uint64_t x = 1; x <= exchangingBit; x++) {
                uint64_t modInv = x*step + x;
                ringProducts[x - 1] = (a * modInv) & RING_SIZE; // ModularConversion(proxy, t); ?
            }
            cout << "stored all ring products..." << endl;
            unsigned char *ptr = proxy->GetBuffer1();
            AddValueToCharArray(ringProducts, &ptr, exchangingBit);
            Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), exchangingBit * 8);

            cout << "sent ring products to helper" << endl;
            // receive fresh share from helper
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), 8); // either the share of modInv or -1 to identify, to continue searching
            ptr = proxy->GetBuffer1();
            uint64_t share = ConvertToLong(&ptr);
            cout << "got fresh share from helper: " << share << endl;
            if (share != -1){
                // the modInv has been found
                cout << "ModularInverse was found: " << share << endl;
                return share;
            }
        }
        return 0;
    }
    else if (proxy->GetPRole() == helper) {
        for (uint16_t step = 0; step < 64; step++) {
            cout << "step " << step << endl;
            Receive(proxy->GetSocketP1(), proxy->GetBuffer1(), exchangingBit * 8);
            Receive(proxy->GetSocketP2(), proxy->GetBuffer2(), exchangingBit * 8);
            unsigned char *ptr1 = proxy->GetBuffer1();
            unsigned char *ptr2 = proxy->GetBuffer2();
            cout << "got ring products from parties..." << endl;

            uint64_t ringProducts_recon[exchangingBit];
            ringProducts_recon[0] = (ConvertToLong(&ptr1) + ConvertToLong(&ptr2)); //modInv = exchangeBit*step + 1
            uint64_t m;
            for(uint64_t i = 1; i < exchangingBit; i++){
                // reconstructed product of a was: exchangeBit * step + i+1
                ringProducts_recon[i] = (ConvertToLong(&ptr1) + ConvertToLong(&ptr2));
                for(uint64_t j = 0; j < i; j++){
                    if(((ringProducts_recon[j] + ringProducts_recon[i]) & RING_SIZE) == 1){
                        //Mod inverse of a is found: i+1 + j+1
                        m = exchangingBit * 2 * step + i + j + 2;
                        cout << "ModularInverse was found: " << m << endl;
                        // SEND fresh share of found modular inverse
                        //reassign buffer because ptr1 and ptr2 were incremented by convert2Long calls.
                        ptr1 = proxy->GetBuffer1();
                        ptr2 = proxy->GetBuffer2();

                        uint64_t tmp = proxy->GenerateRandom();
                        AddValueToCharArray(tmp, &ptr1);
                        AddValueToCharArray(m - tmp, &ptr2);

                        thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), 8);
                        thread thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), 8);
                        thr1.join();
                        thr2.join();
                        cout << "sent fresh share of ModularInverse to parties; m= " << m << endl;
                        return 0;
                    }
                }
            }
            //reassign buffer because ptr1 and ptr2 were incremented by convert2Long calls.
            ptr1 = proxy->GetBuffer1();
            ptr2 = proxy->GetBuffer2();

            uint64_t noValFound = -1;
            AddValueToCharArray(noValFound, &ptr1);
            AddValueToCharArray(noValFound, &ptr2);

            thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), 8);
            thread thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), 8);
            thr1.join();
            thr2.join();
        }
        return 0;
    }
    return -1;
}


// TODO test
/** Compute the vectorized division of a / b where a and b are vectors - not the integer approximation of the result
 *
 * @param proxy : Party instance
 * @param a : vector of dividends
 * @param b : vector of dividers
 * @param size: number of division operations - which is the size of a and b
 * @param first_call : indicates whether the aucDivide call is for the integer part of the division result, i.e. first call.
 * If it is the first call, then there will be the second call of aucDivide for the fractional part of the division result
 * @return vector (a / b)
 */
uint64_t* Divide(Party *const proxy, const uint64_t *a, const uint64_t *b, uint32_t size, int shift = FRACTIONAL_BITS, bool first_call = true) {
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *signs;
        if(first_call) {
            uint64_t *inp1 = new uint64_t[2 * size];
            uint64_t *inp2 = new uint64_t[2 * size];
            for(int i = 0; i < size; i++) {
                inp1[i] = a[i];
                inp1[size + i] = b[i];
                inp2[i] = (uint64_t) 0 - a[i];
                inp2[size + i] = (uint64_t) 0 - b[i];
            }
            signs = MostSignificantBit(proxy, inp1, 2 * size, shift);
            uint64_t *abs_vals = Multiplex(proxy, inp1, inp2, signs, 2 * size, shift);
            a = &abs_vals[0];
            b = &abs_vals[size];

            delete [] inp1;
            delete [] inp2;
        }

        // initialize the variables for quotient and remainder and the Role vector, which is a vector full of the Role value
        uint64_t *Q = new uint64_t[size];
        uint64_t *R = new uint64_t[size];

        // obtain every bit of the dividend
        uint64_t *msb_bits_of_a = new uint64_t[L_BIT * size];

        for(int i = 0; i < size; i++) { // each value
            Q[i] = 0;
            R[i] = 0;
            uint64_t tmp = a[i];
            for(int j = 0; j < L_BIT; j++) { // each bit of the value
                msb_bits_of_a[i * L_BIT + j] = tmp;
                tmp = tmp << 1;
            }
        }
        uint64_t *bits_of_a = MostSignificantBit(proxy, msb_bits_of_a, L_BIT * size, shift, false);

        delete [] msb_bits_of_a;

        // traverse all bits of the dividend
        for (int16_t j = L_BIT - 1; j >= 0; j--) {
            for(int i = 0; i < size; i++) {
                R[i] = R[i] << 1; // shift the remainder
                R[i] += bits_of_a[(i * L_BIT) + (L_BIT - 1 - j)];
            }

            uint64_t *c = Compare(proxy, R, b, size, shift); // compare the current R and divider

            uint64_t *o1 = new uint64_t[2 * size];
            uint64_t *o2 = new uint64_t[2 * size];
            for(int i = 0; i < size; i++) {
                o1[2 * i] = c[i];
                o1[2 * i + 1] = c[i];
                o2[2 * i] = b[i];
                o2[2 * i + 1] = ((uint64_t) proxy->GetPRole()) << j;
            }

            // if the current R is larger than or equal to the divider, subtract the divider from R
            uint64_t *v = Multiply(proxy, o1, o2, 2 * size, shift);
            for(int i = 0; i < size; i++) {
                R[i] = R[i] - v[2 * i];
                Q[i] = Q[i] + v[2 * i + 1];
            }
            delete [] c;
            delete [] o1;
            delete [] o2;
            delete [] v;
        }

        delete [] bits_of_a;

        if(first_call) {
            // determine the selection bits for the signs of the results based on the signs of a's and b's
            // choose the positive result if a < 0 and b < 0, or a >= 0 and b >= 0
            // choose the negative result if a >= 0 and b < 0, or a < 0 and b >= 0
            // This is exactly what XOR does. We mimic XOR arithmetically, i.e. a XOR b = a + b - 2ab
            uint64_t *tmp = Multiply(proxy, signs, &signs[size], size, shift); // for determining the signs of the results
            uint64_t *c = new uint64_t[size]; // for determining the signs of the results - selection bits
            for(int i = 0; i < size; i++) {
                R[i] = R[i] << shift; // prepare the remainder for the second division call
                Q[i] = Q[i] << shift; // prepare the quotient for the final quotient
                c[i] = (signs[i] + signs[i + size]) - 2 * tmp[i]; // for determining the signs of the results
            }
            delete [] tmp;

            uint64_t *neg_Q = new uint64_t[size]; // the negative of the results in case they are the correct ones based on signs
            uint64_t *sec_div = Divide(proxy, R, b, size, shift, false); // second division call for the fractional part of the final quotient
            for(int i = 0; i < size; i++) {
                Q[i] += sec_div[i];
                neg_Q[i] = (uint64_t) 0 - Q[i];
            }

            delete [] sec_div;
            delete [] signs;

            // based on the above analysis, we select the correct version of the final quotient
            uint64_t *div_res = Multiplex(proxy, Q, neg_Q, c, size, shift);
            delete [] c;
            delete [] neg_Q;
            delete [] Q;
            delete [] R;
            delete [] a;
            return div_res;
        }

        delete [] R;

        return Q;
    }
    else if (proxy->GetPRole() == helper) {
        if(first_call) {
            MostSignificantBit(proxy, 0, 2 * size, shift);
            Multiplex(proxy, 0, 0, 0, 2 * size, shift);
        }

        MostSignificantBit(proxy, 0, L_BIT * size, shift, false);

        for (int16_t i = L_BIT - 1; i >= 0; i--) {
            Compare(proxy, 0, 0, size, shift);
            Multiply(proxy, 0, 0, 2 * size, shift);
        }

        if(first_call) {
            Multiply(proxy, 0, 0, size, shift);
            Divide(proxy, 0, 0, size, shift, false);
            Multiplex(proxy, 0, 0, 0, size, shift);
        }
        return NULL;
    }
    return NULL;
}

/** Compute the division of a / b - not the integer approximation of the result
 *
 * @param proxy : Party instance
 * @param a : dividend
 * @param b : divider
 * @param first_call : indicates whether the aucDivide call is for the integer part of the division result, i.e. first call.
 * If it is the first call, then there will be the second call of aucDivide for the fractional part of the division result
 * @return a / b
 */
uint64_t Divide(Party *const proxy, uint64_t a, uint64_t b, int shift = FRACTIONAL_BITS) {
    if (proxy->GetPRole() == helper) {
        Divide(proxy, nullptr, nullptr, 1, shift);
        return 0;
    } else {
        uint64_t* result_array = Divide(proxy, &a, &b, 1, shift);
        uint64_t result = result_array[0];
        delete[] result_array;
        return result;
    }
}

// TODO test
/** Perform division operation, or more specifically normalization operation, of two given inputs. The operation is
 * taken from SecureNN, but it is implemented by using the building blocks of CECILIA. Note that there is an implicit
 * assumption for Normalize to work correctly: the elements of a must be less than the corresponding elements of b.
 *
 * @param proxy
 * @param a: the nominators
 * @param b: the denominators
 * @param size: the number of elements in a and b
 * @return div: uint64_t vector consisting of elementwise division of a/b
 */
uint64_t* Normalise(Party *const proxy, const uint64_t *const a, const uint64_t *const b, uint32_t size, int shift = FRACTIONAL_BITS) {
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *u = new uint64_t[size]; // holds how much needs to be subtracted from the nominator
        uint64_t *div = new uint64_t[size]; // resulting division
        for(int i = 0; i < size; i++) {
            u[i] = 0;
            div[i] = 0;
        }

        // iterate every bit of the fractional part to determine whether they are 1 or 0
        for(int i = 1; i <= shift; i++) {
            // compute the possible remaining of the nominator after subtracting denominator and previously subtracted value
            uint64_t *z = new uint64_t[size];
            for(int j = 0; j < size; j++) {
                z[j] = ((a[j] - u[j]) << i) - b[j];
            }

            uint64_t *msb_z = MostSignificantBit(proxy, z, size, shift);
            delete [] z;

            uint64_t *concat_cont_and_subt = new uint64_t[size * 2];
            uint64_t *twice_msb_z = new uint64_t[size * 2];
            for(int j = 0; j < size; j++) {
                twice_msb_z[j] = (proxy->GetPRole() << shift) - msb_z[j];
                twice_msb_z[j + size] = twice_msb_z[j];
                concat_cont_and_subt[j] = proxy->GetPRole() << (shift - i); // the contribution to the division result
                concat_cont_and_subt[j + size] = Truncate(proxy, b[j], i); // what to subtract from the nominator
            }
            delete [] msb_z;

            // computes possibly what to subtract and what to add & determines if we need to perform those operations
            uint64_t *tmp = Multiply(proxy, twice_msb_z, concat_cont_and_subt, 2 * size, shift);
            delete [] concat_cont_and_subt;
            delete [] twice_msb_z;

            for(int j = 0; j < size; j++) {
                div[j] = div[j] + tmp[j];
                u[j] = u[j] + tmp[j + size];
            }
            delete [] tmp;
        }

        delete [] u;
        return div;
    }
    else if (proxy->GetPRole() == helper) {
        for(int i = 1; i <= shift; i++) {
            MostSignificantBit(proxy, 0, size, shift);
            Multiply(proxy, 0, 0, 2 * size, shift);
        }
    }
    return nullptr;

}


/** Perform several multiplications of matrices of size a_row-by-a_col and a_col-by-b_col stored in a and b.
 *
 * @param a three dimensional matrix
 * @param b three dimensional matrix
 * @param n_matrices number of two-dimensional matrices of @p a and @p b
 * @param a_row number of rows per two-dimensional matrix
 * @param a_col number of columns per two-dimensional matrix of @p a
 * @param b_col number of columns per two-dimensional matrix of @p b
 * @return a matrix of size @p n_matrices by @p a_row by @p b_col
 */
uint64_t*** MatrixMatrixMultiply(
    Party *const proxy,
    const uint64_t *const *const *const a,
    const uint64_t *const *const *const b,
    uint32_t n_matrices,
    uint32_t a_row,
    uint32_t a_col,
    uint32_t b_col,
    int shift= FRACTIONAL_BITS
) {
    int p_role = proxy->GetPRole();
    if (p_role == proxy1 || p_role == proxy2) {
        cout <<"içeride proxy başladı" << endl;
        // templates for E and F matrices
        uint64_t ***E = new uint64_t**[n_matrices];
        uint64_t ***F = new uint64_t**[n_matrices];

        // receive the shares of A, B and C matrices
        if (p_role == proxy2) {
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(),
                    n_matrices * a_row * b_col * 8);
        }
        unsigned char *ptr = proxy->GetBuffer1();

        uint64_t ***mt2 = new uint64_t**[n_matrices];

        for(int i = 0; i < n_matrices; i++) {
            // <X-A>_1
            E[i] = new uint64_t *[a_row];
            for(int j = 0; j < a_row; j++) {
                E[i][j] = new uint64_t[a_col];
                for(int k = 0; k < a_col; k++) {
                    E[i][j][k] = a[i][j][k] - proxy->GenerateCommonRandom2();
                }
            }

            // <Y-B>_1
            F[i] = new uint64_t *[a_col];
            for(int j = 0; j < a_col; j++) {
                F[i][j] = new uint64_t[b_col];
                for(int k = 0; k < b_col; k++) {
                    F[i][j][k] = b[i][j][k] - proxy->GenerateCommonRandom2();
                }
            }

            if (p_role == proxy1) { // proxy1 generates its share of C locally using the common random generator with the helper
                mt2[i] = new uint64_t *[a_row];
                for (int j = 0; j < a_row; j++) {
                    mt2[i][j] = new uint64_t[b_col];
                    for (int k = 0; k < b_col; k++) {
                        mt2[i][j][k] = proxy->GenerateCommonRandom2();
                    }
                }
            }
            else { // proxy2 gets its share of C from the helper
                ConvertTo2dArray(&ptr, mt2[i], a_row, b_col);
            }
        }

        // reconstruct E = X-A and F = Y-B
        uint64_t ***recE = Reconstruct(proxy, E, n_matrices, a_row, a_col);
        uint64_t ***recF = Reconstruct(proxy, F, n_matrices, a_col, b_col);

        for(int n = 0; n < n_matrices; n++) {
            for(int r = 0; r < a_row; r++) {
                delete [] E[n][r];
            }
            delete [] E[n];

            for(int r = 0; r < a_col; r++) {
                delete [] F[n][r];
            }
            delete [] F[n];
        }
        delete [] E;
        delete [] F;

        // compute -role * E * F + <X>_role * F + E * <Y>_role + <C>_role
        uint64_t ***z = new uint64_t**[n_matrices];
        uint64_t ***tmpEF;
        uint64_t ***tmpXF;
        uint64_t ***tmpEY;
        // the part below takes the majority of the time - it can be addressed for optimization
        tmpEF = LocalMatrixMatrixMultiply(recE, recF, n_matrices, a_row, a_col, b_col, 0);
        tmpXF = LocalMatrixMatrixMultiply(a, recF, n_matrices, a_row, a_col, b_col, 0);
        tmpEY = LocalMatrixMatrixMultiply(recE, b, n_matrices, a_row, a_col, b_col, 0);

        for(int n = 0; n < n_matrices; n++) {
            z[n] = new uint64_t * [a_row];
            for(int i = 0; i < a_row; i++) {
                z[n][i] = new uint64_t[b_col];
                for(int j = 0; j < b_col; j++) {
                    z[n][i][j] = (p_role * -1) * tmpEF[n][i][j] + tmpXF[n][i][j] + tmpEY[n][i][j] + mt2[n][i][j];
                    z[n][i][j] = Truncate(proxy, z[n][i][j], shift);
                }
            }
        }

        for(int n = 0; n < n_matrices; n++) {
            for(int r = 0; r < a_row; r++) {
                delete [] recE[n][r];
                delete [] tmpEF[n][r];
                delete [] tmpXF[n][r];
                delete [] tmpEY[n][r];
                delete [] mt2[n][r];
            }
            delete [] recE[n];
            delete [] tmpEF[n];
            delete [] tmpXF[n];
            delete [] tmpEY[n];
            delete [] mt2[n];

            for(int r = 0; r < a_col; r++) {
                delete [] recF[n][r];
            }
            delete [] recF[n];
        }
        delete [] recE;
        delete [] recF;
        delete [] tmpEF;
        delete [] tmpXF;
        delete [] tmpEY;
        delete [] mt2;

        return z;
    }
    else if(p_role == helper) {
        cout <<"içeride helper başladı" << endl;
        unsigned char *ptr_out = proxy->GetBuffer1();

        // temporary matrices to hold the current A, B and C matrices
        uint64_t **tmpA = new uint64_t*[a_row];
        uint64_t **tmpA1 = new uint64_t*[a_row];
        uint64_t **tmpA2 = new uint64_t*[a_row];
        uint64_t **tmpB = new uint64_t*[a_col];
        uint64_t **tmpB1 = new uint64_t*[a_col];
        uint64_t **tmpB2 = new uint64_t*[a_col];
        uint64_t **tmpC1 = new uint64_t*[a_row];
        uint64_t **tmpC;

        for(int i = 0; i < a_row; i++) {
            tmpA[i] = new uint64_t[a_col];
            tmpA1[i] = new uint64_t[a_col];
            tmpA2[i] = new uint64_t[a_col];
            tmpC1[i] = new uint64_t[b_col];
        }

        for(int i = 0; i < a_col; i++) {
            tmpB[i] = new uint64_t[b_col];
            tmpB1[i] = new uint64_t[b_col];
            tmpB2[i] = new uint64_t[b_col];
        }

        uint64_t tmp; // to hold the generated random values

        // matrix generations
        for(uint32_t n = 0; n < n_matrices; n++) {
            // generation of A and its shares
            for(uint32_t i = 0; i < a_row; i++) {
                for(uint32_t j = 0; j < a_col; j++) {
                    tmpA1[i][j] = proxy->GenerateCommonRandom();
                    tmpA2[i][j] = proxy->GenerateCommonRandom2();
                    tmpA[i][j] = tmpA1[i][j] + tmpA2[i][j];
                }
            }

            // generation of B and its shares
            for(uint32_t i = 0; i < a_col; i++) {
                for(uint32_t j = 0; j < b_col; j++) {
                    tmpB1[i][j] = proxy->GenerateCommonRandom();
                    tmpB2[i][j] = proxy->GenerateCommonRandom2();
                    tmpB[i][j] = tmpB1[i][j] + tmpB2[i][j];
                }
            }

            // calculation of A * B
            tmpC = LocalMatrixMatrixMultiply(tmpA, tmpB, a_row, a_col, b_col, 0); // why shift=0?

            // generation of shares of C
            for(uint32_t i = 0; i < a_row; i++) {
                for(uint32_t j = 0; j < b_col; j++) {
                    tmp = proxy->GenerateCommonRandom();
                    AddValueToCharArray(tmpC[i][j] - tmp, &ptr_out);
                }
            }

            for(int i = 0; i < a_row; i++) {
                delete [] tmpC[i];
            }
            delete [] tmpC;
        }

        for(int i = 0; i < a_row; i++) {
            delete [] tmpA[i];
            delete [] tmpA1[i];
            delete [] tmpA2[i];
            delete [] tmpC1[i];
        }
        delete [] tmpA;
        delete [] tmpA1;
        delete [] tmpA2;
        delete [] tmpC1;
        for(int i = 0; i < a_col; i++) {
            delete [] tmpB[i];
            delete [] tmpB1[i];
            delete [] tmpB2[i];
        }
        delete [] tmpB;
        delete [] tmpB1;
        delete [] tmpB2;

        // send these matrices to proxy1 and proxy2, respectively
        Send(proxy->GetSocketP2(), proxy->GetBuffer1(), n_matrices * a_row * b_col * 8);

        return nullptr;
    }
    else {
        return nullptr;
    }
}

/** Perform multiplication of matrices a and b.
 * The function assumes that the number of columns of a equals to the number of rows of b.
 *
 * @param a two dimensional matrix
 * @param b two dimensional matrix
 * @param a_row number of rows of @p a and @p b
 * @param a_col number of columns of @p a
 * @param b_col number of columns of @p b
 * @return a matrix of size @p a_row by @p b_col
 */
uint64_t** MatrixMatrixMultiply(
    Party *const proxy,
    const uint64_t *const *const a,
    const uint64_t *const *const b,
    uint32_t a_row,
    uint32_t a_col,
    uint32_t b_col,
    int shift = FRACTIONAL_BITS
) {
    if (proxy->GetPRole() == helper) {
        cout << "helper başladı" << endl;
        MatrixMatrixMultiply(proxy, nullptr, nullptr, 1, a_row, a_col, b_col, shift);
        return nullptr;
    } else {
        cout << "proxy başladı" << endl;
        uint64_t ***result_array = MatrixMatrixMultiply(proxy, &a, &b, 1, a_row, a_col, b_col, shift);
        uint64_t **result = result_array[0];
        delete[] result_array;
        return result;
    }
}

/** Perform n_matrices multiplications of matrices of size a_row-by-a_col and vectors of size a_col stored in a and
 * b, respectively.
 *
 * @param a three dimensional matrix
 * @param b two dimensional matrix
 * @param n_matrices number of matrices in @p a / vectors in @p b
 * @param a_row number of rows of @p a
 * @param a_col number of columns of @p a / size of @p b
 * @return a two-dimensional matrix of size @p n_matrices by @p a_row
 */
uint64_t** MatrixVectorMultiply(
    Party *const proxy,
    const uint64_t *const *const *const a,
    const uint64_t *const *const b,
    uint32_t n_matrices,
    uint32_t a_row,
    uint32_t a_col,
    int shift = FRACTIONAL_BITS
) {
    int p_role = proxy->GetPRole();
    if (p_role == proxy1 || p_role == proxy2) {
        // reformat the given vectors into matrices, ensuring that one of their dimensions is 1
        uint64_t ***mat_b = new uint64_t**[n_matrices];
        for(uint32_t i = 0; i < n_matrices; i++) {
            mat_b[i] = new uint64_t *[a_col];
            for(uint32_t j = 0; j < a_col; j++) {
                mat_b[i][j] = new uint64_t[1];
                mat_b[i][j][0] = b[i][j];
            }
        }

        uint64_t ***result_array = MatrixMatrixMultiply(proxy, a, mat_b, n_matrices, a_row, a_col, 1);

        for(uint32_t i = 0; i < n_matrices; i++) {
            for(uint32_t j = 0; j < a_col; j++) {
                delete [] mat_b[i][j];
            }
            delete [] mat_b[i];
        }
        delete [] mat_b;

        uint64_t **result = new uint64_t *[n_matrices];
        for(uint64_t i = 0; i < n_matrices; i++) {
            result[i] = new uint64_t[a_row];
            for(uint64_t j = 0; j < a_row; j++) {
                result[i][j] = result_array[i][j][0];
            }
        }
        return result;
    }
    else if(p_role == helper) {
        MatrixMatrixMultiply(proxy, NULL, NULL, n_matrices, a_row, a_col, 1, shift);
        return NULL;
    }
    else {
        return nullptr;
    }
}

/** Perform multiplication of matrix a and vector b. The function assumes that the number of columns of a is equal to
 * the length of b.
 * @param a two dimensional matrix
 * @param b vector
 * @param a_row number of rows of @p a
 * @param a_col number of columns of @p a / size of @p b
 * @return a vector of size @p a_row
 */
uint64_t* MatrixVectorMultiply(
    Party *const proxy,
    const uint64_t *const *const a,
    const uint64_t *const b,
    uint32_t a_row,
    uint32_t a_col,
    int shift = FRACTIONAL_BITS
) {
    if (proxy->GetPRole() == helper) {
        MatrixVectorMultiply(proxy, nullptr, nullptr, 1, a_row, a_col, shift);
        return nullptr;
    } else {
        uint64_t** result_array = MatrixVectorMultiply(proxy, &a, &b, 1, a_row, a_col, shift);
        uint64_t* result = result_array[0];
        delete[] result_array;
        return result;
    }
}


uint64_t Sigmoid2(Party *const proxy, uint64_t z_share){ //having z as a double is more intuitive for me but the convention seems to be othervise.
    if (proxy->GetPRole() == helper) { //at the writing of this code, I have no direct use for the helper. This may be inherently flawed.
        return 0;
    }
    else{

        double lower_limit = -2;
        double upper_limit = 2;

        uint64_t low_share = proxy->CreateShare(lower_limit);
        uint64_t upper_share = proxy->CreateShare(upper_limit); //bunlar gerekmiyordu aslında...

        proxy->SendBytes(coreCompare); //haber vermek lazım
        uint64_t l_share = Compare(proxy, z_share, low_share); //l için gizli pay

        proxy->SendBytes(coreCompare); //haber vermek lazım
        uint64_t h_share = Compare(proxy, z_share, upper_share); //h için gizli pay

        //double l = ConvertToDouble(Reconstruct(proxy,l_share));
        //double h = ConvertToDouble(Reconstruct(proxy,h_share));

        proxy->SendBytes(coreMultiply); //bunları unutmak beni öldürecak
        uint64_t sigmoid_term1 = Multiply(proxy,(h_share - l_share), z_share); // (h-l)*z
        
        proxy->SendBytes(coreMultiply);
        uint64_t sigmoid_term2_5z = Multiply(proxy, z_share, ConvertToUint64(5)); //bu z=0.5 için nasıl 5 gelir lan?!
        cout << "5z: " << ConvertToDouble(Reconstruct(proxy,sigmoid_term2_5z)) << endl;
        proxy->SendBytes(coreDivide); //bölcem birazdan. local bölme yok sanırım
        uint64_t sigmoid_term2 = Divide(proxy,sigmoid_term2_5z, ConvertToUint64(4)); // 5x/4 için fazla uzun değil mi? BURADA DÜZGÜNÜZ
        cout << "bölümden sonra 5z/4: " << ConvertToDouble(Reconstruct(proxy, sigmoid_term2)) << endl;
        uint64_t addendum = 0;
        if (proxy->GetPRole() == proxy1) {
            addendum = ConvertToUint64(0.5);
        }
        sigmoid_term2 += addendum; //ikinci terim 5x/4+0.5, yukarıyı uzatmak istemedim.
        cout << "5z/4+0.5: " << ConvertToDouble(Reconstruct(proxy, sigmoid_term2)) << endl;


        proxy->SendBytes(coreMultiply);
        uint64_t one_minus_h = -h_share; //bu zıkkımın da kerreten yapılması gerekiyordu.
        if (proxy->GetPRole() == proxy1) {
            one_minus_h += ConvertToUint64(1);
        }
        cout << "1-h: " << ConvertToDouble(Reconstruct(proxy, one_minus_h)) << endl;
        uint64_t sigmoid_term3 = Multiply(proxy,(one_minus_h),l_share); //l(1-h)
        cout << "h: " << ConvertToDouble(Reconstruct(proxy,h_share)) << endl;
        cout << "l: " << ConvertToDouble(Reconstruct(proxy,l_share)) << endl;
        proxy->SendBytes(coreMultiply);
        uint64_t sigmoid_shared = sigmoid_term1 + h_share + Multiply(proxy,sigmoid_term2, sigmoid_term3);

        double sigmoid = ConvertToDouble(Reconstruct(proxy, sigmoid_shared));
        cout << "term1: " << ConvertToDouble(Reconstruct(proxy,sigmoid_term1)) << endl;
        cout << "term2: " << ConvertToDouble(Reconstruct(proxy,sigmoid_term2)) << endl;
        cout << "term3: " << ConvertToDouble(Reconstruct(proxy,sigmoid_term3)) << endl;
        cout <<"for "<<ConvertToDouble(Reconstruct(proxy,z_share))<< " sigmoid is: "<< sigmoid << endl;
        return sigmoid_shared;
        }
}


uint64_t* oldSigmoid(Party *const proxy, uint64_t *z_shares, int size){
    if(proxy->GetPRole() == helper){ //sırayla yapacağımız işlemlerin haberi verilmeli. 
        Compare(proxy, nullptr, nullptr, size);
        Compare(proxy, nullptr, nullptr, size);
        Multiply(proxy, nullptr, nullptr, size);
        Multiply(proxy, nullptr, nullptr, size);
        Divide(proxy, nullptr, nullptr, size);
        Multiply(proxy, nullptr, nullptr, size);
        Multiply(proxy, nullptr, nullptr, size);
        Multiply(proxy, nullptr, nullptr, size);
        return nullptr;
    }
    else{ 
        uint32_t param[1] = {(uint32_t)size}; //bunu "SendBytes" demeye devam ettiğim için kullanıyorum.

        double lower_limit[size];
        double higher_limit[size];
        double fives[size];
        double fours[size];
        double ones[size];
        for(int i=0; i<size; i++){
            lower_limit[i] = -2; //çünkü alt limit -2 :p
            higher_limit[i] = 2; //bu salaklığı yapmamın sebebi Compare fonksiyonunun iki array gerektirdiğini zannetmem.
            fives[i] = 5; //çünkü her elemanı 5'le çarpıcam.
            fours[i] = 4; //this is getting old....
            ones[i] = 1; // ben numpy'ı özledim...
        }

        uint64_t *lower_limit_shares = proxy->CreateShare(lower_limit, size);//-2 and 2 are quirks of the approx. alg. Kerreten yapılmaları gerek.
        uint64_t *higher_limit_shares = proxy->CreateShare(higher_limit, size);
        uint64_t *five_shares = proxy->CreateShare(fives, size); //bunları alt alta yazayım.

        uint64_t *l_shares = Compare(proxy, z_shares, lower_limit_shares, size); //Compare(x,y) gives x>y 
        uint64_t *h_shares = Compare(proxy, z_shares, higher_limit_shares, size);//h'ler için kıyas
        uint64_t *four_shares = proxy->CreateShare(fours, size);
        uint64_t *one_shares = proxy->CreateShare(ones, size);


        /*  Let's remember that we cannot have partially defined functions, so one must devine a unified formula.
        *   My solution is the monstosity seen in the funciton description.
        *   I have separated the terms of the formula for clarity and easier debugging.
        *   Now begins the calculations of the arrays of those separated terms.
        */

        uint64_t h_minus_l_shares[size];
        for(int i=0; i<size; i++){
            h_minus_l_shares[i] = h_shares[i] - l_shares[i]; //I need (h-l) and I can do that locally.
        }
        uint64_t *sigmoid_term1 = Multiply(proxy, h_minus_l_shares, z_shares, size);// terms in the from of (h-l)*z
        uint64_t *sigmoid_term2_5z = Multiply(proxy, z_shares, five_shares, size); //bunu adım adım buluyorum çünkü tek satırda çarpıp bölemem

        uint64_t *sigmoid_term2 = Divide(proxy, sigmoid_term2_5z, four_shares, size); //5z/4.   maksut 5z/4 + 0.5

        
        
        uint64_t addendum = ConvertToUint64(0.25); //ikisi de 0.25 eklesin ki 0.5 eklenmiş olsun
        double minus_ones[size];
        for(int i=0; i<size; i++){ //bunun daha zekice yöntemleri var ama önemli diil.
            sigmoid_term2[i] += addendum;
            minus_ones[i] = -1; // -h'leri elde etmek için biraz dolanbaç
        }

        uint64_t *minus_one_shares = proxy->CreateShare(minus_ones, size);
        uint64_t *minus_h_shares = Multiply(proxy, h_shares, minus_one_shares, size);

        
        
        uint64_t *one_minus_h_shares = Add(proxy, minus_h_shares, one_shares, size);
        

        uint64_t *sigmoid_term3 = Multiply(proxy, l_shares, one_minus_h_shares, size); // l(1-h) veriyor, inş...

        uint64_t *sigmoid_term2_times_3 = Multiply(proxy, sigmoid_term2, sigmoid_term3, size); //çarpmam lazım bunları

        
        uint64_t *half_step = Add(proxy, sigmoid_term1, h_shares, size); //(h-l)z+h yani. Bunun hesabını haber vermem gerekmiyor, hehe
        uint64_t *final_sigmoid = Add(proxy, half_step, sigmoid_term2_times_3, size);



        return final_sigmoid;
    }
}

uint64_t oldSigmoid(Party *const proxy, uint64_t z_share){ //tek eleman için çalışan hâli.
    uint64_t z_share_array[1] = {z_share};
    uint64_t *result_of_vectorised = oldSigmoid(proxy, z_share_array, 1);
    uint64_t result_share = result_of_vectorised[0];
    return result_share;
}


/** An approximation of the sigmoid function. UNDER CONSTRUCTION
 *  Operates on singular values.
 *  The approximation function is 0 for z<-2, 1 for z>2, and x/4+0.5 for the values in betwee.
 * @param proxy
 * @param z_shares array of z values 
 * @param size size of the z_shares array. Bunu niçin z_shares'in kendisinden çıkarmadığımızdan emin değilim, sanırım daha temiz bir tatbik oluyor.
 * @return (h-l)z + h + (5z/4 + 0.5) * l * (1-h). see the comments for detail.
 */

// uint64_t Sigmoid(Party *const proxy, uint64_t* z_shares, int size){ 
    
//     if (proxy->GetPRole() == helper) { //at the writing of this code, I have no direct use for the helper. This may be inherently flawed.
        
//         Compare(proxy, 0, 0); //first compare to create l_share
//         Compare(proxy, 0, 0); //to crate h_share

//         return 0;
//     }

//     else{ //if Role == proxy

//         uint64_t low_share = proxy->CreateShare(-2.0);
//         uint64_t upper_share = proxy->CreateShare(2.0); //-2 and 2 are quirks of the approx. alg. Kerreten yapılmaları gerek.

//         proxy->SendBytes(coreCompare);
//         uint64_t l_share = Compare(proxy, z_shares, low_share, size); //this serves as a boolean, which I use in the calculation of the func.

//         proxy->SendBytes(coreCompare); //I guess we have to notify the helper each time.
//         uint64_t h_share = Compare(proxy, z_shares, upper_share); //see two lines above.

//         /*  Let's remember that we cannot have partially defined functions, so one must devine a unified formula.
//         *   My solution is the monstosity seen in the funciton description.
//         *   I have  separated the terms of the formula for clarity and easier debugging.
//         */
//         proxy->SendBytes(coreMultiply); //bunları unutmak beni öldürecak
//         uint64_t sigmoid_term1 = Multiply(proxy,(h_share - l_share), z_shares); // (h-l)*z

//         proxy->SendBytes(coreMultiply);
//         uint64_t sigmoid_term2_5z = Multiply(proxy, z_share, ConvertToUint64(5));
//         proxy->SendBytes(coreDivide); //I don't think we have local division even for scalers.
//         uint64_t sigmoid_term2 = Divide(proxy, sigmoid_term2_5z, ConvertToUint64(4)); 

//         uint64_t addendum = 0; //0.25 de olabilür
//         if (proxy->GetPRole() == proxy1) {
//             addendum = ConvertToUint64(0.5);
//         }
//         sigmoid_term2 += addendum; //deuxième terme est 5z/4+0.5 

//         proxy->SendBytes(coreMultiply);
//         uint64_t one_minus_h = -h_share; //bu zıkkımın da kerreten yapılması gerekiyordu.
//         if (proxy->GetPRole() == proxy1) {
//             one_minus_h += ConvertToUint64(1);
//         }
//         uint64_t sigmoid_term3 = Multiply(proxy,one_minus_h,l_share); //l(1-h)

//         proxy->SendBytes(coreMultiply); //We arrive at the final form:
//         uint64_t sigmoid_shared = sigmoid_term1 + h_share + Multiply(proxy,sigmoid_term2, sigmoid_term3); //(h-l)z + h + (5z/4 + 0.5) * l * (1-h)

//         return sigmoid_shared; 
//         }
// }
double generateLaplaceNoise(double scale) {
    cout <<"scale is: "<< scale << endl;
    //if helper falan yazmadım başına çünkü sadece helper kullanacak 
   
    // onların generateRandom common fonksiyonunu kullanmadım çünkü bize double lazım!!
    // bu kısıma birlikte bakalım!!!
    
    std::random_device rd; //bu seed oluşturmak için
    std::mt19937 gen(rd()); //bu o seed'i kullanarak random bir sayı olışturuyor
    

    std::uniform_real_distribution<> dis(0.0, 1.0);//!!!!!!bu sıkıntı aslında  çünkü uint64 hali ile yapıyoruz işlemi
    
   
    double u = dis(gen) - 0.5;//bu da o random sayıyı kullanıyor
    //auto n= proxy->GenerateCommonRandom();//onlarınkini kullandım burada random sayı oluşturmak için
    //ama ufak bir sorun var bu double değil de uint64 döndürüyor biz bunu nerede çağıracaksak ona göre değiştirmeliyiz
    
    //uint64_t u = dis(n) - 0.5;
    // Apply the inverse CDF of the Laplace distribution
    double noise = scale * ((u >= 0) ? -std::log(1 - 2*u) : std::log(1 + 2*u));
    cout << "noise: "<<noise << endl;
    return noise;
}

double* generateLaplaceNoise(double scale, int size) {
    std::random_device rd; //bu seed oluşturmak için
    std::mt19937 gen(rd()); //bu o seed'i kullanarak random bir sayı olışturuyor
    std::uniform_real_distribution<> dis(0.0, 1.0);//!!!!!!bu sıkıntı aslında  çünkü uint64 hali ile yapıyoruz işlemi
    
    double *noises = new double[size];
    for(int i=0; i<size; ++i){
        double u = dis(gen) - 0.5;//bu da o random sayıyı kullanıyor
        double noise = scale * ((u >= 0) ? -std::log(1 - 2*u) : std::log(1 + 2*u));
        noises[i] = noise;
    }
    return noises;
}

uint64_t *add_noise(Party *const proxy, uint64_t *numbers, uint32_t sz, uint32_t scale_32bit){
    const int length = int(sz);
    double scale = ConvertToDouble(uint64_t(scale_32bit));

    if (proxy->GetPRole() == helper){ //helper gürültünün imalinden mesul
        
        //bu kısımda gürültü imal edilmeli ve onun gizli payları oluşturulmalı.
        double *noise = new double[length];

        noise = generateLaplaceNoise(scale, length);
        
        uint64_t *noise_shares = ConvertToUint64(noise, length);

        uint64_t proxy1_share;
        uint64_t proxy2_share;
        
        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2(); //HELPER'IN buffer'lara pointer'lar
        //cout << "pointerlar falan tamam" << endl;
        for(int i=0; i<length; i++){ //bu optimal olmayabilir.
            //cout << "for'a girdik" << endl;
            uint64_t tempShare = proxy->GenerateRandom();  //elimizde GenerateRandom da var, ikisi de aynı çalışıyor gibi.
            //cout << "rastgele üretim" << endl;
            proxy1_share = tempShare; //ilki rastgele bir sayı olacak.
            proxy2_share = noise_shares[i] - tempShare; //ikincisi de o rastgele ile bizim göndermek istediğimizin farkı.

            AddValueToCharArray(proxy1_share, &ptr1); //Helper proxy1 için olan verisini kendi 1. buffer'ına ve
            AddValueToCharArray(proxy2_share, &ptr2); //proxy2'nin verisini kendi 2. buffer'ına yüklüyor
        }        
        //cout << "for looptan çıktık" << endl;
        thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), length * 8);
        thread thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), length * 8);        

        thr1.join();
        thr2.join();

        delete[] noise;
        return nullptr;
    }
    else{
        //cout << "add_noise proxy başla" << endl;
        uint64_t *noise_shares = new uint64_t[length];
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), 8*length); 
        //cout << "add_noise proxy receive tamam" << endl;

        unsigned char *ptr = proxy->GetBuffer1(); //KENDİ birinci buffer'ı (birincinin buffer'ı değil!!)     

        for(int i=0; i<length; ++i){
            noise_shares[i] = ConvertToLong(&ptr);
        }
        //cout << "add_noise proxy çevirmeler tamam" << endl;
        uint64_t *noiseful_shares = Add(proxy, numbers, noise_shares, length);


        delete[] noise_shares;
        return noiseful_shares;
    }
}


uint64_t *add_noise_old(Party *const proxy, uint64_t *numbers, uint64_t size){
    if (proxy->GetPRole() == helper){
        double noise[size]; //bunun içini noise ile dolduracağız, şimdilik 100 falan koyuyorum.

        for(int i=0; i<size; i++){
            noise[i] = 100;
        }
        uint64_t *noise_shares = proxy->CreateShare(noise, size); //bunları proxy'lere göndereceğiz.

        unsigned char *ptr_out2 = proxy->GetBuffer2(); //tamamen çalıyorum!!
        for(int j=0; j<size; j++){
            AddValueToCharArray(noise_shares[j], &ptr_out2); //sanırım burada noise_shares'in içindekileri ptr'ın gösterdiği yere koyuyoruz
        }
        Send(proxy->GetSocketP2(), proxy->GetBuffer2(), size*8); //neden *8 bilmiyorum ama 2'nin buffer'a yolluyoruz.
        //ben "new" demediğim için delete lazım değil :p
        return 0;
    }
    else {
        uint64_t noises[size];

        if(proxy->GetPRole() == proxy2){ //proxy->GetPRole() == proxy2
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), size*8); //birbirlerinkini kullandıklarını varsaymak makul (?)
            unsigned char *ptr = proxy->GetBuffer1();

            for (int i = 0; i < size; ++i) {
                noises[i] = ConvertToLong(&ptr);
            }
            
        }
        else{
            cout << "proxy 1 babayı alır" <<endl;
        }
        return noises; //şimdilik direkt döndürsün gürültüyü.
    }
}


uint64_t* fonk(Party *const proxy,const uint64_t const x, uint32_t size ){
    if(proxy->GetPRole() == helper){
        uint64_t numbers1[size];
        for(int i=0; i<size; i++){
            numbers1[i] = uint64_t(1); //proxy1'e 1, 2'ye de 2 göndermek istiyorum sadece.
        }
        uint64_t numbers2[size];
        for(int i=0; i<size; i++){
            numbers2[i] = uint64_t(2); //bu proxy2'ye gidecek.
        }

        unsigned char *ptr_out1 = proxy->GetBuffer1();
        unsigned char *ptr_out2 = proxy->GetBuffer2(); //buffer'lara ulaşmak için iki pointer.


        uint64_t mittenda[size]; //0,1,2,3,4
        for(int i=0; i<size; i++){
            uint64_t mittendum = ConvertToUint64(i);
            mittenda[i] = mittendum; //bunları gönderip geri okuyacağım.
        }

        for(int i=0; i<size; i++){
            AddValueToCharArray(mittenda[i], &ptr_out1); //buffer1'e 1ler yerleştiriyorum (?) ŞU AN 10,20,...,50 YOLLUYORUM.
        }
        for(int i=0; i<size; i++){
            AddValueToCharArray(mittenda[i], &ptr_out2); //buffer2'ye 2ler yerleştiriyorum (?) ŞU AN İKİSİNE DE AYNI ŞEYİ GÖNDERMEYİ DENİYORUM.
        }
        Send(proxy->GetSocketP1(),proxy->GetBuffer1(), size *8); 
        Send(proxy->GetSocketP2(),proxy->GetBuffer2(), size *8); //niçin *8 olduğunu anlamıyorum şu an, patlarsak buradan patlarız.

        for(int i=0; i<size; i++){
            cout << mittenda[i] << endl;
        }
        return 0;
    }
    else{ //proxy 1 ve 2 ise
        //uint64_t numbers_received[size]; //helper'dan okuyacakları sayıları barındıran array.
        uint64_t *numbers_received = new uint64_t[size];
        uint32_t param[1] = {(uint32_t)size}; 
        if(proxy->GetPRole() == proxy1){
            cout << "proxy1 here!" << endl;
            //evvela proxy1 girdilerini alsın.
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), size*8); //bu yanlış olabilir..
            unsigned char *ptr = proxy->GetBuffer1();
            
            for (int i = 0; i < size; ++i) {
                numbers_received[i] = ConvertToLong(&ptr); //proxy1 kendince okuyor.
                cout << "bir tane uzattım, o da şu: " << numbers_received[i] << endl;
            } 
        }
        else if(proxy->GetPRole() == proxy2){ //I find this easier to read
            cout << "proxy2 here!" << endl;
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer2(), size*8);
            unsigned char *ptr = proxy->GetBuffer2();

            for (int i = 0; i < size; ++i) {
                numbers_received[i] = ConvertToLong(&ptr); //proxy2 kendince okuyor
                cout << "bir tane uzattım, o da şu: " << numbers_received[i] << endl;
            }
        }   //şimdi okumalar bitti, işlemlere geçiyoruz:
            //yapacağım işlem okuduğum sayıları elde mevcut bir başka array ile çarpmaktan ibaret olacak.

            double tens[size];
            for(int i=0; i<size; i++){
                tens[i] = 10;
            }
            uint64_t *ten_shares =  proxy->CreateShare(tens, size);
            proxy->SendBytes(coreVectorisedMultiply, param, 1);
            uint64_t *result_shares = Multiply(proxy, numbers_received, ten_shares, size);



          return numbers_received;
    }
}


// uint64_t* fonk(Party *const proxy,const uint64_t const x, uint32_t size ){
//     if (proxy->GetPRole() == helper) {
//         uint64_t *c1 = new uint64_t[size];
        
//         for(int i=0; i<size; i++){ //öyle birden sınıra kadar sayılar atıyorum şu an sadece
//             c1[i] = uint64_t(i);
//         }

//         cout << "jenerasyon tamam" << endl;

//         unsigned char *ptr_out2 = proxy->GetBuffer2();
//         for (int j = 0; j < size; j++) {
//             AddValueToCharArray(c1[j], &ptr_out2);
//         }

//         Send( proxy->GetSocketP2(), proxy->GetBuffer2(), size * 8);

//         delete[] c1;
//         cout << "helper tamam" << endl;
//         if (DEBUG_FLAG >= 1)
//             cout << "Returning from PMNF_MUL...\n************************************************************" << endl;
//         return 0;
//     }
//     else{
//         if(proxy->GetPRole() == proxy2){
//             cout << "proxy2 gördü" << endl;

//             Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), size * 8);
//             cout << "recieve ok" << endl;
//             unsigned char *ptr = proxy->GetBuffer1();
//             cout << "buffer ok" << endl;
//             uint64_t numbers[size];

//             for (int i = 0; i < size; ++i) {
//                 numbers[i] = ConvertToLong(&ptr);
//                 cout << "bir tane uzattım, o da şu: " << numbers[i] << endl;
//             }
//             return numbers;
//         }
//         else{
//             cout << "proxy1 gördü" << endl;
//             return nullptr;
//         }

//     }
// }



uint64_t *Min(Party *const proxy, const uint64_t *const x, const uint64_t *const y, uint32_t size) { 
    if(proxy->GetPRole() == helper){
        MostSignificantBit(proxy, nullptr, size);
        return nullptr;
    }
    else{
        uint64_t *diffs = new uint64_t[size];
        for (int i = 0; i < size; i++) {
            diffs[i] = x[i] - y[i];  //if (x-y) < 0 --> (x<y)
        }
        uint64_t* msb_shares = MostSignificantBit(proxy, diffs, size);
        uint64_t *msbs = Reconstruct(proxy, msb_shares, size); //0'sa zaten 0
        
        uint64_t *min_shares = new uint64_t[size];

        for (int i = 0; i < size; i++) {
            if(msbs[i]){ //0 değilse eksi sayı imiş, yani ki x<y imiş yani ki x küçük olanmış
                min_shares[i] = x[i];
            }
            else{
                min_shares[i] = y[i];
            }
        }
        return min_shares;
    }
}

uint64_t *Max(Party *const proxy, const uint64_t *const x, const uint64_t *const y, uint32_t size) { 
    if(proxy->GetPRole() == helper){
        MostSignificantBit(proxy, nullptr, size);
        return nullptr;
    }
    else{
        uint64_t *diffs = new uint64_t[size];
        for (int i = 0; i < size; i++) {
            diffs[i] = x[i] - y[i];  //if (x-y) < 0 --> (x<y)
        }
        uint64_t* msb_shares = MostSignificantBit(proxy, diffs, size);
        uint64_t *msbs = Reconstruct(proxy, msb_shares, size); //0'sa zaten 0
        
        uint64_t *max_shares = new uint64_t[size];

        for (int i = 0; i < size; i++) {
            if(msbs[i]){ //0 değilse eksi sayı imiş, yani ki x<y imiş yani ki x küçük olanmış
                max_shares[i] = y[i];
            }
            else{
                max_shares[i] = x[i];
            }
        }
        return max_shares;
    }
}

uint64_t* awesome_sig(Party *const proxy, uint64_t* z_shares, uint32_t size){
    if(proxy->GetPRole() == helper){
        cout << "selamlars" << endl;
        Divide(proxy, nullptr, nullptr, size);
        cout << "böldük" << endl;
        Min(proxy, nullptr, nullptr, size);
        Max(proxy, nullptr, nullptr, size);
        cout << "iş benden çıktı usta" << endl;
        return nullptr;
    }
    else{
        cout << "sig start" << endl;
        double fours[size];
        double halves[size];
        double ones[size];
        double zeros[size];
        for(int i=0; i<size; i++){
            fours[i] = 4;
            halves[i] = 0.5;
            ones[i] = 1;
            zeros[i] = 0;
        }
        uint64_t *four_shares = proxy->CreateShare(fours, size);
        uint64_t *half_shares = proxy->CreateShare(halves, size);
        uint64_t *one_shares = proxy->CreateShare(ones, size);
        uint64_t *zero_shares = proxy->CreateShare(zeros, size);
        
        uint64_t *z_div_4 = Divide(proxy, z_shares, four_shares, size);
        uint64_t *z_div_4_plus_half = Add(proxy, z_div_4, half_shares, size); //I can do this locally. BU 1 TAM EKLİYOR OLABİLİR ŞU AN

        uint64_t *min_term_shares = Min(proxy, z_div_4_plus_half, one_shares, size);
        uint64_t *max_term_shares = Max(proxy, min_term_shares, zero_shares, size);
        return max_term_shares; //which is the sigmoid 
    }
}






 uint64_t *add_noise_new(Party *const proxy, uint64_t *numbers, double scale, int size){
    if (proxy->GetPRole() == helper){ //helper gürültünün imalinden mesul
        cout << "noise helper" << endl;
        //bu kısımda gürültü imal edilmeli ve onun gizli payları oluşturulmalı.
        double *noise = new double[size];
        //zaten helper'ın içinde kullandığımız için helper'a haber vermeye gerek yok sanki
        //noise boyutu değişebilir iki boyutlu array de gerekebilir ona bir daha bakalım

        //----------------------//
        //Burada gerekli fonkisyonu kullanarak laplace dağılımına göre elemanlar ürettiğimiz zaman gürültü ekleme işini halletmiş oluyoruz!
        for(int i=0; i<size; i++){    
            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            noise[i] = generateLaplaceNoise(scale);
            cout << noise[i] << endl;
        }
        //Ben şimdilik öyle kafadan 100 ekliyorum :p
        //----------------------//
        uint64_t *noise_shares = ConvertToUint64(noise, size);

        uint64_t proxy1_share;
        uint64_t proxy2_share;

        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2(); //HELPER'IN buffer'lara pointer'lar

        for(int i=0; i<size; i++){ //bu optimal olmayabilir.
            uint64_t tempShare = proxy->GenerateCommonRandom();  //elimizde GenerateRandom da var, ikisi de aynı çalışıyor gibi.
            proxy1_share = tempShare; //ilki rastgele bir sayı olacak.
            proxy2_share = noise_shares[i] - tempShare; //ikincisi de o rastgele ile bizim göndermek istediğimizin farkı.

            AddValueToCharArray(proxy1_share, &ptr1); //Helper proxy1 için olan verisini kendi 1. buffer'ına ve
            AddValueToCharArray(proxy2_share, &ptr2); //proxy2'nin verisini kendi 2. buffer'ına yüklüyor
        }        
        
        Send(proxy->GetSocketP1(), proxy->GetBuffer1(), 8*size); //asıl gönderme işlemi.
        Send(proxy->GetSocketP2(), proxy->GetBuffer2(), 8*size); //thread kullanmayı sevmiyorum. o yüzden kullanmadım :p

        delete[] noise;
        return nullptr;
    }
    else{
        cout << "noise proxy" << endl;
        uint64_t *noise_shares = new uint64_t[size];
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), 8*size); 

        unsigned char *ptr = proxy->GetBuffer1(); //KENDİ birinci buffer'ı (birincinin buffer'ı değil!!)     

        for(int i=0; i<size; ++i){
            noise_shares[i] = ConvertToLong(&ptr);
        }
        uint64_t *noiseful_shares = Add(proxy, numbers, noise_shares, size);

        //delete[] noise_shares;
        //return noiseful_shares;
        return noise_shares;
    }
}






/**
 * @param X: girdiler işte, bunun tipinden emin değilim henüz. iki boyutlu bir array olması laızm
 * @param r: row sayısı
 * @param n: feature sayısı
 * @param epsilon: bunu double koyabiliriz yav
 */
uint64_t* GradientDescent_old(Party *const proxy, uint64_t **X, uint64_t *y, uint64_t* theta_shares, int const r, int const n, double const learning_rate, double epsilon, int iterations){ //sanırım burada vektör kullanmak zorunda kalacağım.
    uint32_t param_r[1] = {(uint32_t)r}; //BU DAHA SONRA SİLİNECEK.
    uint32_t param_n[1] = {(uint32_t)n}; //BU DAHA SONRA SİLİNECEK.
    
    //bazı 2d array'lerin transpozunu almak gerekecek.
    auto T = [](uint64_t **matrix, int r, int n){
        uint64_t** transpose = new uint64_t*[n];
        for(int i=0; i<n; i++){
            transpose[i] = new uint64_t[r];} //memory allocation baby!!
        for(int i=0; i<r; ++i){
            for(int j=0; j<n; ++j){
                transpose[j][i] = matrix[i][j];}}
        return transpose;
    };

    if(proxy->GetPRole() == helper){
        //MatrixMatrixMultiply(proxy, nullptr, nullptr, 0, 0, 0);
        return nullptr;
    }
    else{
        double minus_ones[r];     // şekli doğru
        double halves[r];         // şekli doğru
        double eights[n];         // r idi, n'ye tebdil edili.
        double twos[n];           // n'ye tebdil
        double rs[n];             // n'ye tebdil
        double learning_rates[n]; // n'ye tebdil
        double minus_ones_n[n];
        for(int i=0; i<r; i++){
            minus_ones[i] = -1;
            halves[i] = 0.5;
            //eights[i] = 8; DÜZELTİLDİ.
            //twos[i] = 2;
            //rs[i] = r;
            //learning_rates[i] = learning_rate;
        }
        for(int i=0; i<n; i++){
            eights[i] = 8;
            twos[i] = 2;
            rs[i] = r;
            learning_rates[i] = learning_rate;
            minus_ones_n[i] = -1;
        }
        uint64_t *minus_one_shares = proxy->CreateShare(minus_ones, r);
        uint64_t *half_shares = proxy->CreateShare(halves, r);
        uint64_t *eight_shares_n = proxy->CreateShare(eights, n);                 //R'DEN N'YE TEBDİL EDİLDİ.
        uint64_t *two_shares_n = proxy->CreateShare(twos, n);                     //R'DEN N'YE TEBDİL EDİLDİ.
        uint64_t *r_shares_n = proxy->CreateShare(rs, n);                         //R'DEN N'YE TEBDİL EDİLDİ.
        uint64_t *learning_rate_shares_n = proxy->CreateShare(learning_rates, n); //R'DEN N'YE TEBDİL EDİLDİ.
        uint64_t *minus_one_shares_n = proxy->CreateShare(minus_ones, n);         //N olarak eklendi, yoğudu

        double d = n-1;
        double sensitivity = (3 * d) + 0.25 * (n * n); // calculates the sensitivity with the number of features
        uint64_t sens_share = proxy->CreateShare(sensitivity); //müstakil bir sayı ve defaatle kullanılacak

        double scale = sensitivity / epsilon * iterations;  // calculate beta, lokal hesaplanabilir bu, gizli bir yanı yok ki.

        uint32_t noise_param[2] = {uint32_t(scale), uint32_t(n)};

        double *noise1 = new double[r]; //BUNLARIN ŞEKLİNİN NE OLMSAI GEREKTİĞİNİ ŞU AN KESTİREMİYORUM, İMPLEMANTASYONA BAĞLI.
        double *noise2 = new double[r]; //bunlar şimdilik boş kalacaklar, silmeyi unutmamak lazım.

        uint64_t *curiousity = new uint64_t[n];
        for(int i=0; i<iterations; i++){
            cout << "for loop'a girdik" << endl;
            //term1 = np.matmul(X.T, (0.5 - y)) + noise1
            uint64_t **X_transpose = T(X, r, n); //X.T tamam
            proxy->SendBytes(coreVectorisedMultiply, param_r, 1);
            uint64_t *minus_y_shares = Multiply(proxy, y, minus_one_shares, r);
            uint64_t *half_minus_y_shares = Add(proxy, minus_y_shares, half_shares, r);
            cout << "0.5 -y işi tamam" << endl;
            
            uint32_t param_vec_mult[2] = {uint32_t(n),uint32_t(r)};
            proxy->SendBytes(coreMatrixVectorMultiply, param_vec_mult, 2);
            uint64_t *term1_shares = MatrixVectorMultiply(proxy,X_transpose, half_minus_y_shares, n, r);
            

            //term1'e gürültü eklemeyi deniyoruz.
            proxy->SendBytes(lgAddNoiseNew, noise_param, 2);
            uint64_t *term1_shares_noisefull = add_noise_new(proxy, term1_shares, scale, n);
            cout << "scale is: "<< scale << endl;


            uint64_t *term1_for_debug = Reconstruct(proxy, term1_shares, n);
            cout << "term1------------------------------" << endl;
            cout << ConvertToDouble(term1_for_debug[0]) << endl;
            cout << "term1------------------------------" << endl;

            //uint64_t *term1_shares = multiplyMatrixVector(X_transpose, half_minus_y_shares, n, r); //r ve n ters çünkü X transpozun şekli bi' garip
            //term1 -> (n,1) şeklinde
            //term2 = np.matmul(X.T, X) / 8 + noise2 
            cout << "Xt * (0.5-y) kısmı tamam" << endl;
            //uint64_t **X_transpose_X = multiplyMatrices(T(X, r, n), X, n, r, n); //bu n ve r'ler İNANILMAZ KARIŞTI.

            uint32_t param_mult[3] = {uint32_t(n),uint32_t(r),uint32_t(n)};
            proxy->SendBytes(coreMatrixMatrixMultiply, param_mult, 3);
            uint64_t **X_transpose_X = MatrixMatrixMultiply(proxy, T(X,r,n), X, n, r, n);

            cout << "XtX çarpımı tamam" << endl;
            // Xt * X --> (n,n) şeklinde
            //PS. X_transpose'u yeniden kullanmaya çalıştığımda yanlış sonuç alıyordum, o yüzden burada yeniden hesapladım
            //Dilin nasıl çalıştığını daha iyi takip eden biri bunu optimize etmekte zorlanmayacaktır.

            //X.T*X'in her elemanını 8'e bölmek lazım:
            //satır satır vektör muamelesi yapmak kabil.
            //uint64_t **term2_shares = new uint64_t*[r]; //BU N UZUNLUĞUNDA OLMALI!!!
            uint64_t **term2_shares = new uint64_t*[n]; //BU N UZUNLUĞUNDA OLMALI!!!
            for(int j=0; j<n; j++){ //X.T*X'in her satırı için R'YE KADARDI N'YE ÇEVİRDİM.
                proxy->SendBytes(coreVectorisedDivide, param_n, 1);
                term2_shares[j] = Divide(proxy, X_transpose_X[j], eight_shares_n, n); //o satırı 8'e böl N'YE DEĞİŞTİRİLDİ
            }
            uint64_t *theta_reconst = Reconstruct(proxy, theta_shares, n);
            cout << "THEATA------------------------------THEATA" << endl;
            cout << ConvertToDouble(theta_reconst[0]) << endl;
            cout << ConvertToDouble(theta_reconst[1]) << endl;
            cout << ConvertToDouble(theta_reconst[2]) << endl;
            cout << ConvertToDouble(theta_reconst[3]) << endl;
            cout << ConvertToDouble(theta_reconst[4]) << endl;
            cout << "THEATA------------------------------THEATA" << endl;

            cout << "XtX'lerin bölünmesi zıkkımı tamam" << endl;
            //term2_shares = Add(term2_shares, noise2_shares); DAHA GÜRÜLTÜLERİ GETİRMEDİĞİMİZ İÇİN BIRAKIYORUM BÖYLE.
            

            uint32_t param_vec_mult_for_term2_theta[2] = {uint32_t(n),uint32_t(n)};
            //gradient = term1 + 2 * np.matmul(term2, theta) ki burada theta n elemanlı bir vektör
            proxy->SendBytes(coreMatrixVectorMultiply, param_vec_mult_for_term2_theta, 2);
            uint64_t *term2_times_theta = MatrixVectorMultiply(proxy, term2_shares, theta_shares, n, n);
            //uint64_t *term2_times_theta = multiplyMatrixVector(term2_shares, theta_shares, n, n); //benim yazdığım lambda fonksiyonu !!!!!!!!!!

            uint64_t *term2_with_theta_for_debug = Reconstruct(proxy, term2_times_theta, n);
            cout << "------------------------------term2*theta" << endl;
            cout << ConvertToDouble(term2_with_theta_for_debug[0]) << endl;
            cout << ConvertToDouble(term2_with_theta_for_debug[1]) << endl;
            cout << ConvertToDouble(term2_with_theta_for_debug[2]) << endl;
            cout << ConvertToDouble(term2_with_theta_for_debug[3]) << endl;
            cout << "------------------------------term2*theta" << endl;


            proxy->SendBytes(coreVectorisedMultiply, param_n, 1);
            uint64_t *double_term2_times_theta = Multiply(proxy, term2_times_theta, two_shares_n, n); //2 * np.matmul(term2, theta) N'YE TEBDİL
            uint64_t *gradient_shares = Add(proxy, term1_shares, double_term2_times_theta, n); //N'YE TEBDİL
            cout << "vektörize çarpımlar !!" << endl;
            //theta -= learning_rate * (gradient / m)
            proxy->SendBytes(coreVectorisedDivide, param_n, 1);
            uint64_t *gradient_div_r_shares = Divide(proxy,gradient_shares, r_shares_n, n); // gradient / r N'YE TEBDİL
            cout << "gradient/n !!!!" << endl;
            proxy->SendBytes(coreVectorisedMultiply, param_n, 1);
            uint64_t *learn_rate_times_grad_div_r_shares = Multiply(proxy, learning_rate_shares_n, gradient_div_r_shares, n); //N'YE TEBDİL
            cout << "final" << endl;
            proxy->SendBytes(coreVectorisedMultiply, param_n, 1); //ya bu çok saçma ama çok güvenilir bir yöntem çıkarma yapmak için :/
            uint64_t *substract_from_theta = Multiply(proxy, learn_rate_times_grad_div_r_shares, minus_one_shares_n, n); //N'YE TEBDİL
            theta_shares = Add(proxy, theta_shares, substract_from_theta, n); //BU N'LER İLE İLGİLİ BİR HATA VAR!! KESİNLİKLE VAR!!
            
            
            /**
             * bildiğim fakat henüz düzeltmediğim bir hata var:
             * bazı vektörize işlemler satır sayısına göre iken bazıları da (özellikle theta'yı ilgilendirenler) mecburen sütun sayısı
             * ile ilgili.
             * ben hep "size" kullanmaya alışkın olduğum için bu yukarıdaki son satır hariç hep r dedim
             * nerelerde n olması gerektiği çok önemli bir konu fakat kare matrislerde denerken sorun çıkmıyor :P
             */
            
            curiousity = theta_shares; //bunu istediğim çıktığı döndürmek istediğim zaman, debug için kullandım.
        }

        delete[] noise1;
        delete[] noise2;
        
        return theta_shares;
    }
}




uint64_t* Sigmoid(Party *const proxy, uint64_t* z,uint32_t size,int shift = FRACTIONAL_BITS){
     
    if(proxy ->GetPRole()==helper){
        cout << "sig'e girdik" << endl;
        Compare(proxy,nullptr,nullptr,size,shift);
        Compare(proxy,nullptr,nullptr,size,shift);
        Multiply(proxy,nullptr,nullptr,size,shift);
        Multiply(proxy,nullptr,nullptr,size,shift);
        Multiply(proxy,nullptr,nullptr,size,shift);
        cout << "çarplar bitti" << endl;
        return nullptr;
       
    }
    else{
        cout << "proxy sig'e girdi" << endl;
        uint64_t l=proxy->CreateShare(-2.0);
        uint64_t u=proxy->CreateShare(2.0);
        uint64_t *lv = new uint64_t[size];
        uint64_t *uv = new uint64_t[size];
        uint64_t *ov = new uint64_t[size];
        uint64_t *fv = new uint64_t[size];
        uint64_t *hv = new uint64_t[size];

        uint64_t f=proxy->CreateShare(4.0);
        uint64_t o=proxy->CreateShare(1.0);
        uint64_t h=proxy->CreateShare(0.5);
        
        double *quarters = new double[size];

        for(int i=0;i<size;i++){
         lv[i]=l;
         uv[i]=u;
         ov[i]=o;
         fv[i]=f;
         hv[i]=h;
         
            quarters[i] = 0.25;
        }
        uint64_t *quarter_shares = proxy->CreateShare(quarters, size);
          
        uint64_t* l1=Compare(proxy,lv,z,size);
        uint64_t* u1=Compare(proxy,uv,z,size);

        uint64_t *term1= new uint64_t[size];
        for(int i=0;i<size;i++){
            term1 [i]= (o-l1[i]);
        }
       
        //uint64_t* z_div_4=Divide(proxy,z,fv,size,shift);
        uint64_t *z_div_4 = Multiply(proxy, z, quarter_shares, size);

        uint64_t *term2= new uint64_t[size];
         for(int i=0;i<size;i++){
            term2[i]=(z_div_4[i]- hv[i]);
        }
        
        uint64_t* term3=Multiply(proxy,term2,u1,size);
        
          for(int i=0;i<size;i++){
            term3[i]+=ov[i];
          }
      
        uint64_t* final= Multiply(proxy,term1,term3,size);
        
        return final;



    }
}

// X : feature matrix
  // sonuçları vektör olarak döndürüyor
 uint64_t* Predict(Party *const proxy, uint64_t** X, uint64_t* theta, uint32_t row, uint32_t column,int shift = FRACTIONAL_BITS){
    if(proxy ->GetPRole()==helper){

        for(int i=0;i<row;i++){
        DotProduct(proxy,nullptr,nullptr,column,shift);
       }
        Sigmoid(proxy,nullptr,row);
        Compare(proxy,nullptr,nullptr,row,shift);
       return 0;
       
    }else{

        double *h= new double[row];
        uint64_t *result = new uint64_t[row];
        for(int i=0;i<row;i++){
            uint64_t r= DotProduct(proxy,X[i],theta,column,shift);
            result[i]=r;
            h[i]=0.5;

        }
        uint64_t* hs=proxy-> CreateShare(h,row);

        return Compare(proxy,Sigmoid(proxy,result,row),hs,row);// büyük eşittir yoktu belki 1-(s<0.5) şeklinde yazabiliriz
    }

 }


//**************************************************


uint64_t **add_noise(Party *const proxy, uint64_t **numbers, uint32_t size1_32bit, uint32_t size2_32bit, uint32_t scale_32bit){
    double scale = ConvertToDouble(uint64_t(scale_32bit));
    const int size1 = int(size1_32bit);
    const int size2 = int(size2_32bit);

    if (proxy->GetPRole() == helper){ 
        double** noise = new double*[size1];
        for (int i = 0; i < size1; ++i) {
            noise[i] = generateLaplaceNoise(scale, size2); //satır satır doldurmuş olmalı
            //noise[i] = new double[size2];
        }
        cout << ",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,," << endl;
        cout << noise[0][0] << endl;
        cout << noise[0][1] << endl;
        cout << noise[0][2] << endl;
        cout << ",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,," << endl;
        // for(int i=0; i<size1; i++){   
        //     for (int j = 0; j < size2; j++)
        //     {
        //         noise[i][j] = generateLaplaceNoise(scale);//scale parametre olarak alınmalı
        //     }
        // }
            cout << "dönüştürmeden önce" << endl;
            uint64_t **noise_shares = ConvertToUint64(noise, size1,size2);
            cout << "dönüştürmeden sonra" << endl;

            uint64_t proxy1_share;
            uint64_t proxy2_share;

            unsigned char *ptr1 = proxy->GetBuffer1();
            unsigned char *ptr2 = proxy->GetBuffer2(); //HELPER'IN buffer'lara pointer'lar

            for(int i=0; i<size1; i++){ //bu optimal olmayabilir.
                for(int j=0; j<size2;j++){
                
                uint64_t tempShare = proxy->GenerateCommonRandom();  //elimizde GenerateRandom da var, ikisi de aynı çalışıyor gibi.
                proxy1_share = tempShare; //ilki rastgele bir sayı olacak.
                proxy2_share = noise_shares[i][j] - tempShare; //ikincisi de o rastgele ile bizim göndermek istediğimizin farkı.

                AddValueToCharArray(proxy1_share, &ptr1); //Helper proxy1 için olan verisini kendi 1. buffer'ına ve
                AddValueToCharArray(proxy2_share, &ptr2); //proxy2'nin verisini kendi 2. buffer'ına yüklüyor
            }      
            
            //Send(proxy->GetSocketP1(), proxy->GetBuffer1(), 8*size2); //asıl gönderme işlemi.
            //Send(proxy->GetSocketP2(), proxy->GetBuffer2(), 8*size2); //thread kullanmayı sevmiyorum. o yüzden kullanmadım :p
            }
            Send(proxy->GetSocketP1(), proxy->GetBuffer1(), 8*size2*size1); //asıl gönderme işlemi.
            Send(proxy->GetSocketP2(), proxy->GetBuffer2(), 8*size2*size1); //thread kullanmayı sev
            
            thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(),8*size2*size1);
            thread thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(),8*size2*size1);
            thr1.join();
            thr2.join();

            cout << "helper işini bitirdi." << endl;

            delete[]noise;
            return nullptr;
    }
    else{
        
        uint64_t** noise_shares = new uint64_t *[size1];
        for (int i = 0; i < size1; ++i) {
            noise_shares[i] = new uint64_t [size2];
        }
        cout << "receive başlıyor" << endl;
        for(int i=0; i<size1; ++i){ //size1 satır için, size2 tane elemanı oku.
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), 8*size2); 
        }
        cout << "receive bitiyor." << endl;

        unsigned char *ptr = proxy->GetBuffer1(); //KENDİ birinci buffer'ı (birincinin buffer'ı değil!!)     
        int counter = 0;
        for(int i=0; i<size1; i++){
            for(int j=0;j<size2;j++){
                cout << counter << endl;
                noise_shares[i][j] = ConvertToLong(&ptr);
                counter += 1;
            }
            
        }
        uint64_t **noiseful_shares = new uint64_t*[size1];
        for (int i = 0; i < size1; ++i) {
            noiseful_shares[i] = new uint64_t [size2];
        }
        cout << "noisful tayini tamam" << endl;

        for(int i=0; i<size1; i++){
            for(int j=0;j<size2;j++){
                noiseful_shares[i][j] = numbers[i][j]+noise_shares[i][j];
            }
        }
        for (int i = 0; i < size1; ++i) {
            delete[] noise_shares[i];
        }
        delete[] noise_shares;
        return noiseful_shares;
    }
}





uint64_t** invertMatrix(Party *const proxy,uint64_t** A, int n) {
    if (proxy->GetPRole() == helper) {
        // Divide(proxy,0,0);
        // Multiply(proxy,0,0);
        return nullptr;


    }else{
        // Create an augmented matrix by appending the identity matrix
        uint64_t** augmentedMatrix = new uint64_t*[n];
        uint64_t one=proxy->CreateShare(1.0);
        for (int i = 0; i < n; i++) {
            //cout<<"i:"<<i<<endl;
            augmentedMatrix[i] = new uint64_t[2 * n];
            for (int j = 0; j < n; j++) {
                augmentedMatrix[i][j] = A[i][j];
            }
            augmentedMatrix[i][i + n] = one;
        }
        int counter = 0;
        for (int i = 0; i < n; i++) {
            cout << "invertion counter: "<<counter<<endl; counter += 1;
            uint64_t pivot = augmentedMatrix[i][i];

            for (int j = 0; j < 2 * n; j++) {
                proxy->SendBytes(coreDivide); 
                uint64_t div=Divide(proxy,augmentedMatrix[i][j],pivot);
                augmentedMatrix[i][j]=div;
            }
            
            for (int k = 0; k < n; k++) {
                if (k != i) {
                    uint64_t factor = augmentedMatrix[k][i];
                    for (int j = 0; j < 2 * n; j++) {
                        proxy->SendBytes(coreMultiply);
                        augmentedMatrix[k][j] -= Multiply(proxy,factor,augmentedMatrix[i][j]);
                    }
                }
            }
        }
        uint64_t** inverseMatrix = new uint64_t*[n];
        for (int i = 0; i < n; i++) {
            inverseMatrix[i] = new uint64_t[n];
            for (int j = 0; j < n; j++) {
                inverseMatrix[i][j] = augmentedMatrix[i][j + n];
            }
        }

        for (int i = 0; i < n; i++) {
            delete[] augmentedMatrix[i];
        }
        delete[] augmentedMatrix;

        return inverseMatrix;
    }
}




//bazı 2d array'lerin transpozunu almak gerekecek.
    auto T = [](uint64_t **matrix, int r, int n){
        uint64_t** transpose = new uint64_t*[n];
        for(int i=0; i<n; i++){
            transpose[i] = new uint64_t[r];} //memory allocation baby!!
        for(int i=0; i<r; ++i){
            for(int j=0; j<n; ++j){
                transpose[j][i] = matrix[i][j];}}
        return transpose;
    };

uint64_t* GradientDescent_Turkmen(Party *const proxy, uint64_t **X, uint64_t *y, uint64_t *theta_shares, int const r, int const n, int iterations){ //lazımdan fazlasını yazmadım.
    if(proxy->GetPRole() == helper){
        cout << "Türkmen koşuyor" << endl;
        return 0;
    }
    else{

        double zeros[n];
        double minus_quarters[n];
        double minus_ones[n];
        for(int i=0; i<n; ++i){
            zeros[i] = 0;
            minus_quarters[i] = -0.25;
            minus_ones[i] = -1;
        }

        //TÜRKMEN:  beta_old = np.zeros(n)
        uint64_t *beta_old_shares = proxy->CreateShare(zeros,n); //notasyonu doğrudan takip edeceğim, çünkü sıkıldım :p
        //TÜRKMEN:  beta=np.zeros(n)   --> yenisini oluşturmak laızm.
        uint64_t *beta_shares = proxy->CreateShare(zeros,n); 

        //TÜRKMEN:  H=(-0.25)*np.matmul(X.T,X)
        uint64_t **X_transpose = T(X,r,n); //bir kere daha kullanılacak, isim veriyorum.

        uint32_t param_mult[3] = {uint32_t(n),uint32_t(r),uint32_t(n)};
        proxy->SendBytes(coreMatrixMatrixMultiply, param_mult, 3);
        uint64_t **X_transpose_X = MatrixMatrixMultiply(proxy, X_transpose, X, n, r, n); //(n,n) şeklinde olduğu malum
        cout << "X_T_T" << ConvertToDouble(Reconstruct(proxy, X_transpose_X[0][0])) << endl;
        cout << "X_T_T" << ConvertToDouble(Reconstruct(proxy, X_transpose_X[0][1])) << endl;
        cout << "X_T_T" << ConvertToDouble(Reconstruct(proxy, X_transpose_X[0][2])) << endl;
        uint64_t *minus_quarter_shares = proxy->CreateShare(minus_quarters, n); 

        uint64_t **H_shares = new uint64_t*[n]; //r rows HAYIR, bunun da şekli (n,n)!!!
        for(int i=0; i<n; ++i){
            H_shares[i] = new uint64_t[n];
        }

        uint64_t quadrantarium_share = proxy->CreateShare(-0.25);
        for(int i=0; i<n; ++i){
            for(int j=0; j<n; ++j){
                proxy->SendBytes(coreMultiply);
                H_shares[i][j] =  Multiply(proxy,X_transpose_X[i][j], quadrantarium_share);
            }
        }
        cout << "ilk bölüme kadar test ettik zaten" << endl;

        //buraya kadar bir test alalım.
        //buraya kadar çalıştığını BİLİYORUM.

        //TÜRKMEN:  H_inv=np.linalg.inv(H)

        uint64_t **H_inv_shares = invertMatrix(proxy, H_shares, n);

        cout << "Hinv " << ConvertToDouble(Reconstruct(proxy, H_inv_shares[0][0])) << endl;
        cout << "Hinv " << ConvertToDouble(Reconstruct(proxy, H_inv_shares[0][1])) << endl;
        cout << "Hinv " << ConvertToDouble(Reconstruct(proxy, H_inv_shares[0][2])) << endl;

        /*for i in range(iterations):
            beta=beta_old
            s=sigmoid(np.matmul(X,beta))
            beta_temp=y-s
            beta_grad=np.matmul(X.T,beta_temp)
            beta_old=beta-np.matmul(H_inv,beta_grad) */
        
        
        uint32_t param_vec_mult[2] = {uint32_t(r),uint32_t(n)}; //lazım olacak. X İLE ÇARPARKEN
        for(int i=0; i<iterations; ++i){
            beta_shares = beta_old_shares;
            //np.matmul(X,beta)
        
            proxy->SendBytes(coreMatrixVectorMultiply, param_vec_mult, 2);
            uint64_t *X_times_beta_shares = MatrixVectorMultiply(proxy, X , beta_shares, r, n);

            //s = sigmoid(..)
            uint32_t param_sig[1] = {uint32_t(r)};
            proxy->SendBytes(lgSigmoid, param_sig, 1);
            uint64_t *sigmoid_shares = Sigmoid(proxy, X_times_beta_shares, r);



            cout <<"sig of: "<< ConvertToDouble(Reconstruct(proxy,X_times_beta_shares[0]))<< " is " << ConvertToDouble(Reconstruct(proxy,sigmoid_shares[0])) << "\n";
            cout <<"sig of: "<< ConvertToDouble(Reconstruct(proxy,X_times_beta_shares[1]))<< " is " << ConvertToDouble(Reconstruct(proxy,sigmoid_shares[1])) << "\n";
            cout <<"sig of: "<< ConvertToDouble(Reconstruct(proxy,X_times_beta_shares[2]))<< " is " << ConvertToDouble(Reconstruct(proxy,sigmoid_shares[2])) << "\n";


            //buraya kadar çalıştığını biliom

            //beta_temp = y - s

            uint64_t *beta_temp_shares = new uint64_t[r]; //bunu silmek gerekiyor mu bilmiyorum (bellek bakımından)
            for(int i=0; i<r; ++i){
                beta_temp_shares[i] = y[i] - sigmoid_shares[i]; //bunu lokal yapmak lazım. Çıkarma için metot tanımlanmamış, gerek de yok.
            }
            cout << "sig "<< ConvertToDouble(Reconstruct(proxy, sigmoid_shares[0])) << " y " << ConvertToDouble(Reconstruct(proxy, y[0])) << " beta_temp " << ConvertToDouble(Reconstruct(proxy, beta_temp_shares[0]))<< endl;
            cout << "sig "<< ConvertToDouble(Reconstruct(proxy, sigmoid_shares[1])) << " y " << ConvertToDouble(Reconstruct(proxy, y[1])) << " beta_temp " << ConvertToDouble(Reconstruct(proxy, beta_temp_shares[1]))<< endl;
            cout << "sig "<< ConvertToDouble(Reconstruct(proxy, sigmoid_shares[2])) << " y " << ConvertToDouble(Reconstruct(proxy, y[2])) << " beta_temp " << ConvertToDouble(Reconstruct(proxy, beta_temp_shares[2]))<< endl;

            // beta_grad=np.matmul(X.T,beta_temp)
            uint32_t param_vec_mult_XT[2] = {uint32_t(n),uint32_t(r)};
            proxy->SendBytes(coreMatrixVectorMultiply, param_vec_mult_XT, 2);
            uint64_t *beta_grad_shares = MatrixVectorMultiply(proxy, X_transpose, beta_temp_shares, n, r); //yukarıdaki çarpımın tersi.
            //Bunun şekli (n,1)
            cout << "beta grad " << ConvertToDouble(Reconstruct(proxy, beta_grad_shares[0])) << endl;
            cout << "beta grad " << ConvertToDouble(Reconstruct(proxy, beta_grad_shares[1])) << endl;
            cout << "beta grad " << ConvertToDouble(Reconstruct(proxy, beta_grad_shares[2])) << endl;

            // beta_old=beta-np.matmul(H_inv,beta_grad)

            uint32_t param_vec_mult_Hinv_beta_grad[2] = {uint32_t(n),uint32_t(n)}; //şekiller değiştikçe yeni param yazıyorum :p
            proxy->SendBytes(coreMatrixVectorMultiply, param_vec_mult_Hinv_beta_grad, 2);
            uint64_t *H_inv_times_beta_grad_shares = MatrixVectorMultiply(proxy, H_inv_shares, beta_grad_shares, n, n);

            cout << "Hinv * beta grad " << ConvertToDouble(Reconstruct(proxy, H_inv_times_beta_grad_shares[0])) << endl;
            cout << "Hinv * beta grad " << ConvertToDouble(Reconstruct(proxy, H_inv_times_beta_grad_shares[1])) << endl;
            cout << "Hinv * beta grad " << ConvertToDouble(Reconstruct(proxy, H_inv_times_beta_grad_shares[2])) << endl;


            uint64_t *minus_one_shares = proxy->CreateShare(minus_ones, n);
            uint32_t param_vecmult_n[] = {uint32_t(n)};
            proxy->SendBytes(coreVectorisedMultiply, param_vecmult_n, 1);
            uint64_t *minus_h_inv_times_beta_grad_shares = Multiply(proxy, H_inv_times_beta_grad_shares, minus_one_shares, n);

            beta_old_shares = Add(proxy, beta_shares, minus_h_inv_times_beta_grad_shares, n);

        }


        //ÖNCE BELLEK TEMİZLİĞİ YAPMAK LAZIM:
        for(int i=0; i<n; i++){
            delete []H_shares[i]; //bu doğru olmayabilir.
        }
        return beta_shares; //son iterasyonu boşlamış olmuyor muyuz bunu döndürürken?
    }
}








uint64_t* GradientDescent(Party *const proxy, uint64_t **X, uint64_t *y, uint64_t* theta_shares, int const r, int const n, double const learning_rate, double epsilon, int iterations){ //sanırım burada vektör kullanmak zorunda kalacağım.
    uint32_t param_r[1] = {(uint32_t)r}; 
    uint32_t param_n[1] = {(uint32_t)n}; // n ve r gerektiren metotlar için param'lar.
    
    //bazı 2d array'lerin transpozunu almak gerekecek.
    /* auto T = [](uint64_t **matrix, int r, int n){
        uint64_t** transpose = new uint64_t*[n];
        for(int i=0; i<n; i++){
            transpose[i] = new uint64_t[r];} //memory allocation baby!!
        for(int i=0; i<r; ++i){
            for(int j=0; j<n; ++j){
                transpose[j][i] = matrix[i][j];}}
        return transpose;
    }; */

    if(proxy->GetPRole() == helper){
        //Bunu şu an doğru kullanmıyoruz. Halledeceğiz inş.
        return nullptr;
    }
    else{
        double minus_ones[r];     // 
        double halves[r];         // 
        double eights[n];         // 
        double twos[n];           // 
        double rs[n];             // 
        double learning_rates[n]; // 
        double minus_ones_n[n];   // 
        for(int i=0; i<r; i++){
            minus_ones[i] = -1;
            halves[i] = 0.5;
        }
        for(int i=0; i<n; i++){
            eights[i] = 8;
            twos[i] = 2;
            rs[i] = r;
            learning_rates[i] = learning_rate;
            minus_ones_n[i] = -1;
        }
        uint64_t *minus_one_shares = proxy->CreateShare(minus_ones, r);
        uint64_t *half_shares = proxy->CreateShare(halves, r);
        uint64_t *eight_shares_n = proxy->CreateShare(eights, n);                 // Bazı hesaplar için gizli paylar.
        uint64_t *two_shares_n = proxy->CreateShare(twos, n);                     // Bunları hesaplamanın daha iyi bir yolu olabilir.
        uint64_t *r_shares_n = proxy->CreateShare(rs, n);                         //
        uint64_t *learning_rate_shares_n = proxy->CreateShare(learning_rates, n); //
        uint64_t *minus_one_shares_n = proxy->CreateShare(minus_ones, n);         // 

        double d = n-1;
        double sensitivity = (3 * d) + 0.25 * (n * n); // calculates the sensitivity with the number of features
        double scale_double = sensitivity / epsilon * iterations;  // calculate beta, lokal hesaplanabilir bu, gizli bir yanı yok ki.
        
        uint32_t scale = uint32_t(ConvertToUint64(scale_double));
        //add_noise metodu scale'i uint32_t formunda bekliyor. O yüzden elimizle typecast ediyoruz.
        //Bunu buffer-socket işlerini bozduktan sonra denedik.
       
        for(int i=0; i<iterations; i++){
            //cout << "for loop'a girdik" << endl; KOD BOYUNCA BUNUN GİBİ DEBUG PRINTLERİ MEVCUT.

            //term1 = np.matmul(X.T, (0.5 - y)) + noise1        hesaplanacak
            uint64_t **X_transpose = T(X, r, n); //X.T tamam

            proxy->SendBytes(coreVectorisedMultiply, param_r, 1);
            uint64_t *minus_y_shares = Multiply(proxy, y, minus_one_shares, r); //bu çıkarmanın çok uzun bir yolu.

            uint64_t *half_minus_y_shares = Add(proxy, minus_y_shares, half_shares, r);
            //cout << "0.5 -y işi tamam" << endl;
            
            uint32_t param_vec_mult[2] = {uint32_t(n),uint32_t(r)};
            proxy->SendBytes(coreMatrixVectorMultiply, param_vec_mult, 2);
            uint64_t *term1_shares = MatrixVectorMultiply(proxy,X_transpose, half_minus_y_shares, n, r);

            //************************************************** Noise1 ekleniyor
            //term1 -> (n,1) şeklinde
            uint64_t *term1_before_noise = Reconstruct(proxy, term1_shares, n); //bu aşağıdaki debug çıktısı için

            uint32_t param_scale[2] = {uint32_t(n), scale}; //scale'i zaten tyepcast ettim.
            
            proxy->SendBytes(lgAddNoise, param_scale, 2);
            term1_shares=add_noise(proxy,term1_shares,uint32_t(n), scale); //burada n'yi typecast etmenin bir gereği olmasa gerek
            
           
            uint64_t *term1_for_debug = Reconstruct(proxy, term1_shares, n);
            // cout << "FOR TERM1-----------------------------" << endl << endl << endl;
            // cout << "before noise: "<< ConvertToDouble(term1_before_noise[0]) << " after noise: " << ConvertToDouble(term1_for_debug[0]) << endl;
            // cout << "before noise: "<< ConvertToDouble(term1_before_noise[1]) << " after noise: " << ConvertToDouble(term1_for_debug[1]) << endl;
            // cout << "before noise: "<< ConvertToDouble(term1_before_noise[2]) << " after noise: " << ConvertToDouble(term1_for_debug[2]) << endl;
            // cout << "FOR TERM1------------------------------" << endl << endl << endl;
            // *****************************************************

            
            //term2 = np.matmul(X.T, X) / 8 + noise2 

            uint32_t param_mult[3] = {uint32_t(n),uint32_t(r),uint32_t(n)};
            proxy->SendBytes(coreMatrixMatrixMultiply, param_mult, 3);
            uint64_t **X_transpose_X = MatrixMatrixMultiply(proxy, T(X,r,n), X, n, r, n);

            double *XTX_debug = ConvertToDouble(Reconstruct(proxy, X_transpose_X[0], n),n);
            // cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
            // cout << XTX_debug[0] << endl;
            // cout << XTX_debug[1] << endl;
            // cout << XTX_debug[2] << endl;
            // cout << XTX_debug[3] << endl;
            // cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl; //bu gayet çirkin bir debug. ama debug ist debug


            
            // Xt * X --> (n,n) şeklinde
            //X.T*X'in her elemanını 8'e bölmek lazım:
            //satır satır vektör muamelesi yapmak kabil.
            uint64_t **term2_shares = new uint64_t*[n]; 
            for(int j=0; j<n; j++){ 
                proxy->SendBytes(coreVectorisedDivide, param_n, 1); //bunu vektörize yapmak mümkün olmalı
                term2_shares[j] = Divide(proxy, X_transpose_X[j], eight_shares_n, n); 
            }
            //cout << "bölümden çıkıyor" << endl;
            //*********************************** Noise2 ekleniyor
            double *term2_firstline_debug = ConvertToDouble(Reconstruct(proxy, term2_shares[0], n),n); //TERM2'NİN SADECE İLK SATIRI!!!!!!!!
            //yukarıdaki yine aşağıdaki debug çıktısı için.

            for(int satr = 0; satr< n; ++satr){
                proxy->SendBytes(lgAddNoise, param_scale, 2);
                term2_shares[satr] = add_noise(proxy, term2_shares[satr], uint32_t(n), scale);
            }

            // BU AŞAĞIDAKİ ÜÇ SATIR, add_noise FONKSİYONUNU VEKTÖRİZE ETME DENEMEMİZİN BAŞARISIZLIĞI İLE YORUMA ALINMIŞTIR.
            // YUKARIDAKİ FOR LOOP DA AYNI ŞEYİ YAPIYOR, FAKAT 15 DEFA İLETİŞİM KURMA ZORUNLULUĞU VAR.

            // uint32_t param_noise_term2[3] = {uint32_t(n),uint32_t(n), scale}; //again, scale is already typecasted
            // proxy->SendBytes(lgAddNoiseMatrix, param_noise_term2, 3);
            // term2_shares = add_noise(proxy, term2_shares, uint32_t(n),uint32_t(n), scale);



            double *term2_noisy_first_line_debug = ConvertToDouble(Reconstruct(proxy, term2_shares[0], n),n);
            // cout << "TERM2*-*-*-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-" << endl;
            // cout << "before noise: "<<term2_firstline_debug[0] << " after noise: "<< term2_noisy_first_line_debug[0] <<endl;
            // cout << "before noise: "<<term2_firstline_debug[1] << " after noise: "<< term2_noisy_first_line_debug[1] <<endl;
            // cout << "before noise: "<<term2_firstline_debug[2] << " after noise: "<< term2_noisy_first_line_debug[2] <<endl;
            // cout << "before noise: "<<term2_firstline_debug[3] << " after noise: "<< term2_noisy_first_line_debug[3] <<endl;
            // cout << "TERM2*-*-*-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-" << endl;
            
            // ********************************************************************
            
            
            uint64_t *theta_reconst = Reconstruct(proxy, theta_shares, n);
            // cout << "THEATA------------------------------THEATA" << endl;
            // cout << ConvertToDouble(theta_reconst[0]) << endl;
            // cout << ConvertToDouble(theta_reconst[1]) << endl;
            // cout << ConvertToDouble(theta_reconst[2]) << endl;
            // cout << ConvertToDouble(theta_reconst[3]) << endl;
            // cout << ConvertToDouble(theta_reconst[4]) << endl; //thetaların başta sıfır olması lazım.
            // cout << "THEATA------------------------------THEATA" << endl;
            
            //gradient = term1 + 2 * np.matmul(term2, theta) ki burada theta n elemanlı bir vektör
            uint32_t param_vec_mult_for_term2_theta[2] = {uint32_t(n),uint32_t(n)};

            proxy->SendBytes(coreMatrixVectorMultiply, param_vec_mult_for_term2_theta, 2);
            uint64_t *term2_times_theta = MatrixVectorMultiply(proxy, term2_shares, theta_shares, n, n);

            uint64_t *term2_with_theta_for_debug = Reconstruct(proxy, term2_times_theta, n);
            // cout << "------------------------------term2*theta" << endl;
            // cout << ConvertToDouble(term2_with_theta_for_debug[0]) << endl;
            // cout << ConvertToDouble(term2_with_theta_for_debug[1]) << endl;
            // cout << ConvertToDouble(term2_with_theta_for_debug[2]) << endl;
            // cout << ConvertToDouble(term2_with_theta_for_debug[3]) << endl;
            // cout << "------------------------------term2*theta" << endl;

            //theta'yı güncellemek için son hesaplar:
            
            proxy->SendBytes(coreVectorisedMultiply, param_n, 1);
            uint64_t *double_term2_times_theta = Multiply(proxy, term2_times_theta, two_shares_n, n); //2 * np.matmul(term2, theta) 
            uint64_t *gradient_shares = Add(proxy, term1_shares, double_term2_times_theta, n); 
            //cout << "vektörize çarpımlar !!" << endl;
            //theta -= learning_rate * (gradient / m)
            proxy->SendBytes(coreVectorisedDivide, param_n, 1);
            uint64_t *gradient_div_r_shares = Divide(proxy,gradient_shares, r_shares_n, n); // gradient / n
            //cout << "gradient/n !!!!" << endl;
            proxy->SendBytes(coreVectorisedMultiply, param_n, 1);
            uint64_t *learn_rate_times_grad_div_r_shares = Multiply(proxy, learning_rate_shares_n, gradient_div_r_shares, n); 
            //cout << "final" << endl;
            proxy->SendBytes(coreVectorisedMultiply, param_n, 1); //ya bu çok saçma ama çok güvenilir bir yöntem çıkarma yapmak için :/
            uint64_t *substract_from_theta = Multiply(proxy, learn_rate_times_grad_div_r_shares, minus_one_shares_n, n); 
            theta_shares = Add(proxy, theta_shares, substract_from_theta, n);
        }

        // delete[] noise1; muhtemelen silmek lazım.
        // delete[] noise2;
        
        return theta_shares;
    }
}







 

#endif //CORE_H


