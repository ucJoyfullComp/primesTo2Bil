// Calculating all primes from 2 to 2000000000
// optimize for time
//
// Uzi Cohen © 2019
//

// Added code for using all CPUs on machine, adapted from:
//
// Sample of launching C++11 threads that report what CPU they run on.
//
// Eli Bendersky [http://eli.thegreenplace.net]
// This code is in the public domain.
#include <algorithm>
#include <chrono>
#include <mutex>
#include <sched.h>
#include <thread>
#include <vector>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdint>
#include <climits>
#include <cfloat>
#include <cfenv>
#include <cmath>
#include <array>
#include <bitset>

using namespace std;

//#define __128BIT_DATA

#ifdef __128BIT_DATA
typedef __int128 MYINT;
typedef unsigned __int128 MYUINT;
#define BITS_IN_MYUINT ((__CHAR_BIT__)*(__SIZEOF_INT128__))
#define ULL_SHR (7)
#else
typedef long long int MYINT;
typedef unsigned long long int MYUINT;
#define BITS_IN_MYUINT ((__CHAR_BIT__)*(__SIZEOF_LONG_LONG__))
#define ULL_SHR (6)
#endif

#define MAX_NUM (2000000000LL)
#define ULL_MASK (BITS_IN_MYUINT-1)
#define MAX_NUM_SZ (MAX_NUM/BITS_IN_MYUINT)
#define NUM_PRIMES_TO_MASKTAB (12)
#define OFFSET_MASKTAB (BITS_IN_MYUINT)


long int lastprime2Bil_idx;

#define MULTIPLE_THREAD

#define BASE_PRIMES_SZ (5000)
MYUINT base_primes[BASE_PRIMES_SZ];
MYUINT base_mask_primes[NUM_PRIMES_TO_MASKTAB*OFFSET_MASKTAB];

void printOutAll(MYUINT *ints)
{
    stringstream fname;
    
    fname << "primesAllBits.txt";
    
    ofstream fout(fname.str());
    
    for (unsigned long int idx=0; idx<MAX_NUM_SZ; idx++)
        fout << bitset<BITS_IN_MYUINT>(ints[idx]) << endl;
    
    fout.close();
}

void first_primes(void)
{
    vector<unsigned long int> vPrimes, vNums;
    unsigned long int max_prime = sqrt(MAX_NUM);
    
    vPrimes.clear();
    vPrimes.insert(vPrimes.end(),2);
    
    vNums.clear();
    for (uint32_t idx=0; idx<=max_prime; idx++) {
        vNums.insert(vNums.end(),(idx%2)?1:0);
    }
    vNums[0] = 0;
    vNums[1] = 0;
    vNums[2] = 1;
    for (uint32_t idx=3; idx<=max_prime; idx++) {
        if (vNums[idx]) {
            vPrimes.insert(vPrimes.end(),idx);
            for (uint32_t idx_idx=2*idx; idx_idx<=max_prime; idx_idx+=idx) {
                vNums[idx_idx] = 0;
            }
        }
    }
    lastprime2Bil_idx = vPrimes.size() - 1;
    for (long int idx=0; idx <= lastprime2Bil_idx; idx++) base_primes[idx] = vPrimes[idx];
    for (long int idx=lastprime2Bil_idx+1; idx < BASE_PRIMES_SZ; idx++) base_primes[idx] = 0;
}

#ifdef MULTIPLE_THREAD
void sieve_optimized1(unsigned int startd, unsigned int sized)
{
    cout << "Thread id: " << std::this_thread::get_id() << std::endl;
    
//    unsigned int startd = MAX_NUM_SZ/4;
//    unsigned int sized = MAX_NUM_SZ/4;

    unsigned long int prime_idx;
    unsigned long int *primes = new unsigned long int[BASE_PRIMES_SZ];
    MYUINT *bit = new MYUINT[BITS_IN_MYUINT];
    MYUINT *mask_primes = new MYUINT[NUM_PRIMES_TO_MASKTAB*OFFSET_MASKTAB];
    MYUINT *ints = new MYUINT[sized];

    for (int idx=0; idx<BITS_IN_MYUINT; idx++) {
        bit[idx] = ~(((MYUINT)1ULL)<<idx);
    }

    // STEP 1: copying primes table to local copy
    std::copy(&base_primes[0],&base_primes[BASE_PRIMES_SZ],primes);

#ifdef __128BIT_DATA
    std::fill(&ints[0],&ints[sized],((MYUINT)0xAAAAAAAAAAAAAAAAULL<<64)|(MYUINT)0xAAAAAAAAAAAAAAAAULL);
#else
    std::fill(&ints[0],&ints[sized],0xAAAAAAAAAAAAAAAAULL);
#endif
    ints[0] &= ~((MYUINT)0x3ULL);
    ints[0] |=  ((MYUINT)0x4ULL);
    
    // STEP 2: generating mask table for optimization of the sieve.
    std::copy(&base_mask_primes[0],&base_mask_primes[NUM_PRIMES_TO_MASKTAB*OFFSET_MASKTAB],mask_primes);

    long int idx;

    // sieve the primes
    for (prime_idx=1; prime_idx<NUM_PRIMES_TO_MASKTAB; prime_idx++) {
        unsigned long int ints_idx, mask_idx;
        unsigned long int prime;

        // these are the most important optimizations
        // sieve for prime
        prime = primes[prime_idx];
        for (ints_idx=0, mask_idx=startd%prime; ints_idx<sized; ints_idx++) {
            ints[ints_idx] &= mask_primes[prime_idx*OFFSET_MASKTAB+mask_idx];
            ++mask_idx %= prime;
        }
        ints[0] |= ((MYUINT)1ULL) << prime; // set prime as a prime number
    }
     
    for (idx = NUM_PRIMES_TO_MASKTAB; idx <= lastprime2Bil_idx; idx++) {
        MYUINT prime = primes[idx];
        MYUINT lastnum = sized*BITS_IN_MYUINT;
        for (MYUINT num=2*prime; num<lastnum; num+=prime) {
            int ints_idx = num >> ULL_SHR;
            int bit_idx = num & ULL_MASK;
            ints[ints_idx] &= bit[bit_idx];
        }
    }

//    printOutAll(ints);
#ifdef DEBUG
    std::stringstream fname;
    fname << "primesFr" << startd*BITS_IN_MYUINT/1000000 << "MTo" << (startd+sized)*BITS_IN_MYUINT/1000000 << "M.txt";
    std::ofstream fout(fname.str());
    unsigned long int sz = sized*BITS_IN_MYUINT;
    for (unsigned long int num=0; num<sz; num++) {
        unsigned long int idx = num >> ULL_SHR;
        if ((~bit[num & ULL_MASK]) & ints[idx]) 
            fout << num+startd*BITS_IN_MYUINT << endl;
    }
#endif
    
    delete[] primes;
    delete[] bit;
    delete[] mask_primes;
    delete[] ints;
}

void sieve_optimized_other(unsigned int startd, unsigned int sized)
{
    cout << "Thread id: " << std::this_thread::get_id() << std::endl;
    
//    unsigned int startd = MAX_NUM_SZ/4;
//    unsigned int sized = MAX_NUM_SZ/4;

    unsigned long int prime_idx;
    unsigned long int *primes = new unsigned long int[BASE_PRIMES_SZ];
    MYUINT *bit = new MYUINT[BITS_IN_MYUINT];
    MYUINT *mask_primes = new MYUINT[NUM_PRIMES_TO_MASKTAB*OFFSET_MASKTAB];
    MYUINT *ints = new MYUINT[sized];

    for (int idx=0; idx<BITS_IN_MYUINT; idx++) {
        bit[idx] = ~(((MYUINT)1ULL)<<idx);
    }

    // STEP 1: copying primes table to local copy
    std::copy(&base_primes[0],&base_primes[BASE_PRIMES_SZ],primes);

#ifdef __128BIT_DATA
    std::fill(&ints[0],&ints[sized],((MYUINT)0xAAAAAAAAAAAAAAAAULL<<64)|(MYUINT)0xAAAAAAAAAAAAAAAAULL);
#else
    std::fill(&ints[0],&ints[sized],0xAAAAAAAAAAAAAAAAULL);
#endif
    
    // STEP 2: generating mask table for optimization of the sieve.
    std::copy(&base_mask_primes[0],&base_mask_primes[NUM_PRIMES_TO_MASKTAB*OFFSET_MASKTAB],mask_primes);

    long int idx;

    // sieve the primes
    for (prime_idx=1; prime_idx<NUM_PRIMES_TO_MASKTAB; prime_idx++) {
        unsigned long int ints_idx, mask_idx;
        unsigned long int prime;
        
        
        // these are the most important optimizations
        // sieve for prime
        prime = primes[prime_idx];
        for (ints_idx=0, mask_idx=startd%prime; ints_idx<sized; ints_idx++) {
            ints[ints_idx] &= mask_primes[prime_idx*OFFSET_MASKTAB+mask_idx];
            ++mask_idx %= prime;
        }
    }
     
    for (idx = NUM_PRIMES_TO_MASKTAB; idx <= lastprime2Bil_idx; idx++) {
        MYUINT prime = primes[idx];
        MYUINT lastnum = sized*BITS_IN_MYUINT;
        for (MYUINT num=2*prime; num<lastnum; num+=prime) {
            int ints_idx = num >> ULL_SHR;
            int bit_idx = num & ULL_MASK;
            ints[ints_idx] &= bit[bit_idx];
        }
    }
    
#ifdef DEBUG
    std::stringstream fname;
    fname << "primesFr" << startd*BITS_IN_MYUINT/1000000 << "MTo" << (startd+sized)*BITS_IN_MYUINT/1000000 << "M.txt";
    std::ofstream fout(fname.str());
    unsigned long int sz = sized*BITS_IN_MYUINT;
    for (unsigned long int num=0; num<sz; num++) {
        unsigned long int idx = num >> ULL_SHR;
        if ((~bit[num & ULL_MASK]) & ints[idx]) 
            fout << num+startd*BITS_IN_MYUINT << endl;
    }
#endif
    
    delete[] primes;
    delete[] bit;
    delete[] mask_primes;
    delete[] ints;
}

#endif

#ifndef MULTIPLE_THREAD
MYUINT ints[MAX_NUM_SZ];
MYUINT bit[BITS_IN_MYUINT];

void init(void)
{
    long int prime;
    
    for (int idx=0; idx<BITS_IN_MYUINT; idx++) {
        bit[idx] = ~(((MYUINT)1ULL)<<idx);
    }
        
    first_primes();

    // precalculating masks for the lowest value primes (2,3,5,7,11,13)
    // This is the most important optimization
    // as it reduces the number of sieve loops considerably.
    {
        prime =2 ;
        base_mask_primes[0] = ~((MYUINT)0ULL);
        for (long int num=0; num<BITS_IN_MYUINT; num+=prime) {
            base_mask_primes[0] |= ((MYUINT)1ULL) << num;
        }
        base_mask_primes[0] = ~base_mask_primes[0];
    }
    
    for (long int prime_idx=1; prime_idx<NUM_PRIMES_TO_MASKTAB; prime_idx++) {
        long int prime = base_primes[prime_idx];
        for (long int idx=0; idx<prime; idx++)
            base_mask_primes[prime_idx*OFFSET_MASKTAB+idx] = ~((MYUINT)0ULL);
        for (long int num=0; num < prime*BITS_IN_MYUINT; num+=prime) {
            long int idx = num / BITS_IN_MYUINT;
            long int offset = num % BITS_IN_MYUINT;
            base_mask_primes[prime_idx*OFFSET_MASKTAB+idx] &= ~(((MYUINT)1ULL) << offset);
        }
    }

#ifdef __128BIT_DATA
    std::fill(&ints[0],&ints[MAX_NUM_SZ],((MYUINT)0xAAAAAAAAAAAAAAAAULL<<64)|(MYUINT)0xAAAAAAAAAAAAAAAAULL);
#else
    std::fill(&ints[0],&ints[MAX_NUM_SZ],0xAAAAAAAAAAAAAAAAULL);
#endif
}

void sieve_erat(void)
{
    ints[0] &= ~((MYUINT)0x3ULL);
    ints[0] |=  ((MYUINT)0x4ULL);
    // these are the most important optimizations
    for (MYINT prime_idx=1; prime_idx<NUM_PRIMES_TO_MASKTAB; prime_idx++) {
        unsigned long int ints_idx, mask_idx;
        unsigned long int prime;
        
        
        // these are the most important optimizations
        // sieve for prime
        prime = base_primes[prime_idx];
        for (ints_idx=0, mask_idx=0; ints_idx<MAX_NUM_SZ; ints_idx++) {
            ints[ints_idx] &= base_mask_primes[prime_idx*OFFSET_MASKTAB+mask_idx++];
            mask_idx %= prime;
        }
        ints[0] |= ((MYUINT)1ULL) << prime; // set prime as a prime number
    }
    
    for (MYINT idx = NUM_PRIMES_TO_MASKTAB; idx <= lastprime2Bil_idx; idx++) {
        MYUINT pr = base_primes[idx];
        for (MYUINT num=2*pr; num<MAX_NUM; num+=pr) {
            int ints_idx = num >> ULL_SHR;
            int bit_idx = num & ULL_MASK;
            ints[ints_idx] &= bit[bit_idx];
        }
    }
}

int main(int argc, char **argv)
{
    init();
    sieve_erat();
#ifdef DEBUG
    for (long int num=0; num<2000000000; num++) {
        long int idx = num >> ULL_SHR;
        if ((~bit[num & ULL_MASK]) & ints[idx]) 
            cout << num << endl;
    }
#endif
	return 0;
}
#else
void init(void)
{
    first_primes();

    // this is an abberation since 64 (number of bits in MYINT) is divisible by
    // the prime 2
#ifdef __128BIT_DATA
    std::fill(base_mask_primes,&base_mask_primes[NUM_PRIMES_TO_MASKTAB*OFFSET_MASKTAB],((MYUINT)0xFFFFFFFFFFFFFFFFULL<<64)|(MYUINT)0xFFFFFFFFFFFFFFFFULL);
#else
    std::fill(base_mask_primes,&base_mask_primes[NUM_PRIMES_TO_MASKTAB*OFFSET_MASKTAB],0xFFFFFFFFFFFFFFFFULL);
#endif
    
    for (unsigned long int prime_idx=0; prime_idx<BITS_IN_MYUINT; prime_idx+=2)
      base_mask_primes[0] |= ((MYUINT)1ULL) << prime_idx;
    base_mask_primes[0]=~base_mask_primes[0];
    
    for (unsigned long int prime_idx=1; prime_idx<NUM_PRIMES_TO_MASKTAB; prime_idx++) {
        unsigned long int prime=base_primes[prime_idx];
        for (unsigned long int num=0; num < prime*BITS_IN_MYUINT; num+=prime) {
            unsigned long int idx = num / BITS_IN_MYUINT;
            unsigned long int offset = num % BITS_IN_MYUINT;
            base_mask_primes[prime_idx*OFFSET_MASKTAB+idx] &= ~(((MYUINT)1ULL) << offset);
        }
    }
}

//
int main(int argc, const char** argv) {
    auto start = std::chrono::system_clock::now();
    unsigned num_cpus = std::thread::hardware_concurrency();

    std::cout << "Launching " << num_cpus << " threads\n";

    init();
    
    // A mutex ensures orderly access to std::cout from multiple threads.
    std::mutex iomutex;
    std::vector<std::thread> threads(num_cpus-1);

    if (num_cpus > 1)
        for (unsigned thread_idx=0; thread_idx<(num_cpus-1); thread_idx++)
            threads[thread_idx] = std::thread(sieve_optimized_other, (thread_idx+1)*MAX_NUM_SZ/num_cpus, MAX_NUM_SZ/num_cpus);
    sieve_optimized1(0, MAX_NUM_SZ/num_cpus);

    if (num_cpus > 1) {
        for(auto& t : threads) {
            t.join();
        }
    }
#ifndef DEBUG
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Calculation time: " << elapsed.count() << " mili-seconds" << endl;
#endif
    return 0;
}
#endif
