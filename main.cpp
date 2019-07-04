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
#include <climits>
#include <cfloat>
#include <cfenv>
#include <cmath>
#include <array>
#include <bitset>

using namespace std;

#define MAX_NUM (2000000000LL)
#define PAGES (4)
#define MAX_PAGE (MAX_NUM/PAGES)
#define BITS_IN_ULL ((__CHAR_BIT__)*(__SIZEOF_LONG_LONG__))
#define MAX_NUM_SZ (MAX_NUM/BITS_IN_ULL)
#define NUM_PRIMES_TO_MASKTAB (12)
#define OFFSET_MASKTAB (64)


long int lastprime2Bil_idx;

#define MULTIPLE_THREAD

#define BASE_PRIMES_SZ (5000)
unsigned long long int base_primes[BASE_PRIMES_SZ];
unsigned long long int base_mask_primes[NUM_PRIMES_TO_MASKTAB*OFFSET_MASKTAB];

void printOutAll(unsigned long long int *ints)
{
    stringstream fname;
    
    fname << "primesAllBits.txt";
    
    ofstream fout(fname.str());
    
    for (unsigned long int idx=0; idx<MAX_NUM_SZ; idx++)
        fout << bitset<64>(ints[idx]) << endl;
    
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
    unsigned long long int *bit = new unsigned long long int[BITS_IN_ULL];
    unsigned long long int *mask_primes = new unsigned long long int[NUM_PRIMES_TO_MASKTAB*OFFSET_MASKTAB];
    unsigned long long int *ints = new unsigned long long int[sized];

    for (int idx=0; idx<BITS_IN_ULL; idx++) {
        bit[idx] = ~(1ULL<<idx);
    }

    // STEP 1: copying primes table to local copy
    std::copy(&base_primes[0],&base_primes[BASE_PRIMES_SZ],primes);

    std::fill(&ints[0],&ints[sized],0xAAAAAAAAAAAAAAAAULL);
    ints[0] &= ~0x3ULL;
    ints[0] |=  0x4ULL;
    
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
        ints[0] |= 1ULL<<prime; // set prime as a prime number
    }
     
    for (idx = NUM_PRIMES_TO_MASKTAB; idx <= lastprime2Bil_idx; idx++) {
        unsigned long long int prime = primes[idx];
        unsigned long long int lastnum = sized*BITS_IN_ULL;
        for (unsigned long long int num=2*prime; num<lastnum; num+=prime) {
            register int ints_idx = num>>6;
            register int bit_idx = num&63;
            ints[ints_idx] &= bit[bit_idx];
        }
    }

//    printOutAll(ints);
#ifdef DEBUG
    std::stringstream fname;
    fname << "primesFr" << startd*BITS_IN_ULL/1000000 << "MTo" << (startd+sized)*BITS_IN_ULL/1000000 << "M.txt";
    std::ofstream fout(fname.str());
    unsigned long int sz = sized*BITS_IN_ULL;
    for (unsigned long int num=0; num<sz; num++) {
        unsigned long int idx = num>>6;
        if ((~bit[num&63]) & ints[idx]) 
            fout << num+startd*BITS_IN_ULL << endl;
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
    unsigned long long int *bit = new unsigned long long int[BITS_IN_ULL];
    unsigned long long int *mask_primes = new unsigned long long int[NUM_PRIMES_TO_MASKTAB*OFFSET_MASKTAB];
    unsigned long long int *ints = new unsigned long long int[sized];

    for (int idx=0; idx<BITS_IN_ULL; idx++) {
        bit[idx] = ~(1ULL<<idx);
    }

    // STEP 1: copying primes table to local copy
    std::copy(&base_primes[0],&base_primes[BASE_PRIMES_SZ],primes);

    std::fill(&ints[0],&ints[sized],0xAAAAAAAAAAAAAAAAULL);
    
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
        unsigned long long int prime = primes[idx];
        unsigned long long int lastnum = (startd+sized)*BITS_IN_ULL;
        for (unsigned long long int num=((startd*BITS_IN_ULL+prime-1)/prime)*prime; num<lastnum; num+=prime) {
            register int ints_idx = num>>6;
            register int bit_idx = num&63;
            ints[ints_idx-startd] &= bit[bit_idx];
        }
    }
    
#ifdef DEBUG
    std::stringstream fname;
    fname << "primesFr" << startd*BITS_IN_ULL/1000000 << "MTo" << (startd+sized)*BITS_IN_ULL/1000000 << "M.txt";
    std::ofstream fout(fname.str());
    unsigned long int sz = sized*BITS_IN_ULL;
    for (unsigned long int num=0; num<sz; num++) {
        unsigned long int idx = num>>6;
        if ((~bit[num&63]) & ints[idx]) 
            fout << num+startd*BITS_IN_ULL << endl;
    }
#endif
    
    delete[] primes;
    delete[] bit;
    delete[] mask_primes;
    delete[] ints;
}

#endif

#ifndef MULTIPLE_THREAD
unsigned long long int ints[MAX_NUM_SZ];
unsigned long long int bit[BITS_IN_ULL];

void init(void)
{
    long int prime;
    
    for (int idx=0; idx<BITS_IN_ULL; idx++) {
        bit[idx] = ~(1ULL<<idx);
    }
        
    first_primes();

    // precalculating masks for the lowest value primes (2,3,5,7,11,13)
    // This is the most important optimization
    // as it reduces the number of sieve loops considerably.
    {
        prime =2 ;
        base_mask_primes[0] = ~0ULL;
        for (long int num=0; num<BITS_IN_ULL; num+=prime) {
            base_mask_primes[0] |= 1ULL<<num;
        }
        base_mask_primes[0] = ~base_mask_primes[0];
    }
    
    for (long int prime_idx=1; prime_idx<NUM_PRIMES_TO_MASKTAB; prime_idx++) {
        long int prime = base_primes[prime_idx];
        for (long int idx=0; idx<prime; idx++)
            base_mask_primes[prime_idx*OFFSET_MASKTAB+idx] = ~0ULL;
        for (long int num=0; num < prime*BITS_IN_ULL; num+=prime) {
            long int idx = num / BITS_IN_ULL;
            long int offset = num % BITS_IN_ULL;
            base_mask_primes[prime_idx*OFFSET_MASKTAB+idx] &= ~(1ULL<<offset);
        }
    }

    std::fill(&ints[0],&ints[MAX_NUM_SZ],0xAAAAAAAAAAAAAAAAULL);
}

void sieve_erat(void)
{
    ints[0] &= ~0x3ULL;
    ints[0] |=  0x4ULL;
    // these are the most important optimizations
    for (long long int prime_idx=1; prime_idx<NUM_PRIMES_TO_MASKTAB; prime_idx++) {
        unsigned long int ints_idx, mask_idx;
        unsigned long int prime;
        
        
        // these are the most important optimizations
        // sieve for prime
        prime = base_primes[prime_idx];
        for (ints_idx=0, mask_idx=0; ints_idx<MAX_NUM_SZ; ints_idx++) {
            ints[ints_idx] &= base_mask_primes[prime_idx*OFFSET_MASKTAB+mask_idx++];
            mask_idx %= prime;
        }
        ints[0] |= 1ULL<<prime; // set prime as a prime number
    }
    
    for (long long int idx = NUM_PRIMES_TO_MASKTAB; idx <= lastprime2Bil_idx; idx++) {
        unsigned long long int pr = base_primes[idx];
        for (unsigned long long int num=2*pr; num<MAX_NUM; num+=pr) {
            register int ints_idx = num>>6;
            register int bit_idx = num&63;
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
        long int idx = num>>6;
        if ((~bit[num&63]) & ints[idx]) 
            cout << num << endl;
    }
#endif
	return 0;
}
#else
void init(void)
{
    first_primes();

    // this is an abberation since 64 (number of bits in long long int) is divisible by
    // the prime 2
    std::fill(base_mask_primes,&base_mask_primes[NUM_PRIMES_TO_MASKTAB*OFFSET_MASKTAB],0xFFFFFFFFFFFFFFFFULL);
    
    for (unsigned long int prime_idx=0; prime_idx<BITS_IN_ULL; prime_idx+=2)
      base_mask_primes[0] |= 1ULL<<prime_idx;
    base_mask_primes[0]=~base_mask_primes[0];
    
    for (unsigned long int prime_idx=1; prime_idx<NUM_PRIMES_TO_MASKTAB; prime_idx++) {
        unsigned long int prime=base_primes[prime_idx];
        for (unsigned long int num=0; num < prime*BITS_IN_ULL; num+=prime) {
            unsigned long int idx = num / BITS_IN_ULL;
            unsigned long int offset = num % BITS_IN_ULL;
            base_mask_primes[prime_idx*OFFSET_MASKTAB+idx] &= ~(1ULL<<offset);
        }
    }
}

//
int main(int argc, const char** argv) {
    auto start = std::chrono::system_clock::now();
#ifndef DEBUG
    unsigned num_cpus = std::thread::hardware_concurrency();
#endif
    
    num_cpus = 4;
    
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
