# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 13:41:29 2021

@author: trevo
"""

import sys, os, marshal, timeit, numpy as np
 
class Prime:
    # Initialize loop parameters, c=max search range, c=search chunk size,
    # debug_en=debug enable, dbname=database file name
    def __init__(self, n=1000000, c=None, debug_en=False, dbname='primes.dat'):
        self.n        = n
        self.c        = n if c is None else c
        self.dbname  = dbname
        self.debug_en = debug_en
        if os.path.exists(self.dbname) and (self.n > self.c):
            print "Removing existing database file"
            os.remove(self.dbname)
        # We already know the first prime (2)
        self._nprimes = 1
    
    # search will find all primes in a given range (self.n)
    # in one large memory stored array
    def search(self):
        # Max range to search for primes
        SQRTN   = int(np.sqrt(self.n)+1)
        # Generate NDARRAY of bools, assume all are Prime=True to start
        # This array represents odd numbers > 2
        isprime = np.full((self.n-1)/2, True,dtype=bool)
        
        # Start with 3 and search from there
        k = 3
        while k < SQRTN:
            # Mark all factors of k False starting at k**2
            if isprime[k/2-1]: isprime[k**2/2-1::k] = False
            # Index by 2 since we know evens are not prime
            k += 2
        # All True's are prime - convert to numeral and prepend '2'
        self.p = np.concatenate(([2], np.where(isprime)[0]*2 + 3))
        self._nprimes = self.p.size
    
    # _search_p marks factors of known primes
    def _search_p(self, isprime, SQRTL, first, last, arrlen):
        with open(self.dbname, 'rb') as dbfile:
            try:
                while True:
                    # Increment through the database file until we run
                    # out of known primes or exceed search range
                    pload = np.frombuffer(marshal.load(dbfile), dtype=np.int64)
                    for k in pload:
                        if (k<SQRTL):
                            ksq     = k**2
                            # Find first factor using k**2 unless
                            # not in current array range
                            start   = k if ksq < first else ksq
                            isprime = self._update(isprime, k, start, first, last, arrlen)
                        else:
                            return isprime
            except:
                # Run out of known primes
                pass
        return isprime
    
    # _search_k marks factors of new primes
    def _search_k(self, isprime, SQRTL, k, first, last, arrlen):
        # Initiate bool pointer
        ip = 0
        while k < SQRTL:
            if isprime[ip]:
                # Current k**2 must fall within this range or
                # the loop would have exited, start at k**2
                isprime = self._update(isprime, k, k**2, first, last, arrlen)
            # Increment bool pointer by 1
            ip += 1
            # Index value under test by 2 since evens are not prime
            k  += 2
        return isprime, k
    
    # _update marks all factors of current prime False
    def _update(self, isprime, k, start, first, last, arrlen):
        # Solve modulus of value represented by first index
        # and the starting value (k or k**2)
        m   = first%start
        # Odd modulus represents an even factor so find next factor
        ix  = m if m==0 else (start - m)/2 if m&1 else (2*start - m)/2
        # Mark factors of k starting at first index
        if ix < arrlen: isprime[ix::k] = False
        return isprime
    
    # search_c will find all primes in a given range (self.n)
    # by breaking that range into chunks (size.c) and storing 
    # previously found primes (primes.p) on disk
    def search_c(self):
        # Max range to search for primes in big loop
        SQRTN   = int(np.sqrt(self.n)+1)
        # Initiatie loop variables
        k      = 3
        first  = k
        while k < SQRTN or first < self.n:
            # Solve for length of the current array; can only  
            # be smaller than self.c for the last chunk
            arrlen  = min(self.c, (self.n-k+2)/2)
            # Value represented by last index in current iteration
            last    = k + 2*(arrlen - 1)
            # Generate NDARRAY of bools, assume all are Prime=True to start
            # This array represents odd numbers only
            isprime = np.full(arrlen, True,dtype=bool)
            # Max range to search for primes in small loop
            SQRTL   = int(np.sqrt(last))+1
            # Print debug information so we can keep track of search status
            if self.debug_en:
                print("first=%s, last=%s, k=%s, "
                      "SQRTL=%s, SQRTN=%s, nprimes=%s"
                      %(first, last, k, SQRTL, SQRTN, self._nprimes))
            # If any primes have been found, mark their factors False
            if os.path.exists(self.dbname):
                isprime    = self._search_p(isprime, SQRTL, first, last, arrlen)
            # Mark factors of new primes False
            isprime, k = self._search_k(isprime, SQRTL, k, first, last, arrlen)
            # Efficiently serialize (marshal) prime numbers to disk
            with open(self.dbname, 'ab') as dbfile:
                pdump = np.where(isprime)[0]*2 + first
                marshal.dump(pdump, dbfile)
            # Keep track of how many prime numbers have been found
            self._nprimes += pdump.size
            # Value represented by first index in next iteration
            first += 2*arrlen
            # First prime to use when marking factors of new primes
            k = first if pdump.size==0 or first>pdump[-1] else pdump[-1]

def run_search(n=None, c=None, debug_en=None, dbname=None):
    if n is None:               n = int(1e6)
    if c is None:               c = n
    if debug_en is None: debug_en = True
    if dbname is None:   dbname   = 'primes.dat'
    
    prime = Prime(n, c, debug_en, dbname)
    if prime.debug_en:
        print prime.n, prime.c, prime.debug_en, prime.dbname
    
    if prime.n > prime.c:
        # Run chunked search if chunks are smaller than range
        if prime.debug_en:
            print("Storing primes in .dat file")
        prime.search_c()
    else:
        # Run single array search otherwise
        if prime.debug_en:
            print("Completing all processing in memory")
        prime.search()
    return prime

# run file as a script
if __name__=='__main__':
    # Call code by name with at least one input argument
    # Chunk size cap based on 2GB of int64 primes
    # Set second argument to "max" to use max chunk size
    c, debug_en, dbname = [None]*3
    nargin = len(sys.argv)
    print "Running ",sys.argv[0]
    if nargin > 1:
        print "arg1 =",sys.argv[1]
        n = int(sys.argv[1])
    else:
        sys.exit("At least one input argument required for __name__==\'main\'")
    if nargin > 2:
        print "arg2 =",sys.argv[2]
        if sys.argv[2] == "max":
            c = 268435360
        else:
            c = min(268435360, int(sys.argv[2]))
    if nargin > 3:
        print "arg3 =",sys.argv[3]
        debug_en = bool(sys.argv[3])
    if nargin > 4:
        print "arg4 =",sys.argv[4]
        dbname = str(sys.argv[4]).split('.')[0]
        if ".dat" not in dbname:
            dbname += ".dat"
        print "Database file name is",dbname
    
    t0 = timeit.default_timer()
    prime = run_search(n, c, debug_en, dbname)
    t1 = timeit.default_timer()
    td = t1 - t0
    if prime.debug_en:
        print("n = %s, found %s primes in  %s seconds"
              %(prime.n, prime._nprimes, td))
