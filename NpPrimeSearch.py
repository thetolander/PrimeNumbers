# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 13:41:29 2021

@author: trevo
"""

import numpy as np
import time
 
class Prime:
    def __init__(self, n=1000000):
        self.n = n
    
    def search(self):
        # Max range to search
        SQRTN        = int(np.sqrt(self.n)+1)
        # Generate NDARRAY of bools, assume all are Prime=True to start
        # This array represents odd numbers > 2
        isprime      = np.full((self.n-1)/2, True,dtype=bool)
        
        # Start with 3 and search from there
        k = 3
        while k < SQRTN:
            # Mark all factors of k False starting at k**2
            if isprime[k/2-1]: isprime[k**2/2-1::k] = False
            # Index by 2 since we know evens are not prime
            k += 2
        # All True's are prime - convert to numeral and prepend '2'
        self.p = np.concatenate(([2], np.where(isprime)[0]*2 + 3))

# run file as a script
if __name__=='__main__':
    n = int(1e10)
    
    t0 = time.time()
    prime = Prime(n)
    prime.search()
    t1 = time.time() - t0
    print("n = %s, found %s primes in  %s seconds"%(prime.n, prime.p.size, t1))
