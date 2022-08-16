You can compile it with the command:

g++ zuc.cpp zuc.h BIT.h main.cpp -pthread -O3 -std=c++11 -march=native -fopenmp

================================================================================
You can change the configuration in main.cpp.

---If you want to test the bias in the keystream, please define _TEST_KEYSTREAM_.

---If you want to test the bias in the LFSR, please define _TEST_INTERNAL_STATE_.

---They both cannot be defined at the same time.

---You can modify the value of threadNum according to your machines.

===============================================================================
You are required to input something.

E.g., if you want to test the 30-round attack on zuc-v2 (the bias is 2^{-19.2}),

please first input "2".

Then, please input "x 15 0" if the key varies in each sample, where x is an integer you input.

Please make sure that 2^x * threadNum > 2^42 to detect a valid biased linear relation for 30-round zuc-v2 

================================================================================
In my experiment, due to the limits of the machine, threadNum = 110 is used.

Then, I input "2".

Next, I input "35 15 0".
