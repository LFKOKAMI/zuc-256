#include "zuc.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <thread>
#include <random>
#include <ctime>
#include <cstdio>
#include <string>
#include <fstream>
#define _TEST_KEYSTREAM_
//#define _TEST_INTERNAL_STATE_
using namespace std;

const u64 threadNum = 4;

/*
---If you want to test the bias in the keystream, please define _TEST_KEYSTREAM_.
---If you want to test the bias in the LFSR, please define _TEST_INTERNAL_STATE_.
---The both cannot be defined at the same time.

---You can modify the value of threadNum according to your machines.
*/

int testSh = 0, testMask = 0;
u32 targetMask = 0;

//scale: 2^size * 2^size
void WHT(L64* arr, L64 len,u64 &bias,bool &sum,u64 &mask, bool isTarget) {
	L64 h = 1;
	L64 x = 0, y = 0;
	while (h < len) {
		for (u64 i = 0; i < len; ) {
			for (u64 j = i; j < i+h; j++) {
				x = arr[j];
				y = arr[j + h];
				arr[j] = x + y;
				arr[j + h] = x - y;
			}
			i = i + 2 * h;
		}
		h = 2 * h;
	}

	//find the best biased linear relation
	L64 zero = 0;
	bool sign = 0;
	bias = 0;
	mask = 0;
	for (u64 i = 1; i < len; i++) {
		sign = 0;
		if (arr[i] < 0) {
			arr[i] = zero - arr[i];
			sign = 1;//there are more 1 (-1)^1=-1 (there are many -1).
		}

		if (arr[i] == arr[0]) {//sometimes we use a 16-bit table to store a 15-bit value
			continue;
		}

		arr[i] = arr[i] / 2;
		//cout <<i<< " : " <<hex << arr[i] << endl;

		if (bias < arr[i]) {
			bias = arr[i];
			mask = i;
			sum = sign;
		}

		//output the bias of the tested mask
		if (isTarget && i==testMask) {
			//cout << "Test the bias for the masks in the paper";
			//cout<<" --, mask = " << hex << targetMask;
			//cout<< ", sum = " << sign;
			//cout << ", log2(bias) = -" << log2((arr[0] / arr[i])) << endl;
		}
	}
}

void findBestLinearMaskThreadWHT(L64*** table,u64 testTimes,int tabNum,u32 tabSize) {
	L64** finalTab;
	finalTab = new L64 * [tabNum];
	for (int i = 0; i < tabNum; i++) {
		finalTab[i] = new L64[tabSize];
	}
	int end = tabSize;
	u64 maxBias = 0, maxMask = 0, bias = 0, mask = 0;
	bool sum = 0, tSum = 0;;
	int maxSh = 0;
	bool isTarget = false;
	for (int sh = 0; sh < tabNum; sh++) {
		for (int j = 0; j < end; j++) {
			finalTab[sh][j] = 0;
			for (int k = 0; k < threadNum; k++) {
				finalTab[sh][j] = finalTab[sh][j] + table[k][sh][j];
			}
			//cout << hex << j << ": " <<dec<< finalTab[sh][j] << endl;
		}
		//cout << endl;
		if (sh == testSh) {
			isTarget = true;
		}
		else {
			isTarget = false;
		}

		WHT(finalTab[sh], tabSize, bias, tSum, mask,isTarget);
		if (maxBias < bias) {
			maxBias = bias;
			sum = tSum;
			maxMask = mask;
			maxSh = sh;
		}
	}
	
	cout << "The best one among the samples -- ";

#ifdef _TEST_KEYSTREAM_
	if (tabSize == 0x100) {
		cout << "mask = " << hex << (maxMask << (8 * maxSh)) << " ";
		cout << ", sum = " << sum << " ";
		cout << ", log2(bias) = -" << log2((testTimes * threadNum / maxBias)) << endl;
	}
	else if (tabSize == 0x10000) {
		if (maxSh < 4) {
			cout << "mask = " << hex << (maxMask << (8 * maxSh)) << " ";
			cout << ", sum = " << sum << " ";
			cout << ", log2(bias) = -" << log2((testTimes * threadNum / maxBias)) << endl;
		}
		else {
			u32 finalMask = (maxMask & 0xff) << 24;
			finalMask = finalMask | ((maxMask>>8) & 0xff);
			cout << "mask = " << hex << finalMask << " ";
			cout << ", sum = " << sum << " ";
			cout << ", log2(bias) = -" << log2((testTimes * threadNum / maxBias)) << endl;
		}
	}
#endif // _TEST_KEYSTREAM_

#ifdef _TEST_INTERNAL_STATE_
	if (tabSize == 0x100) {
		int shift = 0;
		if (maxSh > 0 ) {
			shift = 7 + 8 * (maxSh - 1);
		}
		cout << "mask = " << hex << (maxMask << shift) << " ";
		cout << ", sum = " << sum << " ";
		cout << ", log2(bias) = -" << log2((testTimes * threadNum / maxBias)) << endl;
	}
	else if (tabSize == 0x10000) {
		int shift = 0;
		if (maxSh > 0) {
			shift = 7 + 8 * (maxSh - 1);
		}
		if (maxSh < 4) {
			cout << "mask = " << hex << (maxMask << shift) << " ";
			cout << ", sum = " << sum << " ";
			cout << ", log(bias) = -" << log2((testTimes * threadNum / maxBias)) << endl;
		}
		else {
			u32 finalMask = (maxMask & 0xff) << 23;
			finalMask = finalMask | ((maxMask >> 8) & 0x7f);
			cout << "mask = " << hex << finalMask << " ";
			cout << ", sum = " << sum << " ";
			cout << ", log(bias) = -" << log2((testTimes * threadNum / maxBias)) << endl;
		}
	}
#endif

	//cout << "sh:" << maxSh << " mask:" << hex << maxMask;
	//cout << " sum:" << sum << " rate:" << (testTimes * threadNum / maxBias) << endl;

	for (int i = 0; i < tabNum; i++) {
		delete[]finalTab[i];
	}
	delete[]finalTab;
}

void assignIV(u32 v[], int start, u64 value) {
	for (int i = 0; i < 8; i++) {
		v[start + i] = (value >> (8 * i)) & 0xff;
	}
}

void computeDiffAndMasks(string signedDiff, int size, u32& diff, u32& orMask, u32& andMask) {
	int len = size-1;
	diff = 0;

	andMask = 0x7fffffff;
	orMask = 0x0;

	u32 mask = 0;
	for (int i = 0; i < size; i++) {
		if (signedDiff[i] == 'n') {
			//cout << len - i << ": n" << endl;
			diff = (diff + (1 << (len - i))) % p;
			mask = 1 << (len - i);
			mask = ~mask;
			mask = mask & 0x7fffffff;
			andMask = andMask & mask;
		}
		else if (signedDiff[i] == 'u') {
			//cout << len - i << ": u" << endl;
			diff = (p - (1 << (len - i)) + diff) % p;
			orMask = orMask | (1 << (len - i));
		}
	}
}

void IVKeyBiasAttackWithWHT(u64 testTimes, int tNumber, int isKeyFixed, L64** table, u32 k[], int rounds,int tabNum,u32 tabSize) {
	u32 diff[16], andMask[16], orMask[16];
	for (int i = 0; i < 16; i++) {
		diff[i] = 0;
		andMask[i] = 0x7fffffff;
		orMask[i] = 0x0;
	}

	string signedDiff[11] = {
		"000nn0nn0000000nn0nn0nn00000n0n",
		"0000u0000000000nnnnn0nn000000n0",
		"0000n00n0000000uu00u000nn00n0nn",
		"0000000u000000000nn0n00nn000n00",
		"0n0000n00000000uuuuuuuu00n000u0",
		"nu00000000nuuuuuu0000000u000000",
		"00n000000000000nuu0000000u00000",
		"000000000000nu00000000000000000",
		"000unnn00000n0n000unnn000000000",
		"000000000000000000000uunnnnnn00",
		"0000000000000000000000u0unnnnn0",
	};

	u32 R2[3] = { 0xc99de9d6, 0xb7b8cf96,0xfaf5498c };

	for (int i = 0; i < 11; i++) {
		computeDiffAndMasks(signedDiff[i], 31, diff[i], orMask[i], andMask[i]);
	}

	//u64 table[tabNum][tabSize];//distribution
	for (int i = 0; i < tabNum; i++) {
		for (int j = 0; j < tabSize; j++) {
			table[i][j] = 0;
		}
	}

	ZUC zuc;

	u32 output = 0, outputNew = 0, mod = 0;
	u32 s[16], st[16], stt[16];
	u32 r[2], rt[2], rtt[16];
	u32 v[25];
	u64 one = 1;

	std::mt19937_64 mt;            // メルセンヌ・ツイスタの32ビット版
	std::random_device rnd;     // 非決定的な乱数生成器
	mt.seed(rnd() + tNumber);

	u32 key[32];
	for (int i = 0; i < 32; i++) {
		key[i] = k[i];
	}
	zuc.correctKey7(key, R2);
	u64 randValue = 0;

	u32 c[2] = { 0,0 };
	u64 cnt = 0;
	u64 unit = one << 24;
	u64 cnt3 = 0;
	for (u64 i = 0; i < testTimes; i++) {
		//generate random iv
		assignIV(v, 0, mt());
		assignIV(v, 8, mt());
		assignIV(v, 16, mt());
		v[24] = mt() & 0x3f;
		for (int i = 17; i < 25; i++) {
			v[i] = v[i] & 0x3f;
		}
		//s8L:000u nnn0 0000 0000=v[3]||v[11]
		v[3] = v[3] & 0xf1;
		v[3] = v[3] | 0x10;

		//assign key
		if (isKeyFixed == 1) {
			for (int j = 0; j < 4; j++) {
				randValue = mt();
				for (int y = 0; y < 8; y++) {
					key[8 * j + y] = (randValue >> (8 * y)) & 0xff;
				}
			}
			zuc.correctKey7(key, R2);
		}

		//correct iv and key
		if (zuc.fulfillIVAutoBestAttack(R2, key, v) == 0) {
			cout << "wrong" << endl;
		}
		//load key and iv
		zuc.load(key, v, s);
		//mask 
		for (int i = 0; i < 11; i++) {
			s[i] = s[i] & andMask[i];
			s[i] = s[i] | orMask[i];
		}
		r[0] = 0;
		r[1] = 0;

#ifdef _TEST_KEYSTREAM_
		zuc.initialization(s, r, rounds - 2);
		st[15] = zuc.workMode();
#endif //_TEST_KEYSTREAM_

#ifdef _TEST_INTERNAL_STATE_
		zuc.initialization(s, r, rounds);
		st[15] = zuc.getState(15);
#endif // _TEST_INTERNAL_STATE_

		//add difference
		for (int i = 0; i < 11; i++) {
			s[i] = (s[i] + diff[i]) % p;
		}
		r[0] = 0;
		r[1] = 0;

#ifdef _TEST_KEYSTREAM_
		zuc.initialization(s, r, rounds - 2);
		stt[15] = zuc.workMode();
		mod = st[15] ^ stt[15];
		if (tabSize == 0x100) {
			table[0][mod & 0xff]++;
			table[1][(mod >> 8) & 0xff]++;
			table[2][(mod >> 16) & 0xff]++;
			table[3][(mod >> 24) & 0xff]++;
		}
		else if (tabSize == 0x10000) {
			table[0][mod & 0xffff]++;
			table[1][(mod >> 8) & 0xffff]++;
			table[2][(mod >> 16) & 0xffff]++;
			u32 index = (mod >> 24) & 0xff;
			index = ((mod & 0xff) << 8) | index;//index=[7:0]||[31:24]
			table[3][index]++;
		}
#endif //_TEST_KEYSTREAM_

#ifdef _TEST_INTERNAL_STATE_
		zuc.initialization(s, r, rounds);
		stt[15] = zuc.getState(15);
		mod = (p - st[15] + stt[15]) % p;
		if (tabSize == 0x100) {
			table[0][mod & 0x7f]++;
			table[1][(mod >> 7) & 0xff]++;
			table[2][(mod >> 15) & 0xff]++;
			table[3][(mod >> 23) & 0xff]++;
		}
		else if (tabSize == 0x10000) {
			table[0][mod & 0x7fff]++;
			table[1][(mod >> 7) & 0xffff]++;
			table[2][(mod >> 15) & 0xffff]++;
			u32 index = (mod >> 23) & 0xff;
			index = ((mod & 0x7f) << 8) | index;//index=[6:0]||[30:23]
			table[3][index]++;
		}
#endif // _TEST_INTERNAL_STATE_
	}
}

void IVKeyBiasAttackNewSchemeWithWHT(u64 testTimes, int tNumber, int isKeyFixed, L64** table, u32 k[], int rounds, int tabNum, u32 tabSize) {
	u32 diff[16], andMask[16], orMask[16];
	for (int i = 0; i < 16; i++) {
		diff[i] = 0;
		andMask[i] = 0x7fffffff;
		orMask[i] = 0x0;
	}

	string signedDiff[11] = {
		"000000nn0000000u00000nn00000u0n",
		"0000n00u0000000uu0u000u0n00n00n",
		"0n000u000000000nnnn00nn00000n0n",
		"00n00nn000000000000u0000n0nn0n0",
		"0u0nn000000000000n00000000n0000",
		"000000nu0000000uuuuuuuu0000000u",
		"000000000000000nnnnnnn0000000n0",
		"0000000n00000000000000000000000",
		"nnnnnn0n0000000nn0nuuu0uu0000uu",
		"000000000000000000nuuuuuuuuu000",
		"00000000000000000ununnnnnnn0000",
	};

	u32 R2[3] = { 0xa21c991b, 0xcf1106f0,0x32f0e1e3 };

	u32 key[32];
	for (int i = 0; i < 32; i++) {
		key[i] = k[i];
	}


	for (int i = 0; i < 11; i++) {
		computeDiffAndMasks(signedDiff[i], 31, diff[i], orMask[i], andMask[i]);
	}

	ZUC zuc;

	//u64 table[tabNum][tabSize];//distribution
	for (int i = 0; i < tabNum; i++) {
		for (int j = 0; j < tabSize; j++) {
			table[i][j] = 0;
		}
	}

	u32 output = 0, outputNew = 0, mod = 0;
	u32 s[16], st[16], stt[16];
	u32 r[2], rt[2], rtt[16];
	u32 v[16];
	u64 one = 1;

	std::mt19937_64 mt;            // メルセンヌ・ツイスタの32ビット版
	std::random_device rnd;     // 非決定的な乱数生成器
	mt.seed(rnd() + tNumber);
	u64 randValue = 0;

	u32 c[2] = { 0,0 };
	u64 cnt = 0;
	u64 unit = one << 24;
	u64 cnt3 = 0;
	for (u64 i = 0; i < testTimes; i++) {
		//assign iv
		for (int j = 0; j < 2; j++) {
			randValue = mt();
			for (int y = 0; y < 8; y++) {
				v[8 * j + y] = (randValue >> (8 * y)) & 0xff;
			}
		}
		//assign key
		if (isKeyFixed == 1) {
			for (int j = 0; j < 4; j++) {
				randValue = mt();
				for (int y = 0; y < 8; y++) {
					key[8 * j + y] = (randValue >> (8 * y)) & 0xff;
				}
			}
		}

		//S8L: nn0n uuu0 uu00 00uu=v[1]||v[9]
		v[9] = v[9] | 0xc3;
		v[1] = v[1] | 0x0e;
		v[1] = v[1] & 0x2f;
		zuc.correctV0(v, R2);

		//correct iv and key
		if (zuc.fulfillIVAutoBestAttackNew(R2, key, v) == 0) {
			cout << "wrong" << endl;
		}

		//load key and iv
		zuc.loadNewScheme(key, v, s);
		//mask 
		for (int i = 0; i < 11; i++) {
			s[i] = s[i] & andMask[i];
			s[i] = s[i] | orMask[i];
		}
		r[0] = 0;
		r[1] = 0;

#ifdef _TEST_KEYSTREAM_
		zuc.initialization(s, r, rounds - 2);
		st[15] = zuc.workMode();
#endif //_TEST_KEYSTREAM_

#ifdef _TEST_INTERNAL_STATE_
		zuc.initialization(s, r, rounds);
		st[15] = zuc.getState(15);
#endif // _TEST_INTERNAL_STATE_

		//add difference
		for (int i = 0; i < 11; i++) {
			s[i] = (s[i] + diff[i]) % p;
		}
		r[0] = 0;
		r[1] = 0;

#ifdef _TEST_KEYSTREAM_
		zuc.initialization(s, r, rounds - 2);
		stt[15] = zuc.workMode();
		mod = st[15] ^ stt[15];
		if (tabSize == 0x100) {
			table[0][mod & 0xff]++;
			table[1][(mod >> 8) & 0xff]++;
			table[2][(mod >> 16) & 0xff]++;
			table[3][(mod >> 24) & 0xff]++;
		}
		else if (tabSize == 0x10000) {
			table[0][mod & 0xffff]++;
			table[1][(mod >> 8) & 0xffff]++;
			table[2][(mod >> 16) & 0xffff]++;
			u32 index = (mod >> 24) & 0xff;
			index = ((mod & 0xff) << 8) | index;//index=[7:0]||[31:24]
			table[3][index]++;
		}
#endif //_TEST_KEYSTREAM_

#ifdef _TEST_INTERNAL_STATE_
		zuc.initialization(s, r, rounds);
		stt[15] = zuc.getState(15);
		mod = (p - st[15] + stt[15]) % p;
		if (tabSize == 0x100) {
			table[0][mod & 0x7f]++;
			table[1][(mod >> 7) & 0xff]++;
			table[2][(mod >> 15) & 0xff]++;
			table[3][(mod >> 23) & 0xff]++;
		}
		else if (tabSize == 0x10000) {
			table[0][mod & 0x7fff]++;
			table[1][(mod >> 7) & 0xffff]++;
			table[2][(mod >> 15) & 0xffff]++;
			u32 index = (mod >> 23) & 0xff;
			index = ((mod & 0x7f) << 8) | index;//index=[6:0]||[30:23]
			table[3][index]++;
		}
#endif // _TEST_INTERNAL_STATE_
	}
}

void findBiasWHT(int shift, int rounds, int isKeyFixed, int choice) {
	L64*** tabTT;
	tabTT = new L64 * *[threadNum];
	const u32 tabSize = 1 << 16;//using 16-bit masks
	const int tabNum = 4;

	for (int i = 0; i < threadNum; i++) {
		tabTT[i] = new L64 * [tabNum];
		for (int j = 0; j < tabNum; j++) {
			tabTT[i][j] = new L64[tabSize];
			for (int u = 0; u < tabSize; u++) {
				tabTT[i][j][u] = 0;
			}
		}
	}

	u32 key[32];
	//generate random key
	srand(time(NULL));
	for (int i = 0; i < 32; i++) {
		key[i] = rand() & 0xff;
	}

	vector<thread> myThreads;
	u64 one = 1;
	u64 testTimes = one << shift;

	cout << "The total number of samples: 0x" << hex << (testTimes * threadNum) << endl;

	for (int j = 0; j < threadNum; j++) {
		if (choice == 1) {
			myThreads.push_back(thread(IVKeyBiasAttackWithWHT, testTimes, j, isKeyFixed, tabTT[j], key, rounds, tabNum, tabSize));
		}
		else if (choice == 2) {
			myThreads.push_back(thread(IVKeyBiasAttackNewSchemeWithWHT, testTimes, j, isKeyFixed, tabTT[j], key, rounds, tabNum, tabSize));
		}
	}
	for (thread& single : myThreads) {
		single.join();
	}
	findBestLinearMaskThreadWHT(tabTT, testTimes, tabNum, tabSize);

	//delete
	for (int i = 0; i < threadNum; i++) {
		for (int j = 0; j < tabNum; j++) {
			delete[]tabTT[i][j];
		}
		delete[]tabTT[i];
	}
	delete[]tabTT;
}

int main() {
	int choice = 31;
	cout << "zuc or zuc-v2, please input (1/2):";
	cin >> choice;
	int isKeyStream = 0;

#ifdef _TEST_INTERNAL_STATE_
	cout << "Test the bias in LFSR!" << endl;
	isKeyStream = 0;
#endif 

#ifdef _TEST_KEYSTREAM_
	cout << "Test the bias in the keystream word!" << endl;
	isKeyStream = 1;
#endif
	
	cout << hex << "There are in total " << threadNum << " threads!" << endl;

	//The detected best masks
	if (isKeyStream == 1 && choice==1) {
		testSh = 0;
		testMask = 0x80;
		targetMask = 0x80;
	}
	else if (isKeyStream == 1 && choice == 2) {
		testSh = 2;
		testMask = 0x4000;
		targetMask = 0x40000000;
	}
	else if (isKeyStream == 0 && choice == 1) {
		testSh = 0;
		testMask = 0x40;
		targetMask = 0x40;
	}
	else if (isKeyStream == 0 && choice == 2) {
		testSh = 2;
		testMask = 0x4000;
		targetMask = 0x20000000;
	}

	//output tips:
	cout << endl << "Tips:" << endl << endl;
	cout << "---If shift = x, there will be 2^x samples in each thread." << endl;
	cout << "---Therefore, there are in total threadNum * 2^x samples." << endl << endl;

	cout << "---If you want to attack r rounds on LFSR, please set rounds = r - 15" << endl;
	cout << "-E.g. if you want to test 31-round attack, set rounds = 16" << endl;
	cout << "---If you want to attack r rounds on the keystream, set rounds  = r" << endl;
	cout << "-E.g. if you want to test 15-round attack, set rounds = 15" << endl << endl;

	cout << "---If isKeyFixed = 1, key will be the same in all threads and in each sample." << endl;
	cout << "---If isKeyFixed = 0, key will vary in all samples in all threads." << endl << endl;


	int shift = 0, rounds = 15, isKeyFixed = 1;
	cout << "input shift, rounds, isKeyFixed(0/1):";
	cin >> shift >> rounds >> isKeyFixed;

	if (isKeyFixed == 1) {
		cout << "key is fixed for different IVs" << endl;
	}
	else {
		cout << "key also varies in each experiment" << endl;
	}

	for (int i = 0; i < 10; i++) {
		cout << endl << "times : " << dec << i << " -- ";
		findBiasWHT(shift, rounds, isKeyFixed, choice);
	}
}