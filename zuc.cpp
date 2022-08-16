#include "zuc.h"
#include "BIT.h"
#include <iostream>
#include <vector>
#include <ctime>
#include <algorithm>
#include <string>
#define _NEW_SCHEME_
//#define _OLD_SCHEME_

using namespace std;

ZUC::ZUC() {
	for (int i = 0; i < 256; i++) {
		S0[i] = SS0[i];
		S1[i] = SS1[i];
	}
	computeInverseSBox(S0INVER, S1INVER);
	//initialization();
}

void ZUC::initialization() {
	valTable0.clear();
	valTable1.clear();
	valTable0.resize(256 * 256);
	valTable1.resize(256 * 256);
	computeDDT(0, DDT0,valTable0);
	computeDDT(1, DDT1,valTable1);
}

u32 ZUC::inverseL2(u32 t) {
	u32 res = 0;
	//L2^{-1}(x) = x^30 + x^28 + x^24 + x^20 + x^18 + x^16 + x^14 + x^10 + x^8 + x^2 + 1
	res = t ^ rot(t, 30) ^ rot(t, 28)^ rot(t, 24)^ rot(t, 20)^ rot(t, 18)^ rot(t, 16);
	res = res ^ rot(t, 14) ^ rot(t, 10) ^ rot(t, 8) ^ rot(t, 2);

	return res;
}

u32 ZUC::L1(u32 t) {
	return (t ^ rot(t, 2) ^ rot(t, 10) ^ rot(t, 18) ^ rot(t, 24));
}

u32 ZUC::L2(u32 t) {
	return (t ^ rot(t, 8) ^ rot(t, 14) ^ rot(t, 22) ^ rot(t, 30));
}

//s-box
u32 ZUC::SFun(u32 t) {
	u32 u[4];
	u[0] = S1[t & 0xff];
	u[1] = S0[(t >> 8) & 0xff];
	u[2] = S1[(t >> 16) & 0xff];
	u[3] = S0[(t >> 24) & 0xff];
	return ((u[3] << 24) | (u[2] << 16) | (u[1] << 8) | u[0]);
}

//inverse of the Sbox
void ZUC::computeInverseSBox(u32 S0Inv[], u32 S1Inv[]) {
	for (int i = 0; i < 256; i++) {
		S0Inv[S0[i]] = i;
		S1Inv[S1[i]] = i;
	}
}

u32 ZUC::SFunInverse(u32 t,u32 S0Inv[],u32 S1Inv[]) {
	u32 u[4];
	u[0] = S1Inv[t & 0xff];
	u[1] = S0Inv[(t >> 8) & 0xff];
	u[2] = S1Inv[(t >> 16) & 0xff];
	u[3] = S0Inv[(t >> 24) & 0xff];
	return ((u[3] << 24) | (u[2] << 16) | (u[1] << 8) | u[0]);
}

inline bool ZUC::trans(u32 DDT[][256], u32 input, u32 output) {
	if (DDT[input][output] > 0) {
		return 1;
	}
	else {
		return 0;
	}
}

u32 ZUC::linear80(u32 t) {
	u32 u = left(t, 8);
	return (u + t) % p;
}

//(2^8+1)
u32 ZUC::inverse80(u32 t) {
	u32 result = 0;
	for (int j = 0; j < 31; j++) {//compute inverse8 multi target
		if (((INVER8 >> j) & 0x1) == 1) {
			result = (result + left(t, j)) % p;
		}
	}
	return result;
}

u32 ZUC::inverse20(u32 t) {
	return left(t, 11);
}

void ZUC::bitOrg() {
	//x[0]=s[15]H || s[14]L
	//x[1]=s[11]L || s[9]H
	//x[2]=s[7]L || s[5]H
	x[0] = (high(s[15]) << 16) | (low(s[14]));
	x[1] = (low(s[11]) << 16) | (high(s[9]));
	x[2] = (low(s[7]) << 16) | (high(s[5]));
}

void ZUC::lsfr() {
	u32 y = 0;
	y = x[0] ^ r[0];
	y = y + r[1];
	y = y >> 1;
	y = (y + s[0]) % p;
	y = (y + left(s[0], 8)) % p;
	y = (y + left(s[4], 20)) % p;
	y = (y + left(s[10], 21)) % p;
	y = (y + left(s[13], 17)) % p;
	y = (y + left(s[15], 15)) % p;
	if (y == 0) {
		y = p;
	}
	//shift
	for (int i = 0; i < 15; i++) {
		s[i] = s[i + 1];
	}
	//update
	s[15] = y;
}

void ZUC::regUpdate() {
	r[0] = r[0] + x[1];
	r[1] = r[1] ^ x[2];
	//shift
	u32 tr0, tr1;
	tr0 = (r[0] << 16) | (r[1] >> 16);
	tr1 = (r[1] << 16) | (r[0] >> 16);
	//shift
	r[0] = (tr0 ^ rot(tr0, 2) ^ rot(tr0, 10) ^ rot(tr0, 18) ^ rot(tr0, 24));
	r[1] = (tr1 ^ rot(tr1, 8) ^ rot(tr1, 14) ^ rot(tr1, 22) ^ rot(tr1, 30));
	//sbox
	u32 u0[4],u1[4];
	for (int i = 0; i < 4; i++) {
		u0[i] = (r[0] >> (8 * i)) & 0xff;
		u1[i] = (r[1] >> (8 * i)) & 0xff;
	}
	u0[3] = S0[u0[3]];
	u0[2] = S1[u0[2]];
	u0[1] = S0[u0[1]];
	u0[0] = S1[u0[0]];

	u1[3] = S0[u1[3]];
	u1[2] = S1[u1[2]];
	u1[1] = S0[u1[1]];
	u1[0] = S1[u1[0]];

	r[0] = (u0[3] << 24) | (u0[2] << 16) | (u0[1] << 8) | u0[0];
	r[1] = (u1[3] << 24) | (u1[2] << 16) | (u1[1] << 8) | u1[0];
}

void ZUC::initialization(u32 st[], u32 rt[], int rounds) {
	for (int i = 0; i < 16; i++) {
		s[i] = st[i];
	}
	r[0] = rt[0];
	r[1] = rt[1];

	for (int i = 0; i < rounds; i++) {
		bitOrg();
		lsfr();
		regUpdate();
	}
}

u32 ZUC::workMode() {
	bitOrg();
	regUpdate();
	//no feedback
	u32 y = 0;
	y = (y + s[0]) % p;
	y = (y + left(s[0], 8)) % p;
	y = (y + left(s[4], 20)) % p;
	y = (y + left(s[10], 21)) % p;
	y = (y + left(s[13], 17)) % p;
	y = (y + left(s[15], 15)) % p;
	if (y == 0) {
		y = p;
	}
	//shift
	for (int i = 0; i < 15; i++) {
		s[i] = s[i + 1];
	}
	//update
	s[15] = y;

	//keystream (no x[3])
	bitOrg();
	x[3]= (low(s[2]) << 16) | (high(s[0]));
	u32 stream = x[0] ^ r[0];
	stream = stream + r[1];
	stream = stream ^ x[3];
	return stream;
}

void ZUC::load(u32 key[], u32 v[]) {
	s[0] = LOAD(key[0], d[0], key[21], key[16]);
	s[1] = LOAD(key[1], d[1], key[22], key[17]);
	s[2] = LOAD(key[2], d[2], key[23], key[18]);
	s[3] = LOAD(key[3], d[3], key[24], key[19]);
	s[4] = LOAD(key[4], d[4], key[25], key[20]);
	u32 t = d[5]| v[17];
	s[5] = LOAD(v[0], t, key[5], key[26]);

	t = d[6]| v[18];
	s[6] = LOAD(v[1], t, key[6], key[27]);

	t = d[7]| v[19];
	s[7] = LOAD(v[10], t, key[7], v[2]);

	t = d[8] | v[20];
	s[8] = LOAD(key[8], t, v[3], v[11]);

	t = d[9] | v[21];
	s[9] = LOAD(key[9], t, v[12], v[4]);

	t = d[10] | v[22];
	s[10] = LOAD(v[5], t, key[10], key[28]);

	t = d[11] | v[23];
	s[11] = LOAD(key[11], t, v[6], v[13]);

	t = d[12] | v[24];
	s[12] = LOAD(key[12], t, v[7], v[14]);

	s[13] = LOAD(key[13], d[13], v[15], v[8]);

	t = (key[31] >> 4) & 0xf;
	t = d[14] | t;
	s[14] = LOAD(key[14], t, v[16], v[9]);

	t = (key[31]) & 0xf;
	t = d[15] | t;
	s[15] = LOAD(key[15], t, key[30], key[29]);
}

void ZUC::load(u32 key[], u32 v[], u32 st[]) {
	st[0] = LOAD(key[0], d[0], key[21], key[16]);
	st[1] = LOAD(key[1], d[1], key[22], key[17]);
	st[2] = LOAD(key[2], d[2], key[23], key[18]);
	st[3] = LOAD(key[3], d[3], key[24], key[19]);
	st[4] = LOAD(key[4], d[4], key[25], key[20]);
	u32 t = d[5] | v[17];
	st[5] = LOAD(v[0], t, key[5], key[26]);

	t = d[6] | v[18];
	st[6] = LOAD(v[1], t, key[6], key[27]);

	t = d[7] | v[19];
	st[7] = LOAD(v[10], t, key[7], v[2]);

	t = d[8] | v[20];
	st[8] = LOAD(key[8], t, v[3], v[11]);

	t = d[9] | v[21];
	st[9] = LOAD(key[9], t, v[12], v[4]);

	t = d[10] | v[22];
	st[10] = LOAD(v[5], t, key[10], key[28]);

	t = d[11] | v[23];
	st[11] = LOAD(key[11], t, v[6], v[13]);

	t = d[12] | v[24];
	st[12] = LOAD(key[12], t, v[7], v[14]);

	st[13] = LOAD(key[13], d[13], v[15], v[8]);

	t = (key[31] >> 4) & 0xf;
	t = d[14] | t;
	st[14] = LOAD(key[14], t, v[16], v[9]);

	t = (key[31]) & 0xf;
	t = d[15] | t;
	st[15] = LOAD(key[15], t, key[30], key[29]);
}

void ZUC::loadNewScheme(u32 key[], u32 v[], u32 st[]) {
	for (int i = 0; i < 7; i++) {
		st[i] = LOAD(key[i], D[i], key[16 + i], key[24 + i]);
	}
	for (int i = 7; i < 15; i++) {
		st[i] = LOAD(key[i], D[i], v[i-7], v[1 + i]);
	}
	st[15] = LOAD(key[15], D[15], key[23], key[31]);
}

void ZUC::getState(u32 st[],u32 rt[]) {
	for (int i = 0; i < 16; i++) {
		st[i] = s[i];
	}
	rt[0] = r[0];
	rt[1] = r[1];
}

u32 ZUC::getState(int i) {
	return s[i];
}

void ZUC::printState() {
	cout << "s: ";
	for (int i = 15; i >= 0; i--) {
		cout << hex << s[i] << " ";
	}
	cout << " || r: ";
	for (int i = 1; i >= 0; i--) {
		cout << hex << r[i] << " ";
	}
	cout << endl;
}

//compute DDT
void ZUC::computeDDT(int num, u32 DDT[][256], vector<vector<u32> >& valTable) {
	u32 ud, u, t, td;
	for (int i = 0; i < 0x100; i++) {
		for (int j = 0; j < 0x100; j++) {
			DDT[i][j] = 0;
			valTable[(j << 8) | i].clear();
			//compute DDT[i][j]
			for (int t = 0; t < 0x100; t++) {
				if (num == 0) {
					u = S0[t];
					ud = S0[t ^ i];
				}
				else {
					u = S1[t];
					ud = S1[t ^ i];
				}
				if ((u ^ ud) == j) {
					DDT[i][j]++;
					valTable[(j << 8) | i].push_back(t);
				}
			}
		}
	}
	/*for (int i = 0; i < 0x100; i++) {
		for (int j = 0; j < 0x100; j++) {
			if (DDT[i][j] > 0) {
				cout << hex << ((j << 8) | i) << ": " << DDT[i][j] << " , " << valTable[(j << 8) | i].size() << endl;
			}
		}
		cout << endl;
	}*/
}

//the 31-round attack
bool ZUC::fulfillIVAutoBestAttack(u32 R2[], u32 key[], u32 v[]) {//correct IV and key
	//R2[0] = S\circ L2 (S5H || S11L)
	//S5H = v[0] || c[0] || v[17] || k[5][7] (8-1-6-1) \\\\ c[0]<<<6 = d[5] (7-bit)
	//S11L = v[6] || v[13] (8-8)
	u32 inv = SFunInverse(R2[0], S0INVER, S1INVER);
	inv = inverseL2(inv);
	v[13] = inv & 0xff;
	v[6] = (inv >> 8) & 0xff;
	if (BIT(inv, 16) == 1) {
		key[5] = key[5] | 0x80;
	}
	else {
		key[5] = key[5] & 0x7f;
	}
	v[17] = (inv >> 17) & 0x3f;//a 6-bit value
	v[0] = (inv >> 24) & 0xff;
	//cout << hex << ((inv >> 16) & 0xffff) << endl;
	//cout << "inv:" << hex << inv << endl;

	//R2[1] = S\circ L2 ( (R2[0]^X[2])L || (R1[0]+X[1])H )
	u32 U = SFunInverse(R2[1], S0INVER, S1INVER);
	U = inverseL2(U);
	//focus on (R2[0]^X[2])L = UH = R2[0]L ^ S6H
	//S6H = v[1] || c[1] || v[18] || k[6][7] \\\\ c[1]<<<6 = d[6] (7 bit)
	u32 uh = (U >> 16) & 0xffff;
	uh = (uh ^ R2[0]) & 0xffff;
	v[1] = (uh >> 8) & 0xff;
	v[18] = (uh >> 1) & 0x3f;
	if (BIT(uh, 0) == 1) {
		key[6] = key[6] | 0x80;
	}
	else {
		key[6] = key[6] & 0x7f;
	}
	//fulfill UL
	u32 S9H = (key[9] << 8) | (1 << 7) | (v[21] << 1) | (v[12] >> 7);
	u32 S7L = (key[7] << 8) | (v[2]);
	u32 R1 = L1((S9H << 16) | S7L);
	R1 = SFun(R1);
	u32 S10H = (v[5] << 8) | (1 << 7) | (v[22] << 1) | (key[10] >> 7);
	u32 low = (R1 + S10H) & 0xffff;
	u32 full = ((U & 0xffff) << 16) | low;//UL || low
	u32 S12L = full - R1;
	S12L = (S12L >> 16)&0xffff;
	v[7] = S12L >> 8;
	v[14] = S12L & 0xff;

	//third round
	u32 add = (S12L << 16) | (S10H);
	u32 R11 = R1 + add;//the output of the first register after adding X2 in the second round
	u32 R11H = (R11) & 0xffff;//shift to the 16 MSB
	u32 R11L = (v[3] << 8) | v[11];
	R11L=R11L ^ ((R2[0] >> 16) & 0xffff);//S8L^R2[0]H
	u32 R12 = (R11H << 16) | R11L;
	R12 = L1(R12);
	R12 = SFun(R12);
	//compute S7H:IV10 || d7 || IV19 || K7[7]
	u32 IN = SFunInverse(R2[2], S0INVER, S1INVER);
	IN = inverseL2(IN);
	u32 HIGH = (IN >> 16) & 0xffff;
	HIGH = (HIGH ^ R2[1]) & 0xffff;
	v[10] = (HIGH >> 8) & 0xff;
	v[19] = (HIGH >> 1) & 0x3f;
	/*if (BIT(HIGH, 0) == 1) {
		key[7] = key[7] | 0x80;
	}
	else {
		key[7] = key[7] & 0x7f;
	}*/
	if (BIT(HIGH, 0) != BIT(key[7], 7)) {
		return false;
	}
	//consider LOW of IN
	u32 LOW = IN & 0xffff;
	//LOW = ( R12+(S13L,S11H) )H
	u32 S11H = (key[11] << 8) | (1 << 7) | (v[23] << 1) | (v[6] >> 7);
	low = (R12 + S11H) & 0xffff;
	full = (LOW << 16) | low;//UL || low
	u32 S13L = (full - R12);
	S13L = S13L >> 16;
	v[8] = S13L & 0xff;
	v[15] = (S13L >> 8) & 0xff;
	return true;
}

void ZUC::correctKey7(u32 key[], u32 R2[]) {
	u32 IN = SFunInverse(R2[2], S0INVER, S1INVER);
	IN = inverseL2(IN);
	u32 HIGH = (IN >> 16) & 0xffff;
	HIGH = (HIGH ^ R2[1]) & 0xffff;
	if (BIT(HIGH, 0) == 1) {
		key[7] = key[7] | 0x80;
	}
	else {
		key[7] = key[7] & 0x7f;
	}
}

void ZUC::correctV0(u32 v[], u32 R2[]) {
	u32 IN = SFunInverse(R2[2], S0INVER, S1INVER);
	IN = inverseL2(IN);
	u32 HIGH = (IN >> 16) & 0xffff;
	HIGH = (HIGH ^ R2[1]) & 0xffff;
	if (BIT(HIGH, 0) == 1) {
		v[0] = v[0] | 0x80;
	}
	else {
		v[0] = v[0] & 0x7f;
	}
}

//the 30-round attack on the new loading scheme
bool ZUC::fulfillIVAutoBestAttackNew(u32 R2[], u32 key[], u32 v[]) {//correct IV and key
	//R2[0] = S\circ L2 (S5H || S11L)
	//S5H = K5 || D5 || K21[7]
	//S11L = IV4 || IV12
	u32 fInput = SFunInverse(R2[0], S0INVER, S1INVER);
	fInput = inverseL2(fInput);
	v[12] = fInput & 0xff;
	v[4] = (fInput >> 8) & 0xff;
	key[5] = (fInput >> 24) & 0xff;
	if (((fInput >> 16) & 0x1) == 1) {
		key[21] = key[21] | 0x80;
	}
	else {
		key[21] = key[21] & 0x7f;
	}

	//second round
	//R2[1] = S\circ L2 ( (R2[0]^X[2])L || (R1[0]+X[1])H )
	u32 U = SFunInverse(R2[1], S0INVER, S1INVER);
	U = inverseL2(U);
	//focus on (R2[0]^X[2])L = UH = R2[0]L ^ S6H
	//S6H = K6 || D6 || K22[7]
	u32 uh = (U >> 16) & 0xffff;
	uh = (uh ^ R2[0]) & 0xffff;
	key[6] = (uh >> 8) & 0xff;
	if ((uh & 0x1) == 1) {
		key[22] = key[22] | 0x80;
	}
	else {
		key[22] = key[22] & 0x7f;
	}
	//fulfill UL (R1[0]=SL(S9H || S7L)), X[1] = S12L || S10H
	//S9H = K9||D9||IV[2][7], S7L=IV0||IV8
	u32 d9 = 0x31;
	u32 S9H = (key[9] << 8) | (d9 << 1) | ((v[2] >> 7) & 0x1);
	u32 S7L = (v[0] << 8) | v[8];
	u32 R1 = L1((S9H << 16) | S7L);
	R1 = SFun(R1);
	//S10H = K10||D10||IV[3][7], S12L=IV5||IV13
	u32 d10 = 0x18;
	u32 S10H = (key[10] << 8) | (d10 << 1) | ((v[3] >> 7) & 0x1);
	u32 low = (R1 + S10H) & 0xffff;
	u32 full = ((U & 0xffff) << 16) | low;//UL || low
	u32 S12L = full - R1;
	S12L = S12L >> 16;
	v[5] = S12L >> 8;
	v[13] = S12L & 0xff;

	//third round
	u32 add = (S12L << 16) | (S10H);
	u32 R11 = R1 + add;//the output of the first register after adding X2 in the second round
	u32 R11H = (R11) & 0xffff;//shift to the 16 MSB
	u32 R11L = (v[1] << 8) | v[9];
	R11L = R11L ^ ((R2[0] >> 16) & 0xffff);//S8L^R2[0]H
	u32 R12 = (R11H << 16) | R11L;
	R12 = L1(R12);
	R12 = SFun(R12);
	//!!!!(old)compute S7H:IV10 || d7 || IV19 || K7[7]
	//compute S7H:K7 || d7 || V0[7] 
	u32 IN = SFunInverse(R2[2], S0INVER, S1INVER);
	IN = inverseL2(IN);
	u32 HIGH = (IN >> 16) & 0xffff;
	HIGH = (HIGH ^ R2[1]) & 0xffff;
	key[7] = (HIGH >> 8) & 0xff;
	if (BIT(HIGH, 0) != BIT(v[0], 7)) {//cannot modify v[8][7]
		return false;
	}
	//consider LOW of IN
	u32 LOW = IN & 0xffff;
	//LOW = ( R12+(S13L,S11H) )H
	u32 S11H = (key[11] << 8) | (D[11]<<1) | (v[4] >> 7);
	low = (R12 + S11H) & 0xffff;
	full = (LOW << 16) | low;//UL || low
	u32 S13L = (full - R12);
	S13L = S13L >> 16;
	v[14] = S13L & 0xff;
	v[6] = (S13L >> 8) & 0xff;

	/*cout << hex << "S7H:";
	toBinary(HIGH, 16, 0);*/

	return true;
}