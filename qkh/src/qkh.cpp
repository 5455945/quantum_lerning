#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <complex>
//#include <random>

// �����ӻ����˶������Ӽ����еĸ��ż���Ӱ�졷--Ҷ�� ���� ��
// ���Ż����� QKH ģ�͵����ӷ������
// QKH:��������������Harperģ��(quantum kecked Harper)

#pragma warning(disable: 4996)
using namespace std;
typedef complex<double> complexe;

//#define M_PI _Pi
#define M_PI 3.14159265

double K, L, hbar;  // Harper ģ�Ͳ���
int nq, nr;    // nq:���ӼĴ��������ӱ�����Ŀ(��1���������ӱ���),nr=nq-1
long int N, Nr;  // Hilbert �ռ�ά��
int steps;       // ��ʱ��Ƭ�㷨�зֽ�������и���

int p_pair, q_pair;
long int p, p_imp, q, q_imp;
int res1, res2;
long int portes;

double gam;    // ���������
int dissipation_err;    // ��Ǻ�ɢ�����Ƿ����
int Ntray=150; // ���ӹ켣����Ŀ
int channel = 4;   // ѡ���ɢ���ŵ�����ģ��:��ֵ�����ŵ�����λ�����ŵ���
int Iter = 10;   // ����ɲ���

// ��������"num"�Ķ����Ʊ�ʾ��"1"����Ŀ֮��
int binsum(int num) {
    int sum, temp;
    sum = 0;
    while (num != 0) {
        temp = num % 2;
        if (temp == 1) {
            sum++;
        }
        num = num / 2;
    }
    return sum;
}

complexe ampdamp2a(complexe *VE, int k) {
    int j2, b, i, j;
    complexe VEtemp, aux1;

    VEtemp = 0.0;
    i = 0;
    for (j2 = 0; j2 < N / 2; j2++) {
        j = j2 / (i << k) * (1 << (k + 1)) + (1 << k) + j2 - j2 / (1 << k) * (1 <<k);
        i = j - (1 << k);
        if (i == 0) {
            VEtemp = VEtemp + conj(VE[j]) * VE[j];
        }
        else {
            if (channel == 0) {
                b = binsum(j);

                aux1 = conj(complexe(sqrt(1. / float(b)), 0.) * VE[j]);
                aux1 = aux1 * complexe(sqrt(1. / float(b)), 0.) * VE[j];
            }
            else {
                aux1 = conj(VE[j]) * VE[j];
            }
            VEtemp = VEtemp + aux1;
        }
    }
    return VEtemp;
}

// ��ֵ�����ŵ�ģ�ͣ����� L_mu |psi>
void ampdamp2b(complexe *VE, complexe *VE3, int k) {
    int j2, b, i, j;

    for (j2 = 0; j2 < N; j2++) {
        VE3[j2] = 0.0;
    }
    for (j2 = 0; j2 < N / 2; j2++) {
        j = j2 / (1 << k) * (1 << (k + 1)) + (1 << k) + j2 - j2 / (1 << k) * (1 << k);
        i = j - (1 << k);
        if (i == 0) {
            VE3[i] = VE[j];
        }
        else if (channel == 0) {
            b = binsum(j);
            VE3[i] = complexe(sqrt(1. / float(b)), 0.) * VE[j];
        }
        else {
            VE3[i] = VE[j];
        }
    }
}

// ȥ�����ŵ�ģ�ͣ����� L_mu |psi>
void depolarizaion2b(complexe *VE, complexe *VE3, int k, int j) {
    int i, i2, j2;
    complexe aux1;

    // ���ط�ת����
    if (j == 0) {
        for (i = 0; i < N / 2; i++) {
            i2 = (i / (1 << k)) * (1 << (k + 1)) + i - (i / (1 << k)) * (1 <<k);
            j2 = i2 + (1 << k);
            VE3[i2] = VE[j2];

            j2 = i2;
            i2 = i2 + (1 << k);
            VE3[i2] = VE[j2];
        }
    }

    // ��λ��ת����
    if (j == 2) {
        for (i = 0; i < N / 2; i++) {
            i2 = (i / (1 << k)) * (1 << (k + 1)) + i - (i / (1 << k)) * (1 << k);
            VE3[i2] = VE[i2];

            i2 = i2 + (1 << k);
            VE3[i2] = -1.0 * VE[i2];
        }
    }

    // ���غ���λͬʱ��ת����
    if (j == 1) {
        for (i = 0; i < N / 2; i++) {
            i2 = (i / (1 << k)) * (1 << (k + 1)) + i - (i / (1 << k)) * (1 << k);
            j2 = i2 + (1 << k);
            VE3[i2] = complexe(0.0, -1.0) * VE[j2];

            j2 = i2;
            i2 = i2 + (1 << k);
            VE3[i2] = complexe(0.0, 1.0) * VE[j2];
        }
    }
}

// ���ӹ켣�Ӻ���
void distep(complexe *VE) {
    int i, j, jj, k, channelNr;
    int init;
    double normalisation, delta;
    double dp, dpm;
    //double dpv[3 * nq];
    complexe VEtemp, *VE3;
    int bsum;

    double* dpv = (double*)malloc(3 * nq *sizeof(double));
    if (dpv == nullptr) {
        fprintf(stderr, "Memory allocaltion of dpv failed.\n");
        exit(1);
    }

    VE3 = (complexe*)malloc(N * sizeof(complexe));
    if (VE3 == nullptr) {
        fprintf(stderr, "Memory allocaltion of VE3 failed.\n");
        exit(1);
    }

    if (channel <= 1) {
        channelNr = nq;
    }
    else {
        if (channel == 5) {
            channelNr = 3 * nq;
        }
        else {
            channelNr = nq;
        }
    }

    for (jj = 1; jj <= Iter; jj++) {
        if (channel <= 1) {  // ��ֵ�����ŵ�
            for (k = 0; k < nq; k++) {
                VEtemp = ampdamp2a(VE, k);
                dp = real(VEtemp) / (double)Iter;
                dpv[k+1] = gam * dp;
            }
        }
        else {  // ȥ�����ŵ�
            for (k = 0; k < nq; k++) {
                if (channel == 5) {
                    for (j = 0; j <= 2; j++) {
                        dp = 1.0 / (double)Iter;
                        dpv[1 + k*3 +j] = gam * dp;
                    }
                }
                else {
                    j = channel - 2;
                    dp = 1.0 / (float)Iter;
                    dpv[1 + k] = gam * dp;
                }
            }
        }

        dp = 0.0;
        for (k = 1; k <= channelNr; k++) {
            dp = dp + dpv[k];
        }

        //delta = (random() / (double)RAND_MAX);  // ���������delta
        delta = (rand() / (double)RAND_MAX);  // ���������delta

        if (delta < dp) {  // ����ԾǨ
            j = 1;
            dpm = dpv[j];
            while ((delta > dpm) && (j < channelNr)) {
                j = j + 1;
                dpm = dpm + dpv[j];
            }
            j = j - 1;

            if (channel <= 1) {
                ampdamp2b(VE, VE3, j);
            }
            else if (channel == 5) {
                depolarizaion2b(VE, VE3, (j - (j % 3)) / 3, (j % 3));
            }
            else {
                depolarizaion2b(VE, VE3, j, channel - 2);
            }

            normalisation = 0.0;
            for (i = 0; i < N; i++) {
                normalisation += norm(VE3[i]);
            }
            for (i = 0; i < N; i++) {
                VE[i] = VE3[i] / sqrt(normalisation);
            }
        }
        else {  // �� H_eff ���ݻ�
            for (i = 0; i < N; i++) {
                VE3[i] = 0.0;
            }
            if (channel <= 1) {
                VE3[0] = VE[0];
                init = 1;
            }
            else {
                init = 0;
            }

            for (i = init; i < N; i++) {
                if (channel <= 1) {
                    if (channel == 1) {
                        bsum = binsum(i);
                        VE3[i] = (1.0 - (gam / 2.0 / (float)Iter * (float)bsum)) * VE[i];
                    }
                    else {
                        VE3[i] = (1.0 - (gam / 2.0 / (float)Iter)) * VE[i];
                    }
                }
                else {
                    VE3[i] = (1.0 - (gam / 2.0 / (float)Iter * (float)channelNr)) * VE[i];
                }
            }

            normalisation = 0.0;
            for (i = 0; i < N; i++) {
                normalisation += norm(VE3[i]);
            }
            for (i = 0; i < N; i++) {
                VE[i] = VE3[i] / sqrt(normalisation);
            }
        }
    }

    if (VE3 != nullptr) {
        free(VE3);
        VE3 = nullptr;
    }
    if (dpv != nullptr) {
        free(dpv);
        dpv = nullptr;
    }
}

// �Ե�j�����ӱ��ؽ��� Hadamard �����ű任
void Porte_Ai(complexe *n, int j) {
    long int i, masque, complement;
    complexe temp0, temp1;
    complexe a11, a12, a21, a22;

    a11 = 1 / sqrt(2.0);
    a12 = a11;
    a21 = a11;
    a22 = -a11;

    masque = (long int)1 << j;
    for (i = 0; i < N; i++) {
        if (!(i & masque)) {
            complement = i + masque;
            temp0 = n[i];
            temp1 = n[complement];
            n[i] = a11 * temp0 + a12 * temp1;
            n[complement] = a21 * temp0 + a22 * temp1;
        }
    }

    if (dissipation_err == 1) {
        distep(n);  // �����ɢ����
    }

    portes++;
}

// ���Ӹ���Ҷ�任�е��ܿ���λ��
void Porte_Bjk(complexe *n, int j, int k, int signe) {
    long int i, masquej, masquek;
    complexe phase;

    masquej = (long int)1 << j;
    masquek = (long int)1 << k;

    phase = polar(1.0, signe * M_PI / ((long int)1 << (k - j)));
    for (i = 0; i < N; i++) {
        if ((i & masquej) && (i & masquek)) {
            n[i] *= phase;
        }
    }

    if (dissipation_err == 1) {
        distep(n);
    }
    portes++;
}

// ���ƽ�Ϊ theta ���ܿ���λ��
void Control_Phase(complexe *n, int j, int k, double theta) {
    long int i, masquej, masquek;
    complexe phase;

    masquej = (long int)1 << j;
    masquek = (long int)1 << k;
    phase = polar(1.0, theta);
    for (i = 0; i < N; i++) {
        if ((i & masquej) && (i & masquek)) {
            n[i] *= phase;
        }
    }

    if (dissipation_err == 1) {
        distep(n);
    }

    portes++;
}

// �������ӱ��� j �� k
void Porte_Echange(complexe *n, int j, int k) {
    long int i, masquej, masquek, complement;
    complexe temp;

    masquej = (long int) 1 << j;
    masquek = (long int) 1 << k;

    for (i = 0; i < N; i++) {
        if ((!(i & masquej)) && (i & masquek)) {
            complement = i - masquek + masquej;
            temp = n[i];
            n[i] = n[complement];
            n[complement] = temp;
        }
    }
}

// ���Ӹ���Ҷ�任, signe ȡ -1 ��ʾ���Ӹ���Ҷ��任
void QFT(complexe *n, int signe) {
    int x, m;
    for (x = nr - 1; x >= 0; x--) {
        for (m = nr - 1; m > x; m--) {
            Porte_Bjk(n, x, m, signe);
        }
        Porte_Ai(n, x);
    }

    for (x = 0; x < (nr / 2); x++) {
        Porte_Echange(n, x, nr - x - 1);
    }
}

// ��ʱ��Ƭ�ƽ��㷨�е���λ����ת��
void Rotation_z(complexe *n, int j, double coef) {
    long int i, masquej;
    complexe phase;

    masquej = (long int) 1 << j;
    phase = polar(1.0, coef / 2.0);
    for (i = 0; i < N; i++) {
        if (i & masquej) {
            n[i] /= phase;
        }
        else {
            n[i] *= phase;
        }
    }

    if (dissipation_err == 1) {
        distep(n);
    }
    portes++;
}

// ��ʱ��Ƭ�ƽ��㷨�е��ܿ�U��
void Control_U(complexe *n, int j, int debut, long int pp) {
    int k;
    for (k = 0; k <= debut; k++) {
        Control_Phase(n, j, debut - k, pp * M_PI / pow(2.0, (double)k));
    }
}

// ��ʱ��Ƭ�ƽ��㷨ʵ�ֵ����������㷨
void Kick(complexe *n, double coef, long int pp, int a) {
    double alpha;
    int i;
    int control_qubit, debut;

    alpha = -coef /steps;
    control_qubit = nr;
    debut = nr - 1 - a;

    Porte_Ai(n, control_qubit);
    Control_U(n, control_qubit, debut, -pp);
    Porte_Ai(n, control_qubit);

    for (i = 0; i < (steps - 1); i++) {
        Rotation_z(n, control_qubit, alpha / 2.0);
        Porte_Ai(n, control_qubit);
        Control_U(n, control_qubit, debut, 2 * pp);
        Porte_Ai(n, control_qubit);
        Rotation_z(n, control_qubit, alpha);
        Porte_Ai(n, control_qubit);
        Control_U(n, control_qubit, debut, -2 * pp);
        Porte_Ai(n, control_qubit);
        Rotation_z(n, control_qubit, alpha / 2.0);
    }

    Porte_Ai(n, control_qubit);
    Control_U(n, control_qubit, debut, pp);
    Porte_Ai(n, control_qubit);
}

// ���������� Harper ģ���ݻ�
void Evolution(complexe *n) {
    QFT(n, 1);  // ���Ӹ���Ҷ�任
    Kick(n, K / hbar, q_imp,  q_pair); // ��� U_theta ���ݻ�
    QFT(n, -1); // ���Ӹ���Ҷ��任
    Kick(n, L / hbar, p_imp, p_pair); // ��� U_p ���ݻ�
}

// ������̬ nh ��Ϊ������ (p0, q0) �ĸ�˹����
void gauss_add(complexe *nh, double p0, double q0, double sigma2) {
    int p;
    double angle, fak, norm, xx, NN, N2;
    complexe val;

    norm = 1.0 / sqrt(sqrt(2.0 * M_PI * sigma2));
    NN = (double) Nr;
    N2 = NN / 2.0;
    for (p = 0; p < Nr; p++) {
        xx = (double) p - p0;
        if (xx > N2) {
            xx -= NN;
        }
        if (xx < -N2) {
            xx += NN;
        }
        xx = xx * xx / 4.0 / sigma2;
        fak = exp (-xx) * norm;
        angle = -q0 * (double)p;
        val = polar(fak, angle);
        nh[p] = nh[p] + val;
    }
}

// ���㱣���
double fidelity(complexe *n_clean, complexe *n) {
    complexe sum;
    int i;

    sum = 0.0;
    for (i = 0; i < Nr; i++) {
        sum += conj(n_clean[i]) * n[i];
    }
    return pow(abs(sum), 2);
}

int main(int argc, char** argv) {
    complexe *n, *nh, *n_clean;
    long int i, j, iterations;
    int o;
    double itoq, itop;
    double normalisation;
    double *fid;   // �����

    // ����5������: Harperģ���в���K��L,���������,����ģ����ʵ�����ӱ�����Ŀ,�Լ�Harperģ�͵�������
    if (argc < 6) {
        puts("syntax: K L dissipation_rate nr iterations");
        exit(1);
    }
    if (sscanf(argv[1], "%lf", &K) == 0) {
        puts("syntax: K L dissipation_rate nr iterations");
        exit(1);
    }
    if (sscanf(argv[2], "%lf", &L) == 0) {
        puts("syntax: K L dissipation_rate nr iterations");
        exit(1);
    }
    if (sscanf(argv[3], "%lf", &gam) == 0) {
        puts("syntax: K L dissipation_rate nr iterations");
        exit(1);
    }
    if (sscanf(argv[4], "%d", &nr) == 0) {
        puts("syntax: K L dissipation_rate nr iterations");
        exit(1);
    }
    if (sscanf(argv[5], "%ld", &iterations) == 0) {
        puts("syntax: K L dissipation_rate nr iterations");
        exit(1);
    }

    //srandom((unsigned int)time(NULL));
    srand((unsigned int)time(NULL));
    Nr = 1 << nr;
    nq = nr + 1;
    N = 1 << nq;
    q = 1; // ��ռ���λ��������������Ŀ
    p = 1; // ��ռ��ж���������������Ŀ
    hbar = 2.0 * M_PI / (6.0 + (sqrt(5.0) + 1.0) / 2.0);
    //hbar = 2.0 * M_PI * p * q / ((double)Nr);

    // ÿ��������������ж�ʱ��Ƭ����Ŀ
    steps = 100;

    itoq = 2.0 * M_PI * q / ((double)Nr);
    itop = 2.0 * M_PI * p / ((double)Nr);
    p_pair = 0;
    p_imp = p;
    while (!(p_imp & 1)) {
        p_imp >>= 1;
        p_pair++;
    }

    q_pair = 0;
    q_imp = q;
    while (!(q_imp & 1)) {
        q_imp >>= 1;
        q_pair++;
    }

    // �ܿ�����̬�ڴ����
    n = (complexe*)malloc(N*sizeof(complexe));
    if (n == nullptr) {
        fprintf(stderr, "Memory allocation of n failed.\n");
        exit(1);
    }

    // �����ݻ�������̬�ڴ����
    n_clean = (complexe*)malloc(N*sizeof(complexe));
    if (n_clean == nullptr) {
        fprintf(stderr, "Memory allocation of n_clean failed.\n");
        exit(1);
    }

    nh = (complexe*)malloc(Nr*sizeof(complexe));
    if (nh == nullptr) {
        fprintf(stderr, "Memory allocation of nh failed.\n");
        exit(1);
    }

    fid = (double*)malloc(iterations*sizeof(double));
    if (fid == nullptr) {
        fprintf(stderr, "Memory allocations of fid failed.\n");
        exit(1);
    }

    for (int i = 0; i < iterations; i++) {
        fid[i] = 0.0;
    }

    // ����̬��ʼ��Ϊһ��˹����
    gauss_add(nh, Nr/2, Nr/2, sqrt(hbar/2.0));

    // ���ӹ켣ѭ��
    for (o = 1; o <= Ntray; o++) {
        for (i = 0; i < Nr; i++) {
            n[i] = nh[i];
        }
        for (i = Nr; i < N; i++) {
            n[i] = 0.0;
        }
        for (i = 0; i < N; i++) {
            n_clean[i] = n[i];
        }

        // Harperģ�͵��ݻ�, ��iterations��
        for (j = 0; j < iterations; j++) {
            fid[j] = fid[j] + fidelity(n_clean, n);

            // �������Ӱ�������̬�ݻ�
            dissipation_err = 1;
            Evolution(n);
            normalisation = 0.0;  // ��һ��
            for (i = 0; i < Nr; i++) {
                normalisation += norm(n[i]);
            }
            for (i = 0; i < Nr; i++) {
                n[i] = n[i] / sqrt(normalisation);
            }

            // �������Ӱ�������̬�ݻ�
            dissipation_err = 0;
            Evolution(n_clean);
            normalisation = 0.0;
            for (i = 0; i < Nr; i++) {
                normalisation += norm(n_clean[i]);
            }
            for (i = 0; i < Nr; i++) {
                n_clean[i] = n_clean[i] / sqrt(normalisation);
            }
        }
    }

    for (i = 0; i < iterations; i++) {
        fid[i] = fid[i] / (double)Ntray;
        printf("%f\n", fid[i]);
    }

    if (n != nullptr) {
        free(n);
        n = nullptr;
    }
    if (n_clean != nullptr) {
        free(n_clean);
        n_clean = nullptr;
    }
    if (nh != nullptr) {
        free(nh);
        nh = nullptr;
    }
    if (fid != nullptr) {
        free(fid);
        fid = nullptr;
    }

    return 0;
}