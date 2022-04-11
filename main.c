// Adapted from Prof. Koyama's MPF code in his textbook
// Author: Chuanqi Zhu
// Created on: 2022/2/16
// Updated on 2022/04/11

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define NDX 40 //差分計算における計算領域一辺の分割数
#define NDY 40
#define NDZ 40

#define N 2 //考慮する結晶方位の数＋１(MPF0.cppと比較して、この値を大きくしている)

int ndx = NDX;
int ndy = NDY;
int ndz = NDZ;
int ndmx = NDX - 1;
int ndmy = NDY - 1; //計算領域の一辺の差分分割数(差分ブロック数), ND-1を定義
int ndmz = NDZ - 1;
int nm = N - 1, nmm = N - 2; //考慮する結晶方位の数、N-2（考慮する結晶方位の数－１）を定義
double PI = 3.141592;        //π、計算カウント数
double RR = 8.3145;          //ガス定数

double aij[N][N]; //勾配エネルギー係数
double wij[N][N]; //ペナルティー項の係数
double mij[N][N]; //粒界の易動度
double fij[N][N]; //粒界移動の駆動力
int anij[N][N];
double thij[N][N];
double vpij[N][N];
double etaij[N][N];
int phinum;

int i, j, k, l, ii, jj, kk, ll, it; //整数
int ip, im, jp, jm, lp, lm;         //整数
int n1, n2, n3;                     //整数

int istep = 0;
// int n000;		//位置(i,j)において、pが０ではない方位の個数（n00>=n000）
int nstep; //計算カウント数の最大値（計算終了カウント）
int pstep;
double dtime, L, dx;     // L計算領域の一辺の長さ(nm), 差分プロック１辺の長さ(m)
double M0;               //粒界の易動度
double W0;               //ペナルティー項の係数
double A0;               //勾配エネルギー係数
double S0;               //粒界移動の駆動力
double temp;             //温度
double sum1, sum2, sum3; //各種の和の作業変数
double pddtt;            //フェーズフィールドの時間変化率

double gamma0; //粒界エネルギ密度
double delta;  //粒界幅（差分ブロック数にて表現）
double mobi;   //粒界の易動度
double vm0;    //モル体積

double astre;
double th, vp, eta;
double thetax, thetay;
double epsilon0;
double termiikk, termjjkk;

double phidx, phidy, phidz;
double phidxx, phidyy, phidzz;
double phidxy, phidxz, phidyz;
double phiabs;

double xxp, xyp, xzp, yxp, yyp, yzp, zxp, zyp, zzp;

double phidxp, phidyp, phidzp;
double phidxpx, phidypx, phidzpx;
double phidxpy, phidypy, phidzpy;
double phidxpz, phidypz, phidzpz;
double ep, epdx, epdy, epdz;
double term0;
double termx, termx0, termx1, termx0dx, termx1dx;
double termy, termy0, termy1, termy0dy, termy1dy;
double termz, termz0, termz1, termz0dz, termz1dz;

int x11, y11, x1h[10], y1h[10]; //初期核の座標
double t, r0, r;

//******* メインプログラム ******************************************
int main(int argc, char *argv[])
{
    nstep = 501;
    pstep = 100;
    dtime = 5.0;
    temp = 1000.0;
    L = 2000.0;
    vm0 = 7.0e-6;
    delta = 7.0;
    mobi = 1.0;

    dx = L / 40.0 * 1.0e-9;              //差分プロック１辺の長さ(m)
    gamma0 = 0.5 * vm0 / RR / temp / dx; //粒界エネルギ密度（0.5J/m^2）を無次元化
    astre = 0.05;
    A0 = 8.0 * delta * gamma0 / PI / PI; //勾配エネルギー係数[式(4.40)]
    W0 = 4.0 * gamma0 / delta;           //ペナルティー項の係数[式(4.40)]
    M0 = mobi * PI * PI / (8.0 * delta); //粒界の易動度[式(4.40)]
    S0 = 50.0 / RR / temp;               //粒界移動の駆動力

    for (i = 0; i <= nm; i++)
    {
        for (j = 0; j <= nm; j++)
        {
            wij[i][j] = W0;
            aij[i][j] = A0;
            mij[i][j] = M0;
            fij[i][j] = 0.0;
            anij[i][j] = 0;
            thij[i][j] = 0.0;
            vpij[i][j] = 0.0;
            etaij[i][j] = 0.0;
            if ((i == 0) || (j == 0))
            {
                fij[i][j] = S0;
                anij[i][j] = 1;
            }
            if (i < j)
            {
                fij[i][j] = -fij[i][j];
            }
            if (i == j)
            {
                wij[i][j] = 0.0;
                aij[i][j] = 0.0;
                mij[i][j] = 0.0;
                fij[i][j] = 0.0;
                anij[i][j] = 0;
            }
        }
    }
    // thij[1][2] = PI / 4.0;
    // thij[2][1] = PI / 4.0;
    // vpij[1][2] = PI / 4.0;
    // vpij[2][1] = PI / 4.0;
    // etaij[1][2] = PI / 4.0;
    // etaij[2][1] = PI / 4.0;

    double phi[N][NDX][NDY][NDZ], phi2[N][NDX][NDY][NDZ]; //フェーズフィールド、フェーズフィールド補助配列
    int phiIdx[N + 1][NDX][NDY][NDZ];                     //位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の番号
    int phiNum[NDX][NDY][NDZ];

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (l = 0; l <= ndmz; l++)
            {
                if ((i - NDX / 2) * (i - NDX / 2) + (j - NDY / 2) * (j - NDY / 2) + (l - NDZ / 2) * (l - NDZ / 2) < 36)
                {
                    phi[1][i][j][l] = 1.0;
                    phi[0][i][j][l] = 0.0;
                }
                else
                {
                    phi[1][i][j][l] = 0.0;
                    phi[0][i][j][l] = 1.0;
                }
            }
        }
    }
    // }

start:;

    if ((((int)(istep) % pstep) == 0))
    {
        FILE *stream;
        char buffer[30];
        sprintf(buffer, "data/phi/3d%d.vtk", istep);
        stream = fopen(buffer, "a");

        fprintf(stream, "# vtk DataFile Version 1.0\n");
        fprintf(stream, "phi_%d.vtk\n", istep);
        fprintf(stream, "ASCII\n");
        fprintf(stream, "DATASET STRUCTURED_POINTS\n");
        fprintf(stream, "DIMENSIONS %d %d %d\n", NDX, NDY, NDZ);
        fprintf(stream, "ORIGIN 0.0 0.0 0.0\n");
        fprintf(stream, "ASPECT_RATIO 1.0 1.0 1.0\n");
        fprintf(stream, "\n");
        fprintf(stream, "POINT_DATA %d\n", NDX * NDY * NDZ);
        fprintf(stream, "SCALARS scalars float\n");
        fprintf(stream, "LOOKUP_TABLE default\n");

        for (i = 0; i <= ndmx; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (l = 0; l <= ndmz; l++)
                {
                    fprintf(stream, "%e\n", phi[1][i][j][l]);
                }
            }
        }
        fclose(stream);
    }

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (l = 0; l <= ndmz; l++)
            {
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                lp = l + 1;
                lm = l - 1;
                if (i == ndmx)
                {
                    ip = 0;
                }
                if (i == 0)
                {
                    im = ndmx;
                }
                if (j == ndmy)
                {
                    jp = 0;
                }
                if (j == 0)
                {
                    jm = ndmy;
                }
                if (l == ndmz)
                {
                    lp = 0;
                }
                if (l == 0)
                {
                    lm = ndmz;
                }

                //--- 位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の個数---
                phinum = 0;
                for (ii = 0; ii <= nm; ii++)
                {
                    if ((phi[ii][i][j][l] > 0.0) ||
                        ((phi[ii][i][j][l] == 0.0) && (phi[ii][ip][j][l] > 0.0) ||
                         (phi[ii][im][j][l] > 0.0) ||
                         (phi[ii][i][jp][l] > 0.0) ||
                         (phi[ii][i][jm][l] > 0.0) ||
                         (phi[ii][i][j][lp] > 0.0) ||
                         (phi[ii][i][j][lm] > 0.0)))
                    {
                        phinum++;
                        phiIdx[phinum][i][j][l] = ii;
                    }
                }
                phiNum[i][j][l] = phinum;
            }
        }
    }

    // Evolution Equations
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (l = 0; l <= ndmz; l++)
            {
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                lp = l + 1;
                lm = l - 1;
                if (i == ndmx)
                {
                    ip = 0;
                }
                if (i == 0)
                {
                    im = ndmx;
                } //周期的境界条件
                if (j == ndmy)
                {
                    jp = 0;
                }
                if (j == 0)
                {
                    jm = ndmy;
                }
                if (l == ndmz)
                {
                    lp = 0;
                }
                if (l == 0)
                {
                    lm = ndmz;
                }

                for (n1 = 1; n1 <= phiNum[i][j][l]; n1++)
                {
                    ii = phiIdx[n1][i][j][l];
                    pddtt = 0.0;
                    for (n2 = 1; n2 <= phiNum[i][j][l]; n2++)
                    {
                        jj = phiIdx[n2][i][j][l];
                        sum1 = 0.0;
                        for (n3 = 1; n3 <= phiNum[i][j][l]; n3++)
                        {
                            kk = phiIdx[n3][i][j][l];

                            // calculate the interface normal and deirivatives of the phase field
                            phidx = (phi[kk][ip][j][l] - phi[kk][im][j][l]) / 2.0;
                            phidy = (phi[kk][i][jp][l] - phi[kk][i][jm][l]) / 2.0;
                            phidz = (phi[kk][i][j][lp] - phi[kk][i][j][lm]) / 2.0;

                            phidxx = (phi[kk][ip][j][l] + phi[kk][im][j][l] - 2.0 * phi[kk][i][j][l]);
                            phidyy = (phi[kk][i][jp][l] + phi[kk][i][jm][l] - 2.0 * phi[kk][i][j][l]);
                            phidzz = (phi[kk][i][j][lp] + phi[kk][i][j][lm] - 2.0 * phi[kk][i][j][l]);

                            phidxy = (phi[kk][ip][jp][l] + phi[kk][im][jm][l] - phi[kk][im][jp][l] - phi[kk][ip][jm][l]) / 4.0;
                            phidxz = (phi[kk][ip][j][lp] + phi[kk][im][j][lm] - phi[kk][im][j][lp] - phi[kk][ip][j][lm]) / 4.0;
                            phidyz = (phi[kk][i][jp][lp] + phi[kk][i][jm][lm] - phi[kk][i][jm][lp] - phi[kk][i][jp][lm]) / 4.0;

                            phiabs = phidx * phidx + phidy * phidy + phidz * phidz;

                            if (anij[ii][kk] == 1 && phiabs != 0.0)
                            {
                                epsilon0 = sqrt(aij[ii][kk]);

                                th = thij[ii][kk];
                                vp = vpij[ii][kk];
                                eta = etaij[ii][kk];

                                xxp = cos(th) * cos(vp);
                                yxp = sin(th) * cos(vp);
                                zxp = sin(vp);
                                xyp = -sin(th) * cos(eta) - cos(th) * sin(vp) * sin(eta);
                                yyp = cos(th) * cos(eta) - sin(th) * sin(vp) * sin(eta);
                                zyp = cos(vp) * sin(eta);
                                xzp = sin(eta) * sin(th) - cos(eta) * cos(th) * sin(vp);
                                yzp = -sin(eta) * cos(th) - cos(eta) * sin(th) * sin(vp);
                                zzp = cos(eta) * cos(vp);

                                phidxp = phidx * xxp + phidy * yxp + phidz * zxp;
                                phidyp = phidx * xyp + phidy * yyp + phidz * zyp;
                                phidzp = phidx * xzp + phidy * yzp + phidz * zzp;

                                phidxpx = phidxx * xxp + phidxy * yxp + phidxz * zxp;
                                phidypx = phidxx * xyp + phidxy * yyp + phidxz * zyp;
                                phidzpx = phidxx * xzp + phidxy * yzp + phidxz * zzp;

                                phidxpy = phidxy * xxp + phidyy * yxp + phidyz * zxp;
                                phidypy = phidxy * xyp + phidyy * yyp + phidyz * zyp;
                                phidzpy = phidxy * xzp + phidyy * yzp + phidyz * zzp;

                                phidxpz = phidxz * xxp + phidyz * yxp + phidzz * zxp;
                                phidypz = phidxz * xyp + phidyz * yyp + phidzz * zyp;
                                phidzpz = phidxz * xzp + phidyz * yzp + phidzz * zzp;

                                ep = epsilon0 * (1.0 - 3.0 * astre + 4.0 * astre * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 2.0));

                                epdx = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpx + pow(phidyp, 3.0) * phidypx + pow(phidzp, 3.0) * phidzpx) / pow(phiabs, 2.0) - (phidx * phidxx + phidy * phidxy + phidz * phidxz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));
                                epdy = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpy + pow(phidyp, 3.0) * phidypy + pow(phidzp, 3.0) * phidzpy) / pow(phiabs, 2.0) - (phidx * phidxy + phidy * phidyy + phidz * phidyz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));
                                epdz = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpz + pow(phidyp, 3.0) * phidypz + pow(phidzp, 3.0) * phidzpz) / pow(phiabs, 2.0) - (phidx * phidxz + phidy * phidyz + phidz * phidzz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));

                                term0 = 2.0 * ep * epdx * phidx + phidxx * ep * ep + 2.0 * ep * epdy * phidy + phidyy * ep * ep + 2.0 * ep * epdz * phidz + phidzz * ep * ep;

                                termx0 = (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / phiabs;
                                termy0 = (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / phiabs;
                                termz0 = (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / phiabs;

                                termx1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidx / pow(phiabs, 2.0);
                                termy1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidy / pow(phiabs, 2.0);
                                termz1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidz / pow(phiabs, 2.0);

                                termx0dx = (3.0 * pow(phidxp, 2.0) * phidxpx * xxp + 3.0 * pow(phidyp, 2.0) * phidypx * xyp + 3.0 * pow(phidzp, 2.0) * phidzpx * xzp) / phiabs - (2.0 * phidx * phidxx + 2.0 * phidy * phidxy + 2.0 * phidz * phidxz) * (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / pow(phiabs, 2.0);
                                termy0dy = (3.0 * pow(phidxp, 2.0) * phidxpy * yxp + 3.0 * pow(phidyp, 2.0) * phidypy * yyp + 3.0 * pow(phidzp, 2.0) * phidzpy * yzp) / phiabs - (2.0 * phidx * phidxy + 2.0 * phidy * phidyy + 2.0 * phidz * phidyz) * (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / pow(phiabs, 2.0);
                                termz0dz = (3.0 * pow(phidxp, 2.0) * phidxpz * zxp + 3.0 * pow(phidyp, 2.0) * phidypz * zyp + 3.0 * pow(phidzp, 2.0) * phidzpz * zzp) / phiabs - (2.0 * phidx * phidxz + 2.0 * phidy * phidyz + 2.0 * phidz * phidzz) * (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / pow(phiabs, 2.0);

                                termx1dx = ((phidxx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidx * (4.0 * pow(phidxp, 3.0) * phidxpx + 4.0 * pow(phidyp, 3.0) * phidypx + 4.0 * pow(phidzp, 3.0) * phidzpx))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxx + phidy * phidxy + phidz * phidxz) * phidx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);
                                termy1dy = ((phidyy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidy * (4.0 * pow(phidxp, 3.0) * phidxpy + 4.0 * pow(phidyp, 3.0) * phidypy + 4.0 * pow(phidzp, 3.0) * phidzpy))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxy + phidy * phidyy + phidz * phidyz) * phidy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);
                                termz1dz = ((phidzz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidz * (4.0 * pow(phidxp, 3.0) * phidxpz + 4.0 * pow(phidyp, 3.0) * phidypz + 4.0 * pow(phidzp, 3.0) * phidzpz))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxz + phidy * phidyz + phidz * phidzz) * phidz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);

                                termx = 16.0 * epsilon0 * astre * (epdx * (termx0 - termx1) + ep * (termx0dx - termx1dx));
                                termy = 16.0 * epsilon0 * astre * (epdy * (termy0 - termy1) + ep * (termy0dy - termy1dy));
                                termz = 16.0 * epsilon0 * astre * (epdz * (termz0 - termz1) + ep * (termz0dz - termz1dz));

                                termiikk = term0 + termx + termy + termz;
                            }
                            else
                            {
                                termiikk = aij[ii][kk] * (phidxx + phidyy + phidzz);
                            }

                            if (anij[jj][kk] == 1 && phiabs != 0.0)
                            {
                                epsilon0 = sqrt(aij[jj][kk]);

                                th = thij[jj][kk];
                                vp = vpij[jj][kk];
                                eta = etaij[jj][kk];

                                xxp = cos(th) * cos(vp);
                                yxp = sin(th) * cos(vp);
                                zxp = sin(vp);
                                xyp = -sin(th) * cos(eta) - cos(th) * sin(vp) * sin(eta);
                                yyp = cos(th) * cos(eta) - sin(th) * sin(vp) * sin(eta);
                                zyp = cos(vp) * sin(eta);
                                xzp = sin(eta) * sin(th) - cos(eta) * cos(th) * sin(vp);
                                yzp = -sin(eta) * cos(th) - cos(eta) * sin(th) * sin(vp);
                                zzp = cos(eta) * cos(vp);

                                phidxp = phidx * xxp + phidy * yxp + phidz * zxp;
                                phidyp = phidx * xyp + phidy * yyp + phidz * zyp;
                                phidzp = phidx * xzp + phidy * yzp + phidz * zzp;

                                phidxpx = phidxx * xxp + phidxy * yxp + phidxz * zxp;
                                phidypx = phidxx * xyp + phidxy * yyp + phidxz * zyp;
                                phidzpx = phidxx * xzp + phidxy * yzp + phidxz * zzp;

                                phidxpy = phidxy * xxp + phidyy * yxp + phidyz * zxp;
                                phidypy = phidxy * xyp + phidyy * yyp + phidyz * zyp;
                                phidzpy = phidxy * xzp + phidyy * yzp + phidyz * zzp;

                                phidxpz = phidxz * xxp + phidyz * yxp + phidzz * zxp;
                                phidypz = phidxz * xyp + phidyz * yyp + phidzz * zyp;
                                phidzpz = phidxz * xzp + phidyz * yzp + phidzz * zzp;

                                ep = epsilon0 * (1.0 - 3.0 * astre + 4.0 * astre * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 2.0));

                                epdx = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpx + pow(phidyp, 3.0) * phidypx + pow(phidzp, 3.0) * phidzpx) / pow(phiabs, 2.0) - (phidx * phidxx + phidy * phidxy + phidz * phidxz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));
                                epdy = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpy + pow(phidyp, 3.0) * phidypy + pow(phidzp, 3.0) * phidzpy) / pow(phiabs, 2.0) - (phidx * phidxy + phidy * phidyy + phidz * phidyz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));
                                epdz = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpz + pow(phidyp, 3.0) * phidypz + pow(phidzp, 3.0) * phidzpz) / pow(phiabs, 2.0) - (phidx * phidxz + phidy * phidyz + phidz * phidzz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));

                                term0 = 2.0 * ep * epdx * phidx + phidxx * ep * ep + 2.0 * ep * epdy * phidy + phidyy * ep * ep + 2.0 * ep * epdz * phidz + phidzz * ep * ep;

                                termx0 = (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / phiabs;
                                termy0 = (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / phiabs;
                                termz0 = (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / phiabs;

                                termx1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidx / pow(phiabs, 2.0);
                                termy1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidy / pow(phiabs, 2.0);
                                termz1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidz / pow(phiabs, 2.0);

                                termx0dx = (3.0 * pow(phidxp, 2.0) * phidxpx * xxp + 3.0 * pow(phidyp, 2.0) * phidypx * xyp + 3.0 * pow(phidzp, 2.0) * phidzpx * xzp) / phiabs - (2.0 * phidx * phidxx + 2.0 * phidy * phidxy + 2.0 * phidz * phidxz) * (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / pow(phiabs, 2.0);
                                termy0dy = (3.0 * pow(phidxp, 2.0) * phidxpy * yxp + 3.0 * pow(phidyp, 2.0) * phidypy * yyp + 3.0 * pow(phidzp, 2.0) * phidzpy * yzp) / phiabs - (2.0 * phidx * phidxy + 2.0 * phidy * phidyy + 2.0 * phidz * phidyz) * (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / pow(phiabs, 2.0);
                                termz0dz = (3.0 * pow(phidxp, 2.0) * phidxpz * zxp + 3.0 * pow(phidyp, 2.0) * phidypz * zyp + 3.0 * pow(phidzp, 2.0) * phidzpz * zzp) / phiabs - (2.0 * phidx * phidxz + 2.0 * phidy * phidyz + 2.0 * phidz * phidzz) * (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / pow(phiabs, 2.0);

                                termx1dx = ((phidxx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidx * (4.0 * pow(phidxp, 3.0) * phidxpx + 4.0 * pow(phidyp, 3.0) * phidypx + 4.0 * pow(phidzp, 3.0) * phidzpx))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxx + phidy * phidxy + phidz * phidxz) * phidx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);
                                termy1dy = ((phidyy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidy * (4.0 * pow(phidxp, 3.0) * phidxpy + 4.0 * pow(phidyp, 3.0) * phidypy + 4.0 * pow(phidzp, 3.0) * phidzpy))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxy + phidy * phidyy + phidz * phidyz) * phidy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);
                                termz1dz = ((phidzz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidz * (4.0 * pow(phidxp, 3.0) * phidxpz + 4.0 * pow(phidyp, 3.0) * phidypz + 4.0 * pow(phidzp, 3.0) * phidzpz))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxz + phidy * phidyz + phidz * phidzz) * phidz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);

                                termx = 16.0 * epsilon0 * astre * (epdx * (termx0 - termx1) + ep * (termx0dx - termx1dx));
                                termy = 16.0 * epsilon0 * astre * (epdy * (termy0 - termy1) + ep * (termy0dy - termy1dy));
                                termz = 16.0 * epsilon0 * astre * (epdz * (termz0 - termz1) + ep * (termz0dz - termz1dz));

                                termjjkk = term0 + termx + termy + termz;
                            }
                            else
                            {
                                termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);
                            }

                            sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j][l];
                        }
                        pddtt += -2.0 * mij[ii][jj] / (double)phiNum[i][j][l] * (sum1 - 8.0 / PI * fij[ii][jj] * sqrt(phi[ii][i][j][l] * phi[jj][i][j][l]));
                        //フェーズフィールドの発展方程式[式(4.31)]
                    }
                    phi2[ii][i][j][l] = phi[ii][i][j][l] + pddtt * dtime; //フェーズフィールドの時間発展（陽解法）
                    if (phi2[ii][i][j][l] >= 1.0)
                    {
                        phi2[ii][i][j][l] = 1.0;
                    } //フェーズフィールドの変域補正
                    if (phi2[ii][i][j][l] <= 0.0)
                    {
                        phi2[ii][i][j][l] = 0.0;
                    }
                }
            } // j
        }     // i
    }

    for (k = 0; k <= nm; k++)
    {
        for (i = 0; i <= ndmx; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (l = 0; l <= ndmz; l++)
                {
                    phi[k][i][j][l] = phi2[k][i][j][l];
                }
            }
        }
    }

    //
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (l = 0; l <= ndmz; l++)
            {
                sum1 = 0.0;
                for (k = 0; k <= nm; k++)
                {
                    sum1 += phi[k][i][j][l];
                }
                for (k = 0; k <= nm; k++)
                {
                    phi[k][i][j][l] = phi[k][i][j][l] / sum1;
                }
            }
        }
    }

    istep = istep + 1.;
    if (istep < nstep)
    {
        goto start;
    }

end:;
    return 0;
}
