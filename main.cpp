// Adapted from Prof. Koyama's MPF code in his textbook
// Author: Chuanqi Zhu
// Created on: 2022/2/16
// Updated on 2022/04/20

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include "CImg.h" //CImg ライブラリ（描画用）使用のためのヘッダ

using namespace cimg_library;

#define NDX 100 //差分計算における計算領域一辺の分割数
#define NDY 1
#define NDZ 1

#define N 2

int ndmx = NDX - 1;
int ndmy = NDY - 1; //計算領域の一辺の差分分割数(差分ブロック数), ND-1を定義
int ndmz = NDZ - 1;
int nm = N - 1;
double PI = 3.141592; //π、計算カウント数

int phinum;

CImg<unsigned char> phi_fldxy(NDX, NDY, 1, 3);
char outFilePhi_xy[64];

int i, j, k, ni, ii, jj, kk, ll, it; //整数
int ip, im, jp, jm, kp, km;          //整数
int n1, n2, n3;                      //整数

int istep = 0;
// int n000;		//位置(i,j)において、pが０ではない方位の個数（n00>=n000）
int nstep; //計算カウント数の最大値（計算終了カウント）
int pstep;
double dtime, L, dx; // L計算領域の一辺の長さ(nm), 差分プロック１辺の長さ(m)
double M0;           //粒界の易動度
double W0;           //ペナルティー項の係数
double A0;           //勾配エネルギー係数
double S0;           //粒界移動の駆動力
double F0;
double temp0;            //温度
double sum1, sum2, sump; //各種の和の作業変数
double pddtt;            //フェーズフィールドの時間変化率

double Tm, dH;
double sph_s, kap_s, sph_l, kap_l;
double Tg, Tv, Tr, T_left, T_right;

double Cdt_s = kap_s / sph_s;
double Cdt_l = kap_l / sph_l;

double gamma0; //粒界エネルギ密度
double delta;  //粒界幅（差分ブロック数にて表現）
double mobi;   //粒界の易動度
double vm0;    //モル体積

double astre, astrem;
double th, vp, eta;
double thii, vpii, etaii;
double thetax, thetay;
double epsilon0;
double termiikk, termjjkk;
double miijj;

double phidx, phidy, phidz;
double phidxii, phidyii, phidzii;
double phidxx, phidyy, phidzz;
double phidxy, phidxz, phidyz;
double phiabs, phiabsii;

double xxp, xyp, xzp, yxp, yyp, yzp, zxp, zyp, zzp;
double xxpii, xypii, xzpii, yxpii, yypii, yzpii, zxpii, zypii, zzpii;

double phidxp, phidyp, phidzp;
double phidxpii, phidypii, phidzpii;
double phidxpx, phidypx, phidzpx;
double phidxpy, phidypy, phidzpy;
double phidxpz, phidypz, phidzpz;
double ep, epdx, epdy, epdz;
double term0;
double termx, termx0, termx1, termx0dx, termx1dx;
double termy, termy0, termy1, termy0dy, termy1dy;
double termz, termz0, termz1, termz0dz, termz1dz;

double tempip, tempim, Tddtt;
double sumplane, inttemp, intvel;
int hasS, allS;
int intpos, frapass, curpos, prepos, curst, prest;

int x11, y11, x1h[10], y1h[10]; //初期核の座標
double t, r0, r;

double ****phi, ****phi2;
double ***temp, ***temp2;
int ***phiNum, ****phiIdx;
double **aij, **wij, **mij;
double **anij, **thij, **vpij, **etaij;

//******* メインプログラム ******************************************
int main(int argc, char *argv[])
{
    phi = new double ***[N];
    phi2 = new double ***[N];
    for (ni = 0; ni <= nm; ni++)
    {
        phi[ni] = new double **[NDX];
        phi2[ni] = new double **[NDX];
        for (i = 0; i <= ndmx; i++)
        {
            phi[ni][i] = new double *[NDY];
            phi2[ni][i] = new double *[NDY];
            for (j = 0; j <= ndmy; j++)
            {
                phi[ni][i][j] = new double[NDZ];
                phi2[ni][i][j] = new double[NDZ];
            }
        }
    }

    phiIdx = new int ***[N + 1];
    for (ni = 0; ni <= N; ni++)
    {
        phiIdx[ni] = new int **[NDX];
        for (i = 0; i <= ndmx; i++)
        {
            phiIdx[ni][i] = new int *[NDY];
            for (j = 0; j <= ndmy; j++)
            {
                phiIdx[ni][i][j] = new int[NDZ];
            }
        }
    }

    temp = new double **[NDX];
    temp2 = new double **[NDX];
    for (i = 0; i <= ndmx; i++)
    {
        temp[i] = new double *[NDY];
        temp2[i] = new double *[NDY];
        for (j = 0; j <= ndmy; j++)
        {
            temp[i][j] = new double[NDZ];
            temp2[i][j] = new double[NDZ];
        }
    }

    phiNum = new int **[NDX];
    for (i = 0; i <= ndmx; i++)
    {
        phiNum[i] = new int *[NDY];

        for (j = 0; j <= ndmy; j++)
        {
            phiNum[i][j] = new int[NDZ];
        }
    }

    aij = new double *[N];
    wij = new double *[N];
    mij = new double *[N];
    anij = new double *[N];
    thij = new double *[N];
    vpij = new double *[N];
    etaij = new double *[N];
    for (ni = 0; ni <= nm; ni++)
    {
        aij[ni] = new double[N];
        wij[ni] = new double[N];
        mij[ni] = new double[N];
        anij[ni] = new double[N];
        thij[ni] = new double[N];
        vpij[ni] = new double[N];
        etaij[ni] = new double[N];
    }

    nstep = 1000001;
    pstep = 100000;

    dx = 1.0e-5;
    dtime = 1.0e-6;
    delta = 5.0 * dx;
    mobi = 1.68e-9;

    gamma0 = 0.5;
    astre = 0.03;
    astrem = 0.26;

    Tm = 1687.0;
    dH = 4.122e9;
    temp0 = 1686.95;

    Tg = 8.0e3;
    Tv = 1.5e-4;
    Tr = Tg * Tv;
    T_left = 1683.4 - NDX / 4 * dx * Tg;
    T_right = T_left + Tg * NDX * dx;

    sph_s = 2.29e6;
    kap_s = 22.0;
    sph_l = 2.53e6;
    kap_l = 54.0;
    Cdt_s = kap_s / sph_s;
    Cdt_l = kap_l / sph_l;

    A0 = 8.0 * delta * gamma0 / PI / PI; //勾配エネルギー係数[式(4.40)]
    W0 = 4.0 * gamma0 / delta;           //ペナルティー項の係数[式(4.40)]
    M0 = mobi * PI * PI / (8.0 * delta); //粒界の易動度[式(4.40)]

    for (i = 0; i <= nm; i++)
    {
        for (j = 0; j <= nm; j++)
        {
            wij[i][j] = W0;
            aij[i][j] = A0;
            mij[i][j] = M0;
            anij[i][j] = 0;
            thij[i][j] = 0.0;
            vpij[i][j] = 0.0;
            etaij[i][j] = 0.0;
            if ((i == 0) || (j == 0))
            {
                anij[i][j] = 0;
            }
            if (i == j)
            {
                wij[i][j] = 0.0;
                aij[i][j] = 0.0;
                mij[i][j] = 0.0;
                anij[i][j] = 0;
            }
        }
    }
    // thij[1][0] = PI / 4.0;
    // thij[0][1] = PI / 4.0;
    // vpij[1][0] = PI / 8.0;
    // vpij[0][1] = PI / 8.0;
    // vpij[2][0] = PI / 4.0;
    // vpij[0][2] = PI / 4.0;
    // etaij[1][0] = PI / 4.0;
    // etaij[0][1] = PI / 4.0;

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                if (i < NDX / 4)
                // if (((i - NDX / 2) * (i - NDX / 2) + (j - NDY / 2) * (j - NDY / 2) + (k - NDZ / 2) * (k - NDZ / 2)) < (NDX / 8) * (NDX / 8))
                {
                    phi[1][i][j][k] = 1.0;
                    phi[0][i][j][k] = 0.0;
                }
                else
                {
                    phi[1][i][j][k] = 0.0;
                    phi[0][i][j][k] = 1.0;
                }
                temp[i][j][k] = T_left + +i * dx * Tg;
            }
        }
    }
    // }

start:;

    if ((((int)(istep) % pstep) == 0))
    {
        std::cout << "interface position is: " << intpos << std::endl;
        std::cout << "interface temperature is: " << inttemp << std::endl;
        std::cout << "interface veocity is: " << intvel << std::endl;
        // ****** XY *******
        cimg_forXY(phi_fldxy, x, y)
        {
            phi_fldxy(x, y, 0) = 255. * (phi[1][x][y][NDZ / 2]); // red
            phi_fldxy(x, y, 1) = 255. * (phi[1][x][y][NDZ / 2]); // green
            phi_fldxy(x, y, 2) = 255. * (phi[1][x][y][NDZ / 2]); // blue
        }
        sprintf(outFilePhi_xy, "figures/phi/2dxy%d.png", istep);
        phi_fldxy.save_jpeg(outFilePhi_xy);

        FILE *stream; //ストリームのポインタ設定
        char buffer[30];
        sprintf(buffer, "data/phi/1d%d.csv", istep);
        stream = fopen(buffer, "a"); //書き込む先のファイルを追記方式でオープン

        for (i = 0; i <= ndmx; i++)
        {
            fprintf(stream, "%e   ", phi[1][i][0][0]);
            fprintf(stream, "\n");
        }
        fclose(stream); //ファイルをクローズ

        FILE *streamt; //ストリームのポインタ設定
        char buffert[30];
        sprintf(buffert, "data/temp/1d%d.csv", istep);
        streamt = fopen(buffert, "a"); //書き込む先のファイルを追記方式でオープン

        for (i = 0; i <= ndmx; i++)
        {
            fprintf(streamt, "%e   ", (temp[i][0][0] - Tm));
            fprintf(streamt, "\n");
        }
        fclose(streamt); //ファイルをクローズ

        // FILE *stream;
        // char buffer[30];
        // sprintf(buffer, "data/phi/3d%d.vtk", istep);
        // stream = fopen(buffer, "a");

        // fprintf(stream, "# vtk DataFile Version 1.0\n");
        // fprintf(stream, "phi_%d.vtk\n", istep);
        // fprintf(stream, "ASCII\n");
        // fprintf(stream, "DATASET STRUCTURED_POINTS\n");
        // fprintf(stream, "DIMENSIONS %d %d %d\n", NDX, NDY, NDZ);
        // fprintf(stream, "ORIGIN 0.0 0.0 0.0\n");
        // fprintf(stream, "ASPECT_RATIO 1.0 1.0 1.0\n");
        // fprintf(stream, "\n");
        // fprintf(stream, "POINT_DATA %d\n", NDX * NDY * NDZ);
        // fprintf(stream, "SCALARS scalars float\n");
        // fprintf(stream, "LOOKUP_TABLE default\n");

        // for (k = 0; k <= ndmz; k++)
        // {
        //     for (j = 0; j <= ndmy; j++)
        //     {
        //         for (i = 0; i <= ndmx; i++)
        //         {
        //             // sump = 0.0;
        //             // for (ii = 0; ii <= nm; ii++)
        //             // {
        //             //     sump += phi[ii][i][j][k] * phi[ii][i][j][k];
        //             // }
        //             fprintf(stream, "%e\n", phi[1][i][j][k]);
        //         }
        //     }
        // }
        // fclose(stream);
    }

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                kp = k + 1;
                km = k - 1;
                if (i == ndmx)
                {
                    ip = ndmx;
                }
                if (i == 0)
                {
                    im = 0;
                }
                if (j == ndmy)
                {
                    jp = 0;
                }
                if (j == 0)
                {
                    jm = ndmy;
                }
                if (k == ndmz)
                {
                    kp = 0;
                }
                if (k == 0)
                {
                    km = ndmz;
                }

                //--- 位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の個数---
                phinum = 0;
                for (ii = 0; ii <= nm; ii++)
                {
                    if ((phi[ii][i][j][k] > 0.0) ||
                        ((phi[ii][i][j][k] == 0.0) && (phi[ii][ip][j][k] > 0.0) ||
                         (phi[ii][im][j][k] > 0.0) ||
                         (phi[ii][i][jp][k] > 0.0) ||
                         (phi[ii][i][jm][k] > 0.0) ||
                         (phi[ii][i][j][kp] > 0.0) ||
                         (phi[ii][i][j][km] > 0.0)))
                    {
                        phinum++;
                        phiIdx[phinum][i][j][k] = ii;
                    }
                }
                phiNum[i][j][k] = phinum;
            }
        }
    }

    // Evolution Equations
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                kp = k + 1;
                km = k - 1;
                if (i == ndmx)
                {
                    ip = ndmx;
                }
                if (i == 0)
                {
                    im = 0;
                }
                if (j == ndmy)
                {
                    jp = 0;
                }
                if (j == 0)
                {
                    jm = ndmy;
                }
                if (k == ndmz)
                {
                    kp = 0;
                }
                if (k == 0)
                {
                    km = ndmz;
                }

                for (n1 = 1; n1 <= phiNum[i][j][k]; n1++)
                {
                    ii = phiIdx[n1][i][j][k];

                    phidxii = (phi[ii][ip][j][k] - phi[ii][im][j][k]) / 2.0 / dx;
                    phidyii = (phi[ii][i][jp][k] - phi[ii][i][jm][k]) / 2.0 / dx;
                    phidzii = (phi[ii][i][j][kp] - phi[ii][i][j][km]) / 2.0 / dx;
                    phiabsii = phidxii * phidxii + phidyii * phidyii + phidzii * phidzii;

                    pddtt = 0.0;
                    for (n2 = 1; n2 <= phiNum[i][j][k]; n2++)
                    {
                        jj = phiIdx[n2][i][j][k];
                        sum1 = 0.0;
                        for (n3 = 1; n3 <= phiNum[i][j][k]; n3++)
                        {
                            kk = phiIdx[n3][i][j][k];

                            // calculate the interface normal and deirivatives of the phase field
                            phidx = (phi[kk][ip][j][k] - phi[kk][im][j][k]) / 2.0 / dx;
                            phidy = (phi[kk][i][jp][k] - phi[kk][i][jm][k]) / 2.0 / dx;
                            phidz = (phi[kk][i][j][kp] - phi[kk][i][j][km]) / 2.0 / dx;

                            phidxx = (phi[kk][ip][j][k] + phi[kk][im][j][k] - 2.0 * phi[kk][i][j][k]) / dx / dx;
                            phidyy = (phi[kk][i][jp][k] + phi[kk][i][jm][k] - 2.0 * phi[kk][i][j][k]) / dx / dx;
                            phidzz = (phi[kk][i][j][kp] + phi[kk][i][j][km] - 2.0 * phi[kk][i][j][k]) / dx / dx;

                            phidxy = (phi[kk][ip][jp][k] + phi[kk][im][jm][k] - phi[kk][im][jp][k] - phi[kk][ip][jm][k]) / 4.0 / dx / dx;
                            phidxz = (phi[kk][ip][j][kp] + phi[kk][im][j][km] - phi[kk][im][j][kp] - phi[kk][ip][j][km]) / 4.0 / dx / dx;
                            phidyz = (phi[kk][i][jp][kp] + phi[kk][i][jm][km] - phi[kk][i][jm][kp] - phi[kk][i][jp][km]) / 4.0 / dx / dx;

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
                            sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j][k];
                        }

                        thii = thij[ii][jj];
                        vpii = vpij[ii][jj];
                        etaii = etaij[ii][jj];

                        xxpii = cos(thii) * cos(vpii);
                        yxpii = sin(thii) * cos(vpii);
                        zxpii = sin(vpii);
                        xypii = -sin(thii) * cos(etaii) - cos(thii) * sin(vpii) * sin(etaii);
                        yypii = cos(thii) * cos(etaii) - sin(thii) * sin(vpii) * sin(etaii);
                        zypii = cos(vpii) * sin(etaii);
                        xzpii = sin(etaii) * sin(thii) - cos(etaii) * cos(thii) * sin(vpii);
                        yzpii = -sin(etaii) * cos(thii) - cos(etaii) * sin(thii) * sin(vpii);
                        zzpii = cos(etaii) * cos(vpii);

                        phidxpii = phidxii * xxpii + phidyii * yxpii + phidzii * zxpii;
                        phidypii = phidxii * xypii + phidyii * yypii + phidzii * zypii;
                        phidzpii = phidxii * xzpii + phidyii * yzpii + phidzii * zzpii;

                        if (anij[ii][jj] == 1 && phiabsii != 0.0)
                        {
                            miijj = mij[ii][jj] * (1.0 - 3.0 * astrem + 4.0 * astrem * (pow(phidxpii, 4.0) + pow(phidypii, 4.0) + pow(phidzpii, 4.0)) / pow(phiabsii, 2.0));
                        }
                        else
                        {
                            miijj = mij[ii][jj];
                        }
                        if (ii == 1 && jj == 0)
                        {
                            F0 = -(temp[i][j][k] - Tm) * dH / Tm;
                        }
                        else if (ii == 0 && jj == 1)
                        {
                            F0 = (temp[i][j][k] - Tm) * dH / Tm;
                        }
                        pddtt += -2.0 * miijj / (double)phiNum[i][j][k] * (sum1 - 8.0 / PI * F0 * sqrt(phi[ii][i][j][k] * phi[jj][i][j][k]));
                        //フェーズフィールドの発展方程式[式(4.31)]
                    }
                    phi2[ii][i][j][k] = phi[ii][i][j][k] + pddtt * dtime; //フェーズフィールドの時間発展（陽解法）
                    if (phi2[ii][i][j][k] >= 1.0)
                    {
                        phi2[ii][i][j][k] = 1.0;
                    } //フェーズフィールドの変域補正
                    if (phi2[ii][i][j][k] <= 0.0)
                    {
                        phi2[ii][i][j][k] = 0.0;
                    }
                    // termperature increase from release of latent heat
                    if (ii == 1)
                    {
                        temp[i][j][k] += pddtt * dtime * 30.0 * pow(phi[1][i][j][k], 2.0) * pow((1.0 - phi[1][i][j][k]), 2.0) * dH / sph_s;
                    }
                }
            } // j
        }     // i
    }

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                for (ii = 0; ii <= nm; ii++)
                {
                    phi[ii][i][j][k] = phi2[ii][i][j][k];
                }
            }
        }
    }

    //
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                sum1 = 0.0;
                for (ii = 0; ii <= nm; ii++)
                {
                    sum1 += phi[ii][i][j][k];
                }
                for (ii = 0; ii <= nm; ii++)
                {
                    phi[ii][i][j][k] = phi[ii][i][j][k] / sum1;
                }
            }
        }
    }

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                kp = k + 1;
                km = k - 1;
                // right boundary
                if (i == ndmx)
                {
                    tempip = T_right;
                }
                else
                {
                    tempip = temp[ip][j][k];
                }
                // left boundary
                if (i == 0)
                {
                    tempim = T_left;
                }
                else
                {
                    tempim = temp[im][j][k];
                }
                if (j == ndmy)
                {
                    jp = 0;
                }
                if (j == 0)
                {
                    jm = ndmy;
                }
                if (k == ndmz)
                {
                    kp = 0;
                }
                if (k == 0)
                {
                    km = ndmz;
                }
                Tddtt = (Cdt_l * phi[0][i][j][k] + Cdt_s * phi[1][i][j][k]) * (tempip + tempim + temp[i][jp][k] + temp[i][jm][k] + temp[i][j][kp] + temp[i][j][km] - 6.0 * temp[i][j][k]) / dx / dx;
                temp2[i][j][k] = temp[i][j][k] + Tddtt * dtime;
            }
        }
    }

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                temp[i][j][k] = temp2[i][j][k];
            }
        }
    }

    T_left -= Tr * dtime;
    T_right -= Tr * dtime;

    sumplane = 0.0;
    hasS = 0;
    allS = 0;
    // check if the bottom is solid
    for (j = 0; j <= ndmy; j++)
    {
        for (k = 0; k <= ndmz; k++)
        {
            if (phi[0][0][j][k] != 0.0)
            {
                hasS = 0;
            }
            sumplane += phi[0][0][j][k];
            if ((sumplane == 0.0) && (j == ndmy) && (k == ndmz))
            {
                hasS = 1;
            }
        }
    }
    // search interface front
    if (hasS == 1)
    {
        allS = 1;
        for (i = 0; i <= ndmx; i++)
        {
            if (allS == 0)
            {
                intpos = i - 1;
                inttemp = temp[intpos][0][0];
                break;
            }
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    if (phi[0][i][j][k] > 0.0)
                    {
                        allS = 0;
                        break;
                    }
                }
                if (allS == 0)
                {
                    break;
                }
            }
        }
    }

    if (intpos > (NDX / 4))
    {
        frapass += 1;
        curpos = frapass;
        curst = istep;
        intvel = double(curpos - prepos) * dx / double(curst - prest) / dtime;
        prepos = curpos;
        prest = curst;
        for (i = 0; i <= ndmx - 1; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    for (ii = 0; ii <= nm; ii++)
                    {
                        phi[ii][i][j][k] = phi[ii][i + 1][j][k];
                    }
                    temp[i][j][k] = temp[i + 1][j][k];
                }
            }
        }
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                for (ii = 0; ii <= nm; ii++)
                {
                    phi[ii][ndmx][j][k] = phi[ii][ndmx - 1][j][k];
                }
                temp[ndmx][j][k] = temp[ndmx - 1][j][k] + Tg * dx;
            }
        }
        T_left += Tg * dx;
        T_right += Tg * dx;
    }

    istep = istep + 1.;
    if (istep < nstep)
    {
        goto start;
    }

end:;
    return 0;
}
