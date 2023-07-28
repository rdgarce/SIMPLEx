#include "SIMPLEx_test.h"

// Test 1 definition
#define PLtest1_m 4
#define PLtest1_n 6
double PLtest1_AT[PLtest1_n][PLtest1_m] = {
                                            {-5, 0, 3, 2},
                                            {4, 1, 4, 1},
                                            {1, 0, 0, 0},
                                            {0, 1, 0, 0},
                                            {0, 0, 1, 0},
                                            {0, 0, 0, 1}
                                        };
double PLtest1_b[PLtest1_m] = {20, 6.5, 36, 16};
uint16_t PLtest1_nbIDx[PLtest1_n - PLtest1_m] ={0, 1};
uint16_t PLtest1_bIDx[PLtest1_m] ={2, 3, 4, 5};
double PLtest1_Zcnb[PLtest1_n - PLtest1_m] = {-11, -10};
double PLtest1_Zcb[PLtest1_m] = {0, 0, 0, 0};

#define PLtest1_Wcnb NULL
#define PLtest1_Wcb NULL
#define PLtest1_Wb 0

uint16_t PLtest1_bIDx_opt[PLtest1_m] ={2, 3, 1, 0};
double PLtest1_b_opt[PLtest1_m] = {28.8, 1.7, 4.8, 5.6};
#define PLtest1_Zopt -109.6

// Test 2 definition
#define PLtest2_m 4
#define PLtest2_n 8
double PLtest2_AT[PLtest2_n][PLtest2_m] = {
                                            {1,0,2,4},
                                            {0,1,5,2},
                                            {1,0,0,0},
                                            {0,1,0,0},
                                            {0,0,-1,0},
                                            {0,0,0,-1},
                                            {0,0,1,0},
                                            {0,0,0,1},
                                        };
double PLtest2_b[PLtest2_m] = {6,4,10,8};
uint16_t PLtest2_nbIDx[PLtest2_n - PLtest2_m] = {0,1,4,5};
uint16_t PLtest2_bIDx[PLtest2_m] = {2,3,6,7};
double PLtest2_Zcnb[PLtest2_n - PLtest2_m] = {-1,-1,0,0};
double PLtest2_Zcb[PLtest2_m] = {0,0,0,0};

double PLtest2_Wcnb[PLtest2_n - PLtest2_m] = {-6,-7,1,1};
double PLtest2_Wcb[PLtest2_m] = {0,0,0,0};
#define PLtest2_Wb -18

uint16_t PLtest2_bIDx_opt[PLtest2_m] ={5, 4, 1, 0};
double PLtest2_b_opt[PLtest2_m] = {24, 22, 4, 6};
#define PLtest2_Zopt -10

// Test 3 definition
#define PLtest3_m 2
#define PLtest3_n 4
double PLtest3_AT[PLtest3_n][PLtest3_m] = {
                                            {1, 1},
                                            {2, 1.0/3.0},
                                            {1, 0},
                                            {0, 1},
                                        };
double PLtest3_b[PLtest3_m] = {4, 2};
uint16_t PLtest3_nbIDx[PLtest3_n - PLtest3_m] ={0, 1};
uint16_t PLtest3_bIDx[PLtest3_m] ={2, 3};
double PLtest3_Zcnb[PLtest3_n - PLtest3_m] = {-5, -4};
double PLtest3_Zcb[PLtest3_m] = {0, 0};

#define PLtest3_Wcnb NULL
#define PLtest3_Wcb NULL
#define PLtest3_Wb 0

uint16_t PLtest3_bIDx_opt[PLtest3_m] ={1, 0};
double PLtest3_b_opt[PLtest3_m] = {6.0/5.0, 8.0/5.0};
#define PLtest3_Zopt -64.0/5.0

// Test 4 definition
#define PLtest4_m 2
#define PLtest4_n 6
double PLtest4_AT[PLtest4_n][PLtest4_m] = {
                                            {4, -5},
                                            {3, -3},
                                            {0, 1},
                                            {2, 0},
                                            {1, 0},
                                            {0, 1},
                                        };
double PLtest4_b[PLtest4_m] = {2, 1};
uint16_t PLtest4_nbIDx[PLtest4_n - PLtest4_m] ={0, 1, 2, 3};
uint16_t PLtest4_bIDx[PLtest4_m] ={4, 5};
double PLtest4_Zcnb[PLtest4_n - PLtest4_m] = {7, 2, -5, -1};
double PLtest4_Zcb[PLtest4_m] = {0, 0};

#define PLtest4_Wcnb NULL
#define PLtest4_Wcb NULL
#define PLtest4_Wb 0

uint16_t PLtest4_bIDx_opt[PLtest4_m] ={0, 2};
double PLtest4_b_opt[PLtest4_m] = {1.0/2.0, 7.0/2.0};
#define PLtest4_Zopt -14

// Test 5 definition
#define PLtest5_m 3
#define PLtest5_n 7
double PLtest5_AT[PLtest5_n][PLtest5_m] = {
                                            {1.0/4.0, 1.0/2.0, 0},
                                            {-8, -12, 0},
                                            {-1, -1.0/2.0, 1},
                                            {9, 3, 0},
                                            {1, 0, 0},
                                            {0, 1, 0},
                                            {0, 0, 1},
                                        };
double PLtest5_b[PLtest5_m] = {0, 0, 1};
uint16_t PLtest5_nbIDx[PLtest5_n - PLtest5_m] ={0, 1, 2, 3};
uint16_t PLtest5_bIDx[PLtest5_m] ={4, 5, 6};
double PLtest5_Zcnb[PLtest5_n - PLtest5_m] = {-3.0/4.0, 20, -1.0/2.0, 6};
double PLtest5_Zcb[PLtest5_m] = {0, 0, 0};

#define PLtest5_Wcnb NULL
#define PLtest5_Wcb NULL
#define PLtest5_Wb 0

uint16_t PLtest5_bIDx_opt[PLtest5_m] ={2, 4, 0};
double PLtest5_b_opt[PLtest5_m] = {1, 3.0/4.0, 1};
#define PLtest5_Zopt -1.25

// Test 6 definition
#define PLtest6_m 1
#define PLtest6_n 2
double PLtest6_AT[PLtest6_n][PLtest6_m] = {
                                            {49},
                                            {1},
                                        };
double PLtest6_b[PLtest6_m] = {57};
uint16_t PLtest6_nbIDx[PLtest6_n - PLtest6_m] ={0};
uint16_t PLtest6_bIDx[PLtest6_m] ={1};
double PLtest6_Zcnb[PLtest6_n - PLtest6_m] = {-83};
double PLtest6_Zcb[PLtest6_m] = {0};

double PLtest6_Wcnb[PLtest6_n - PLtest6_m] = {-49};
double PLtest6_Wcb[PLtest6_m] = {0};
#define PLtest6_Wb -57

uint16_t PLtest6_bIDx_opt[PLtest6_m] = {0};
double PLtest6_b_opt[PLtest6_m] = {1.16326530612244894};
#define PLtest6_Zopt -96.55102041

// Test 7 definition (Test not passed)
#define PLtest7_m 2
#define PLtest7_n 7
double PLtest7_AT[PLtest7_n][PLtest7_m] = {
                                            {3, 9},
                                            {10, 1},
                                            {1, 6},
                                            {1, 6},
                                            {4, 11},
                                            {1, 0},
                                            {0, 1},
                                        };
double PLtest7_b[PLtest7_m] = {7, 9};
uint16_t PLtest7_nbIDx[PLtest7_n - PLtest7_m] ={0, 1, 2, 3, 4};
uint16_t PLtest7_bIDx[PLtest7_m] ={5, 6};
double PLtest7_Zcnb[PLtest7_n - PLtest7_m] = {-8, -2, -6, -6, -11};
double PLtest7_Zcb[PLtest7_m] = {0, 0};

#define PLtest7_Wcnb NULL
#define PLtest7_Wcb NULL
#define PLtest7_Wb 0

uint16_t PLtest7_bIDx_opt[PLtest7_m] = {0, 1}; // Not yet known
double PLtest7_b_opt[PLtest7_m] = {1,1};       // Not yet known
#define PLtest7_Zopt -9.559322

// Inserimento dei test nell'array
PLtest_t PLtest_array[] = {
                                    {
                                        .AT = (double *)PLtest1_AT,
                                        .b = (double *)PLtest1_b,
                                        .b_opt = (double *)PLtest1_b_opt,
                                        .bIDx = PLtest1_bIDx,
                                        .bIDx_opt = PLtest1_bIDx_opt,
                                        .m = PLtest1_m,
                                        .n = PLtest1_n,
                                        .nbIDx = PLtest1_nbIDx,
                                        .Wb = PLtest1_Wb,
                                        .Wcb = (double *)PLtest1_Wcb,
                                        .Wcnb = (double *)PLtest1_Wcnb,
                                        .Zcb = (double *)PLtest1_Zcb,
                                        .Zcnb = (double *)PLtest1_Zcnb,
                                        .Zopt_opt = PLtest1_Zopt,
                                        .ret = 0
                                    },
                                    {
                                        .AT = (double *)PLtest2_AT,
                                        .b = (double *)PLtest2_b,
                                        .b_opt = (double *)PLtest2_b_opt,
                                        .bIDx = PLtest2_bIDx,
                                        .bIDx_opt = PLtest2_bIDx_opt,
                                        .m = PLtest2_m,
                                        .n = PLtest2_n,
                                        .nbIDx = PLtest2_nbIDx,
                                        .Wb = PLtest2_Wb,
                                        .Wcb = (double *)PLtest2_Wcb,
                                        .Wcnb = (double *)PLtest2_Wcnb,
                                        .Zcb = (double *)PLtest2_Zcb,
                                        .Zcnb = (double *)PLtest2_Zcnb,
                                        .Zopt_opt = PLtest2_Zopt,
                                        .ret = 0
                                    },
                                    {
                                        .AT = (double *)PLtest3_AT,
                                        .b = (double *)PLtest3_b,
                                        .b_opt = (double *)PLtest3_b_opt,
                                        .bIDx = PLtest3_bIDx,
                                        .bIDx_opt = PLtest3_bIDx_opt,
                                        .m = PLtest3_m,
                                        .n = PLtest3_n,
                                        .nbIDx = PLtest3_nbIDx,
                                        .Wb = PLtest3_Wb,
                                        .Wcb = (double *)PLtest3_Wcb,
                                        .Wcnb = (double *)PLtest3_Wcnb,
                                        .Zcb = (double *)PLtest3_Zcb,
                                        .Zcnb = (double *)PLtest3_Zcnb,
                                        .Zopt_opt = PLtest3_Zopt,
                                        .ret = 0
                                    },
                                    {
                                        .AT = (double *)PLtest4_AT,
                                        .b = (double *)PLtest4_b,
                                        .b_opt = (double *)PLtest4_b_opt,
                                        .bIDx = PLtest4_bIDx,
                                        .bIDx_opt = PLtest4_bIDx_opt,
                                        .m = PLtest4_m,
                                        .n = PLtest4_n,
                                        .nbIDx = PLtest4_nbIDx,
                                        .Wb = PLtest4_Wb,
                                        .Wcb = (double *)PLtest4_Wcb,
                                        .Wcnb = (double *)PLtest4_Wcnb,
                                        .Zcb = (double *)PLtest4_Zcb,
                                        .Zcnb = (double *)PLtest4_Zcnb,
                                        .Zopt_opt = PLtest4_Zopt,
                                        .ret = 0
                                    },
                                    {
                                        .AT = (double *)PLtest5_AT,
                                        .b = (double *)PLtest5_b,
                                        .b_opt = (double *)PLtest5_b_opt,
                                        .bIDx = PLtest5_bIDx,
                                        .bIDx_opt = PLtest5_bIDx_opt,
                                        .m = PLtest5_m,
                                        .n = PLtest5_n,
                                        .nbIDx = PLtest5_nbIDx,
                                        .Wb = PLtest5_Wb,
                                        .Wcb = (double *)PLtest5_Wcb,
                                        .Wcnb = (double *)PLtest5_Wcnb,
                                        .Zcb = (double *)PLtest5_Zcb,
                                        .Zcnb = (double *)PLtest5_Zcnb,
                                        .Zopt_opt = PLtest5_Zopt,
                                        .ret = 0
                                    },
                                    {
                                        .AT = (double *)PLtest6_AT,
                                        .b = (double *)PLtest6_b,
                                        .b_opt = (double *)PLtest6_b_opt,
                                        .bIDx = PLtest6_bIDx,
                                        .bIDx_opt = PLtest6_bIDx_opt,
                                        .m = PLtest6_m,
                                        .n = PLtest6_n,
                                        .nbIDx = PLtest6_nbIDx,
                                        .Wb = PLtest6_Wb,
                                        .Wcb = (double *)PLtest6_Wcb,
                                        .Wcnb = (double *)PLtest6_Wcnb,
                                        .Zcb = (double *)PLtest6_Zcb,
                                        .Zcnb = (double *)PLtest6_Zcnb,
                                        .Zopt_opt = PLtest6_Zopt,
                                        .ret = 0
                                    }
                                };

#define N_TESTS sizeof(PLtest_array)/sizeof(PLtest_array[0])