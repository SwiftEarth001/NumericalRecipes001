#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

// the mantissa and exponent sizes for
//   the define statements below are
//   for presumed decimal storage;
#define MANTISSA_SIZE 5
// #define EXPONENT_SIZE 4
#define BASE 10

const static unsigned long int mantissa_bound = pow((BASE), (MANTISSA_SIZE));


//---------------------------------------------------------
// 1. floating point storage
//---------------------------------------------------------

// float base 10 storage;
// the mantissa and exponent sizes for
//   the bit field below are for
//   binary storage of decimal numbers;
// normalized mantissa;
// unsigned unbiased exponent;
typedef struct FLOAT_A {
    unsigned int mantissa : 31;
    unsigned int exponent : 8;
    unsigned int sign : 1;
} FLOAT_A;

// single precision float base 2;
// normalized mantissa with implicit leading bit;
// exponent bias -127;
typedef union FLOAT_CAST {
  float f;
  struct {
    unsigned int mantissa : 23;  // mantissa + (unsigned int)(0x1000000)
    unsigned int exponent : 8;  // exponent - 127
    unsigned int sign : 1;
  } parts;
} float_cast;

// double precision float base 2;
// normalized mantissa with implicit leading bit;
// exponent bias -1023;
typedef union DOUBLE_CAST {
  double df;
  struct {
    unsigned long long int mantissa : 52;  // mantissa + (unsigned long long int)(0x20000000000000)
    unsigned long long int exponent : 11;  // exponent - 1023
    unsigned long long int sign : 1;
  } parts;
} double_cast;


//---------------------------------------------------------
// 2. function headers and inline functions
//---------------------------------------------------------

/*
 * @function    translate
 * @discussion  stores binary single precision floating-point
 *              value x as a FLOAT_A "decimal" precision
 *              floating-point
 * @return      pointer to a FLOAT_A
*/
FLOAT_A* translate(float x);

/*
 * @function    addition
 * @discussion  adds two FLOAT_A values and returns a
 *              sum FLOAT_A
 * @return      pointer to a FLOAT_A
*/
FLOAT_A* addition(FLOAT_A &a1, FLOAT_A &a2);

/*
 * @function    accumulation
 * @discussion  translates an array of binary single-
 *              precision floats to FLOAT_A, and
 *              returns their cumulative sum
 * @return      pointer to a FLOAT_A
*/
FLOAT_A* accumulation(float x[], int size);

inline auto print_FLOAT_A(const FLOAT_A &flx) {
    printf("sign = %u\n", flx.sign);
    printf("exponent = %d\n", (signed int) ((signed char) flx.exponent));
    printf("mantissa = %u\n", flx.mantissa);   
}

inline auto FLOAT_A_to_double(const FLOAT_A &flx) -> double {
    double_cast fd;
    double x = (double)(flx.mantissa) * pow((BASE), (signed int)((signed char)flx.exponent));
    x = (flx.sign == 1) ? -x : x;
    fd.df = x;
    return (fd.df);
}

inline auto aggregate_floats(const float* fv, int size) -> float {
    float sum = fv[0];
    for (int j=1; j < size; j++) {
        sum = sum + fv[j];
    }
    return (sum);
}

inline auto round_to_nearest(int r, unsigned int &m, unsigned int &e) {
    m = (r >= 5) ? m+1 : m;
    if (m == mantissa_bound) {
        m = m / (BASE);
        e++;
    }
}


//---------------------------------------------------------
// 3. function definitions
//---------------------------------------------------------

FLOAT_A* translate(float x) {
    FLOAT_A* flx = new FLOAT_A;
    float_cast x1;
    x1.f = x;
    flx->sign = x1.parts.sign;

    float a = abs(x);
    unsigned int m, e=0;
    int remainder;

    if (a > mantissa_bound) {
        m = (unsigned int) a;
        while ( m >= mantissa_bound ) {
            remainder = m % (BASE);
            m = m / (BASE);
            e++;
        }
        round_to_nearest(remainder, m, e);
    } else {
        while ( a < mantissa_bound / (BASE) ) {
            a = a * (BASE);
            e--;
        }
        m = (unsigned int) a;
        remainder = (unsigned int) ((BASE) * (a - m));
        round_to_nearest(remainder, m, e);
    }

    flx->mantissa = m;
    flx->exponent = e;
    return (flx);
}

FLOAT_A* addition(FLOAT_A &a1, FLOAT_A &a2) {
    unsigned int e1 = a1.exponent;
    unsigned long long int m1 = a1.mantissa;
    
    unsigned int e2 = a2.exponent;
    unsigned long long int m2 = a2.mantissa;
    
    signed char shift = e1 - e2;
    unsigned int e;
    if (shift >= 0) {
        m1 = m1 * (unsigned long long int)(pow((BASE), shift));
        e = e2;
    } else {
        m2 = m2 * (unsigned long long int)(pow((BASE), -shift));
        e = e1;
    }
    unsigned long long int m3;
    if ((a1.sign == 0) && (a2.sign == 0)) {
        m3 = m1 + m2;
    } else if ((a1.sign == 1) && (a2.sign == 0)) {
        m3 = -m1 + m2;
    } else if ((a1.sign == 0) && (a2.sign == 1)) {
        m3 = m1 - m2;
    } else if ((a1.sign == 1) && (a2.sign == 1)) {
        m3 = -m1 - m2;
    }

    unsigned int m, sign = ((signed long long int)(m3) < 0) ? 1 : 0;
    m3 = abs((signed long long int)(m3));
    int remainder;
    if (m3 > mantissa_bound) {
        while ( m3 >= mantissa_bound ) {
            remainder = m3 % (BASE);
            m3 = m3 / (BASE);
            e++;
        }
        m = (unsigned int) m3;
        round_to_nearest(remainder, m, e);
    } else {
        while ( m3 < mantissa_bound / (BASE) ) {
            m3 = m3 * (BASE);
            e--;
        }
        m = (unsigned int) m3;
    }

    FLOAT_A* flx = new FLOAT_A;
    flx->mantissa = m;
    flx->exponent = e;
    flx->sign = sign;

    return (flx);
}

FLOAT_A* accumulation(float x[], int size) {
    FLOAT_A* s;
    FLOAT_A* sum = translate(x[0]);
    for (int k = 1; k < size; k++) {
        s = translate(x[k]);
        sum = addition(*sum, *s);
    }
    return (sum);
}


//---------------------------------------------------------
// 4. main method
//---------------------------------------------------------

int main(int argc, char *argv[])
{
    std::string floats_file;
    int numrows, numcols;
    if( argc == 4 ) { 
        floats_file = argv[1];
        numrows = std::stoi(argv[2]);
        numcols = std::stoi(argv[3]);
    } else {
        std::cout << "Usage is \" ./floating_point.exe FloatsFile numrows numcols \" ";
        std::cout << "where FloatsFile is a space delimited tabular array of floats.\n";
        std::cout << "Using floats.txt as FloatsFile...\n" << std::endl;
        floats_file = "floats.txt";
        numrows = 100;
        numcols = 100;
    }
    std::ifstream inFile(floats_file);
    std::string line;
    std::vector< std::vector<float> > all_floats;
    while ( getline(inFile, line) ) {
        std::istringstream is(line);
        all_floats.push_back(
            std::vector<float>( std::istream_iterator<float>(is),
                                std::istream_iterator<float>() ) 
                            );
    }
    inFile.close();
    float** randfloats = new float*[numrows];
    for (int i=0; i < numrows; i++) {
        randfloats[i] = new float[numcols];
        for (int j=0; j < numcols; j++) {
            randfloats[i][j] = all_floats[i][j];
        }
    }
    all_floats.clear();

    //------------------------------
    // main::test01
    //------------------------------

    std::ofstream test1file;
    test1file.open("test1_results.txt");
    double* test1doubles = new double[numcols];
    // we use the first row of randfloats for
    //   our test here;
    // we figured there could be some precision
    //   loss in converting from FLOAT_A to
    //   double, but there appears to be none;
    // float32 converts precisely to double, or
    //   float64; 
    // it turns out that when BASE and MANTISSA
    //   are larger than float capacity, our
    //   converted FLOAT_A's-to-doubles equal
    //   the converted float32's-to-doubles,
    //   so that their difference is zero;
    for (int j=0; j < numcols; j++) {
        FLOAT_A* x = translate(randfloats[0][j]);
        double d1 = FLOAT_A_to_double(*x);
        double d2 = (double)(randfloats[0][j]);
        test1doubles[j] = (d1 - d2) / d2;
        test1file << test1doubles[j] << std::endl;
    }
    test1file.close();

    //------------------------------
    // main::test02
    //------------------------------

    std::ofstream test2file;
    test2file.open("test2_results.txt");
    double* test2doubles = new double[numcols];
    // we use the first two rows of randfloats
    //   for our test here;
    for (int j=0; j < numcols; j++) {
        FLOAT_A* x = translate(randfloats[0][j]);
        FLOAT_A* y = translate(randfloats[1][j]);
        FLOAT_A* sum = addition(*x, *y);
        double d1 = FLOAT_A_to_double(*sum);
        double d2 = (double)(randfloats[0][j] + randfloats[1][j]);
        test2doubles[j] = (d1 - d2) / d2;
        test2file << test2doubles[j] << std::endl;
    }
    test2file.close();

    //------------------------------
    // main::test03
    //------------------------------

    std::ofstream test3file;
    test3file.open("test3_results.txt");
    double* test3doubles = new double[numrows];
    // here we aggregate each row of randfloats;
    for (int i=0; i < numrows; i++) {
        FLOAT_A* fxsum = accumulation(randfloats[i], numcols);
        float fsum = aggregate_floats(randfloats[i], numcols);
        double d1 = FLOAT_A_to_double(*fxsum);
        double d2 = (double)(fsum);
        test3doubles[i] = (d1 - d2) / d2;
        test3file << test3doubles[i] << std::endl;
    }
    test3file.close();
}