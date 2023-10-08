#include <iostream>

typedef struct bitfield
{
    signed int   a:3;
    unsigned int b:13;
    unsigned int c:1;
} bitfield;


int main()
{
    std::cout << "Hello World" << std::endl;
    std::cout << sizeof(bitfield) << std::endl;
}