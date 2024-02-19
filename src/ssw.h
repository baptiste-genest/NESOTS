#ifndef SSW_H
#define SSW_H

#include "utils.h"

#include <execution>
#include <iostream>

struct array_element{
    scalar value;
    bool belong_to_b;
    int index_a;
    int index_b;
};

using array_elements = std::vector<array_element>;

array_elements merge_sorted_arrays(const scalars& a,const scalars& b);

//REQUIRES A AND B TO BE SORTED
std::tuple<ints,scalars,array_elements> cumulative_function_difference(const scalars& a,const scalars& b);

//https://www.geeksforgeeks.org/program-to-find-weighted-median-of-a-given-array/
int weightedMedian(const ints& arr,const scalars& W);

int computeOptimalCut2(const scalars& a,const scalars& b);

scalar transport1DW1(const labeled<scalar>& a,const labeled<scalar>& b,std::pair<int,int> cuts);

#endif // SSW_H
