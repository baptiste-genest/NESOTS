#include "ssw.h"

int weightedMedian(const ints& arr,const scalars& W)
{

    // Store pr of arr[i] and W[i]
    std::vector<std::tuple<int, scalar,int>> pr;

    for(int index = 0;
        index < arr.size();
        index++)
        pr.push_back({arr[index],
                  W[index],index});

    // Sort the list of pr w.r.t.
    // to their arr[] values
    sort(pr.begin(), pr.end());

    scalar sums = 0;

    // If N is odd
    if (arr.size() % 2 != 0)
    {

        // Traverse the set pr
        // from left to right
        for(int i = 0;i<pr.size();i++)
        {
            sums += std::get<1>(pr[i]);
            if (sums > 0.5)
                return std::get<2>(pr[i]);
        }
    }

    // If N is even
    else
    {

        // For lower median traverse
        // the set pr from left
        int i = 0;
        for(auto element : pr)
        {

            // Update sums
            sums += std::get<1>(element);

            // When sum >= 0.5
            if (sums >= 0.5)
            {
                return std::get<2>(element);
                break;
            }
            i++;
        }

        // For upper median traverse
        // the set pr from right
        sums = 0;
        for(int index = pr.size() - 1;
            index >= 0;
            index--)
        {
            scalar weight = std::get<1>(pr[index]);

            // Update sums
            sums += weight;

            // When sum >= 0.5
            if (sums >= 0.5)
            {
                return std::get<2>(pr[index]);
            }
        }
    }
    std::cerr << "[weighted median] : fail -> " << sums <<std::endl;
    assert(false);
    return -1;
}

std::tuple<ints, scalars, array_elements> cumulative_function_difference(const scalars &a, const scalars &b)
{
    auto C = merge_sorted_arrays(a,b);
    ints y(C.size());scalars d(C.size());
    int la = 0,lb = 0;
    for (int i = 0;i<C.size();i++) {
            for (;a[la] < C[i].value && la < a.size();la++){};
            for (;b[lb] < C[i].value && lb < b.size();lb++){};
        y[i] = la-lb;
        if (i)
            d[i] = C[i].value-C[i-1].value;
        else{
            auto x= std::abs(C[0].value-C.back().value);
            d[0] = std::min(x,1-x);
        }
    }
    return {y,d,C};
}

array_elements merge_sorted_arrays(const scalars &arr1, const scalars &arr2)
{
    auto n = arr1.size();
    auto m = n;
    int i = 0, j = 0; // pointers
    array_elements Union; // Uninon vector
    while (i < n && j < m) {
        if (arr1[i] <= arr2[j]) // Case 1 and 2
            Union.push_back({arr1[i],false,i++,j});
        else // case 3
            Union.push_back({arr2[j],true,i,j++});
    }
    while (i < n) // IF any element left in arr1
        Union.push_back({arr1[i],false,i++,j});
    while (j < m) // If any elements left in arr2
        Union.push_back({arr2[j],true,i,j++});
    if (!std::is_sorted(Union.begin(),Union.end(),[](const array_element& a,const array_element& b){
                return a.value < b.value;
})){
        std::cerr << "NON SORTED UNION" << std::endl;
        assert(false);
    }
    return Union;
}

int computeOptimalCut2(const scalars &a, const scalars &b)
{
    auto&& [y,d,C] = cumulative_function_difference(a,b);
    int cut_id = weightedMedian(y,d);
    return (C[cut_id].index_b + a.size()-C[cut_id].index_a)%b.size();
}
