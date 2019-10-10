/**
 *  @file sort.cpp
 *     Definitions for the sort routines.
 */

#ifdef WIN32
#pragma warning(disable:4786)
#endif

#include "sort.h"

using namespace std;

//----------------------------------------------------------------------------------------------------------------------------------
namespace mdpUtil
{

//==================================================================================================================================
// sort (x,y) pairs by x
void heapsort(std::vector<double>& x, std::vector<int>& y)
{
    int n = x.size();
    if (n < 2) {
        return;
    }
    double rra;
    int rrb;
    int ll = n/2;
    int iret = n-1;

    while (1 > 0) {
        if (ll > 0) {
            ll--;
            rra = x[ll];
            rrb = y[ll];
        } else {
            rra = x[iret];
            rrb = y[iret];
            x[iret] = x[0];
            y[iret] = y[0];
            iret--;
            if (iret == 0) {
                x[0] = rra;
                y[0] = rrb;
                return;
            }
        }

        int i = ll;
        int j = ll + ll + 1;

        while (j <= iret) {
            if (j < iret) {
                if (x[j] < x[j+1]) {
                    j++;
                }
            }
            if (rra < x[j]) {
                x[i] = x[j];
                y[i] = y[j];
                i = j;
                j = j + j + 1;
            } else {
                j = iret + 1;
            }
        }
        x[i] = rra;
        y[i] = rrb;
    }
}
//==================================================================================================================================
void heapsort(std::vector<double>& x, std::vector<double>& y)
{
    int n = x.size();
    if (n < 2) {
        return;
    }
    double rra;
    double rrb;
    int ll = n/2;
    int iret = n-1;

    while (1 > 0) {
        if (ll > 0) {
            ll--;
            rra = x[ll];
            rrb = y[ll];
        } else {
            rra = x[iret];
            rrb = y[iret];
            x[iret] = x[0];
            y[iret] = y[0];
            iret--;
            if (iret == 0) {
                x[0] = rra;
                y[0] = rrb;
                return;
            }
        }

        int i = ll;
        int j = ll + ll + 1;

        while (j <= iret) {
            if (j < iret) {
                if (x[j] < x[j+1]) {
                    j++;
                }
            }
            if (rra < x[j]) {
                x[i] = x[j];
                y[i] = y[j];
                i = j;
                j = j + j + 1;
            } else {
                j = iret + 1;
            }
        }
        x[i] = rra;
        y[i] = rrb;
    }
}
//==================================================================================================================================
template <typename T >
void heapsortT(std::vector<double>& x, std::vector< T >& y)
{
    int n = x.size();
    if (n < 2) {
        return;
    }
    double rra;
    T rrb;
    int ll = n/2;
    int iret = n-1;

    while (1 > 0) {
        if (ll > 0) {
            ll--;
            rra = x[ll];
            rrb = y[ll];
        } else {
            rra = x[iret];
            rrb = y[iret];
            x[iret] = x[0];
            y[iret] = y[0];
            iret--;
            if (iret == 0) {
                x[0] = rra;
                y[0] = rrb;
                return;
            }
        }

        int i = ll;
        int j = ll + ll + 1;

        while (j <= iret) {
            if (j < iret) {
                if (x[j] < x[j+1]) {
                    j++;
                }
            }
            if (rra < x[j]) {
                x[i] = x[j];
                y[i] = y[j];
                i = j;
                j = j + j + 1;
            } else {
                j = iret + 1;
            }
        }
        x[i] = rra;
        y[i] = rrb;
    }
}
//==================================================================================================================================
template void heapsortT<int>(std::vector<double>& x, std::vector< int >& y);
template void heapsortT<size_t>(std::vector<double>& x, std::vector< size_t >& y);
template void heapsortT<double>(std::vector<double>& x, std::vector< double >& y);


}
//----------------------------------------------------------------------------------------------------------------------------------
