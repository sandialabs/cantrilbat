


#include "ApplBase_print.h"

#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;
namespace ca_ab
{
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void print_char(const char c, const int nTimes)
{
    for (int i = 0; i < nTimes; i++) {
        cout << c;
    }
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void print_bool(const bool b)
{
    if (b) {
        cout << "yes";
    } else {
        cout << "no ";
    }
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void pr_if(const int i, const int w)
{
    cout.width(w);
    cout << i;
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void pr_sf(const string s, const int w)
{
    int sz = s.size();
    if (sz < w) {
        int num = w - sz;
        for (int i = 0; i < num; i++) {
            cout << " ";
        }
    }
    cout << s;
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void pr_sf_lj(const string s, const int w, const int crop)
{
    int sz = s.size();
    if (crop && sz > w) {
        const char* pos = s.c_str();
        for (int i = 0; i < w; i++, pos++) {
            cout << *pos;
        }
    } else {
        cout << s;
    }
    if (sz < w) {
        int num = w - sz;
        for (int i = 0; i < num; i++) {
            cout << " ";
        }
    }
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
 * prf_df():
 *   print with a fixed precision and width in fixed format.
 *      d = value
 *      w = width
 *      p = precision
 */
void pr_df(const double d, const int w, const int p)
{
    int wmax = w;
    if (d < 0.0) {
        wmax--;
    }
    double dlnd10 = 0.0;
    double dabs = fabs(d);
    if (dabs >=10.0) {
        dlnd10 = log10(dabs);
    }
    int idlnd10 = (int) dlnd10;
    int pmax = wmax - 2 - idlnd10;
    int puse = p;
    if (puse > pmax) {
        puse = pmax;
    }
    printf("%*.*f", w, puse, d);
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/**
 * pr_dfp()
 *   Print with a fixed precision and a variable width in a
 *   fixed format (i.e., non scientific notation)
 */

void pr_dfp(const double d, const int p)
{
    int pp = cout.precision(p);
    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout << d;
    pp = cout.precision(6);
    cout.setf(ios_base::fmtflags(0), ios_base::floatfield);
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void pr_de(const double d, const int w, const int p)
{
    cout.setf(ios_base::scientific | ios_base::uppercase);
    cout.width(w);
    int pp = cout.precision(p);
    cout <<  d;
    pp = cout.precision(pp);
    cout.unsetf(ios_base::scientific);
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void pr_dg(const double d, const int w, const int p)
{
    if (d == 0.0) {
        pr_df(d, w, p);
    } else if (fabs(d) < 0.1) {
        pr_de(d, w, p);
    } else if (fabs(d) > 1000.) {
        pr_de(d, w, p);
    } else {
        pr_df(d, w, p);
    }
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void dnt(const int i)
{
    if (i == 0) {
        /* Section headings */
        print_char(' ', 5);
    } else if (i == 1) {
        /* informative text below section headings */
        print_char(' ', 8);
    } else if (i == 2) {
        /* indentation of tables themselves */
        print_char(' ', 15);
    } else if (i == 3) {
        /* indentation for txt pertaining to 1 */
        print_char(' ', 18);
    } else if (i == 4) {
        /* indentation for small tables */
        print_char(' ', 10);
    } else {
        print_char(' ', 8);
    }
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void print_map(const map<string,double>& m, const string& prefix)
{
    if (prefix.size() > 0) {
        cout << prefix;
    }
    map<string,double>::const_iterator it;
    for (it = m.begin(); it != m.end(); it++) {
        if (it != m.begin()) {
            cout << " ";
        }
        cout << "(";
        cout << it->first << ", " << it->second << ")";
    }
}
}
