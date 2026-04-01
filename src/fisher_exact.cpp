#include <cmath>
#include <iostream>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////
//// reimplementation from
//// https://github.com/painyeph/FishersExactTest
////////////////////////////////////////////////////////////////////////////////////////

extern "C" double exact_fisher(int a, int ab, int ac, int abcd, int alternative){
    if (a < 0 || ab < 0 || ac < 0 || abcd < 0) return 1.;
    int a_min = max(0, ab + ac - abcd);
    int a_max = min(ab, ac);
    if (a_min == a_max) return 0.;
    double p0 = lgamma(ab+1) + lgamma(ac+1) + lgamma(abcd-ac+1) + lgamma(abcd-ab+1) - lgamma(abcd+1);
    double pa = lgamma(a+1) + lgamma(ab-a+1) + lgamma(ac-a+1) + lgamma(abcd-ab-ac+a+1);
    double st = 1.;
    if (alternative == 0){ // two-sided
        if (ab * ac < a * abcd){
            for (int i = min(a-1, int(round(ab * ac / abcd))); i >= a_min; --i){
                double pi = lgamma(i+1) + lgamma(ab-i+1) + lgamma(ac-i+1) + lgamma(abcd-ab-ac+i+1);
                if (pi < pa) continue;
                double st_new = st + exp(pa - pi);
                if (st_new == st) break;
                st = st_new;
            }
            for (int i = a + 1; i <= a_max; ++i){
                double pi = lgamma(i+1) + lgamma(ab-i+1) + lgamma(ac-i+1) + lgamma(abcd-ab-ac+i+1);
                double st_new = st + exp(pa - pi);
                if (st_new == st) break;
                st = st_new;
            }
        }
        else {
            for (int i = a - 1; a >= a_min; --i){
                double pi = lgamma(i+1) + lgamma(ab-i+1) + lgamma(ac-i+1) + lgamma(abcd-ab-ac+i+1);
                double st_new = st + exp(pa - pi);
                if (st_new == st) break;
                st = st_new;
            }
            for (int i = max(a+1, int(round(ab * ac / abcd))); i <= a_max; ++i){
                double pi = lgamma(i+1) + lgamma(ab-i+1) + lgamma(ac-i+1) + lgamma(abcd-ab-ac+i+1);
                if (pi < pa) continue;
                double st_new = st + exp(pa - pi);
                if (st_new == st) break;
                st = st_new;
            }
        }
        return exp(-max(0., pa - p0 - log(st)));
    }
    else if (alternative == 1){ // left-tailed
        if (ab * ac < a * abcd){
            double sr = 0.;
            for (int i = a + 1; i <= a_max; ++i){
                double sr_new = sr + exp(pa - lgamma(i+1) - lgamma(ab-i+1) - lgamma(ac-i+1) - lgamma(abcd-ab-ac+i+1));
                if (sr_new == sr) break;
                sr = sr_new;
            }
            return 1. - max(0., exp(p0 - pa) * sr);
        }
        else {
            double sl = 1.;
            for (int i = a - 1; i >= a_min; --i){
                double sl_new = sl + exp(pa - lgamma(i+1) - lgamma(ab-i+1) - lgamma(ac-i+1) - lgamma(abcd-ab-ac+i+1));
                if (sl_new == sl) break;
                sl = sl_new;
            }
            return exp(-max(0., pa - p0 - log(sl)));
        }
    }
    // right-tailed
    if (ab * ac > a * abcd){
        double sl = 0.;
        for (int i = a - 1; i >= a_min; --i){
            double sl_new = sl + exp(pa - lgamma(i+1) - lgamma(ab-i+1) - lgamma(ac-i+1) - lgamma(abcd-ab-ac+i+1));
            if (sl_new == sl) break;
            sl = sl_new;
        }
        return 1. - max(0., exp(p0 - pa) * sl);
    }
    else {
        double sr = 1.;
        for (int i = a + 1; i <= a_max; ++i){
            double sr_new = sr + exp(pa - lgamma(i+1) - lgamma(ab-i+1) - lgamma(ac-i+1) - lgamma(abcd-ab-ac+i+1));
            if (sr_new == sr) break;
            sr = sr_new;
        }
        return exp(-max(0., pa - p0 - log(sr)));
    }
}

