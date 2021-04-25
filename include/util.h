


#pragma once


namespace gadg {
    // gcd of more than two multivectors
    int gcd(const int a, const int b) {
        if (a == 0)
            return b;
        return gcd(b % a, a);
    }

    // Function to find gcd of 4 values
    int gcd(const int w, const int x, const int y, const int z) {
        int res_gcd = w;
        res_gcd = gcd(x, res_gcd);
        if (res_gcd == 1) {
            return 1;
        }
        res_gcd = gcd(y, res_gcd);
        if (res_gcd == 1) {
            return 1;
        }
        return gcd(z, res_gcd);
    }
}
