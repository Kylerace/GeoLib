
pub fn factorial(n: anytype) @TypeOf(n) {
    var result: @TypeOf(n) = 1;
    var k: @TypeOf(n) = n;
    while (k > 1) {
        result *= k;
        k -= 1;
    }
    //for (2..n+1) |k| {
    //    result *= k;
    //}
    return result;
}

pub fn cosh(num: anytype) @TypeOf(num) {
    return 0.5 * (@exp(num) + @exp(-num));
}

pub fn sinh(num: anytype) @TypeOf(num) {
    return 0.5 * (@exp(num) - @exp(-num));
}