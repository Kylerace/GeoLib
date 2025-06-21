
// in this file a subset is a bitfield where bit i == 1 iff the subset contains element i
// said element is almost always a blade in the context of our comptime implementation of geometric algebra

/// counts the number of 1 bits in the v bitfield
pub fn count_bits(v: anytype) @TypeOf(v) {
    var val = v;
    var count: @TypeOf(v) = 0;
    //var parity: bool = false;
    while (val != 0) {
        //parity = !parity;
        count += 1;
        val &= val - 1;
    }
    return count;
}

/// counts the number of flips the result of a standard geometric multiplication operation (e.g. gp, op, ip, rp) will have
/// given the two input blades left and rhs. 
pub fn count_flips(left: anytype, rhs: anytype) @TypeOf(left) {
    var lhs = left >> 1;
    var flips: @TypeOf(left) = 0;
    while (lhs != 0) {
        flips += count_bits(lhs & rhs);
        lhs >>= 1;
    }
    //flips += count_bits((left >> (1 + i)) + rhs)
    return flips;
}

/// returns whether the number of flips for left and rhs is odd or even (in integer form)
pub fn flip_parity(left: anytype, rhs: @TypeOf(left)) @TypeOf(left) {
    return count_flips(left, rhs) & 1;
}

pub fn flips_permuted_parity(flip_field: anytype, lhs: anytype, rhs: @TypeOf(lhs)) @TypeOf(lhs) {

    return count_flips(lhs,rhs) ^ (flip_field >> lhs) ^ (flip_field >> rhs) ^ (flip_field >> (lhs ^ rhs)) & 1;
}

pub fn flips_permuted_neg_parity(neg_mask: anytype, flip_field: @TypeOf(neg_mask), lhs: anytype, rhs: @TypeOf(lhs)) @TypeOf(lhs) {
    return (count_bits(neg_mask & (lhs & rhs)) ^ count_flips(lhs,rhs) ^ (flip_field >> lhs) ^ (flip_field >> rhs) ^ (flip_field >> (lhs ^ rhs))) & 1;
}

pub fn unsigned_to_signed_parity(unsigned_parity: anytype) isize {
    return (-2 * @as(isize, unsigned_parity & 1) + 1);
}

/// returns true iff tjhe given subset includes the given blade e.g. 101 contains 0 and 11 but not 10
pub fn subset_includes_blade(subset: anytype, blade: anytype) bool {
    return (1 << blade) & subset != 0;
}