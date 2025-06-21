const std = @import("std");
const debug = std.debug;
const comptimePrint = std.fmt.comptimePrint;

const bit_utils = @import("utils/bit_utils.zig");

// this file handles our representation of an algebra as a whole without anything to do with how we hold data within that algebra
// at runtime


pub fn BladeMult(T: type) type {
    return struct {
        ///1 where the ith 1 blade is included, 0 elsewhere
        /// count_bits(blade)th form
        blade: usize,
        value: T,
        left_value: T = 0,
        right_value: T = 0,
    };
}

pub fn ProdTerm(T: type) type {
    return struct { 
        mult: T, 
        left_blade: usize, 
        right_blade: usize,
        left_value: T,
        right_value: T, 


        pub fn format(
            slf2: @This(),
            comptime _: []const u8,
            _: std.fmt.FormatOptions,
            writer2: anytype
        ) !void {
            _ = try writer2.print("{d:.1} from {d:.1} {b} * {d:.1} {b}", .{slf2.mult, slf2.left_value, slf2.left_blade, slf2.right_value, slf2.right_blade});
        }

    };
}

pub fn UnaryProdTerm(T: type) type {
    return struct {
        mult: T,
        source: usize
    };
}

///contains all information necessary and/or convenient for GAlgebraOp
pub const GAlgebraSignature = struct {
    ///num positive-squaring basis 1blades
    p: usize,
    ///num negative-squaring basis 1blades
    n: usize,
    ///num zero-squaring basis 1blades
    z: usize,
    ///num basis 1blades
    d: usize,
    basis_size: usize,

    pseudoscalar: usize,

    ///holds a 1 in bit index i (0 indexed) if that blade squares to a positive number
    pos_mask: usize,
    ///holds a 1 in bit index i (0 indexed) if that blade squares to a negative number
    neg_mask: usize,
    ///holds a 1 in bit index i (0 indexed) if that blade squares to 
    zero_mask: usize,
    /// holds a 1 in the ith bit index (0-indexed) if blade i in our desired basis has a negative parity relative to the "canonical" blade composed of strictly increasing 1blade factors from left to right
    /// for example, in a basis of 1blades {x,y,z}, we can define the pseudoscalar in our basis to be yxz, which is flipped relative to the canonical binary basis of xyz
    /// because the canonical ordering has the 1st 1blade x in the first place, the 2nd 1blade y in the 2nd, etc. so we would have a 1 in index 111 binary 7 decimal of the bitfield (index 0 is the scalar blade)
    /// 10000000 -> I = e213, e132, or e321, 10000000 >> 111 = 1
    /// 76543210 bit indices
    /// 1000 -> blade e12 = binary 11 is flipped (e21)
    /// (this >> a blade) & 1 == 1 iff that blade is flipped, and is 0 otherwise
    flipped_parities: usize,
    /// i dont think this is actually necessary for any functions in GAlgebraOp
    /// but just in case
    dual: bool,
};

/// handles most operations in geometric algebra 
pub fn GAlgebraOp(T: type) type {
    return enum {
        pub const Self: type = @This();
        pub const BladeTerm: type = BladeMult(T);

        GeometricProduct,
        OuterProduct,
        InnerProduct,
        RegressiveProduct,
        SandwichProduct,
        ScalarMul,
        //TODO: figure out how to handle ops with arbitrary input and output sizes,
        //  with outputs being possibly nullable
        //  e.g. Add would be a map from (BladeTerm, BladeTerm) -> (BladeTerm, ?BladeTerm)
        //  note: is this yak shaving?
        Add,
        Sub,

        //unary
        Dual,
        Undual,
        Reverse,
        Involution,

        ///values here are the k in n choose k (n = basis_len)
        pub fn n_ary(self: @This()) comptime_int {
            return switch (self) {
                .GeometricProduct => 2,
                .OuterProduct => 2,
                .InnerProduct => 2,
                .RegressiveProduct => 2,
                .SandwichProduct => 2,
                .ScalarMul => 2,
                .Dual => 1,
                .Undual => 1,
                .Reverse => 1,
                .Involution => 1,
                .Add => 2,
                .Sub => 2,
            };
        }

        pub fn Codomain(self: @This()) type {
            return switch(self) {
                .GeometricProduct => BladeTerm,
                .OuterProduct => BladeTerm,
                .InnerProduct => BladeTerm,
                .RegressiveProduct => BladeTerm,
                .SandwichProduct => BladeTerm,
                .ScalarMul => BladeTerm,
                .Dual => BladeTerm,
                .Undual => BladeTerm,
                .Reverse => BladeTerm,
                .Involution => BladeTerm,
                .Add => struct{BladeTerm, ?BladeTerm},
                .Sub => struct{BladeTerm, ?BladeTerm},
            };
        }

        pub fn max_output_size(self: @This()) comptime_int {
            return switch(self) {
                .GeometricProduct => 1,
                .OuterProduct => 1,
                .InnerProduct => 1,
                .RegressiveProduct => 1,
                .SandwichProduct => 1,
                .ScalarMul => 1,
                .Dual => 1,
                .Undual => 1,
                .Reverse => 1,
                .Involution => 1,
                .Add => 2,
                .Sub => 2,
            };
        }

        pub fn BladeOp(self: Self, sig: GAlgebraSignature) type {
            return switch (self) {
                
                .GeometricProduct => struct {
                    pub fn eval(left: BladeTerm, right: BladeTerm) Codomain(.GeometricProduct) {
                        //if(@abs(left.value) < 0.00001 or @abs(right.value) < 0.00001) {
                        //    return BladeTerm{.value = 0, .blade = 0};
                        //}
                        const zero_mask = sig.zero_mask;
                        const neg_mask = sig.neg_mask;
                        const flip_field = sig.flipped_parities;

                        const lhs = left.blade;
                        const rhs = right.blade;
                        const squares = lhs & rhs;
                        const blade: usize = lhs ^ rhs;
                        if (squares & zero_mask != 0) {
                            return BladeTerm{.value = 0, .blade = 0 };
                        }
                        //const neg_powers_mod_2: isize = @intCast(((squares & neg_mask) ^ count_flips(lhs, rhs)) & 1);
                        //const value: isize = -2 * neg_powers_mod_2 + 1;
                        const value: isize = bit_utils.unsigned_to_signed_parity(bit_utils.flips_permuted_neg_parity(neg_mask, flip_field, lhs, rhs));
                        const ret_val = left.value * right.value * @as(f64, @floatFromInt(value));
                        //if(@abs(left.value) > 0.0001 and @abs(right.value) > 0.0001) {
                        //    @compileLog(comptimePrint("GeometricProduct.eval(): {d:.1} {b} * {d:.1} {b} = (value {}-> retval {d:.1}){b} with neg_mask {b} and flip_field {b} and zero_mask {b} squares {b}", .{left.value, left.blade, right.value, right.blade, value, ret_val, blade, neg_mask, flip_field, zero_mask, squares}));
                        //}
                        if(@round(ret_val) != ret_val) {
                            @compileLog(comptimePrint("\n\t\tgp lhs val {} * rhs val {} * parity {} = noninteger val {}", .{left.value, right.value, @as(f64, @floatFromInt(value)), ret_val}));
                        }
                        if(@abs(ret_val) > 10.0) {
                            @compileLog(comptimePrint("\n\t------ret_val from .GeometricProduct is greater than 10: {}", ret_val));
                        }
                        //@compileLog(comptimePrint("\nret_val = {}", .{ret_val}));
                        return BladeTerm{ .value = ret_val, .blade = blade };
                    }
                },
                .OuterProduct => struct {
                    pub fn eval(left: BladeTerm, right: BladeTerm) Codomain(.OuterProduct) {
                        //if(@abs(left.value) < 0.00001 or @abs(right.value) < 0.00001) {
                        //    return BladeTerm{.value = 0, .blade = 0};
                        //}
                        const lhs = left.blade;
                        const rhs = right.blade;
                        const squares = lhs & rhs;
                        const blade = lhs ^ rhs;
                        if (squares != 0) {
                            return .{ .value = 0, .blade = blade };
                        }
                        //return .{ .value = if (count_flips(lhs, rhs) & 1) -1 else 1, .blade = blade };
                        return .{.value = bit_utils.unsigned_to_signed_parity(bit_utils.flips_permuted_parity(sig.flipped_parities, lhs, rhs)), .blade = blade};
                    }
                },
                .InnerProduct => struct {
                    pub fn eval(left: BladeTerm, right: BladeTerm) Codomain(.InnerProduct) {
                        //if(@abs(left.value) < 0.00001 or @abs(right.value) < 0.00001) {
                        //    return BladeTerm{.value = 0, .blade = 0};
                        //}
                        const neg_mask: usize = sig.neg_mask;
                        const zero_mask = sig.zero_mask;
                        const flip_field = sig.flipped_parities;

                        const lhs = left.blade;
                        const rhs = right.blade;
                        const squares = lhs & rhs;
                        const blade = lhs ^ rhs;
                        if (squares & zero_mask != 0) {
                            return .{.value = 0, .blade = blade };
                        }
                        const ul = lhs & blade;
                        const ur = rhs & blade;

                        //const neg_powers_mod_2: isize = @intCast(((squares & neg_mask) ^ count_flips(lhs, rhs)) & 1);
                        //const value: isize = -2 * neg_powers_mod_2 + 1;
                        const value: isize = bit_utils.unsigned_to_signed_parity(bit_utils.flips_permuted_neg_parity(neg_mask, flip_field, lhs, rhs));
                        return .{ .value = if (ul == 0 or ur == 0) left.value * right.value * value else 0, .blade = blade };
                    }
                },

                .RegressiveProduct => struct {
                    pub fn eval(left: BladeTerm, right: BladeTerm) Codomain(.RegressiveProduct) {
                        //if(@abs(left.value) < 0.00001 or @abs(right.value) < 0.00001) {
                        //    return BladeTerm{.value = 0, .blade = 0};
                        //}
                        const neg_mask = sig.neg_mask;
                        const flip_field = sig.flipped_parities;
                        const pseudoscalar = sig.pseudoscalar;

                        const lhs = left.blade ^ pseudoscalar;
                        const rhs = right.blade ^ pseudoscalar;
                        const squares = lhs & rhs;
                        const blade = lhs ^ rhs ^ pseudoscalar;
                        if (squares != 0) {
                            return .{ .value = 0, .blade = blade };
                        }
                            
                        //const value = -2 * (((squares & neg_mask) ^ count_flips(lhs, rhs)) & 1) + 1;
                        
                        return .{ .value = bit_utils.unsigned_to_signed_parity(bit_utils.flips_permuted_neg_parity(neg_mask, flip_field, lhs, rhs)), .blade = blade };
                        //const dual = Operation.Dual.BladeOp();
                        //const undual = Operation.Undual.BladeOp();
                        //const outer = Operation.OuterProduct.BladeOp();
                        //return undual.eval(outer.eval(dual.eval(left), dual.eval(right)));
                    }
                },

                //execution is something like first * first * second * coeff
                .SandwichProduct => struct {
                    pub fn eval(left: BladeTerm, right: BladeTerm) Codomain(.SandwichProduct) {
                        //if(@abs(left.value) < 0.00001 or @abs(right.value) < 0.00001) {
                        //    return BladeTerm{.value = 0, .blade = 0};
                        //}
                        //const bread = left.blade;
                        //const meat = right.blade;
                        const gp = GAlgebraOp.GeometricProduct.BladeOp();
                        const reverse = GAlgebraOp.Reverse.BladeOp();
                        return gp.eval(gp.eval(left, right), reverse.eval(left));
                    }
                },


                //TODO: idk how im supporting this if left isnt known at compile time
                .ScalarMul => struct {
                    pub fn eval(left: f64, right: BladeTerm) Codomain(.ScalarMul) {
                        //const bread = left.blade;
                        //const meat = right.blade;
                        return BladeTerm{.value = left * right.value, .blade = right.blade};
                    }
                },

                .Dual => struct {
                    pub fn eval(us: BladeTerm) Codomain(.Dual) {
                        //const squares = us.blade & sig.pseudoscalar;
                        return .{ .value = us.value * bit_utils.unsigned_to_signed_parity(bit_utils.flips_permuted_neg_parity(sig.neg_mask, sig.flipped_parities, us, sig.pseudoscalar)), .blade = us.blade ^ sig.pseudoscalar };
                    }
                },

                .Undual => struct {
                    pub fn eval(us: BladeTerm) Codomain(.Undual) {
                        //var pseudosquarelar = BladeTerm{.value = 1, .blade = pseudoscalar_blade};
                        //const op = Operation.GeometricProduct.BladeOp();
                        //pseudosquarelar = op.eval(pseudosquarelar, pseudosquarelar);
                        //return op.eval(pseudosquarelar, us);
                        const op = @This().Dual.BladeOp();
                        return op.eval(op.eval(op.eval(us)));
                    }
                },

                .Reverse => struct {
                    pub fn eval(us: BladeTerm) Codomain(.Reverse) {
                        const blade: usize = @truncate(us.blade);
                        const one: usize = 0b1;
                        const mult: f64 = comptime blk: {
                            var x: usize = bit_utils.count_bits(blade) & ~one;
                            x >>= 1;
                            if (x & 1 == 0) {
                                break :blk 1.0;
                            }
                            break :blk -1.0;
                        };
                        return BladeTerm{.blade = blade, .value = mult};
                    }
                },

                .Involution => struct {
                    pub fn eval(us: BladeTerm) Codomain(.Involution) {
                        const blade: usize = @truncate(us.blade);
                        return .{.blade = blade, .value = if(blade & 1) -1.0 else 1.0};
                    }
                },

                .Add => struct {
                    pub fn eval(left: BladeTerm, right: BladeTerm) Codomain(.Add) {
                        if(left.blade == right.blade) {
                            return .{BladeTerm{.blade = left.blade, .value = left.value + right.value}, null};
                        }
                        return .{left, right};
                    }
                },
                
                .Sub => struct {
                    pub fn eval(left: BladeTerm, right: BladeTerm) Codomain(.Sub) {
                        if(left.blade == right.blade) {
                            return .{BladeTerm{.blade = left.blade, .value = left.value - right.value}, null};
                        }
                        return .{left, BladeTerm{.blade = right.blade, .value = -right.value}};
                    }
                }
            };
        }

    };
}

const BasisMetric = enum {
    positive,
    negative,
    zero,

    pub fn square(self: @This()) isize {
        return switch(self) {
            .positive => 1,
            .negative => -1,
            .zero => 0
        };
    }
    pub fn from_square(val: isize) @This() {
        if(val < 0) {
            return .negative;
        }
        if(val == 0) {
            return .zero;
        }
        return .positive;
    } 
};

const BladeFloat = struct {
    blade: usize,
    val: f64,
};

///null val = runtime known, nonnull = comptime known
const RuntimeOrComptimeBlade = struct {
    blade: usize,
    val: ?f64,
};

//const pga3_aliases = .{.{name: "Motor", runtime_subset: "s,ix,iy,iz,xy,yz,zx"}, 
//                       .{name: "Point", runtime_subset: "izy,ixz,iyx,xyz", specific_aliases: "izy=px,ixz=py,iyx=pz,xyz=pi"},
//                       .{name: "NormPoint", runtime_subset: "izy,ixz,iyx", comp_subset: "xyz=1", specific_aliases: "izy=px,ixz=py,iyx=pz,xyz=pi"},
//                       .{name: "Line", runtime_subset: "ix,iy,iz,xy,yz,zx"},
//                       .{name: "IdealLine", runtime_subset: "ix,iy,iz"}}
///named subset of a multivector with basis vectors in a specified order and specified runtime and comptime subsets
pub const GAlgebraAlias = struct {
    const Self: type = @This();
    name: ?[]const u8,
    ///blade names, possibly out of order from the canonical ordering
    ordering: []const u8,
    //specific_aliases: ?[]const struct{usize, []const u8},
    pub fn format(
        self: @This(),
        comptime _: []const u8,
        _: std.fmt.FormatOptions,
        writer: anytype
    ) !void {
        _=try writer.print("GAlgebraAlias name: {s}, ordering: {s}", .{
            self.name orelse "",
            self.ordering
        });
    }

};

/// an alias except with additional info created from other arguments in GAlgebraInfo()
pub const GAlgebraCompiledAlias = struct {
    const Self: type = @This();

    name: ?[]const u8,
    /// original ordering string from the inputted alias
    ordering: []const u8,
    /// the parsed blades in binary representatin (cannot be flipped relative to the algebras native ordering right now)
    values: []const usize,

    pub fn format(
        self: @This(),
        comptime _: []const u8,
        _: std.fmt.FormatOptions,
        writer: anytype
    ) !void {
        _=try writer.print("GAlgebraAlias name: {s}, ordering: {s}, values: {any}", .{
            self.name orelse "",
            self.ordering,
            self.values
        });
    }
    
};

/// for each char in the input at index i, looks for its match in the ordering string at index j and sets 
/// output[i] = j
/// comptime so that i can return an array without an allocator (which wouldnt work at comptime)
pub fn chars_to_binary_ordering_indices(comptime len: usize, comptime input: []const u8, comptime ordering: []const u8) [len]usize {
    var output: [len]usize = undefined;
    if(input.len > ordering.len or input.len != len) {
        debug.panic("non matching lens between len {}, input {s}, and ordering {s}", .{len, input, ordering});
    }
    for(0..input.len) |i| {
        const char: u8 = input[i];
        var found: bool = false;
        for(0..ordering.len) |j| {
            const ordered_char: u8 = ordering[j];
            if(char != ordered_char) {
                continue;
            }
            found = true;
            output[i] = j;
        }
        if(found == false) {
            debug.panic("could not find char {s} from input {s} inside of ordering {s}!", .{&.{char}, input, ordering});
        }
    }
    return output;
}

/// assumes all elements are distinct
/// also assumes T is an integer type, signed or unsigned
/// sorts in ascending order.
/// used to determine the number of flips in a basis blade characterized by a string instead of a bitfield
pub fn cycle_sort_ret_writes(T: anytype, slice: []T) usize {
    var cycle_elements: usize = 0;
    var cycles: usize = 0;
    for (0..slice.len) |cycle_start| {
        var item = slice[cycle_start];
        var pos = cycle_start;
        //everything behind us is definitely sorted, some things ahead of us may be sorted
        for (slice[cycle_start + 1 .. slice.len]) |ahead| {
            if (ahead < item) {
                pos += 1;
            }
        }

        if (pos != cycle_start) {
            cycles += 1;
            //cycle_elements += 1;
        }

        //if correct and real pos's are different, then there is a cycle of some length that includes item and the item_at_correct_pos's pos
        //and item_at_correct_pos's_correct_pos and so on until item_at_correct_pos's_correct_pos's..._correct_pos = pos
        //all cycles are disjoint, and the number of writes is the sum of the length of each cycle - number of cycles
        while (pos != cycle_start) {
            //writes += 1;
            pos = cycle_start;
            for (slice[cycle_start + 1 .. slice.len]) |ahead| {
                if (ahead < item) {
                    pos += 1;
                }
            }

            cycle_elements += 1;

            const swapped_item = slice[pos];
            slice[pos] = item;
            item = swapped_item;
        }
    }
    return cycle_elements - cycles;
}

///badly named, its actually a basis blade
pub const BasisVector = struct {
    name: []const u8,
    blade: usize,
    square: BasisMetric,

    pub fn format(
        vector: @This(),
        comptime _: []const u8,
        _: std.fmt.FormatOptions,
        writer: anytype
    ) !void {
        _=try writer.print("BasisVector {{{s}, {b}, {c}}}", .{
            vector.name,
            vector.blade,
            switch(vector.square) {
                .positive => '+',
                .negative => '-',
                .zero => '0'
            }
        });
    }

    pub fn equals(self: BasisVector, other: BasisVector) bool {
        if(!std.mem.eql(u8, self.name, other.name)) {
            return false;
        }
        if(self.blade != other.blade) {
            return false;
        }
        if(self.square != other.square) {
            return false;
        }
        return true;
    }

    pub fn flips(comptime self: @This(), comptime ordering: []const u8) isize {
        const name = self.name;
        if(name.len <= 1) {
            return 1;
        }
        var name_indices = chars_to_binary_ordering_indices(name.len, name, ordering);
        return cycle_sort_ret_writes(usize, &name_indices);
    }
};

pub fn BasisVectorMult(T: type) type {
    return struct {
        basis: BasisVector,
        value: T,
    };
}

pub fn BasisVectorMultOpt(T: type) type {
    return struct {
        basis: BasisVector,
        value: ?T,
    };
}

//fn GAlgebraInfo(dual = true, basis = "0i,+x,+y,+z", ordering="", aliases=pga3_aliases)
//holds info that 2 different MVecSubset / impl types use to distinguish between subsets of their own space vs subsets of other spaces
//1. parse 1-form basis, create p,n,z, d, and masks
//2. generate 2^d basis
//  how to deal with reordering and general ops?
//  cayley_table() doesnt need the full array it just needs subsets + and the square masks, 
//  then when we produce the scatter masks we check what actual index the given opr blade belongs to
//  for each multivector subset operand.
//  so we need a "canonical" ordering like the original Algebra impl + a way to convert from that
//  to the correct orderings
//
//3. generate reordered masks
//4. have a function that goes from canonical blade -> reordered index for that blade
//5. have a function that goes from reordered index -> canonical blade
pub fn GAlgebraInfo(comptime _dual: bool, comptime _basis: []const u8, comptime _ordering: ?[]const u8, comptime _aliases: []const GAlgebraAlias) type {
    var _p: usize = 0;
    var _z: usize = 0;
    var _n: usize = 0;
    var _d: usize = 0;

    var _p_mask: usize = 0;
    var _n_mask: usize = 0;
    var _z_mask: usize = 0;

    //ith blade is 1 if we actually use a flipped parity version of that blade as a basis
    var _flipped_parities: usize = 0;

    var parsed_1basis: []const BasisVector = &.{};
    const scalar = BasisVector{.name = "s",.blade = 0,.square = BasisMetric.positive};
    parsed_1basis = parsed_1basis ++ .{scalar};

    const BasisParseStates = enum {
        ///expects a 0,-, or +, if the next character isnt one of these it explodes
        ParseMetric,
        ///expects an ascii character (non ,)
        /// if it finds one, adds that char to the name of curr_basis_blade and moves to ParseNthBasisChar
        /// otherwise, creates a compile error
        ParseFirstBasisChar,
        ///expects either an ascii character (non ,) or a ,
        /// if it finds an ascii letter it adds it to the name of the current basis blade 
        /// if it finds a , it sets parsed_basis to parsed_basis ++ curr_basis_blade and creates a new basis vector to add to
        /// if the for loop ends at any other stage a compile error is thrown
        ParseNthBasisChar,

    };
    var curr_parse_state: BasisParseStates = .ParseMetric;
    var curr_basis_blade: BasisVector = undefined;
    var curr_basis_i: usize = 0;
    for(0.._basis.len+1) |i| {
        switch(curr_parse_state) {
            .ParseMetric => {
                const curr_char = _basis[i];
                if(curr_char == ' ') {
                    continue;
                }
                if(curr_char == '+') {
                    curr_basis_blade.square = .positive;
                    _p += 1;
                    _p_mask |= (1 << curr_basis_i);
                } else if(curr_char == '-') {
                    curr_basis_blade.square = .negative;
                    _n += 1;
                    _n_mask |= (1 << curr_basis_i);
                } else if(curr_char == '0') {
                    curr_basis_blade.square = .zero;
                    _z += 1;
                    _z_mask |= (1 << curr_basis_i);
                } else {
                    @compileError(comptimePrint("bad char found in basis parsing of GAlgebraInfo! basis: {s}", .{_basis}));
                }
                _d += 1;
                curr_basis_i += 1;
                curr_parse_state = .ParseFirstBasisChar;
            },
            .ParseFirstBasisChar => {
                const curr_char = _basis[i];
                if(curr_char == ',') {
                    @compileError(comptimePrint("GAlgebraInfo error: skipped basis vector name (expected non , but got a ,)", .{}));
                }
                //TODO: make sure curr_char is a letter
                //TODO: allow for e_123 basis blades
                curr_basis_blade.name = &.{curr_char};
                curr_basis_blade.blade = (1 << (curr_basis_i - 1));
                curr_parse_state = .ParseNthBasisChar;
            },
            .ParseNthBasisChar => {
                if(i == _basis.len or _basis[i] == ',') {
                    //const _curr_basis_blade = curr_basis_blade;
                    parsed_1basis = parsed_1basis ++ .{curr_basis_blade};
                    curr_basis_blade = undefined;
                    curr_parse_state = .ParseMetric;
                    continue;
                }
                const curr_char = _basis[i];
                if((curr_char < 'A' or curr_char > 'z') or (curr_char > 'Z' and curr_char < 'a')) {
                    @compileError(comptimePrint("non letter char passed into basis string of GAlgebraInfo! _basis: {s}", .{_basis})); 
                }
                curr_basis_blade.name = curr_basis_blade.name ++ .{curr_char};
            }
        }
    }

    const basis_size: usize = std.math.pow(usize, 2, _d);

    var parsed_basis: [basis_size]BasisVector = undefined;
    parsed_basis[0] = parsed_1basis[0];

    @setEvalBranchQuota(1000000);

    for(1..basis_size) |blade| {
        const blade_size: usize = bit_utils.count_bits(blade);
        if(blade_size < 2) {
            for(0.._d) |i| {
                if(blade & (1 << i) == 0) {
                    continue;
                }
                parsed_basis[blade] = parsed_1basis[i+1];
                break;
            }
            continue;
        }

        var blade_name: []const u8 = "";
        var blade_metric = BasisMetric.positive;
        
        for(0.._d) |i| {
            if(blade & (1 << i) == 0) {
                continue;
            }
            blade_name = blade_name ++ .{parsed_1basis[i+1].name[0]};
            const factor_metric: BasisMetric = parsed_1basis[i+1].square;
            blade_metric = BasisMetric.from_square(factor_metric.square() * blade_metric.square());
        }
        blade_metric = BasisMetric.from_square(blade_metric.square() * (-2 * @as(isize, bit_utils.flip_parity(blade, blade)) + 1));
        const new_blade: BasisVector = .{.blade = blade, .name = blade_name, .square = blade_metric};
        parsed_basis[blade] = new_blade;
    }

    //const pseudoscalar_name: []const u8 = parsed_basis[basis_size-1].name;
    //association of [char in _basis] : [index of that char in all blades in the binary blade ordering scheme]
    //also includes scalar basis char 's' as 0
    const binary_ordered_chars: [_d+1]u8 = blk: {
        var ret: [_d+1]u8 = .{'!'} ** (_d + 1);
        for(0..parsed_1basis.len) |i| {
            ret[i] = parsed_1basis[i].name[0];
        }
        break :blk ret;
    };

    // the index of each basis blade is their binary blade rep, but names are flipped if the user flips them
    var _canonical_basis_with_flips: [basis_size]BasisVector = parsed_basis;
    var ordered_basis: [basis_size]BasisVector = parsed_basis;
    //ordering = "s,xy,yz,zx,x,xyz,y,z"
    const OrderingParseStates = enum {
        ParseFirstBasisChar,
        ParseNthBasisChar,
    };
    curr_basis_blade = undefined;
    curr_basis_i = 0;
    var curr_order_parse_state: OrderingParseStates = .ParseFirstBasisChar;
    if(_ordering != null and _ordering.?.len > 0) {
        const ordering = _ordering orelse unreachable;
        //i need to go through the ordering list and parse out every full blade
        //for each one, find the blade matching every character in its name in parsed_basis
        //if there are switched chars, find how many flips that represents and adjust the square

        for(0..ordering.len+1) |i| {
            switch(curr_order_parse_state) {
                //step 1: parse out the next blade, move to step 2 when we encounter an end (, or i == ordering.len)
                //i cannot be at the end here
                .ParseFirstBasisChar => {
                    if(i == ordering.len) {
                        @compileError(comptimePrint("something went wrong with parsing ordering, basis: {any} ordering: {any}", .{_basis, ordering}));
                    }
                    const curr_char = ordering[i];
                    if(curr_char == ' ') {
                        continue;
                    }
                    if(curr_char == ',') {
                        @compileError(comptimePrint("ordering string had a double ,,! ordering: {any}", .{ordering}));
                    }
                    if((curr_char < 'A' or curr_char > 'z') or (curr_char > 'Z' and curr_char < 'a')) {
                        @compileError(comptimePrint("non letter char passed to the order string! ordering: {any}", ordering));
                    }

                    curr_basis_blade = .{.blade = undefined, .name = "", .square = .positive};
                    curr_basis_blade.name = curr_basis_blade.name ++ .{ordering[i]};
                    //at this point it is legal for the current name to end, when it does (parsenth will detect this)
                    //we will find out what blade it belongs to
                    curr_order_parse_state = .ParseNthBasisChar;
                },
                .ParseNthBasisChar => {
                    if(i == ordering.len or ordering[i] == ',') {
                        //we have the full name, now find the blade, and how many flips our name has rel. to the canonical name
                        var blade: usize = 0;
                        const name: []const u8 = curr_basis_blade.name;
                        if(name.len == 0) {
                            @compileError(comptimePrint("empty name found when parsing ordering string! ordering: {any}", .{ordering}));
                        }

                        var name_indices: [name.len]usize = chars_to_binary_ordering_indices(name.len, name, &binary_ordered_chars);
                        //create the blade
                        for(0..name.len) |name_i| {
                            const blade_idx: usize = name_indices[name_i];
                            if(blade_idx == 0) {
                                if(name.len > 1) {
                                    @compileError(comptimePrint("somehow scalar char s got into the ordering string for some blade, name: {s} ordering {s}", .{name, ordering}));
                                }
                                break;
                            }
                            blade |= (1 << (blade_idx - 1)); //s is idx 0
                        }
                        const flips: usize = cycle_sort_ret_writes(usize, &name_indices);
                        const parity: isize = bit_utils.unsigned_to_signed_parity(flips & 1);
                        if(parity < 0) {
                            //const old_flipped = _flipped_parities;
                            _flipped_parities |= (1 << blade);
                            //@compileLog(comptimePrint("flipped blade {s} found to be binary blade {b} with {} flips, _flipped_parities {b} -> {b}, ordering {s}", .{name, blade, flips, old_flipped, _flipped_parities, ordering}));
                        }
                        curr_basis_blade.blade = blade;
                        curr_basis_blade.square = parsed_basis[blade].square;
                        
                        curr_order_parse_state = .ParseFirstBasisChar;

                        ordered_basis[curr_basis_i] = curr_basis_blade;
                        _canonical_basis_with_flips[blade] = curr_basis_blade;
                        curr_basis_i += 1;

                    } else { //parse the next character of the name
                        const curr_char: u8 = ordering[i];
                        if((curr_char < 'A' or curr_char > 'z') or (curr_char > 'Z' and curr_char < 'a')) {
                            @compileError(comptimePrint("non letter char passed to the order string! ordering: {any}", ordering));
                        }
                        if(curr_char == 's' and curr_basis_blade.name.len != 0) {
                            @compileError(comptimePrint("s added to nonzero len blade name! ordering: {s}", .{ordering}));
                        }
                        var in_1basis: bool = false;
                        for(0..parsed_1basis.len) |j| {
                            if(parsed_1basis[j].name[0] == curr_char) {
                                in_1basis = true;
                                break;
                            }
                        }
                        if(in_1basis == false) {//TODO: add to other state
                            @compileError(comptimePrint("char in ordering string blade not in 1 basis! char {s} ordering {s} parsed_1basis: {any}", .{&.{curr_char}, ordering, parsed_1basis}));
                        }
                        curr_basis_blade.name = curr_basis_blade.name ++ .{curr_char};
                    }
                }
            }
        }
        if(curr_basis_i != basis_size) {
            @compileError(comptimePrint("non null ordering string doesnt include all basis elements! ordering {s}, number parsed: {}, ordered_basis: {any}", .{ordering, curr_basis_i, ordered_basis[0..curr_basis_i]}));
        }

    }

    var _compiled_aliases: [_aliases.len]GAlgebraCompiledAlias = undefined;
    var _aliases_by_subset: [basis_size]?GAlgebraCompiledAlias = undefined;
    for(0..basis_size) |i| {
        _aliases_by_subset[i] = null;
    }

    for(0.._aliases.len) |i| {
        const existing: GAlgebraAlias = _aliases[i];
        //var ally: GAlgebraCompiledAlias = .{.name = existing.name, .order_string = existing.ordering, .ordering = undefined};
        var compiled_values: []const usize = &.{};
        var alias_subset: usize = 0;
        
        curr_basis_blade = undefined;
        curr_basis_i = 0;
        curr_order_parse_state = .ParseFirstBasisChar;
        if(existing.ordering.len > 0) {
            const ordering = existing.ordering;

            for(0..ordering.len+1) |j| {
                switch(curr_order_parse_state) {
                    //step 1: parse out the next blade, move to step 2 when we encounter an end (, or i == ordering.len)
                    //i cannot be at the end here
                    .ParseFirstBasisChar => {
                        if(j == ordering.len) {
                            @compileError(comptimePrint("something went wrong with parsing ordering, basis: {any} ordering: {any}", .{_basis, ordering}));
                        }
                        const curr_char = ordering[j];
                        if(curr_char == ' ') {
                            continue;
                        }
                        if(curr_char == ',') {
                            @compileError(comptimePrint("ordering string had a double ,,! ordering: {any}", .{ordering}));
                        }
                        if((curr_char < 'A' or curr_char > 'z') or (curr_char > 'Z' and curr_char < 'a')) {
                            @compileError(comptimePrint("non letter char passed to the order string! ordering: {any}", ordering));
                        }

                        curr_basis_blade = .{.blade = 0, .name = "", .square = .positive};
                        curr_basis_blade.name = curr_basis_blade.name ++ .{ordering[j]};
                        //at this point it is legal for the current name to end, when it does (parsenth will detect this)
                        //we will find out what blade it belongs to
                        curr_order_parse_state = .ParseNthBasisChar;
                    },
                    .ParseNthBasisChar => {
                        if(j == ordering.len or ordering[j] == ',') {
                            //we have the full name, now find the blade. right now we ignore any flips relative to the algebras actual blade defined from _ordering, this will be changed in the future
                            var blade: usize = 0;
                            const name: []const u8 = curr_basis_blade.name;
                            if(name.len == 0) {
                                @compileError(comptimePrint("empty name found when parsing ordering string! ordering: {any}", .{ordering}));
                            }

                            const name_indices: [name.len]usize = chars_to_binary_ordering_indices(name.len, name, &binary_ordered_chars);
                            //create the blade
                            for(0..name.len) |name_j| {
                                const blade_idx: usize = name_indices[name_j];
                                if(blade_idx == 0) {
                                    if(name.len > 1) {
                                        @compileError(comptimePrint("somehow scalar char s got into the ordering string for some blade for a passed in alias, name: {s} ordering {s}", .{name, ordering}));
                                    }
                                    break;
                                }
                                blade |= (1 << (blade_idx - 1)); //s is idx 0
                            }
                            
                            curr_order_parse_state = .ParseFirstBasisChar;

                            alias_subset |= (1 >> blade);
                            const _blade = blade;
                            compiled_values = compiled_values ++ .{_blade};

                            //@compileLog(comptimePrint("parsing alias i {} j {} char: {c}: {any}. parsed compiled values into {any} with subset {b}. blade name {s}", .{i, j, if (j < ordering.len) ordering[j] else 0, existing, compiled_values, blade, curr_basis_blade.name}));

                            curr_basis_i += 1;

                        } else { //parse the next character of the name
                            const curr_char: u8 = ordering[j];
                            if((curr_char < 'A' or curr_char > 'z') or (curr_char > 'Z' and curr_char < 'a')) {
                                @compileError(comptimePrint("non letter char passed to the order string! ordering: {any}", ordering));
                            }
                            if(curr_char == 's' and curr_basis_blade.name.len != 0) {
                                @compileError(comptimePrint("s added to nonzero len blade name! ordering: {s}", .{ordering}));
                            }
                            var in_1basis: bool = false;
                            for(0..parsed_1basis.len) |k| {
                                if(parsed_1basis[k].name[0] == curr_char) {
                                    in_1basis = true;
                                    break;
                                }
                            }
                            if(in_1basis == false) {//TODO: add to other state
                                @compileError(comptimePrint("char in ordering string blade not in 1 basis! char {s} ordering {s} parsed_1basis: {any}", .{&.{curr_char}, ordering, parsed_1basis}));
                            }
                            curr_basis_blade.name = curr_basis_blade.name ++ .{curr_char};
                        }
                    }
                }
            }
        }
        const calias_subset: usize = alias_subset;
        const ccompiled_values: []const usize = compiled_values;

        const compiled_alias: GAlgebraCompiledAlias = .{.name = existing.name, .ordering = existing.ordering, .values = ccompiled_values};

        _compiled_aliases[i] = compiled_alias;
        _aliases_by_subset[calias_subset] = compiled_alias;
    }

    const __compiled_aliases = _compiled_aliases;
    const __aliases_by_subset = _aliases_by_subset;

    //@compileLog(comptimePrint("_basis: {s}", .{_basis}));
    //@compileLog(comptimePrint("parsed 1 basis: {s}", .{parsed_1basis}));
    //@compileLog(comptimePrint("parsed basis: {s}", .{parsed_basis}));

    //if(_ordering) |ordering| {
    //    @compileLog(comptimePrint("_ordering {s}", .{ordering}));
    //} else {
    //    @compileLog(comptimePrint("_ordering: {?}", .{_ordering}));
    //}
    //@compileLog(comptimePrint("ordered basis: {any}", .{ordered_basis}));
    //@compileLog(comptimePrint("canonical basis with flips: {any}", .{_canonical_basis_with_flips}));

    const sig: GAlgebraSignature = .{
        .p = _p,
        .n = _n,
        .z = _z,
        .d = _d,
        .basis_size = basis_size,
        .pseudoscalar = basis_size - 1,
        .pos_mask = _p_mask,
        .neg_mask = _n_mask,
        .zero_mask = _z_mask,
        .flipped_parities = _flipped_parities,
        .dual = _dual,
    };

    const _parsed_1basis = parsed_1basis;
    const _parsed_basis = parsed_basis;
    const _ordered_basis = ordered_basis;
    const __canonical_basis_with_flips = _canonical_basis_with_flips;
    const ret = struct {
        const Self = @This();
        /// so that users of this can recreate it in type functions using us as a parameter without 
        /// passing in every arg manually or using us as just being type as : type
        /// so that the lang server can see every field in the struct (i hate this part of zig)
        pub const GAlgebraInfoArgs = struct {
            pub const __dual: bool = _dual;
            pub const __basis: []const u8 = _basis;
            pub const __ordering: ?[]const u8 = _ordering;
            pub const __aliases: []const GAlgebraAlias = _aliases;
        };

        ///i wish zig had better ways of constraining types to those returned by a type fn
        pub const IsGAlgebraInfo: bool = true;

        pub const signature: GAlgebraSignature = sig;
        pub const basis_len = basis_size;
        pub const dual = _dual;
        pub const basis_1_blades = _parsed_1basis;
        /// canonical binary ordering, blades cannot be flipped
        /// TODO: do i even want this one
        pub const canonical_basis = _parsed_basis;
        /// custom ordering and blades can be flipped
        pub const ordered_basis_blades = _ordered_basis;
        /// binary ordering but blades can be flipped
        pub const canonical_basis_with_flips = __canonical_basis_with_flips;
        pub const aliases = __compiled_aliases;
        pub const aliases_by_subset = __aliases_by_subset;

        //cayley space of n-ary blade function = 
        //  operand blades -> result blade: rank n Tensor(BladeTerm, &.{basis_size, basis_size,...}),
        //  result blade mults -> rank 1 Tensor(T, &.{basis_size}),
        //  result blade -> operand blades: rank n Tensor(usize, )
        //BladeMult = {value: T, blade: usize}
        //ProdTerm = {mult: T, opr1: usize, opr2: usize, ..., oprn-ary: usize}
        pub fn BinaryCayleyReturn(T: type) type {
            return struct { 
                operands_to_results: [basis_size][basis_size]BladeMult(T), 
                ///1st index = result blade, 2nd index is 0..inv_terms[result blade] independent terms that add to this result blade
                /// ProdTerm = mult{T, opr1 blade, opr2 blade, ..., oprn-ary blade}
                results_to_terms: [basis_size][basis_size]ProdTerm(T), 
                results_to_num_terms: [basis_size]usize,

                pub fn format(
                    slf: @This(),
                    comptime _: []const u8,
                    _: std.fmt.FormatOptions,
                    writer: anytype
                ) !void {
                    _=try writer.print("BinaryCayleyReturn (basis_size {}) {{", .{
                        basis_size
                    });
                    const OperandsDisplayType = struct {
                        lhs: usize,
                        rhs: usize,
                        blade: usize,
                        value: T,
                        left_value: T,
                        right_value: T,

                        pub fn format(
                            slf2: @This(),
                            comptime _: []const u8,
                            _: std.fmt.FormatOptions,
                            writer2: anytype
                        ) !void {
                            _ = try writer2.print("{d:.1} {b} * {d:.1} {b} = {d:.2} {b}", .{slf2.left_value, slf2.lhs, slf2.right_value, slf2.rhs, slf2.value, slf2.blade});
                        }
                    };
                    var operands_arr: [basis_size * basis_size]OperandsDisplayType = undefined;
                    const ResultsToTermsDisplayType = struct {
                        result_blade: usize, 
                        terms: []const ProdTerm(T),
                        
                        pub fn format(
                            slf2: @This(),
                            comptime _: []const u8,
                            _: std.fmt.FormatOptions,
                            writer2: anytype
                        ) !void {
                            _ = try writer2.print("{b} = {any}", .{slf2.result_blade, slf2.terms});
                        }
                    };
                    var results_arr: [basis_size]ResultsToTermsDisplayType = undefined;
                    
                    for(0..basis_size) |i| {
                        for(0..basis_size) |j| {
                            const idx = i + basis_size * j;
                            operands_arr[idx] = .{
                                .lhs = i, 
                                .rhs = j, 
                                .blade = slf.operands_to_results[i][j].blade, 
                                .value = slf.operands_to_results[i][j].value,
                                .left_value = slf.operands_to_results[i][j].left_value,
                                .right_value = slf.operands_to_results[i][j].right_value,
                                
                                };
                        }

                        results_arr[i] = .{.result_blade = i, .terms = slf.results_to_terms[i][0..slf.results_to_num_terms[i]]};
                    }
                    

                    _=try writer.print("operands_to_results ({}): {any}, results_to_terms: {any}, results_to_num_terms: {any} }}", .{
                        slf.operands_to_results.len, operands_arr, results_arr, slf.results_to_num_terms
                    });
                }
            };
        }

        //step 1: create the masked cayley space
        //if this is generic to the transformation we need to tell it how to do each operand blade n-tuple
        pub fn binary_cayley_table(T: type, op: anytype, lhs_subset: usize, rhs_subset: usize) BinaryCayleyReturn(T) {
            @setEvalBranchQuota(10000000);
            //inputs to outputs under op = the image, unsimplified
            // nth axis here is the nth operand acted upon by op. the value stored here is the value returned by op given those operands
            //       lhs blade    rhs blade  = ret blade
            var ret: [basis_size][basis_size]BladeMult(T) = .{.{.{.blade = 0, .value = 0.0, .left_value = 0, .right_value = 0}} ** basis_size} ** basis_size;
            //outputs to inputs under op = the pre image
            // the first axis is always the output, the second axis contains up to basis_len values containing what operand blades contributed to them
            //       result blade =  terms.. mult, lhs blade, rhs blade
            var inv: [basis_size][basis_size]ProdTerm(T) = .{.{.{.mult = 0, .left_blade = 0, .right_blade = 0, .left_value = 0, .right_value = 0}} ** basis_size} ** basis_size;

            var inv_indices: [basis_size]usize = .{0} ** basis_size;

            for (0..basis_size) |lhs| {
                const lhs_value: T = if((1 << lhs) & lhs_subset != 0) 1 else 0;
                for (0..basis_size) |rhs| {

                    const rhs_value: T = if((1 << rhs) & rhs_subset != 0) 1 else 0;

                    var res = op.eval(
                        .{ .value = lhs_value, .blade = lhs}, 
                        .{ .value = rhs_value, .blade = rhs}
                    );
                    res.left_value = lhs_value;
                    res.right_value = rhs_value;
                    ret[lhs][rhs] = res;

                    if(res.value == 0 and ((1 << lhs) & lhs_subset) != 0 and ((1 << rhs) & rhs_subset) != 0) {
                        //@compileLog(comptimePrint("binary_cayley (EARLY RET res.value == 0): {d} {b} in subset {b} * {d} {b} in subset {b} = {d} {b}", .{lhs_value, lhs, lhs_subset, rhs_value, rhs, rhs_subset, res.value, res.blade}));
                    }
                    if (res.value == 0 or ((1 << lhs) & lhs_subset) == 0 or ((1 << rhs) & rhs_subset) == 0) {
                        //@compileLog(comptimePrint("binary_cayley (EARLY RET LAST): {d} {b} in subset {b} * {d} {b} in subset {b} = {d} {b}", .{lhs_value, lhs, lhs_subset, rhs_value, rhs, rhs_subset, res.value, res.blade}));
                        continue;
                    }
                    //@compileLog(comptimePrint("binary_cayley_table(lhs {b}, rhs {b}) = {?}", .{lhs, rhs, res}));
                    //inv[result term blade][0..n] = .{magnitude coeff, lhs operand, rhs operand}
                    // where 0..inv_indices[res.blade] are the filled in indices
                    inv[res.blade][inv_indices[res.blade]] = .{ 
                        .mult = res.value, 
                        .left_blade = lhs, 
                        .right_blade = rhs,
                        .left_value = lhs_value,
                        .right_value = rhs_value,
                    };
                    //@compileLog(comptimePrint("binary_cayley (INVERSE SET term {}): {d} {b} in subset {b} * {d} {b} in subset {b} = {d} {b}", .{inv_indices[res.blade], lhs_value, lhs, lhs_subset, rhs_value, rhs, rhs_subset, res.value, res.blade}));
                    inv_indices[res.blade] += 1;
                }
            }
            //@compileLog(comptimePrint("ops_to_results {any},,,,,,,,,,,, results_to_terms {any}", .{ret, inv}));
            return .{ .operands_to_results = ret, .results_to_terms = inv, .results_to_num_terms = inv_indices };
        }


        pub fn UnaryCayleyReturn(T: type) type {
            return struct {
                operands_to_results: [basis_size]BladeMult(T),
                results_to_operands: [basis_size]UnaryProdTerm(T)
            };
        }

        pub fn unary_cayley_table(T: type, op: anytype, source_subset: usize) UnaryCayleyReturn(T) {

            var ret: [basis_size]BladeMult(T) = .{.{.blade = 0, .value = 0}} ** basis_size;
            var inv: [basis_size]UnaryProdTerm(T) = .{.{.mult = 0, .source = 0}} ** basis_size;

            for (0..basis_size) |blade| {
                //const rhs_curr_idx = count_bits(rhs_subset & ((0b1 << rhs) - 1));
                const value: T = if((1 << blade) & source_subset != 0) 1 else 0;
                if(value == 0) {
                    continue;
                }

                const res = op.eval(.{ .value = value, .blade = blade });
                ret[blade] = res;
                if (res.value == 0) {
                    continue;
                }
                //inv[result term blade][0..n] = .{magnitude coeff, lhs operand, rhs operand}
                // where 0..n are the filled in indices, n = inv_indices[]
                inv[res.blade] = .{ .mult = res.value, .source = blade };
            }

            return .{ .operands_to_results = ret, .results_to_operands = inv};
        }


    };

    const _ret = ret;
    return _ret;
}