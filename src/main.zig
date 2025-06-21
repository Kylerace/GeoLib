const std = @import("std");
const debug = std.debug;
const builtin = @import("builtin");
const meta = std.meta;
const comptimePrint = std.fmt.comptimePrint;

//const mecha = @import("mecha");

//to investigate an alternative structure of library design where "logical" types representing
//an abstract geometric algebra are decoupled from "data / implementation" types that hold runtime data


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

fn cosh(num: anytype) @TypeOf(num) {
    return 0.5 * (@exp(num) + @exp(-num));
}

fn sinh(num: anytype) @TypeOf(num) {
    return 0.5 * (@exp(num) - @exp(-num));
}

//TODO: do i want these in the Blade type?
fn count_bits(v: anytype) @TypeOf(v) {
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

fn count_flips(left: anytype, rhs: anytype) @TypeOf(left) {
    var lhs = left >> 1;
    var flips: @TypeOf(left) = 0;
    while (lhs != 0) {
        flips += count_bits(lhs & rhs);
        lhs >>= 1;
    }
    //flips += count_bits((left >> (1 + i)) + rhs)
    return flips;
}

fn flip_parity(left: anytype, rhs: @TypeOf(left)) @TypeOf(left) {
    return count_flips(left, rhs) & 1;
}

fn flips_permuted_parity(flip_field: anytype, lhs: anytype, rhs: @TypeOf(lhs)) @TypeOf(lhs) {

    return count_flips(lhs,rhs) ^ (flip_field >> lhs) ^ (flip_field >> rhs) ^ (flip_field >> (lhs ^ rhs)) & 1;
}

fn flips_permuted_neg_parity(neg_mask: anytype, flip_field: @TypeOf(neg_mask), lhs: anytype, rhs: @TypeOf(lhs)) @TypeOf(lhs) {
    return (count_bits(neg_mask & (lhs & rhs)) ^ count_flips(lhs,rhs) ^ (flip_field >> lhs) ^ (flip_field >> rhs) ^ (flip_field >> (lhs ^ rhs))) & 1;
}

fn unsigned_to_signed_parity(unsigned_parity: anytype) isize {
    return (-2 * @as(isize, unsigned_parity & 1) + 1);
}

inline fn is_comptime(val: anytype) bool {
    return @typeInfo(@TypeOf(.{val})).Struct.fields[0].is_comptime;
}

///returns the type of either x if its not a pointer or x's pointed type if it is, and whether it is a ptr or not
pub fn UnwrapPtrType(comptime x: type) struct {type, bool} {
    const x_info = @typeInfo(x);
    const x_is_ptr = x_info == .Pointer;
    const x_type = blk: {
        if(x_is_ptr) {
            break :blk x_info.Pointer.child;
        }
        break :blk x;
    };
    return .{x_type, x_is_ptr};
}

///handles logical operations given an algebraic signature expressed as pos,negative, and zero squaring basis bitfields (must add to 2^n)
/// only needs to know the masks of the algebra, + the included 
//pub fn GAlgebraOps(comptime zero_mask: usize, comptime pos_mask: usize, comptime neg_mask: usize) type {
//const GAlgebraOps = enum {
//    const full_subset: usize = pos_mask | neg_mask | zero_mask;
// 
//    const basis_len: usize = count_bits(full_subset);
//    if(std.math.log2_int(usize, basis_len) != std.math.log2_int_ceil(usize, basis_len)) {
//        @compileError(comptimePrint("non power of 2 number of 1s in mask bitfields passed to GAlgebraOps! p {b} {} n {b} {} z {b} {}", .{pos_mask,pos_mask,neg_mask,neg_mask,zero_mask,zero_mask}));
//    }
//    const d: usize = std.math.log2_int(usize, basis_len);
//    const pseudoscalar: usize = blk: {
//        var val: usize = 0;
//        for(0..d) |i| {
//            val |= (1 << i);
//        }
//        break :blk val;
//    };
// 
//    const ret = struct {
//        const Self = @This();
// 
//        pub const BladeTerm = struct {
//            ///1 where the ith 1 blade is included, 0 elsewhere
//            /// count_bits(blade)th form
//            blade: usize,
//            value: f64
//        };
// 
//    };
//    return ret;
//}

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
const GAlgebraSignature = struct {
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
                        const value: isize = unsigned_to_signed_parity(flips_permuted_neg_parity(neg_mask, flip_field, lhs, rhs));
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
                        return .{.value = unsigned_to_signed_parity(flips_permuted_parity(sig.flipped_parities, lhs, rhs)), .blade = blade};
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
                        const value: isize = unsigned_to_signed_parity(flips_permuted_neg_parity(neg_mask, flip_field, lhs, rhs));
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
                        
                        return .{ .value = unsigned_to_signed_parity(flips_permuted_neg_parity(neg_mask, flip_field, lhs, rhs)), .blade = blade };
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
                        return .{ .value = us.value * unsigned_to_signed_parity(flips_permuted_neg_parity(sig.neg_mask, sig.flipped_parities, us, sig.pseudoscalar)), .blade = us.blade ^ sig.pseudoscalar };
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
                            var x: usize = count_bits(blade) & ~one;
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
const GAlgebraAlias = struct {
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
const GAlgebraCompiledAlias = struct {
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
const BasisVector = struct {
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
        const blade_size: usize = count_bits(blade);
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
        blade_metric = BasisMetric.from_square(blade_metric.square() * (-2 * @as(isize, flip_parity(blade, blade)) + 1));
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
                        const parity: isize = unsigned_to_signed_parity(flips & 1);
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

/// this is what MVecArgs is processed into within MVecSubset
const MVecProperties = struct {
    values: []const usize,
    name: ?[]const u8,
    compiled_alias: ?GAlgebraCompiledAlias,
};

/// the type given to MVecSubset to inform what properties it has within the given algebra.
/// values and subset are mutually exclusive
/// name is compatible with anything except alias
pub const MVecArgs = struct {
    const Self: type = @This();
    /// if true, the MVecSubset attempts to find an alias in the algebra matching whatever properties have been defined.
    /// if false, we will only be set to an alias if one of the alias-specific properties have been set (e.g. alias, raw_alias, or alias_name)
    find_alias: bool = true,

    /// array of blade usize fields specifying what blades exist in this subset in what exact order
    values: ?[]const usize = null,
    /// usize bitfield where bit i == 1 iff this MVec contains the ith blade of its algebra (in the canonical binary ordering of blades). doesnt allow for reordered indices of course.
    /// if find_alias == true, then this will find the first alias in the algebra with the same subset
    subset: ?usize = null,
    /// if set, searches through the algebra for an alias matching the given name and sets the MVec's properties to that if found, erroring if not
    alias_name: ?[]const u8 = null,
    /// if set, searches through the algebra for the corresponding compiled alias for the given raw alias and sets the MVec's properties to that if found, erroring if not
    alias: ?GAlgebraAlias = null,
    /// if set, sets properties of this mvec to those specified by the given alias. 
    compiled_alias: ?GAlgebraCompiledAlias = null,
    /// the name the mvec subset will be printed with. overrides any alias if BOTH set. note that if this is not set and there is no alias then the MVec's name will be "Multivector"
    name: ?[]const u8 = null,

    /// checks if the passed in compiled alias is actually from the given algebra
    fn try_find_alias_from_compiled_alias(comptime _alg_info: type, alias: GAlgebraCompiledAlias) ?GAlgebraCompiledAlias {
        const info_args: type = _alg_info.GAlgebraInfoArgs;
        const alg_info: type = GAlgebraInfo(info_args.__dual, info_args.__basis, info_args.__ordering, info_args.__aliases);

        for(0..alg_info.aliases.len) |i| {
            const alias_i: GAlgebraCompiledAlias = alg_info.aliases[i];
            if(alias_i.name != null and alias.name != null) {
                if(std.mem.eql(u8, alias_i.name orelse unreachable, alias.name orelse unreachable) == false) {
                    continue;
                }
            } else if ((alias_i.name == null) != (alias.name == null)) {
                continue;
            }

            if(std.mem.eql(u8, alias_i.ordering, alias.ordering) == false) {
                continue;
            }

            if(std.mem.eql(usize, alias_i.values, alias.values)) {
                return alias_i;
            }
        }
    }

    fn try_find_alias_from_alias_name(comptime _alg_info: type, alias_name: []const u8) ?GAlgebraCompiledAlias {
        const info_args: type = _alg_info.GAlgebraInfoArgs;
        const alg_info: type = GAlgebraInfo(info_args.__dual, info_args.__basis, info_args.__ordering, info_args.__aliases);

        for(0..alg_info.aliases.len) |i| {
            const alias_i: GAlgebraCompiledAlias = alg_info.aliases[i];
            if(alias_i.name) |alias_i_name| {
                if(std.mem.eql(u8, alias_i_name, alias_name)) {
                    return alias_i;
                }
            } 

        }
        return null;

    }

    fn try_find_alias_from_raw_alias(comptime _alg_info: type, alias: GAlgebraAlias) ?GAlgebraCompiledAlias {
        const info_args: type = _alg_info.GAlgebraInfoArgs;
        const alg_info: type = GAlgebraInfo(info_args.__dual, info_args.__basis, info_args.__ordering, info_args.__aliases);

        for(0..alg_info.aliases.len) |i| {
            const alias_i: GAlgebraCompiledAlias = alg_info.aliases[i];
            if(alias_i.name and alias.name) {
                if(std.mem.eql(u8, alias_i.name orelse unreachable, alias.name orelse unreachable) == false) {
                    continue;
                }
            } else if ((alias_i.name == null) != (alias.name == null)) {
                continue;
            }

            if(std.mem.eql(u8, alias_i.ordering, alias.ordering)) {
                return alias_i;
            }
        }
        return null;

    }

    fn try_find_alias_from_name(comptime _alg_info: type, name: []const u8) ?GAlgebraCompiledAlias {
        const info_args: type = _alg_info.GAlgebraInfoArgs;
        const alg_info: type = GAlgebraInfo(info_args.__dual, info_args.__basis, info_args.__ordering, info_args.__aliases);

        for(0..alg_info.aliases.len) |i| {
            const alias_i: GAlgebraCompiledAlias = alg_info.aliases[i];
            if(alias_i.name) |alias_name| {
                if(std.mem.eql(u8, name, alias_name)) {
                    return alias_i;
                }
            }
        }
        return null;
    }

    fn try_find_alias_from_subset(comptime _alg_info: type, subset: usize) ?GAlgebraCompiledAlias {
        const info_args: type = _alg_info.GAlgebraInfoArgs;
        const alg_info: type = GAlgebraInfo(info_args.__dual, info_args.__basis, info_args.__ordering, info_args.__aliases);

        return alg_info.aliases_by_subset[subset];
    }

    fn try_find_alias_from_values(comptime _alg_info: type, values: []const usize) ?GAlgebraCompiledAlias {
        const info_args: type = _alg_info.GAlgebraInfoArgs;
        const alg_info: type = GAlgebraInfo(info_args.__dual, info_args.__basis, info_args.__ordering, info_args.__aliases);
        
        for(0..alg_info.aliases.len) |i| {
            const alias_i: GAlgebraCompiledAlias = alg_info.aliases[i];
            const alias_values: []const usize = alias_i.values;
            if(std.mem.eql(usize, values, alias_values)) {
                return alias_i;
            }
        }
        return null;
    }

    pub fn set_ret_from_alias(ret: *MVecProperties, alias: GAlgebraCompiledAlias) void {
        ret.name = alias.name;
        ret.compiled_alias = alias;
        ret.values = alias.values;
    }

    pub fn finalize(comptime self: Self, comptime _alg_info: anytype) MVecProperties {
        const info_args: type = _alg_info.GAlgebraInfoArgs;
        const alg_info: type = GAlgebraInfo(info_args.__dual, info_args.__basis, info_args.__ordering, info_args.__aliases);

        var ret: MVecProperties = .{.values = &.{}, .name = null, .compiled_alias =  null};
        const LatestSet = enum {
            const Us = @This();

            none,
            values,
            subset,
            alias_name,
            alias,
            compiled_alias,

            pub fn lower_precedence_than(us: Us, other: Us) bool {
                return @intFromEnum(us) < @intFromEnum(other);
            }

            pub fn format(
                enm: @This(),
                comptime _: []const u8,
                _: std.fmt.FormatOptions,
                writer: anytype
            ) !void {
                _=try writer.print("{s}", .{
                    @tagName(enm)
                });
            }
        };
        var latest_set: LatestSet = .none;

        if(self.compiled_alias) |compiled_alias| {
            latest_set = .compiled_alias;
            //alias specific property thus doesnt care about find_alias
            const returned_alias: ?GAlgebraCompiledAlias = try_find_alias_from_compiled_alias(alg_info, compiled_alias);
            if(returned_alias) |confirmed_alias| {
                set_ret_from_alias(&ret, confirmed_alias);
            } else {
                @compileError(comptimePrint("MVecArgs could not be finalized because the compiled_alias was set to a GAlgebraCompiledAlias that didn't match anythign in the algebra! alias: {any}", .{compiled_alias}));
            }
        }

        if(self.alias) |alias| {
            if(LatestSet.lower_precedence_than(.alias, latest_set)) {
                @compileError(comptimePrint("MVecArgs had alias set when {} was already set! These are mutually exclusive. self: {any}", .{latest_set, self}));
            }
            latest_set = .alias;
            //alias specific property thus doesnt care about find_alias
            const returned_alias: ?GAlgebraCompiledAlias = try_find_alias_from_raw_alias(_alg_info, alias);
            if(returned_alias) |confirmed_alias| {
                set_ret_from_alias(&ret, confirmed_alias);
            } else {
                @compileError(comptimePrint("MVecArgs could not be finalized because alias didn't match any alias in the algebra! alias: {any}", .{alias}));
            }
        }

        if(self.alias_name) |alias_name| {
            if(LatestSet.lower_precedence_than(.alias_name, latest_set)) {
                @compileError(comptimePrint("MVecArgs had alias_name set when {} was already set! These are mutually exclusive. self: {any}", .{latest_set, self}));
            }
            latest_set = .alias_name;
            //alias specific property thus doesnt care about find_alias
            const returned_alias: ?GAlgebraCompiledAlias = try_find_alias_from_alias_name(_alg_info, alias_name);
            if(returned_alias) |confirmed_alias| {
                set_ret_from_alias(&ret, confirmed_alias);
            } else {
                @compileError(comptimePrint("MVecArgs could not be finalized because alias_name didn't match any alias in the algebra! name: {s}", .{alias_name}));
            }
        }

        if(self.values) |values| {
            if(LatestSet.lower_precedence_than(.values, latest_set)) {
                @compileError(comptimePrint("MVecArgs had values set when {} was already set! These are mutually exclusive. self: {any}", .{latest_set, self}));
            }
            latest_set = .values;
            if(self.find_alias) {
                const returned_alias: ?GAlgebraCompiledAlias = try_find_alias_from_values(_alg_info, values);
                if(returned_alias) |confirmed_alias| {
                    set_ret_from_alias(&ret, confirmed_alias);
                } else {
                    @compileError(comptimePrint("MVecArgs could not be finalized because values didn't match any alias in the algebra! values: {any}", .{values}));
                }
            } else {
                ret.values = values;
            }
        }

        if(self.subset) |subset| {
            if(LatestSet.lower_precedence_than(.subset, latest_set)) {
                @compileError(comptimePrint("MVecArgs had subset set when {} was already set! These are mutually exclusive. self: {any}", .{latest_set, self}));
            }
            latest_set = .subset;
            if(self.find_alias) {
                const returned_alias: ?GAlgebraCompiledAlias = try_find_alias_from_subset(_alg_info, subset);
                if(returned_alias) |confirmed_alias| {
                    set_ret_from_alias(&ret, confirmed_alias);
                } else {
                    @compileError(comptimePrint("MVecArgs could not be finalized because subset didn't match any alias in the algebra! subset: {b}", .{subset}));
                }
            } else {
                var values: []const usize = &.{};
                //canonical (binary) ordering
                for(0..alg_info.basis_len) |blade| {
                    if(subset & (1 << blade) != 0) {
                        values = values ++ .{blade};
                    }
                }

                ret.values = values;
            }

        }

        if(self.name) |name| {
            ret.name = name;
        }

        return ret;
    }
};

// const Point: type = MVecSubset(alg, f64, .{.alias_name = "Point"});

/// impl type - dense subset of the full multivector of the algebra, with runtime and comptime known parts
/// constraint: dont allow flipped basis blades relative to the algebras basis blades - thats dumb
///     if the user tries that then just create a compile error
/// DO: allow reordered terms, its trivial at the mask creation stage to reoder them as needed for an op
/// 
/// _alg_info: the GAlgebraInfo object telling us the properties of our containing algebra e.g. the number of basis vectors
/// T: the datatype of our array elements. most likely a float
/// values: a slice of blades (1 on bit i iff this subset contains the ith blade of the given algebra). blades can be defined out of order from the canonical 
///  ordering of the algebra, which is the binary ordering counting up from 0 to (num basis vectors ** 2) - 1
pub fn MVecSubset(_alg_info: type, T: type, comptime args: MVecArgs) type {


    const info_args: type = _alg_info.GAlgebraInfoArgs;
    const alg_info: type = GAlgebraInfo(info_args.__dual, info_args.__basis, info_args.__ordering, info_args.__aliases);

    const properties: MVecProperties = args.finalize(_alg_info);
    const _values: []const usize = properties.values;


    const basis_size: usize = alg_info.signature.basis_size;

    var _basis_values: [_values.len]BasisVector = undefined;
    //var _basis_orders: [_values.len]struct {usize, BasisVector} = undefined;
    //we need an array of type binary blade -> real idx
    var _canonical_blade_to_idx: [basis_size]?usize = .{null} ** basis_size;

    var _subset: usize = 0;

    for(0.._values.len) |i| {
        const blade: usize = _values[i];

        if(blade >= alg_info.canonical_basis_with_flips.len) {
            @compileLog(comptimePrint("MVecSubset i: {} blade {b} alg_info.canonical_basis_with_flips: {any} _values.len {} _values {any}", .{i, blade, alg_info.canonical_basis_with_flips, _values.len, _values}));
        }
        const basis: BasisVector = alg_info.canonical_basis_with_flips[blade];
        //for(0..basis_size) |j| {
        //    const ordered_blade = alg_info.ordered_basis_blades[j];
        //    if(!basis.equals(ordered_blade)) {
        //        //@compileLog(comptimePrint("{} != {}", .{basis, ordered_blade}));
        //        continue;
        //    }
        //    _basis_orders[i] = .{j, basis}_values;
        //    _basis_values[i] = basis;
        //}
        _basis_values[i] = basis;

        //@compileLog(comptimePrint("i: {}, blade: {}, basis: {any}, _basis_values[i]: {}", .{i, blade, basis, _basis_values[i]}));
        //_basis_values = _basis_values ++ .{basis};
        _subset |= (1 << blade);
    }

    //@compileLog(comptimePrint("1: _values: {any}, _basis_values: {any}", .{_values, _basis_values}));

    //insertion sort
    //for(1.._values.len) |i| {
    //    //@compileLog(comptimePrint("_basis_orders: {s}, alg_info.ordered_basis_blades: {s}", .{_basis_orders, alg_info.ordered_basis_blades}));
    //    const elem: struct {usize, BasisVector} = _basis_orders[i];
    //    const key: usize = elem.@"0";
    //    var j: isize = i-1;
    //    while(j >= 0 and _basis_orders[j].@"0" > key) {
    //        _basis_orders[j + 1] = _basis_orders[j];
    //        j -= 1;
    //    }
    //    _basis_orders[j + 1] = elem;
    //}

    //@compileLog(comptimePrint("2: _values: {any}, _basis_orders: {any}, _basis_values: {any}", .{_values, _basis_orders, _basis_values}));

    for(0.._values.len) |i| {
        //_basis_values[i] = _basis_orders[i].@"1";
        _canonical_blade_to_idx[_basis_values[i].blade] = i;
    }

    //@compileLog(comptimePrint("3: _values: {any}, _basis_values: {any}", .{_values, _basis_values}));
    //const canonical_basis_with_flips: [basis_size]BasisVector = alg_info.canonical_basis_with_flips;
 
    const __name: []const u8 = blk: {
        if(properties.name) |name| {
            break :blk name;
        }
        break :blk "Multivector";
    };
    const __subset = _subset;
    const __basis_values = _basis_values;
    const __canonical_blade_to_idx = _canonical_blade_to_idx;
    const ret = struct {
        const Self = @This();
        pub const Algebra = alg_info;
        pub const Scalar = T;

        pub const subset = __subset;

        pub const basis_values: []const BasisVector = &__basis_values;
        /// partial function of sig binary blade -> idx in our terms array of that blade.
        /// null element means we don't contain the ith basis blade.
        pub const canonical_blade_to_idx: [basis_size]?usize = __canonical_blade_to_idx;
        pub const values_arg: []const usize = _values;

        pub const num_terms = count_bits(subset);
        pub const name: []const u8 = __name;

        ///runtime int needed to fit a single blade
        pub const ud: type = std.math.IntFittingRange(0, num_terms - 1);
        ///runtime int needed to fit our maximum subset
        pub const ux: type = std.math.IntFittingRange(0, std.math.pow(u128, 2, num_terms) - 1);

        pub const BladeEnum: type = blk: {
            var fields: [num_terms]std.builtin.Type.EnumField = undefined;
            var decls = [_]std.builtin.Type.Declaration{};
            for(0..num_terms) |i| {
                const basis_vector: BasisVector = basis_values[i];
                const field = &fields[i];
                field.* = .{.name = basis_vector.name, .value = i};
            }
            break :blk @Type(.{.Enum = .{
                .tag_type = usize,
                .fields = &fields,
                .decls = &decls,
                .is_exhaustive = true,
            }});
        };
        ///takes a blade enum and returns its index in the terms array
        pub fn blade_enum_to_index(blade: BladeEnum) usize {
            for(0..num_terms) |i| {
                const basis_vector: BasisVector = basis_values[i];
                const basis_name: []const u8 = basis_vector.name;
                const blade_name: [:0]const u8 = @tagName(blade);
                var all_same: bool = true;
                for(0..basis_name.len) |j| {
                    if(blade_name[j] == 0 or (j == basis_name.len - 1 and blade_name[j + 1] != 0)) {
                        all_same = false;
                        break;
                    }
                    if(basis_name[j] != blade_name[j]) {
                        all_same = false;
                        break;
                    }
                }
                if(all_same == true) {
                    return i;
                }
            }
        }

        pub const TermsType: type = @Vector(num_terms, T);
        pub const TermsArrayType: type = [num_terms]T;
        terms: TermsType,

        //pub fn blade_array_is_correct(comptime vals: [].{BladeEnum, T}) bool {
        //    var last_index: ?usize = null;
        //    for(0..vals.len) |i| {
        //        const blade_enum: BladeEnum = vals[i].@"0";
        //        if(last_index) |last_idx| {
        //            if(last_idx > ) //this is dumb
        //        } else {
        //            last_index = blade_enum_to_index(blade_enum);
        //        }
        //        if(blade_enum_to_index(blade_enum) < i) {
        //            return false;
        //        }
        //    }
        //    return true;
        //}
        // 
        //pub fn init(comptime vals: [].{Self.BladeEnum, T}) Self {
        //    if(vals.len > )
        //}

        pub fn init_raw(values: TermsArrayType) Self {
            return Self{.terms = values};
        }

        pub fn init_all(value: T) Self {
            return Self{.terms = @splat(value)};
        }

        pub fn init_counting_up_from(value: T) Self {
            var ret: Self = init_all(value);
            var val: T = value + 1;
            //var terms: [@typeInfo(@TypeOf(ret.terms)).Vector.len]T = ret.terms;
            for(1..@typeInfo(@TypeOf(ret.terms)).Vector.len) |i| {
                ret.terms[i] = val;
                val += 1; 
            }

            return ret;
        }

        pub const SetValsError = error {
            outOfBounds
        };
        ///TODO: make this able to use raw positions or blade enum specifiers (figure out a way to do that in general)
        pub fn set_vals(self: *Self, indices_values: []const struct {usize, T}) *Self {
            for(0..indices_values.len) |i| {
                const index: usize = indices_values[i].@"0";
                const value: T = indices_values[i].@"1";
                if(index < 0 or index > self.terms.len) {
                    continue;
                    //return error.outOfBounds;
                }
                self.terms[index] = value;
            }

            return self;
        }

        pub fn to_zeros(self: *Self) *Self {
            self.terms = @splat(0);
            return self;
        }

        pub fn to_ones(self: *Self) *Self {
            self.terms = @splat(1);
            return self;
        }

        pub fn to_x(self: *Self, val: T) *Self {
            self.terms = @splat(val);
            return self;
        }

        pub fn format(
            mvec: Self,
            comptime _: []const u8,
            _: std.fmt.FormatOptions,
            writer: anytype
        ) !void {
            _ = try writer.print("{s} {{", .{Self.name});
            for(0..num_terms) |i| {
                _ = try writer.print("{s}: {d:.3}", .{Self.basis_values[i].name, mvec.terms[i]});
                if(i < num_terms - 1) {
                    _ = try writer.print(", ", .{});
                }
            }
            _ = try writer.print("}}", .{});
        }



        pub fn BinaryScatterMasksReturn(comptime len: usize) type {
            return struct {
                lhs_masks: [basis_size][len]i32,
                rhs_masks: [basis_size][len]i32,
                coeff_masks: [basis_size][len]T
            };
        }

        pub fn binary_create_scatter_masks_from_cayley_table_inv(comptime inv: [basis_size][basis_size]ProdTerm(T), comptime inv_indices: [basis_size]usize, comptime lhs_type: type, comptime rhs_type: type, comptime res_type: type, comptime len: usize) BinaryScatterMasksReturn(len) {
            const lhs_subset: usize = lhs_type.subset;
            const rhs_subset: usize = rhs_type.subset;
            const superset: usize = res_type.subset;

            var lhs_masks: [basis_size][len]i32 = .{.{-1} ** len} ** basis_size;
            var rhs_masks: [basis_size][len]i32 = .{.{-1} ** len} ** basis_size;
            var coeff_masks: [basis_size][len]T = .{.{0} ** len} ** basis_size;
            var masks_idx: usize = 0;
            const max_terms = basis_size;

            for(0..max_terms) |ith_term| {
                var remaining_terms: bool = false;
                for (0..basis_size) |result_blade| {

                    if ((superset & (1 << result_blade)) != 0 and inv_indices[result_blade] >= ith_term) {
                        remaining_terms = true;
                        break;
                    }
                }
                if (remaining_terms == false) {
                    break;
                }
                defer masks_idx += 1;

                var lhs_mask: [len]i32 = .{-1} ** len;
                var rhs_mask: [len]i32 = .{-1} ** len;
                var coeff_mask: [len]T = .{0} ** len;

                for (0..basis_size) |result_blade| {
                    if (superset & (1 << result_blade) == 0) {
                        continue;
                    }
                    //result_blade value = sigma(ith_term from 0 to end)((mult)[ith_term])
                    const inv_term = inv[result_blade][ith_term];

                    const coeff = inv_term.mult;
                    const left_operand_blade = inv_term.left_blade;
                    const right_operand_blade = inv_term.right_blade;

                    //const lhs_canon_idx = count_bits(lhs_subset & ((1 << left_operand_blade) -% 1));
                    //const rhs_canon_idx = count_bits(rhs_subset & ((1 << right_operand_blade) -% 1));
                    const _lhs_idx = lhs_type.canonical_blade_to_idx[left_operand_blade];
                    const _rhs_idx = rhs_type.canonical_blade_to_idx[right_operand_blade];
                    if(_lhs_idx == null) {
                        continue;
                        //@compileError(comptimePrint("lhs blade {b} with unordered index {b} not found in subset {b} because its index is null! canonical_blade_to_idx: {any}, left_operand_blade: {b}, right_operand_blade: {b}, res blade: {b}, result subset: {b}, lhs_type: {s}, rhs_type: {s}, res_type: {s}", .{left_operand_blade, lhs_canon_idx, lhs_subset, lhs_type.canonical_blade_to_idx, left_operand_blade, right_operand_blade, result_blade, superset, lhs_type.name, rhs_type.name, res_type.name}));
                    }
                    if(_rhs_idx == null) {
                        continue;
                        //@compileError(comptimePrint("rhs blade {b} with unordered index {b} not found in subset {b} because its index is null! canonical_blade_to_idx: {any}, subset & (1 << rhs) = {}", .{right_operand_blade, rhs_canon_idx, rhs_subset, rhs_type.canonical_blade_to_idx, rhs_subset & (1 << right_operand_blade)}));
                    }
                    const lhs_idx: usize = _lhs_idx.?;
                    const rhs_idx: usize = _rhs_idx.?;

                    //rhs_subset = 10111, 
                    if(lhs_subset & (1 << left_operand_blade) == 0) {
                        @compileError(comptimePrint("lhs_subset {b} doesnt contain blade {b} that's supposed to be an operand from it! lhs {b} rhs {b} res {b}", .{lhs_subset, left_operand_blade, lhs_subset, rhs_subset, superset}));
                    }
                    if(rhs_subset & (1 << right_operand_blade) == 0) {
                        @compileError(comptimePrint("rhs_subset {b} doesnt contain blade {b} that's supposed to be an operand from it! lhs {b} rhs {b} res {b}", .{rhs_subset, right_operand_blade, lhs_subset, rhs_subset, superset}));
                    }

                    const res_canon_idx = count_bits(superset & ((1 << result_blade) -% 1));
                    const _res_idx = res_type.canonical_blade_to_idx[result_blade];
                    if(_res_idx == null) {
                        @compileError(comptimePrint("res blade {b} with unordered index {b} not found in subset {b} because its index is null! canonical_blade_to_idx: {any}", .{result_blade, res_canon_idx, superset, res_type.canonical_blade_to_idx}));
                    }
                    const res_idx: usize = _res_idx.?;
                    //@compileLog(comptimePrint("term {}: lhs {b}/{} of {any} * rhs {b}/{} of {any} = {} res {b} {} of {any}", .{ith_term, left_operand_blade, lhs_idx, lhs_type.canonical_blade_to_idx, right_operand_blade, rhs_idx, rhs_type.canonical_blade_to_idx, coeff, result_blade, res_idx, res_type.canonical_blade_to_idx}));
                    //@compileLog(comptimePrint("flipped parities: {b}", .{lhs_type.Algebra.signature.flipped_parities}));

                    lhs_mask[res_idx] = lhs_idx;
                    rhs_mask[res_idx] = rhs_idx;
                    coeff_mask[res_idx] = coeff;
                }

                @memcpy(lhs_masks[ith_term][0..len], lhs_mask[0..len]);
                @memcpy(rhs_masks[ith_term][0..len], rhs_mask[0..len]);
                @memcpy(coeff_masks[ith_term][0..len], coeff_mask[0..len]);
            }

            const ret_lhs: [basis_size]@Vector(len, i32) = lhs_masks;
            const ret_rhs: [basis_size]@Vector(len, i32) = rhs_masks;
            const ret_coeff: [basis_size]@Vector(len, T) = coeff_masks;

            return .{ .lhs_masks = ret_lhs, .rhs_masks = ret_rhs, .coeff_masks = ret_coeff };
        }

        //shuffle, mul, and add
        pub fn execute_sparse_vector_op(lhs: anytype, rhs: anytype, res: anytype, comptime len: usize, comptime masks: BinaryScatterMasksReturn(len),  comptime DataType: type) void {
            const lhs_masks = masks.lhs_masks;
            const rhs_masks = masks.rhs_masks;
            const coeff_masks = masks.coeff_masks;

            //const res_type = @typeInfo(@TypeOf(res)).Pointer.child;

            //std.debug.print("\n\t pre execute res.terms: {d:<.1}, lhs.terms: {d:<.1}, rhs.terms: {d:<.1}, ", .{res.terms, lhs.terms, rhs.terms});


            const ones: @Vector(len, T) = @splat(1.0);
            const invalid_mask: @Vector(len, i32) = @splat(-1.0);
            //const zeros: @Vector(len, T) = @splat(0.0);

            const pos: @Vector(len, T) = @splat(1);
            const neg: @Vector(len, T) = @splat(-1);

            //decreasing delta further doesnt detect near zeros.
            const delta: @Vector(len, T) = @splat(0.00001);

            inline for (lhs_masks, rhs_masks, coeff_masks) |_lhs_mask, _rhs_mask, _coeff_mask| {
                const lhs_mask: @Vector(_lhs_mask.len, i32) = _lhs_mask;
                const rhs_mask: @Vector(_lhs_mask.len, i32) = _rhs_mask;
                const coeff_mask: @Vector(_lhs_mask.len, T) = _coeff_mask;
                const invalid = comptime blk: {
                    //@compileLog(comptimePrint("lhs_mask type: {}, rhs_mask type: {}, coeff_mask type: {}", .{@TypeOf(lhs_mask), @TypeOf(rhs_mask), @TypeOf(coeff_mask)}));
                    break :blk _lhs_mask.len == 0 or _rhs_mask.len == 0 or @reduce(.And, lhs_mask == invalid_mask) or @reduce(.And, rhs_mask == invalid_mask) or @reduce(.And, @abs(coeff_mask) < delta);
                };
                if (!invalid) {
                    //@compileLog(comptimePrint("\nlhs mask: {d: <3.0}, rhs_mask: {d: <3.0}, coeff_mask: {d: <3.0}", .{ lhs_mask, rhs_mask, coeff_mask }));

                    //left: {3,7,8,9} gp right: {0,1,2,3}:
                    // for whatever reason, we determinehas {0,3,4,5,7,8,9}
                    // first = @shuffle(D, left, ones, {-1,0,-1,-1,1,2,3}) =     {1, l0, 1, 1, l1, l2, l3}
                    // second = @shuffle(D, right, ones, {0,-1,-1,1,-1,-1,-1}) = {r0, 1, 1, r1, 1, 1, 1}
                    // res.terms += {r0, l0, 1, r1, l1, l2, l3}

                    //std.debug.print("\ncoeff_mask == zeros is {}, reduce: {}", .{@as(@Vector(len, T), @splat(@as(f64, 2.0))) * coeff_mask == coeff_mask, @reduce(.And, @as(@Vector(len, T), @splat(@as(f64, 2.0))) * coeff_mask == coeff_mask)});
                    const first = @shuffle(DataType, lhs.terms, ones, lhs_mask);
                    //@compileLog(comptimePrint("execute_sparse_vector_op(): @TypeOf(lhs).subset: {b}, lhs_mask: {any}, @TypeOf(rhs).subset: {b}, rhs_mask: {any}, @TypeOf(res).subset: {b}, len: {}", .{UnwrapPtrType(@TypeOf(lhs)).@"0".subset, lhs_mask, @TypeOf(rhs).subset, rhs_mask, UnwrapPtrType(@TypeOf(res)).@"0".subset, len}));
                    const second = @shuffle(DataType, rhs.terms, ones, rhs_mask);
                    //std.debug.print("\n\tfirst: {d: <.1} second: {d: <.1}, coeff: {}", .{first, second, coeff_mask});
                    if (comptime @reduce(.And, coeff_mask == pos)) {
                        res.terms += first * second;
                    } else if (comptime @reduce(.And, coeff_mask == neg)) {
                        res.terms -= first * second;
                    } else {
                        res.terms += first * second * coeff_mask;
                    }
                    //std.debug.print("\n\tres.terms: {d:<.1}", .{res.terms});
                }
            }
        }

        pub fn unary_inverse_table_to_scatter_mask(comptime inv: [basis_size]UnaryProdTerm(T), comptime SourceType: type, comptime DrainType: type, comptime len: usize) struct { [len]i32, [len]T } {
            const res_subset: usize = DrainType.subset;
            const opr_subset: usize = SourceType.subset;
            comptime var coeff_mask: [len]T = .{0} ** len;
            comptime var opr_mask: [len]i32 = .{-1} ** len;
            const one: usize = 0b1;
            for (0..basis_size) |res_blade| {
                if (res_subset & (one << res_blade) == 0) {
                    continue;
                }
                const term: UnaryProdTerm(T) = inv[res_blade];
                const coeff = term.mult;

                if(opr_subset & (1 << term.source) == 0) {
                    @compileError(comptimePrint("term {any} with blade {b} is not in source's subset {b}", .{term, term.source, opr_subset}));
                }

                const opr_index = comptime count_bits(opr_subset & ((one << term.source) -% 1));
                const res_idx = DrainType.canonical_blade_to_idx[res_blade].?;//comptime count_bits(res_subset & ((one << res_blade) -% 1));

                opr_mask[res_idx] = opr_index;
                coeff_mask[res_idx] = coeff;
                //@compileLog(comptimePrint("unary_inverse_table_to_scatter_mask: opr_subset: {b}, term.source: {b}, one << term.source: {b}, one << term.source -% 1: {b}, & opr_subset: {b}, count_bits(<-): {}", .{opr_subset, term.source, one << term.source, (one << term.source) -% 1, opr_subset & ((one << term.source) -% 1), opr_index}));
            }

            const cf_mask: @Vector(len, T) = coeff_mask;
            const operand_mask: @Vector(len, i32) = opr_mask;
            return comptime .{ operand_mask, cf_mask };
        }

        pub fn unary_execute_sparse_vector_op(opr: anytype, res: anytype, comptime masks: anytype, comptime len: usize, comptime DataType: type) void {
            const opr_mask = masks[0];
            const coeff_mask = masks[1];
            //const ones: @Vector(len, i32) = @splat(1);
            //const invalid_mask: @Vector(len, i32) = @splat(-1);
            //const zeros: @Vector(len, T) = @splat(0);
            const pos: @Vector(len, T) = @splat(1);
            //const neg: @Vector(len, T) = @splat(-1);

            //@compileLog(comptimePrint("unary_execute_sparse_vector_op: opr_mask: {any}, coeff_mask: {any}", .{opr_mask, coeff_mask}));
            const gather = @shuffle(DataType, opr.terms, pos, opr_mask);
            res.terms = gather * coeff_mask;
            //debug.print("\nunary_execute_sparse_vector_op: {} -> {}", .{opr, res});
        }

        //TODO: make this work for variably known / unknown terms
        pub fn nonzero_subset(comptime terms: [basis_size][basis_size]ProdTerm(T), comptime indices: [basis_size]usize) usize {
            _ = terms;
            var ret_subset: usize = 0;
            //const delta: T = 0.00001;
            for(0..basis_size) |blade| {
                const len = indices[blade];
                if(len == 0) {
                    continue;
                }
                //var sum: T = 0;
                //for(0..len) |i| {
                //    sum += terms[blade][i].mult;
                //}
                //if(@abs(sum) < delta) {
                //    continue;
                //}
                ret_subset |= 1 << blade;

            }
            return ret_subset;
        }

        pub fn unary_nonzero_subset(comptime terms: [basis_size]UnaryProdTerm(T)) usize {

            var ret_subset: usize = 0;
            const delta: T = 0.00001;
            for(0..basis_size) |blade| {
                const result = terms[blade];
                if(@abs(result.mult) < delta) {
                    continue;
                }
                ret_subset |= 1 << blade;

            }
            return ret_subset;
        }

        ///finds the first (if any) alias of the given algebra (_alg) that has every blade in subset in SOME order
        pub fn find_alias_matching_runtime_bitfield(comptime _alg: anytype, comptime input_subset: usize) ?GAlgebraCompiledAlias {
            const algargs = _alg.GAlgebraInfoArgs;
            const alg = GAlgebraInfo(algargs.__dual, algargs.__basis, algargs.__ordering, algargs.__aliases);
            for(0..alg.aliases.len) |i| {
                const alias = alg.aliases[i];
                var runtime_contains_all: bool = true;
                //@compileLog(comptimePrint("find_alias_matching_runtime_bitfield {}: alias: {any}", .{i, alias}));
                for(0..alias.values.len) |j| {
                    const alias_blade: usize = alias.values[j];
                    if(input_subset & (1 << alias_blade) == 0) {
                        runtime_contains_all = false;
                        //@compileLog(comptimePrint("find_alias_matching_runtime_bitfield {}: alias: {any} breaking because input subset {b} doesnt contain 1<<{b}={b}", .{i, alias, input_subset, alias_blade, 1 << alias_blade}));
                        break;
                    }
                } 
                if(runtime_contains_all == false) {
                    continue;
                }
                return alias;
            }
             
            return null;
        }

        ///given an alias, creates the corresponding _values argument usable for creating an MVecSubset matching that alias
        pub fn alias_to_values_arg(comptime alias: GAlgebraCompiledAlias) []const usize {
            return alias.values;
        }

        pub fn subset_to_values_arg(comptime sbset: usize) []const usize {

            var ret: []const usize = &.{};
            for(0..@bitSizeOf(usize)) |i| { 
                if(sbset & (1 << i) != 0) {
                    ret = ret ++ .{i};
                }   
            }
            const _ret = ret;
            return _ret;
        }

        pub fn MVecSubsetFromSubset(comptime sbset: usize) type {
            const arg = subset_to_values_arg(sbset);
            return MVecSubset(_alg_info, T, .{.find_alias = false, .values = arg});
        }

        /// tries to match the given subset to an existing alias in the algebra and creates an MVec for that if found, and creates one without an alias if not
        pub fn MatchAlias(comptime ret_subset: usize) type {
            const matching_alias = comptime find_alias_matching_runtime_bitfield(Algebra, ret_subset);
            if(matching_alias) |alias| {
                //@compileLog(comptimePrint("MatchAlias ret_subset: {b}, matching_alias picked: {any}", .{ret_subset, alias}));
                //const arg = alias_to_values_arg(alias);
                return MVecSubset(_alg_info, T, .{.compiled_alias = alias});
            } else {
                const arg = subset_to_values_arg(ret_subset);
                //@compileLog(comptimePrint("MatchAlias ret_subset: {b}, matching_alias NOT picked, arg: {any}", .{ret_subset, arg}));
                //@compileLog(comptimePrint("BinaryRes subset_to_values_arg({b}) = {any}", .{ret_subset, arg}));
                return MVecSubset(_alg_info, T, .{.find_alias = false, .values = arg});
            }
        }

        /// returns the output type of the given operation given the operand and result types.
        /// if the result type is null or void then the output type is determined by the possible blades created by the operation
        /// given the operands. if the result type is a pointer or a type then the output is forced to be that type
        pub fn BinaryRes(comptime opr: GAlgebraOp(T), comptime left_outer_type: type, comptime right_outer_type: type, comptime res_outer_type: type) type {
            const res = res_outer_type;//ResType(res_outer_type);
            //users of users of this function can pass in a pointer, type, or null for the result arg
            //pointer -> use the .child type as the logical type of the result, and fill in the results in that array
            //type -> create a new instance of that type as the result
            //null -> infer the type of the result based on opr

            //for the left and right types they can pass in a pointer or struct instance
            //pointer -> use the .child type as the logical type of that arg
            //struct -> use that structs type as the logical type, and proceed as normal

            //we should be guaranteed that res_outer_type is either .Pointer, .Struct, or .Null (@TypeOf(type) is lossy)

            const res_info = @typeInfo(res);
            //@compileLog(comptimePrint("BinaryRes: res: {s} info: {}", .{@typeName(res), res_info}));
            if(res_info == .Type) {
                @compileError(comptimePrint("unguarded type passed into BinaryRes()! @TypeOf(type) is lossy so if a user passes in a type you have to pass in that type directly as an argument instead of passing in @TypeOf(that type)", .{}));
            }
            if(res_info == .Pointer or res_info == .Struct) {
                return res;
            }

            const left_ret = UnwrapPtrType(left_outer_type);
            const right_ret = UnwrapPtrType(right_outer_type);
            
            const left_type: type = left_ret.@"0";
            const right_type: type = right_ret.@"0";
 
            const table = comptime Algebra.binary_cayley_table(T, opr.BladeOp(Algebra.signature), left_type.subset, right_type.subset);
            //comptime {
            //    var str: [1000]u8 = .{0} ** 1000;
            //    var c = 0;
            //    for(0..basis_size) |b| {
            //        const slc: []u8 = try std.fmt.bufPrint(str[c..str.len], "blade {b} =", .{b});
            //        c += slc.len;
            //        for(0..table.results_to_num_terms[b]) |n| {
            //            const fmt: []const u8 = if(n == 0) " {} from {b} * {b}" else " + {} from {b} * {b}";
            //            const term: ProdTerm(T) = table.results_to_terms[b][n];
            //
            //            const slc2: []u8 = try std.fmt.bufPrint(str[c..str.len], fmt, .{term.mult, term.left_blade, term.right_blade});
            //            c += slc2.len;
            //        }
            //        for(0..3) |_| {
            //            str[c] = ' ';
            //            c += 1;
            //        }
            //    }
            //    @compileLog(str[0..c]);
            //}
            const ret_subset: usize = comptime nonzero_subset(table.results_to_terms, table.results_to_num_terms);
            return MatchAlias(ret_subset);
            //return MVecSubset(ret_subset);
        }

        pub fn UnaryRes(comptime opr: GAlgebraOp(T), source: type, drain: type) type {
            const res_info = @typeInfo(drain);

            if(res_info == .Type) {
                @compileError(comptimePrint("unguarded type passed into UnaryRes()! @TypeOf(type) is lossy so you have to pass in the type directly if a user passes in a type", .{}));
            }
            if(res_info == .Pointer or res_info == .Struct) {
                return drain;
            }

            if(res_info != .Null and res_info != .Void) {
                @compileError(comptimePrint("invalid result argument passed to UnaryRes()!"));
            }

            const source_ret = UnwrapPtrType(source);

            const source_t: type = source_ret.@"0";

            const table = comptime Algebra.unary_cayley_table(T, opr.BladeOp(Algebra.signature), source_t.subset);
            const ret_subset: usize = unary_nonzero_subset(table.results_to_operands);
            
            return MatchAlias(ret_subset);

            //return MVecSubset(ret_subset);

        }

        ///if passed_in_res is a pointer, returns that. otherwise creates a new instance of the true_res type passed in
        pub fn get_result(passed_in_res: anytype, comptime true_res: type) true_res {
            if(@typeInfo(@TypeOf(passed_in_res)) == .Pointer) {
                return passed_in_res;
            }
            var res: true_res = undefined;
            _=res.to_zeros();
            return res;
        }

        /// multiplicative binary operation that produces an inferred multivector type. doesnt work for addition
        pub fn binary_op(comptime opr: GAlgebraOp(T), left: anytype, right: anytype, res: anytype) BinaryRes(opr, @TypeOf(left), @TypeOf(right), if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res)) {
            const Result: type = BinaryRes(opr, @TypeOf(left), @TypeOf(right), if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res));
            var result = get_result(res, Result);

            const Left = @TypeOf(left);
            const Right = @TypeOf(right);

            const left_info = @typeInfo(Left);
            const right_info = @typeInfo(Right);   
            const res_info = @typeInfo(Result);             

            const LeftInner = if(left_info == .Pointer) left_info.Pointer.child else Left;
            const RightInner = if(right_info == .Pointer) right_info.Pointer.child else Right;
            const ResInner = if(res_info == .Pointer) res_info.Pointer.child else Result;
            

            //@compileLog(comptimePrint("gp left: {s}, right: {s},,, left ordering: {any},,, right ordering: {any},,, Result ordering {any}", .{@typeName(LeftInner), @typeName(RightInner), LeftInner.values_arg, RightInner.values_arg, ResInner.values_arg}));

            const left_set: usize = LeftInner.subset;
            const right_set: usize = RightInner.subset;
            const res_set: usize = ResInner.subset;

            const table = comptime Algebra.binary_cayley_table(T, opr.BladeOp(Algebra.signature), left_set, right_set);
            //tell the table inverter the res size, then we can expect it to return an array with comptime known size
            //@compileLog(comptimePrint("table: {}", .{table}));
            const res_size = comptime count_bits(res_set);
            const operation = comptime binary_create_scatter_masks_from_cayley_table_inv(table.results_to_terms, table.results_to_num_terms, LeftInner, RightInner, ResInner, res_size);


            //@compileLog(comptimePrint("\nleft_set: {b}, right_set: {b}, res_set: {b}, table: {}, res_size: {}, operation: {}", .{left_set, right_set, res_set, table, res_size, operation}));

            const left_arg = blk: {
                if(left_info == .Pointer) {
                    break :blk left;
                }
                break :blk &left;
            };

            const res_to_use = blk: {
                if(@typeInfo(Result) == .Pointer) {
                    break :blk result;
                }
                break :blk &result;
            };

            execute_sparse_vector_op(left_arg, right, res_to_use, comptime res_size, comptime operation, T);
            return result;
        }

        /// geometric product
        pub fn gp(left: anytype, right: anytype, res: anytype) BinaryRes(.GeometricProduct, @TypeOf(left), @TypeOf(right), if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res)) {
            return binary_op(GAlgebraOp(T).GeometricProduct, left, right, res);
            //const Result: type = BinaryRes(.GeometricProduct, @TypeOf(left), @TypeOf(right), if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res));
            //var result = get_result(res, Result);
            //
            //const Left = @TypeOf(left);
            //const Right = @TypeOf(right);
            //
            //const left_info = @typeInfo(Left);
            //const right_info = @typeInfo(Right);   
            //const res_info = @typeInfo(Result);             
            //
            //const LeftInner = if(left_info == .Pointer) left_info.Pointer.child else Left;
            //const RightInner = if(right_info == .Pointer) right_info.Pointer.child else Right;
            //const ResInner = if(res_info == .Pointer) res_info.Pointer.child else Result;
            //
            //
            ////@compileLog(comptimePrint("gp left: {s}, right: {s},,, left ordering: {any},,, right ordering: {any},,, Result ordering {any}", .{@typeName(LeftInner), @typeName(RightInner), LeftInner.values_arg, RightInner.values_arg, ResInner.values_arg}));
            //
            //const left_set: usize = LeftInner.subset;
            //const right_set: usize = RightInner.subset;
            //const res_set: usize = ResInner.subset;
            //
            //const table = comptime Algebra.binary_cayley_table(T, GAlgebraOp(T).GeometricProduct.BladeOp(Algebra.signature), left_set, right_set);
            ////tell the table inverter the res size, then we can expect it to return an array with comptime known size
            //const res_size = comptime count_bits(res_set);
            //const operation = comptime binary_create_scatter_masks_from_cayley_table_inv(table.results_to_terms, table.results_to_num_terms, LeftInner, RightInner, ResInner, res_size);
            //
            //
            ////@compileLog(comptimePrint("\nleft_set: {b}, right_set: {b}, res_set: {b}, table: {}, res_size: {}, operation: {}", .{left_set, right_set, res_set, table, res_size, operation}));
            //
            //const left_arg = blk: {
            //    if(left_info == .Pointer) {
            //        break :blk left;
            //    }
            //    break :blk &left;
            //};
            //
            //const res_to_use = blk: {
            //    if(@typeInfo(Result) == .Pointer) {
            //        break :blk result;
            //    }
            //    break :blk &result;
            //};
            //
            //execute_sparse_vector_op(left_arg, right, res_to_use, comptime res_size, comptime operation, T);
            //return result;
        }

        /// outer product
        pub fn op(left: anytype, right: anytype, res: anytype) BinaryRes(.OuterProduct, @TypeOf(left), @TypeOf(right), if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res)) {
            return binary_op(GAlgebraOp(T).OuterProduct, left, right, res);
        }

        /// inner product
        pub fn ip(left: anytype, right: anytype, res: anytype) BinaryRes(.InnerProduct, @TypeOf(left), @TypeOf(right), if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res)) {
            return binary_op(GAlgebraOp(T).InnerProduct, left, right, res);
        }

        /// regressive product
        pub fn rp(left: anytype, right: anytype, res: anytype) BinaryRes(.RegressiveProduct, @TypeOf(left), @TypeOf(right), if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res)) {
            return binary_op(GAlgebraOp(T).RegressiveProduct, left, right, res);
        }

        pub fn Add(comptime lhs: type, comptime rhs: type, comptime res: type) type {
            const res_info = @typeInfo(res);
            if(res_info == .Type or res_info == .Pointer or res_info == .Struct) {
                return res;
            }

            const lr = comptime UnwrapPtrType(lhs);
            const rr = comptime UnwrapPtrType(rhs);

            const LeftInner: type = lr.@"0";
            const RightInner: type = rr.@"0";

            const left_set: usize = LeftInner.runtime_subset;
            const right_set: usize = RightInner.runtime_subset;

            var res_set: usize = 0;

            for(0..basis_size) |blade_i| {
                for(0..basis_size) |blade_j| {
                    if(left_set & (1 << blade_i) != 0) {
                        res_set |= (1 << blade_i);
                    }
                    if(right_set & (1 << blade_j) != 0) {
                        res_set |= (1 << blade_j);
                    }
                }
            }
            return MVecSubsetFromSubset(res_set);
        }


        pub fn add(lhs: anytype, rhs: anytype, res: anytype) Add(@TypeOf(lhs), @TypeOf(rhs), if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res)) { 
            const Result = Add(@TypeOf(lhs), @TypeOf(rhs), if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res));
            var result = get_result(res, Result);

            const left_type: type = UnwrapPtrType(@TypeOf(lhs)).@"0";
            const right_type: type = UnwrapPtrType(@TypeOf(rhs)).@"0";
            //@compileLog(comptimePrint("right type name: {s}, info {} ------------------------------------- rhs name {s} info {}", .{@typeName(right_type), @typeInfo(right_type), @typeName(@TypeOf(rhs)), @typeInfo(@TypeOf(rhs))}));
            const ResInner: type = UnwrapPtrType(Result).@"0";

            const left_set: usize = left_type.subset;
            const right_set: usize = comptime blk: {
                if(@typeInfo(right_type) == .Pointer) {
                    @compileError(comptimePrint("right type info {} rhs type info {}",.{@typeInfo(right_type), @typeInfo(@TypeOf(rhs))}));
                }
                break :blk right_type.subset;
            };
            const res_set: usize = comptime blk: {
                break :blk ResInner.subset;
            };
            const len = comptime count_bits(res_set);

            comptime var _lhs_mask: @Vector(len, i32) = @splat(-1);
            comptime var _rhs_mask: @Vector(len, i32) = @splat(-1);
            comptime {
                for (0..basis_size) |result_blade| {
                    if (res_set & (1 << result_blade) == 0) {
                        continue;
                    }
                    //left_set = 0b01011101, if res_blade = 101 (5)
                    //then if (0b01011101 & (0b1 << 101 = 0b100000) != 0 (true)),
                    //then the index in left for blade 101 is:
                    //count_bits(0b01011101 & (0b100000 -% 1 = 0b011111) = 0b00011101) = 4
                    //which is correct
                    
                    const _left_idx = left_type.canonical_blade_to_idx[result_blade];
                    const _right_idx = right_type.canonical_blade_to_idx[result_blade];
                    if(_left_idx == null) {
                        @compileError(comptimePrint("add(): lhs blade {b} not found in subset {b} because its index is null! canonical_blade_to_idx: {any}", .{result_blade, left_set, left_type.canonical_blade_to_idx}));
                    }
                    if(_right_idx == null) {
                        @compileError(comptimePrint("add(): rhs blade {b} not found in subset {b} because its index is null! canonical_blade_to_idx: {any}", .{result_blade, right_set, right_type.canonical_blade_to_idx}));
                    }
                    const left_idx: usize = _left_idx.?;
                    const right_idx: usize = _right_idx.?;

                    //must exist, and if an operand idx exists this is >= that
                    //assuming res is a superset of lhs.subset_field | rhs.subset_field
                    const _res_idx = ResInner.canonical_blade_to_idx[result_blade];
                    if(_res_idx == null) {
                        @compileError(comptimePrint("add(): res blade {b} not found in subset {b} because its index is null! canonical_blade_to_idx: {any}", .{result_blade, res_set, ResInner.canonical_blade_to_idx}));
                    }
                    const res_idx = _res_idx.?;

                    if (left_set & (1 << result_blade) != 0) {
                        _lhs_mask[res_idx] = left_idx;
                    }
                    if (right_set & (1 << result_blade) != 0) {
                        _rhs_mask[res_idx] = right_idx;
                    }
                }
            }
            const lhs_mask = _lhs_mask;
            const rhs_mask = _rhs_mask;
            const zero: @Vector(len, T) = @splat(0);

            const left = @shuffle(T, lhs.terms, zero, lhs_mask);
            const right = @shuffle(T, rhs.terms, zero, rhs_mask);

            result.terms = left + right;

            return result;
        }

        pub fn sub(lhs: anytype, rhs: anytype, res: anytype) Add(@TypeOf(lhs), @TypeOf(rhs), if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res)) { 
            const Result = Add(@TypeOf(lhs), @TypeOf(rhs), if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res));
            var result = get_result(res, Result);

            const left_type: type = UnwrapPtrType(@TypeOf(lhs)).@"0";
            const right_type: type = UnwrapPtrType(@TypeOf(rhs)).@"0";
            const ResInner: type = UnwrapPtrType(Result).@"0";

            const left_set: usize = comptime blk: {
                break :blk left_type.runtime_subset;
            };
            const right_set: usize = comptime blk: {
                break :blk right_type.runtime_subset;
            };
            const res_set: usize = comptime blk: {
                break :blk ResInner.runtime_subset;
            };
            const len = comptime count_bits(res_set);

            comptime var _lhs_mask: @Vector(len, i32) = @splat(-1);
            comptime var _rhs_mask: @Vector(len, i32) = @splat(-1);
            comptime {
                for (0..basis_size) |result_blade| {
                    if (res_set & (1 << result_blade) == 0) {
                        continue;
                    }
                    //left_set = 0b01011101, if res_blade = 101 (5)
                    //then if (0b01011101 & (0b1 << 101 = 0b100000) != 0 (true)),
                    //then the index in left for blade 101 is:
                    //count_bits(0b01011101 & (0b100000 -% 1 = 0b011111) = 0b00011101) = 4
                    //which is correct
                    
                    const _left_idx = left_type.canonical_blade_to_idx[result_blade];
                    const _right_idx = right_type.canonical_blade_to_idx[result_blade];
                    if(_left_idx == null) {
                        @compileError(comptimePrint("add(): lhs blade {b} not found in subset {b} because its index is null! canonical_blade_to_idx: {any}", .{result_blade, left_set, left_type.canonical_blade_to_idx}));
                    }
                    if(_right_idx == null) {
                        @compileError(comptimePrint("add(): rhs blade {b} not found in subset {b} because its index is null! canonical_blade_to_idx: {any}", .{result_blade, right_set, right_type.canonical_blade_to_idx}));
                    }
                    const left_idx: usize = _left_idx.?;
                    const right_idx: usize = _right_idx.?;

                    //must exist, and if an operand idx exists this is >= that
                    //assuming res is a superset of lhs.subset_field | rhs.subset_field
                    const _res_idx = ResInner.canonical_blade_to_idx[result_blade];
                    if(_res_idx == null) {
                        @compileError(comptimePrint("add(): res blade {b} not found in subset {b} because its index is null! canonical_blade_to_idx: {any}", .{result_blade, res_set, ResInner.canonical_blade_to_idx}));
                    }
                    const res_idx = _res_idx.?;

                    if (left_set & (1 << result_blade) != 0) {
                        _lhs_mask[res_idx] = left_idx;
                    }
                    if (right_set & (1 << result_blade) != 0) {
                        _rhs_mask[res_idx] = right_idx;
                    }
                }
            }
            const lhs_mask = _lhs_mask;
            const rhs_mask = _rhs_mask;
            const zero: @Vector(len, T) = @splat(0);

            const left = @shuffle(T, lhs.terms, zero, lhs_mask);
            const right = @shuffle(T, rhs.terms, zero, rhs_mask);


            result.terms = left - right;

            return result;
        }

        pub fn mul_scalar(self: anytype, scalar: T, res: anytype) BinaryRes(.GeometricProduct, @TypeOf(self), MVecSubset(Algebra, T, .{.find_alias = false, .subset = 1}), if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res)) {
            //@compileLog(comptimePrint("mul_scalar self: {s}, res name: {s}, res type info: {?}", .{@typeName(@TypeOf(self)), @typeName(@TypeOf(res)), @typeInfo(@TypeOf(res))}));
            const Result = BinaryRes(.GeometricProduct, @TypeOf(self), MVecSubset(alg_info, T, .{.find_alias = false, .subset = 1}), if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res));
            //var result = get_result(res, Result);
            var result = blk: {

                const passed_res_info = blk2: {
                    const ti = @typeInfo(@TypeOf(res));
                    if(ti == .Type) {
                        break :blk2 @typeInfo(res);
                    }
                    break :blk2 ti;
                };
                if(passed_res_info == .Pointer) {
                    //const ret: true_res = @as(passed_in_res, true_res);
                    //@compileLog(comptimePrint("passed_res_info child: {?}, true res: {s}", .{passed_res_info.Pointer.child, @typeName(true_res)}));
                    break :blk res;
                }
                var ret: Result = undefined;
                _=ret.to_zeros();
                break :blk ret;
            };
            const result_ref = if(@typeInfo(@TypeOf(result)) == .Pointer) result else &result;
            _mul_scalar(self, scalar, result_ref);
            return result;
        }

        pub fn _mul_scalar(self: anytype, scalar: T, result: anytype) void {
            const Result: type = @TypeOf(result);

            if(Result == Self or Result == *Self) {
                result.terms = self.terms * @as(@TypeOf(self.terms), @splat(scalar));
            } else {
                const ResInner: type = UnwrapPtrType(Result).@"0";
                inline for(0..basis_size) |blade_t| {
                    //const blade_t: ud = @truncate(blade);
                    if(ResInner.subset_field & (1 << blade_t) != 0 and Self.subset_field & (1 << blade_t) != 0) {
                        const us_idx = Self.canonical_blade_to_idx[blade_t].?;
                        const res_idx = ResInner.canonical_blade_to_idx[blade_t].?;
                        result.terms[res_idx] = self.terms[us_idx] * scalar;
                    }
                }
            }
        }


        pub fn MeetRes(left: type, right: type, res: type) type {
            if(Self.Algebra.dual == true) {
                return BinaryRes(.OuterProduct, left, right, res);
            } else {
                return BinaryRes(.RegressiveProduct, left, right, res);
            }
        }

        pub fn meet(left: anytype, right: anytype, res: anytype) MeetRes(@TypeOf(left), @TypeOf(right), if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res)) {
            if(comptime Self.Algebra.dual == true) {
                return Self.op(left,right,res);
            } else {
                return Self.rp(left,right,res);
            }
        }

        pub fn JoinRes(left: type, right: type, res: type) type {
            if(Self.Algebra.dual == true) {
                return BinaryRes(.RegressiveProduct, left, right, res);
            } else {
                return BinaryRes(.OuterProduct, left, right, res);
            }
        }

        pub fn join(left: anytype, right: anytype, res: anytype) JoinRes(@TypeOf(left), @TypeOf(right), if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res)) {
            if(comptime Self.Algebra.dual == true) {
                return Self.rp(left,right,res);
            } else {
                return Self.op(left,right,res);
            }
        }

        ///project if grade is comptime known
        pub fn CompProject(comptime grade: ?usize) type {
            if(grade == null) {
                @compileError(comptimePrint("Project requires either a non null result type or a non null grade!"));
            }
            const tgrade = grade.?;
            var res_subset: usize = 0;

            for(0..basis_size) |blade| {
                if(count_bits(blade) == tgrade) {
                    res_subset |= (1 << blade);
                }
            }
            return MatchAlias(res_subset);
        }

        pub fn Project(grade: ?usize, comptime res: type) type {
            const res_info = @typeInfo(res);
            //@compileLog(comptimePrint("BinaryRes: res: {s} info: {}", .{@typeName(res), res_info}));
            if(res_info == .Type) {
                @compileError(comptimePrint("unguarded type passed into Project()! @TypeOf(type) is lossy so if a user passes in a type you have to pass in that type directly as an argument instead of passing in @TypeOf(that type)", .{}));
            }
            if(res_info == .Pointer or res_info == .Struct) {
                return res;
            }

            //res is null, so infer the correct type based on what grade is
            if(is_comptime(grade)) {
                return CompProject(comptime grade);
            } else {
                @compileError(comptimePrint("Project() requires a comptime known grade arg if res is unknown!", .{}));
            }
        }

        /// project self onto res / the inferred type with all blades matching grade
        pub fn project(self: *Self, grade: ?usize, res: anytype) Project(grade, if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res)) {
            const Result: type = Project(grade, if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res));
            const result = get_result(res, Result);

            inline for (0..Self.Algebra.basis_len) |i| {
                const _res_idx = Result.canonical_blade_to_idx[i];
                if(_res_idx == null) {
                    continue;
                }
                const res_idx = _res_idx.?;

                const _self_idx = Self.canonical_blade_to_idx[i];
                if(_self_idx == null) {
                    continue;
                }
                const self_idx = _self_idx.?;

                if(grade) |g| {
                    if(count_bits(i) == g) {
                        result.terms[res_idx] = self.terms[self_idx];
                    }
                } else {
                    result.terms[res_idx] = self.terms[self_idx];
                }
            }
            return result;
        }


        ///subset invariant
        /// 0 1 10, 11, 100, 101, 110, 111
        /// +,+,-,-,+,+,-,-,+,+,...
        pub fn revert(self: *Self) *Self {
            inline for (0..basis_size) |blade| {
                if (Self.subset & (1 << blade) != 0) {
                    const i = Self.canonical_blade_to_idx[blade].?;
                    const mult: T = comptime blk: {
                        var x: usize = count_bits(blade); //this is correct
                        x >>= 1;
                        if (x & 1 == 0) {
                            break :blk 1.0;
                        }
                        break :blk -1.0;
                    };
                    self.terms[i] = mult * self.terms[i];
                }
            }
            return self;
        }

        pub fn reverse(source: *Self, drain: anytype) UnaryRes(.Reverse, @TypeOf(source), @TypeOf(drain)) {
            const Result = UnaryRes(.Reverse, @TypeOf(source), @TypeOf(drain));
            var result = get_result(drain, Result);

            const SourceInner = UnwrapPtrType(@TypeOf(source)).@"0";
            const DrainInner = UnwrapPtrType(@TypeOf(result)).@"0";

            const source_set: usize = SourceInner.subset;
            const drain_set: usize = DrainInner.subset;

            const table = comptime alg_info.unary_cayley_table(T, GAlgebraOp(T).Reverse.BladeOp(Algebra.signature), source_set);
            const operation = comptime Self.unary_inverse_table_to_scatter_mask(table.results_to_operands, SourceInner, DrainInner, count_bits(drain_set));

            Self.unary_execute_sparse_vector_op(source, &result, comptime operation, count_bits(drain_set), Self.Scalar);
            debug.print("\nreverse: {} -> {}", .{source, result});
            return result;
        }

        pub fn invert(self: *Self) *Self {
            _=self.revert();
            //const mag = self.magnitude_squared();
            //var x = self.copy(Self);
            //_=x.mul_scalar(1.0 / mag, .ptr, .{.ptr = self});
            return self.normalize(.ptr, .{.ptr=self});
        }


        ///changes subsets
        /// if we have a blade at bit index i, our dual must have a blade at bit index basis_len - 1 - i
        /// in other words our dual has our subset but mirror across the middle at minimum
        pub fn dual(source: *Self, drain: anytype) UnaryRes(.Dual, @TypeOf(source), if(@typeInfo(@TypeOf(drain)) == .Type) drain else @TypeOf(drain)) {
            const Result = UnaryRes(.Dual, @TypeOf(source), if(@typeInfo(@TypeOf(drain)) == .Type) drain else @TypeOf(drain));
            var result = get_result(drain, Result);

            const SourceInner = UnwrapPtrType(@TypeOf(source)).@"0";
            const DrainInner = UnwrapPtrType(Result).@"0";

            const source_set: usize = SourceInner.subset_field;
            const drain_set: usize = DrainInner.subset_field;

            const table = comptime alg_info.unary_cayley_table(source_set, GAlgebraOp(T).Dual.BladeOp(), null);
            const operation = comptime Self.unary_inverse_table_to_scatter_mask(table.snd, SourceInner, DrainInner, count_bits(drain_set));

            const res_to_use = blk: {
                if(@typeInfo(Result) == .Pointer) {
                    break :blk result;
                }
                break :blk &result;
            };

            Self.unary_execute_sparse_vector_op(source, res_to_use, comptime operation, count_bits(drain_set), T);
            return result;
        }

        pub fn undual(source: anytype, drain: anytype) UnaryRes(.Undual, @TypeOf(source), if(@typeInfo(@TypeOf(drain)) == .Type) drain else @TypeOf(drain)) {
            const Result = UnaryRes(.Undual, @TypeOf(source), if(@typeInfo(@TypeOf(drain)) == .Type) drain else @TypeOf(drain));
            var result = get_result(drain, Result);

            const SourceInner: type = UnwrapPtrType(@TypeOf(source)).@"0";
            const DrainInner: type = UnwrapPtrType(Result).@"0";

            const source_set: usize = SourceInner.subset_field;
            const drain_set: usize = DrainInner.subset_field;

            const table = comptime alg_info.unary_cayley_table(source_set, GAlgebraOp(T).Undual.BladeOp(), null);
            const operation = comptime Self.unary_inverse_table_to_scatter_mask(table.snd, SourceInner, DrainInner, count_bits(drain_set));

            const res_to_use = blk: {
                if(@typeInfo(Result) == .Pointer) {
                    break :blk result;
                }
                break :blk &result;
            };
            //comperror("drain set: {b}, Result: {s}, DrainInner: {s}", .{drain_set, @typeName(Result), @typeName(DrainInner)});
            Self.unary_execute_sparse_vector_op(source, res_to_use, comptime operation, comptime count_bits(drain_set), T);
            return result;
        }


        pub fn involution(source: *Self, drain: anytype) UnaryRes(.Involution, @TypeOf(source), if(@typeInfo(@TypeOf(drain)) == .Type) drain else @TypeOf(drain)) {
            const Result = UnaryRes(.Involution, @TypeOf(source), if(@typeInfo(@TypeOf(drain)) == .Type) drain else @TypeOf(drain));
            const result = get_result(drain, Result);

            const SourceInner = UnwrapPtrType(@TypeOf(source)).@"0";
            const DrainInner = UnwrapPtrType(@TypeOf(drain)).@"0";

            const source_set: usize = SourceInner.subset_field;
            const drain_set: usize = DrainInner.subset_field;

            const table = comptime alg_info.unary_cayley_table(source_set, GAlgebraOp(T).Dual.BladeOp(Algebra.signature), null);
            const operation = comptime Self.unary_inverse_table_to_scatter_mask(table.snd, SourceInner, DrainInner, count_bits(drain_set));

            Self.unary_execute_sparse_vector_op(source, result, comptime operation, count_bits(drain_set), T);
            return result;
        }

        pub fn magnitude_squared(self: *Self) T {
            @setEvalBranchQuota(100000);
            //var scalar: T = 0;

            var gp_res = MVecSubsetFromSubset((1 << Algebra.basis_len) - 1){.terms = undefined};
            _=(&gp_res).to_zeros();

            //_= self.gp(self, gp_res.zeros());
            //var cop = self.copy(Self);
            var _cop = self.copy(Self);
            const cop = _cop.reverse(null);
            //return @reduce(.Add, (&cop).gp(self, gp_res).terms);
            const ret_mvec = self.gp(cop, &gp_res);
            const ret = ret_mvec.terms[0];
            debug.print("\nmagnitude_squared({}): gp_res = {}, copy = {}, ret mvec = {}, ret = {}", .{self, gp_res, cop, ret_mvec, ret});
            return ret;
        }

        pub fn magnitude(self: *Self) T {
            return @sqrt(@abs(self.magnitude_squared()));
        }

        pub fn normal(self: Self) Self {
            var slf = self; 
            const mag = slf.magnitude();
            const ret = self.mul_scalar(1.0 / mag, Self);
            //debug.print("\nnormal: {} -> {}, mag {}", .{self, ret, mag});
            return ret;
        }

        pub fn normalize(self: *Self) *Self {
            return self.mul_scalar(1.0 / self.magnitude(), self);
        }

        /// divides everything in this mvec by the value of the specified blade if it exists
        pub fn normalize_by_blade(self: Self, blade: usize) Self {
            if(Self.canonical_blade_to_idx[blade] == null) {
                return self;
            }
            const mag = self.terms[Self.canonical_blade_to_idx[blade].?];
            if(mag == 0) {
                return self;
            }
            const ret = self.mul_scalar(1.0 / mag, Self);
            debug.print("\nnormalize_by_blade: self {}, blade {b}, ret {}", .{self, blade, ret});
            return ret;
        }

        //in non degenerate algebras only squares of blades are scalars,
        //in degenerate algebras the only non squares that equal 0 are degenerate blad
        //
        //for bivector B:
        //if B ** 2 = -1 then e^tB = cos(t) + Bsin(t) = rotor
        //if B ** 2 = 0 then e^tB = 1 + tB = motor
        //if B ** 2 = 1 then e^tB = cosh(t) + Bsinh(t) = boost
        //e^x = 1 + x^1/1! + x^2/2! + x^3/3! + ...
        //pub fn exp(self: *Self, coeff: T) NVec {
        //    var res = 
        //}

        pub fn copy(self: *Self, comptime res_type: type) res_type {
            var cop = res_type{.terms = @splat(0)};
            _=cop.to_zeros();
            var copcop = res_type{.terms = @splat(0)};
            _=res_type.add(&copcop, self, &cop);
            //@memcpy(&cop.terms, &self.terms);
            return cop;
        }

        pub fn sandwich(bread: anytype, meat: anytype, res: anytype) BinaryRes(.SandwichProduct, @TypeOf(bread), @TypeOf(meat), if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res)) {
            const Result = BinaryRes(.SandwichProduct, @TypeOf(bread), @TypeOf(meat), if(@typeInfo(@TypeOf(res)) == .Type) res else @TypeOf(res));
            var result = get_result(res, Result);

            const LeftInner = UnwrapPtrType(@TypeOf(bread)).@"0";
            const RightInner = UnwrapPtrType(@TypeOf(meat)).@"0";
            const ResInner = UnwrapPtrType(Result).@"0";
            
            const left_set: ux = LeftInner.subset;
            const right_set: ux = RightInner.subset;
            const res_set: ux = ResInner.subset;
            const table = comptime alg_info.binary_cayley_table(T, GAlgebraOp(T).InnerProduct, left_set, right_set);
            //tell the table inverter the res size, then we can expect it to return an array with comptime known size
            const res_size = comptime count_bits(res_set);
            const operation = comptime Self.inverse_table_and_op_res_to_scatter_masks(table.snd, table.third, left_set, right_set, res_set, res_size);

            const left_arg: *Self = blk: {
                if(@TypeOf(bread) == LeftInner) {
                    break :blk &bread;
                }
                break :blk bread;
            };

            const res_to_use = blk: {
                if(@typeInfo(Result) == .Pointer) {
                    break :blk result;
                }
                break :blk &result;
            };

            Self.execute_sparse_vector_op(left_arg, meat, res_to_use, comptime operation, comptime res_size, T);
            return result;
        }

        fn get_square(self: *Self) Self {
            var x = self.copy(Self);
            var y = self.copy(Self);
            //var neg_one = Self{.terms=@splat(0)};
            return x.gp(&y, null);
        }

        pub fn exp(self: *Self, terms: usize) Self {
            @setEvalBranchQuota(1000000000);
            var ret = Self{.terms = @splat(0)};
            var cop = self.copy(Self);

            const square = self.get_square();
            // simple bivectors square to a scalar, and have closed form equations for exp
            // instead of a taylor series
            if(@reduce(.Add, square.terms) == square.terms[0]) {

                //debug.print("\n\t{} ^ 2 == {}, a scalar.", .{self, square});

                const squared_value = square.terms[0];
                const alpha = @sqrt(@abs(square.terms[0]));

                if(squared_value < 0.0) {
                    _=cop.mul_scalar(@sin(alpha) / alpha, &ret);
                    ret.terms[0] += @cos(alpha);
                    //debug.print("\n\t\texp({}) = cos(alpha) + (sin(alpha)/alpha)*self = {} (square is {d:.4}, alpha {d:.4})", .{self, ret, squared_value, alpha});
                    return ret;
                }
                if(squared_value > 0.0) {
                    _=cop.mul_scalar(sinh(alpha) / alpha, &ret);
                    ret.terms[0] += cosh(alpha);
                    //debug.print("\n\t\texp({}) = cosh(mag) + (sinh(mag)/mag)*self = {} (square is {}, mag {})", .{self, ret, alpha, mag});
                    return ret;
                }
                if(squared_value == 0.0) {
                    cop.terms[0] += 1.0;
                    //ebug.print("\n\t\texp({}) = 1 + self = {} (square is {}, mag {})", .{self, ret, alpha, self.magnitude()});
                    return cop;
                }    
            }
            //debug.print("\n\t{} ^ 2 == {}, not a scalar", .{self, square});

            //1 + x + (x^2)/2! + (x^3)/3! + (x^4)/4! + ...
            ret.terms[0] = 1.0;
            for(1..terms) |t| {
                const denominator: f64 = @floatFromInt(factorial(t));
                var term = self.copy(Self);

                for(1..t) |_| {
                    //var tst = NVec{.terms=@splat(0)};
                    _ = term.gp(&cop, null).add(term, &term);
                }

                //var tst = NVec{.terms=@splat(0)};
                var rhs = term.mul_scalar((1.0/denominator), Self);
                //@compileLog(comptimePrint("rhs type name {s}", .{@typeName(@TypeOf(rhs))}));
                ret = ret.add(&rhs, Self);
            }
            return ret;
        }

        ///returns true if we are a subset of other, or if our only nonzero indices are in other
        pub fn is(self: *Self, other: type) bool {
            if(comptime Self.subset_field & other.subset_field == Self.subset_field) {
                return true;
            }
            const one: ux = 0b1;
            for(0..alg_info.signature.basis_size) |blade| {
                const blade_t: ud = @truncate(blade);
                if(other.subset_field & (one << blade_t) != 0 and Self.subset_field & (one << blade_t) != 0) {
                    if(self.terms[count_bits(Self.subset_field & ((one << blade_t) -% 1))] != 0.0) {
                        return false;
                    }
                }
            }
            return true;
        }

    };
    return ret;
}

//requirements: 
//  need to be able to separate logic from implementation
//  need to be able to specify a dual or non dual algebra
//  implementation types need to know about each other when theyre part of the same algebra (implies a wrapping impl type around MVecSubset)
//  need to be able to name basis vectors whatever i want (x,y,z,e012), but also what those basis vectors square to (0i,+x,+y,+z,-t)
//  need to be able to define a standard ordering
//  need to be able to name specific subsets in the given algebra, each with their own ordering and possibly additional basis blade alias's ("Motor", "s,ix,iy,iz,xy,yz,zx")
//  implementors and consumers of those implementors need to be able to easily iterate through mvec arrays / find what basis blades correspond to an index
//  logic type must have a function equivalent to cayley_table in the original impl, but 1. comptime known vals work, and 2. different orderings of the basis blades work
//  the outer impl type should not be the horizon of what mvecs can operate with, the algebra that was defined could always be a subset of some other algebra with a larger signature

//idea 1:
//NO outer logic type that just takes (dual, p,z,n), but a bunch of functions that impls can use at comptime to create the data transformations they want at runtime
//const BasisBlade = struct {val: f64 = 1.0, metric: BasisMetric, binblade: usize} //should be generic to value?
//
//const GAOperation = enum{GeometricProduct, OuterProduct,...}
//
//(and another for unary ops)
//[for impls] determine the logical type both operands belong to -> create comptime basisblade arrays representing the values of each operand known at compile time not to be 0
//-> pass those to nary_cayley_table() -> get a tuple {result blade @ index -> 0..k [n operand blades] terms}
//-> given that, the impl type (logic types dont know about arrays) creates the scatter masks for lhs and rhs for the <=k terms in the result
//-> impl type calls a componentwise scatter add function given those masks
//if lhs term[i] is e21, and rhs term[j] is e12, but result[res idx] is in e12, then you do result[idx] += op(e21,e12)
//pub fn binary_cayley_table(comptime lhs: []BasisBlade, comptime rhs: []BasisBlade, comptime res_basis: []BasisBlade, comptime operation: Operation)

//const GAAlias = struct {
//  name: ?[]const u8, //= "Multivector", ?
//  runtime_subset: []const BasisVector,
//  comptime_subset: []const BasisVector,
//  specific_aliases: ?[]const struct{BasisVector, []const u8},
//
//
//}
//const pga3_aliases = .{.{name: "Motor", runtime_subset: "s,ix,iy,iz,xy,yz,zx"}, 
//                       .{name: "Point", runtime_subset: "izy,ixz,iyx,xyz", specific_aliases: "izy=x,ixz=y,iyx=z,xyz=i"},
//                       .{name: "NormPoint", runtime_subset: "izy,ixz,iyx", comp_subset: "xyz=1", specific_aliases: "izy=x,ixz=y,iyx=z,xyz=i"},
//                       .{name: "Line", runtime_subset: "ix,iy,iz,xy,yz,zx"},
//                       .{name: "IdealLine", runtime_subset: "ix,iy,iz"}}
// space info type
//fn GAlgebraInfo(dual = true, basis = "0i,+x,+y,+z", ordering="",aliases=pga3_aliases, use_e = false)
//
// actual runtime type. i think i just have to accept that the signature is going to be long
//impl type MVecSubset(comptime T: type, comptime alg: GAlgebraInfo, comptime runtime_subset: []const u8, comptime comptime_subset: []const u8, comptime name: ?[]const u8) type
//      //used to interface with logic fn's
//      const basis_blades: [basis_len]BasisBlade = create_basis_array(...)
//      terms: @Vector(T, our_len)
//
//      pub fn create_scatter_masks(...)
//      pub fn componentwise_multiscatter_add(...) void
//
//      pub fn BinaryResult(lhs: type, rhs: type, res: type, opr: GAOperation) type
//          when inferring res type, do what the old impl did except when inferring the output 
//          and its signature matches an alias'd type, use that alias specifically
//      pub fn get_result(...)
//
//      impl every op
//      pub fn gp(lhs: anytype, rhs: anytype, res: anytype) BinaryResult(...)
//          get basis blades of all
//          var res = get_result(...)
//          const inv = comptime ga_logic.binary_cayley_table(lhs_basis,rhs_basis,res_basis, GAOperation.GeometricProduct)
//          const opr = comptime create_scatter_masks(inv)
//          componentwise_multiscatter_add(&lhs,&rhs,&res,opr);
//          return res
   
// implementation type of a dense subset of a multivector following some GAlgebra
// point: the logical data transformations are dependent only on the metric (incl. dimensionality) + whether the algebra is dual or not
// consideration: what if i want to permute the orders of basis elements?
// consideration: what if i want to name the basis elements?
// consideration: what if i want to select blades based on defined properties? e.g. euclidian blades
// if all of the above are parametric, how do i find the geometric product of any 2 multivectors with different values of all of these, including subset?
//pub fn MVec(comptime alg: type, comptime T: type, )

pub fn main() !void {
    //const sig: [2]BasisVector = {}
    //const galg = GAlgebraFull(2, .{.{.name = &.{'x'}, .square = .positive}, .{.name = &.{'y'}, .square = .positive}});
    //galg.test_func();
    //debug.print("AAAAAAAA", .{});
    const ordering: []const u8 = "s,p,x,y,z,px,py,pz,xy,xz,yz,pxy,pyz,pxz,xyz,pxyz";//"ixy,xy,s,ix,iy,y,x,i";
    const point_alias: GAlgebraAlias = .{.name="Point", .ordering = "pyz,pxz,pxy,xyz"};
    const line_alias: GAlgebraAlias = .{.name = "Line", .ordering = "px,py,pz,yz,xz,xy"};
    const plane_alias: GAlgebraAlias = .{.name = "Plane", .ordering = "p,x,y,z"};
    const galg = GAlgebraInfo(true, "0p,+x,+y,+z", ordering, &.{point_alias, line_alias, plane_alias});
    const MVecT = MVecSubset(galg, f64, .{.find_alias = false, .values = &.{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}}); //s,p,x,y,z
    _ = MVecT;

    //@compileLog(comptimePrint("MVecT subset: {b}", .{MVecT.subset})`);
    //const mv1: MVecT = MVecT.init_raw(.{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0});
    //const mv2: MVecT = MVecT.init_raw(.{3.2, -4.1, 1.2, 7.1, -0.5, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0});
    ////@compileLog(comptimePrint("mv1: {}, mv2: {}", .{mv1, mv2}));
    //const mv3 = mv1.gp(mv2, null);
    //debug.print("{} * {} = \n\t{}\n\t ({} * {} = \n\t\t{})\n", .{mv1, mv2, mv3, mv1.terms, mv2.terms, mv3.terms});

    const Plane = MVecSubset(galg, f64, .{.alias_name = "Plane"});
    const Line = MVecSubset(galg, f64, .{.alias_name = "Line"});
    const Point = MVecSubset(galg, f64, .{.alias_name = "Point"});

    const plane1: Plane = Plane{.terms = .{1, 1, 0, 0}};//Plane.init_counting_up_from(1).normalize_by_blade(1);
    const plane2: Plane = Plane{.terms = .{1, 0, 1, 0}};//Plane.init_counting_up_from(5).normalize_by_blade(1);//.set_vals(&.{.{1, 6}, .{2, 1}, .{3, 9}});
    const plane3: Plane = Plane{.terms = .{1, 0, 0, 1}};//Plane.init_counting_up_from(9).normalize_by_blade(1);

    const line1: Line = plane1.meet(plane2, null);
    const line2: Line = Line.init_counting_up_from(1);

    const pseudoscalar = line1.meet(line2, null);

    const point2: Point = plane1.meet(plane2, null).meet(plane3, null);

    debug.print("\n{} meet {} = {}", .{plane1, plane2, line1});
    debug.print("\n{} meet {} = {}", .{line1, line2, pseudoscalar});
    debug.print("\n{} join {} join {} = {}", .{plane1, plane2, plane3, point2});

    //const mv4 = mv1.op(mv2, null);
    //debug.print("{} ^ {} = \n\t{} \n\t({} ^ {} = \n\t\t{})", .{mv1, mv2, mv4, mv1.terms, mv2.terms, mv4.terms});
    //@compileLog(comptimePrint("{} * {} = {}", .{mv1, mv2, mv3}));
}