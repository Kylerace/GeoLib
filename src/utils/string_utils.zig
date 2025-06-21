const std = @import("std");
const comptimePrint = std.fmt.comptimePrint;

const algebra = @import("../algebra.zig");
const BasisVector = algebra.BasisVector;

pub const BasisParseStates = enum {
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

pub fn parse_1basis(basis: []const u8) []const BasisVector {
    var parsed_1basis: []const BasisVector = &.{.{.name = "s", .blade = 0, .square = .positive}};

    var curr_parse_state: BasisParseStates = .ParseMetric;
    var curr_basis_blade: BasisVector = undefined;
    var curr_basis_i: usize = 0;
    for(0..basis.len+1) |i| {
        switch(curr_parse_state) {
            .ParseMetric => {
                const curr_char = basis[i];
                if(curr_char == ' ') {
                    continue;
                }
                if(curr_char == '+') {
                    curr_basis_blade.square = .positive;
                } else if(curr_char == '-') {
                    curr_basis_blade.square = .negative;
                } else if(curr_char == '0') {
                    curr_basis_blade.square = .zero;
                } else {
                    //@compileError(comptimePrint("bad char found in basis parsing of GAlgebraInfo! basis: {s}", .{basis}));
                    // no square symbol -> square is positive. we need to go back and reread this char though
                    curr_basis_blade.square = .positive;

                    curr_basis_i += 1;
                    curr_parse_state = .ParseNthBasisChar;

                    curr_basis_blade.name = &.{curr_char};
                    curr_basis_blade.blade = (1 << (curr_basis_i - 1));
                    curr_parse_state = .ParseNthBasisChar;
                    continue;
                }
                curr_basis_i += 1;
                curr_parse_state = .ParseFirstBasisChar;
            },
            .ParseFirstBasisChar => {
                const curr_char = basis[i];
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
                if(i == basis.len or basis[i] == ',') {
                    //const _curr_basis_blade = curr_basis_blade;
                    parsed_1basis = parsed_1basis ++ .{curr_basis_blade};
                    curr_basis_blade = undefined;
                    curr_parse_state = .ParseMetric;
                    continue;
                }
                const curr_char = basis[i];
                if((curr_char < 'A' or curr_char > 'z') or (curr_char > 'Z' and curr_char < 'a')) {
                    @compileError(comptimePrint("non letter char passed into basis string of GAlgebraInfo! _basis: {s}", .{basis})); 
                }
                curr_basis_blade.name = curr_basis_blade.name ++ .{curr_char};
            }
        }
    }

    return parsed_1basis;
}

/// given a string directing what 1 blades exist in this algebra, return the full algebra ordering where all blades are arranged in canonical (binary) ordering
pub fn canonical_ordering_string(comptime blades: []const u8) []const u8 {
    const parsed_1basis: []const BasisVector = parse_1basis(blades);
    return canonical_ordering_string_from_parsed_basis(parsed_1basis);
}

pub fn canonical_ordering_string_from_parsed_basis(parsed_basis: []const BasisVector) []const u8 {
    var ret: []const u8 = "";
    for(parsed_basis, 0..) |blade, i| {
        if(i > 0) {
            ret = ret ++ ",";
        }
        ret = ret ++ blade.name;
    }

    return ret;
}