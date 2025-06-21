
# GeoLib
GeoLib is a math library in active development implementing<details><summary> Geometric Algebra </summary>
[Geometric algebra](https://www.wikipedia.com/wiki/geometric_algebra) is similar to linear algebra except basis vectors can be multiplied together to produce higher dimensional vectors called blades. An "algebra" or "algebra signature" is the result of n separate basis vectors which combine into 2^n total graded (the grade of a blade being its dimensionality) basis blades each of grades 0 to n (e.g. an algebra with input basis vectors x and y has the basis blades: scalar (grade 0), x (grade 1), y (grade 1), and xy (grade 2)). An element of a geometric algebra is called a multivector, which much like how a vector is the sum of multiples of its vector spaces basis vectors, a multivector is the sum of multiples of its geometric algebras basis blades. Each algebra with a given signature (p basis blades squaring to +1, q squaring to -1, r squaring to 0) has one full multivector with all 2^(n=p+q+r) elements. Yet the full multivector is rarely useful, oftentimes specific subsets (containing only some of the blades) of the full multivector holds properties we deem useful. For example in [Projective Geometric Algera](https://www.wikipedia.com/wiki/plane-based_geometric_algebra), you can represent any plane in 3d space in any position with only the e0, e1, e2, and e3 basis blades and you can represent any line in 3d space (also in any position) with the e01, e02, e03, e12, e23, and e13 blades. 
</details>using comptime generated multivector subsets. GeoLib is designed to counteract issues with most geometric algebra libraries in that they often require compromises in either performance or developer convenience and productivity. Geometric algebra libraries typically can't provide both arbitrary algebra signatures (which are useful for different applications but change the size and meaning of runtime data needed to be stored and operated on) and strong time and space performance guarantees. GeoLib currently does both and is being actively optimized to do so better. GeoLib allows for:
* Defining an arbitrary geometric algebra (GAlgebraInfo) with any default ordering of basis blades with any permutation of 1-blade factors for each 
* For a given algebra, defining a type (MVecSubset) that contains specific subsets of that algebra's full multivector in any ordering. An MVecSubset only contains at runtime an array representing the values of the basis blades its supposed to have
* Operating on MVecSubset's with geometric algebra operations which find what nonzero elements should exist in the result and output the corresponding MVecSubset by default, or force the value into a specific subset if explicitly directed to do so by the user. 

Basic usecase example for projective geometric algebra (PGA):
```zig

    //define what order the basis blades of the algebra will have unless the multivector subset is aliased (which can override ordering)
    const ordering: []const u8 = "s,p,x,y,z,px,py,pz,xy,xz,yz,pxy,pyz,pxz,xyz,pxyz";

    // define aliases for specific subsets and what order you want their elements to be. 
    // if the result of an operation matches the signature of a defined alias then that alias will be the result type unless the user explicitly passes in another one.
    // PGA has mvec subsets useful for planes, points, lines, directions, motors, etc. which can have any location in 3d euclidian space independent of an origin.
    const point_alias: algebra.GAlgebraAlias = .{.name="Point", .ordering = "pyz,pxz,pxy,xyz"};
    const line_alias: algebra.GAlgebraAlias = .{.name = "Line", .ordering = "px,py,pz,yz,xz,xy"};
    const plane_alias: algebra.GAlgebraAlias = .{.name = "Plane", .ordering = "p,x,y,z"};

    //define the algebra:
    //the first argument is whether its dual or not
    //second is the signature of the basis 1blades. p*p = 0, x*x = +1, if a blade -n was introduced it would square to -1
    //third is the exact ordering that unalias'd multivector subsets have (also defines the parity of each basis blade, e.g. if you define a blade as pzy instead of pyz some operations have a negative sign there instead of positive)
    //fourth is what aliases this algebra has.
    const galg = GAlgebraInfo(true, "0p,+x,+y,+z", ordering, &.{point_alias, line_alias, plane_alias});

    //define the full multivector type of this algebra, alongside what datatype its internal array is holding. .subset takes in a bitfield where bit i = 1 iff the subset contains the ith blade, it allows for defining what elements you want without caring what order those elements are in (which would give it the algebra's default ordering)
    // a more convenient way to define multivectors is being actively worked on. 
    const MVecT = MVecSubset(galg, f64, .{.find_alias = false, .subset = ~@as(usize, 0)});
    _ = MVecT;

    //create zig types for each multivector alias, .alias_name will make it find the corresponding alias
    const Plane = MVecSubset(galg, f64, .{.alias_name = "Plane"});
    const Line = MVecSubset(galg, f64, .{.alias_name = "Line"});
    const Point = MVecSubset(galg, f64, .{.alias_name = "Point"});

    //create actual instances of planes, .terms is the actual array hoding the data. easier ways of initializing and changing arrays based on basis blades is being actively worked on.
    const plane1: Plane = Plane{.terms = .{1, 1, 0, 0}};
    const plane2: Plane = Plane{.terms = .{1, 0, 1, 0}};
    const plane3: Plane = Plane{.terms = .{1, 0, 0, 1}};

    // meet is an operation that in dual algebras (like in PGA) increases the grade of each basis blade in the result. In PGA planes have 1-grade blades, lines have 2-grade blades, and points have 3-grade blades, so the meet of two planes (1 grade elements) produces a line (2 grade elements)
    // the second arguement is an optional type for the result, if its non null the result will be forced to be that type. if null its inferred
    const line1: Line = plane1.meet(plane2, null);
    const line2: Line = Line.init_counting_up_from(1); //first elem = 1, second = 2, ...

    //most lines dont intersect, so their meet will usually not be a geometric object
    const pseudoscalar = line1.meet(line2, null);

    // the meet (intersection) of three orthogonal planes is a point
    const point2: Point = plane1.meet(plane2, null).meet(plane3, null);

    //Plane {p: 1.000, x: 1.000, y: 0.000, z: 0.000} meet Plane {p: 1.000, x: 0.000, y: 1.000, z: 0.000} = Line {px: -1.000, py: 1.000, pz: 0.000, yz: 0.000, xz: 0.000, xy: 1.000}
    debug.print("\n{} meet {} = {}", .{plane1, plane2, line1});

    //Line {px: -1.000, py: 1.000, pz: 0.000, yz: 0.000, xz: 0.000, xy: 1.000} meet Line {px: 1.000, py: 2.000, pz: 3.000, yz: 4.000, xz: 5.000, xy: 6.000} = Multivector {pxyz: -6.000}
    debug.print("\n{} meet {} = {}", .{line1, line2, pseudoscalar});

    //Plane {p: 1.000, x: 1.000, y: 0.000, z: 0.000} join Plane {p: 1.000, x: 0.000, y: 1.000, z: 0.000} join Plane {p: 1.000, x: 0.000, y: 0.000, z: 1.000} = Point {pyz: 1.000, pxz: -1.000, pxy: 1.000, xyz: 1.000}
    debug.print("\n{} join {} join {} = {}", .{plane1, plane2, plane3, point2});
```

# Building
run `zig build run` on zig version 13.0
