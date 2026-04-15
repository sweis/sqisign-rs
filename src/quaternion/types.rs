//! Core quaternion-algebra types, ported from `quaternion.h` / `finit.c`.

use super::intbig::*;
use std::fmt;

// ---------------------------------------------------------------------------
// Integer vectors and matrices

pub type IbzVec2 = [Ibz; 2];
pub type IbzVec4 = [Ibz; 4];
pub type IbzMat2x2 = [[Ibz; 2]; 2];
pub type IbzMat4x4 = [[Ibz; 4]; 4];

/// Hand-rolled rational `[numerator, denominator]` (ported from `ibq_t`).
pub type Ibq = [Ibz; 2];

#[inline]
pub fn ibz_vec_2_init() -> IbzVec2 {
    Default::default()
}
#[inline]
pub fn ibz_vec_4_init() -> IbzVec4 {
    Default::default()
}
#[inline]
pub fn ibz_mat_2x2_init() -> IbzMat2x2 {
    Default::default()
}
#[inline]
pub fn ibz_mat_4x4_init() -> IbzMat4x4 {
    Default::default()
}

// ---------------------------------------------------------------------------
// Quaternion algebra

/// Quaternion algebra ramified at `p ≡ 3 (mod 4)` and ∞.
#[derive(Clone, Default)]
pub struct QuatAlg {
    pub p: Ibz,
}

impl QuatAlg {
    pub fn new(p: &Ibz) -> Self {
        Self { p: p.clone() }
    }
    pub fn from_ui(p: u32) -> Self {
        Self {
            p: super::intbig::ibz_from_i64(i64::from(p)),
        }
    }
}

/// Element of the quaternion algebra in basis `(1, i, j, ij)` with
/// `i² = -1`, `j² = -p`, `ij = k`, represented as `coord/denom`.
#[derive(Clone)]
pub struct QuatAlgElem {
    pub denom: Ibz,
    pub coord: IbzVec4,
}

impl Default for QuatAlgElem {
    fn default() -> Self {
        Self {
            denom: Ibz::from(1),
            coord: ibz_vec_4_init(),
        }
    }
}

/// Full-rank lattice in dimension 4: columns of `basis/denom` are algebra
/// elements in basis `(1,i,j,ij)`.
#[derive(Clone)]
pub struct QuatLattice {
    pub denom: Ibz,
    pub basis: IbzMat4x4,
}

impl Default for QuatLattice {
    fn default() -> Self {
        Self {
            denom: Ibz::from(1),
            basis: ibz_mat_4x4_init(),
        }
    }
}

/// Left ideal of a maximal order. `parent_order` mirrors C's
/// `const quat_lattice_t *` and points into static precomputed data.
#[derive(Clone, Default)]
pub struct QuatLeftIdeal {
    pub lattice: QuatLattice,
    pub norm: Ibz,
    pub parent_order: Option<&'static QuatLattice>,
}

/// Extremal maximal order (precomputed data).
#[derive(Clone, Default, Debug)]
pub struct QuatPExtremalMaximalOrder {
    pub order: QuatLattice,
    pub z: QuatAlgElem,
    pub t: QuatAlgElem,
    pub q: u32,
}

/// Parameters for `quat_represent_integer`.
#[derive(Debug)]
pub struct QuatRepresentIntegerParams {
    pub primality_test_iterations: i32,
    pub order: &'static QuatPExtremalMaximalOrder,
    pub algebra: &'static QuatAlg,
}

// ---------------------------------------------------------------------------
// Debug printing (ported from printer.c)

impl fmt::Debug for QuatAlg {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "quaternion algebra ramified at {} and infinity", self.p)
    }
}

impl fmt::Debug for QuatAlgElem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "denominator: {}\ncoordinates: {} {} {} {}",
            self.denom, self.coord[0], self.coord[1], self.coord[2], self.coord[3]
        )
    }
}

impl fmt::Debug for QuatLattice {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "lattice\ndenominator: {}", self.denom)?;
        write!(f, "basis: ")?;
        for i in 0..4 {
            for j in 0..4 {
                write!(f, "{} ", self.basis[i][j])?;
            }
            if i != 3 {
                write!(f, "\n       ")?;
            }
        }
        Ok(())
    }
}

impl fmt::Debug for QuatLeftIdeal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "left ideal\nnorm: {}", self.norm)?;
        writeln!(f, "{:?}", self.lattice)?;
        match self.parent_order {
            Some(o) => write!(f, "parent order {o:?}"),
            None => write!(f, "Parent order not given!"),
        }
    }
}
