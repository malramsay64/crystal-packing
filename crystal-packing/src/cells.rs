//
// cells.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! A cell list, for finding the shapes which could be in contact without checking every pair.
//!
//! Space is divided into a grid of cells at least as wide as the contact distance, so two shapes in
//! contact are in the same or neighbouring cells. Moving a shape then only needs checking against
//! the shapes in the 9 cells around it, rather than every shape. The grid is periodic, with each
//! position given as a fraction of the extent of the grid in each direction, which for a crystal
//! are the fractional coordinates in the reduced basis of the lattice.
//!
//! A grid of fewer than 3 cells in a direction can't rule anything out, so then the cell list is
//! disabled and every pair needs checking.

/// The smallest number of cells in each direction for the grid to be used
const MIN_CELLS: usize = 3;

/// The smallest number of items for the grid to be used, below which checking every pair is
/// quicker than looking through the neighbouring cells
const MIN_ITEMS: usize = 12;

#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct CellList {
    /// The number of cells in each direction, or zero when disabled
    nx: usize,
    ny: usize,
    /// The items in each cell, as (group, index), stored by row
    cells: Vec<Vec<(u32, u32)>>,
    /// The cell of each item, by group and index
    cell_of: Vec<Vec<usize>>,
}

impl CellList {
    /// The number of cells along a direction of a grid whose width across that direction is
    /// `width`, so each cell is at least `contact` wide, limited to `max`
    pub(crate) fn divisions(width: f64, contact: f64, max: usize) -> usize {
        // A small margin keeps cells wider than the contact distance despite rounding
        let n = libm::floor(width / (contact * (1. + 1e-9)));
        if n.is_finite() && n >= 0. {
            (n as usize).min(max)
        } else {
            0
        }
    }

    /// Whether a grid with this many cells is worth using for this many items, ruling out pairs
    pub(crate) fn useful(nx: usize, ny: usize, items: usize) -> bool {
        items >= MIN_ITEMS && nx >= MIN_CELLS && ny >= MIN_CELLS && nx * ny > 9
    }

    /// Whether there are enough items for a grid to be worth considering, which is checked
    /// before working out the size of the grid
    pub(crate) fn enough_items(items: usize) -> bool {
        items >= MIN_ITEMS
    }

    /// Stop using the grid, until it is rebuilt
    pub(crate) fn disable(&mut self) {
        self.nx = 0;
        self.ny = 0;
    }

    /// Rebuild the grid with the given number of cells for the items at the fractional positions
    /// given by `position(group, index)`, with `sizes[group]` items in each group.
    ///
    /// The grid is disabled when it would have too few cells to rule out any pairs.
    pub(crate) fn rebuild(
        &mut self,
        nx: usize,
        ny: usize,
        sizes: &[usize],
        position: impl Fn(usize, usize) -> (f64, f64),
    ) {
        let enabled = CellList::useful(nx, ny, sizes.iter().sum());
        let (nx, ny) = if enabled { (nx, ny) } else { (0, 0) };
        self.nx = nx;
        self.ny = ny;
        self.cells.resize_with(nx * ny, Vec::new);
        self.cells.truncate(nx * ny);
        for cell in self.cells.iter_mut() {
            cell.clear();
        }
        self.cell_of.resize_with(sizes.len(), Vec::new);
        for (group, &size) in sizes.iter().enumerate() {
            self.cell_of[group].clear();
            self.cell_of[group].resize(size, 0);
        }
        if !enabled {
            return;
        }
        for (group, &size) in sizes.iter().enumerate() {
            for index in 0..size {
                let (x, y) = position(group, index);
                let cell = self.cell(x, y);
                self.cells[cell].push((group as u32, index as u32));
                self.cell_of[group][index] = cell;
            }
        }
    }

    pub(crate) fn is_enabled(&self) -> bool {
        self.nx > 0
    }

    /// Whether the current grid can be kept when the largest numbers of cells which are still at
    /// least the contact distance wide become `nx` and `ny`.
    ///
    /// The cells must stay wide enough, and the grid is rebuilt when it could have twice as many
    /// cells in a direction, so it doesn't stay coarse after an expansion. A disabled grid is kept
    /// while a grid still wouldn't be useful.
    pub(crate) fn still_fits(&self, nx: usize, ny: usize, items: usize) -> bool {
        if self.is_enabled() {
            self.nx <= nx && self.ny <= ny && 2 * self.nx > nx && 2 * self.ny > ny
        } else {
            !CellList::useful(nx, ny, items)
        }
    }

    /// The cell containing a fractional position, wrapped into the grid
    fn cell(&self, x: f64, y: f64) -> usize {
        self.index(y, self.ny) * self.nx + self.index(x, self.nx)
    }

    /// The index of a fractional coordinate wrapped into n cells
    fn index(&self, f: f64, n: usize) -> usize {
        // The integer below f, by truncation corrected for negative values, which is exact and
        // much quicker than a call to floor
        let t = f as i64 as f64;
        let floor = if t > f { t - 1. } else { t };
        ((f - floor) * n as f64).min(n as f64 - 1.) as usize
    }

    /// Move an item to the cell of its new fractional position
    pub(crate) fn update(&mut self, group: usize, index: usize, x: f64, y: f64) {
        if !self.is_enabled() {
            return;
        }
        let new = self.cell(x, y);
        let old = self.cell_of[group][index];
        if new == old {
            return;
        }
        let item = (group as u32, index as u32);
        let cell = &mut self.cells[old];
        let position = cell
            .iter()
            .position(|&i| i == item)
            .expect("Every item is in its cell");
        cell.swap_remove(position);
        self.cells[new].push(item);
        self.cell_of[group][index] = new;
    }

    /// Call `f` with every item in the cells neighbouring the cell of an item, including its own,
    /// as (group, index), stopping early when `f` returns true
    pub(crate) fn any_neighbour(
        &self,
        group: usize,
        index: usize,
        mut f: impl FnMut(usize, usize) -> bool,
    ) -> bool {
        let cell = self.cell_of[group][index];
        let (cx, cy) = (cell % self.nx, cell / self.nx);
        for dy in [self.ny - 1, 0, 1] {
            let y = (cy + dy) % self.ny;
            for dx in [self.nx - 1, 0, 1] {
                let x = (cx + dx) % self.nx;
                for &(g, i) in self.cells[y * self.nx + x].iter() {
                    if f(g as usize, i as usize) {
                        return true;
                    }
                }
            }
        }
        false
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::random::Rng;

    #[test]
    fn divisions() {
        assert_eq!(CellList::divisions(10., 2., 100), 4);
        assert_eq!(CellList::divisions(10., 2., 3), 3);
        // Exactly divisible widths lose a cell to the margin, keeping cells wider than contact
        assert_eq!(CellList::divisions(8., 2., 100), 3);
        assert_eq!(CellList::divisions(f64::INFINITY, 2., 100), 0);
    }

    #[test]
    fn small_grids_are_disabled() {
        let mut cells = CellList::default();
        let position = |_, i: usize| (i as f64 / 20., 0.5);
        cells.rebuild(3, 3, &[20], position);
        assert!(!cells.is_enabled());
        cells.rebuild(2, 8, &[20], position);
        assert!(!cells.is_enabled());
        // Too few items to be worth it
        cells.rebuild(4, 4, &[4], position);
        assert!(!cells.is_enabled());
        cells.rebuild(4, 4, &[20], position);
        assert!(cells.is_enabled());
    }

    /// The items found near each item are exactly those within a cell width, periodically
    #[test]
    fn neighbours_match_brute_force() {
        let mut rng = Rng::seed_from_u64(1);
        let (nx, ny) = (5, 7);
        let mut positions: Vec<Vec<(f64, f64)>> = vec![vec![], vec![]];
        for group in positions.iter_mut() {
            for _ in 0..40 {
                group.push((3. * rng.uniform() - 1., rng.uniform()));
            }
        }
        let mut cells = CellList::default();
        let sizes: Vec<usize> = positions.iter().map(|g| g.len()).collect();
        cells.rebuild(nx, ny, &sizes, |g, i| positions[g][i]);
        // Move some items, keeping the grid up to date
        for _ in 0..30 {
            let (g, i) = (rng.index(2), rng.index(40));
            positions[g][i] = (rng.uniform(), 2. * rng.uniform());
            cells.update(g, i, positions[g][i].0, positions[g][i].1);
        }
        let frac = |f: f64| f - f.floor();
        // The periodic distance between two positions, in cell widths
        let apart = |a: (f64, f64), b: (f64, f64)| {
            let d = |p: f64, q: f64, n: usize| {
                let d = (frac(p) - frac(q)).abs() * n as f64;
                d.min(n as f64 - d)
            };
            d(a.0, b.0, nx).max(d(a.1, b.1, ny))
        };
        for g in 0..2 {
            for i in 0..40 {
                let mut found = vec![];
                cells.any_neighbour(g, i, |h, j| {
                    found.push((h, j));
                    false
                });
                for h in 0..2 {
                    for j in 0..40 {
                        let near = apart(positions[g][i], positions[h][j]) < 1.;
                        // Everything within a cell width is found
                        assert!(!near || found.contains(&(h, j)), "{} {} {} {}", g, i, h, j);
                    }
                }
                // Everything found is within two cell widths, so nothing is found twice
                assert!(found
                    .iter()
                    .all(|&(h, j)| apart(positions[g][i], positions[h][j]) < 2.));
                let mut unique = found.clone();
                unique.sort();
                unique.dedup();
                assert_eq!(unique.len(), found.len());
            }
        }
    }
}
