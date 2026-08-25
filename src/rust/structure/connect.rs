use itertools::Itertools;
use numpy::ndarray::Array1;
use numpy::ndarray::ArrayView2;
use numpy::{IntoPyArray, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions;
use pyo3::prelude::*;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet, VecDeque};

use crate::structure::bonds::{Bond, BondList, BondType};
use crate::util::check_signals_periodically;

/// Internal rust function for :func:`biotite.structure.connect.connect_via_residue_names`.
///
/// Parameters
/// ----------
/// res_names
///     The residue name of each atom.
/// atom_names
///     The atom name of each atom.
/// residue_starts
///     The start index of each residue (including exclusive stop index).
/// bond_dict
///     The bond dictionary mapping connected atom name pairs to their bond type if they form a
///     covalent bond.
///
/// Returns
/// -------
/// bonds
///     The bonds between the atoms.
#[pyfunction]
pub fn connect_via_residue_names<'py>(
    _py: Python<'py>,
    res_names: Vec<String>,
    atom_names: Vec<String>,
    residue_starts: PyReadonlyArray1<'py, i64>,
    bond_dict: HashMap<String, Vec<(String, String, BondType)>>,
) -> PyResult<BondList> {
    let mut bond_list = BondList::empty(atom_names.len());
    for (start_i, stop_i) in residue_starts.as_array().iter().tuple_windows() {
        let start_i = *start_i as usize;
        let stop_i = *stop_i as usize;
        let bonds_res = match bond_dict.get(&res_names[start_i]) {
            Some(bonds_res) => bonds_res,
            None => continue,
        };
        let atom_names_in_res = &atom_names[start_i..stop_i];
        for (atom_name1, atom_name2, bond_type) in bonds_res.iter() {
            let atom_indices_1 = find(atom_names_in_res, atom_name1);
            let atom_indices_2 = find(atom_names_in_res, atom_name2);
            // In rare cases the same atom name may appear multiple times
            // (e.g. in altlocs)
            // -> create all possible bond combinations
            for atom_index_i in atom_indices_1.iter() {
                for atom_index_j in atom_indices_2.iter() {
                    unsafe {
                        bond_list.add(Bond::new(
                            (start_i + atom_index_i) as isize,
                            (start_i + atom_index_j) as isize,
                            *bond_type,
                            atom_names.len(),
                        )?)?;
                    }
                }
            }
        }
    }
    bond_list.remove_redundant_bonds();
    Ok(bond_list)
}

/// Internal rust function for :func:`biotite.structure.connect.connect_via_distances`.
///
/// Parameters
/// ----------
/// coord
///     The coordinates of each atom (shape: Nx3, dtype: float32).
/// radii
///     The covalent radius of each atom.
///     NaN for atoms with unknown radius, which consequently never form bonds.
/// residue_starts
///     The start index of each residue (including exclusive stop index).
/// tolerance
///     The tolerance added to the sum of two covalent radii to obtain the maximum
///     bond distance.
/// default_bond_type
///     The bond type to assign to all created bonds.
///
/// Returns
/// -------
/// bonds
///     The bonds between the atoms.
#[pyfunction]
pub fn connect_via_distances<'py>(
    _py: Python<'py>,
    coord: PyReadonlyArray2<'py, f32>,
    radii: PyReadonlyArray1<'py, f32>,
    residue_starts: PyReadonlyArray1<'py, i64>,
    tolerance: f32,
    default_bond_type: BondType,
) -> PyResult<BondList> {
    let coord = coord.as_array();
    let radii = radii.as_array();
    let residue_starts = residue_starts.as_array();
    let n_atoms = radii.len();

    let mut bond_list = BondList::empty(n_atoms);

    for (start, stop) in residue_starts.iter().tuple_windows() {
        let start = *start as usize;
        let stop = *stop as usize;

        for atom_i in start..stop {
            for atom_j in start..atom_i {
                let max_distance = radii[atom_i] + radii[atom_j] + tolerance;
                // Compare squared distances to avoid expensive square roots
                // Both comparands are NaN for atoms with unknown radius or missing
                // coordinates, which makes the comparison false
                if squared_distance(&coord, atom_i, atom_j) <= max_distance * max_distance {
                    unsafe {
                        // atom_j < atom_i is guaranteed by the loop range
                        bond_list.add(Bond {
                            atom1: atom_j,
                            atom2: atom_i,
                            bond_type: default_bond_type,
                        })?;
                    }
                }
            }
        }
    }

    Ok(bond_list)
}

/// Internal rust function for :func:`biotite.structure.connect._connect_inter_residue`.
///
/// Parameters
/// ----------
/// atom_names
///     The atom name of each atom.
/// residue_starts
///     The start index of each residue (including exclusive stop index).
/// link_types
///     The link type of each residue. None if the residue is not found in the CCD.
/// is_disconnected
///     Whether each residue is disconnected from the next residue due to chain or
///     residue discontinuity.
///     Length is ``len(residue_starts) - 2``.
#[pyfunction]
pub fn connect_inter_residue<'py>(
    _py: Python<'py>,
    atom_names: Vec<String>,
    residue_starts: PyReadonlyArray1<'py, i64>,
    link_types: Vec<Option<String>>,
    is_disconnected: PyReadonlyArray1<'py, bool>,
) -> PyResult<BondList> {
    let residue_starts = residue_starts.as_array();
    let is_disconnected = is_disconnected.as_array();
    let n_residues = residue_starts.len() - 1;

    let mut bond_list = BondList::empty(atom_names.len());

    if n_residues < 2 {
        return Ok(bond_list);
    }
    for i in 0..n_residues - 1 {
        // Check if residues are disconnected
        if is_disconnected[i] {
            continue;
        }

        let curr_start = residue_starts[i] as usize;
        let curr_stop = residue_starts[i + 1] as usize;
        let next_start = residue_starts[i + 1] as usize;
        let next_stop = residue_starts[i + 2] as usize;

        // Get link type for both residues
        let (curr_link, next_link) = match (&link_types[i], &link_types[i + 1]) {
            (Some(curr), Some(next)) => (curr.as_str(), next.as_str()),
            _ => continue,
        };

        let (curr_connect_name, next_connect_name) =
            if is_peptide_link(curr_link) && is_peptide_link(next_link) {
                ("C", "N")
            } else if is_nucleic_link(curr_link) && is_nucleic_link(next_link) {
                ("O3'", "P")
            } else {
                // Create no bond if the connection types of consecutive
                // residues are not compatible
                continue;
            };

        // Find first connector atom in each residue
        let curr_connect_idx = atom_names[curr_start..curr_stop]
            .iter()
            .position(|name| name == curr_connect_name);
        let next_connect_idx = atom_names[next_start..next_stop]
            .iter()
            .position(|name| name == next_connect_name);

        if let (Some(curr_idx), Some(next_idx)) = (curr_connect_idx, next_connect_idx) {
            unsafe {
                // curr_start < next_start is guaranteed
                // by iterating consecutive residues
                bond_list.add(Bond {
                    atom1: curr_start + curr_idx,
                    atom2: next_start + next_idx,
                    bond_type: BondType::Single,
                })?;
            }
        }
    }

    Ok(bond_list)
}

/// Bond orders higher than a triple bond are not assigned.
const MAX_BOND_ORDER: u8 = 3;

/// Internal rust function for :func:`biotite.structure.connect.infer_bond_types`.
///
/// Implements the bond order assignment algorithm from Kim & Kim
/// (https://doi.org/10.1002/bkcs.10334):
/// *Valence states*, i.e. combinations of the allowed valences of all atoms, are
/// enumerated lazily via best-first search.
/// In deviation from the reference, which enumerates purely in decreasing order of
/// total valence, states are ordered by increasing total absolute formal charge first
/// and decreasing total valence second:
/// This finds the chemically preferred state much faster (the state space grows
/// exponentially with the number of multivalent atoms) and prefers resonance
/// structures with fewer formal charges.
/// For each valence state multiple bonds are greedily assigned between adjacent
/// unsaturated atoms, restarting from different atoms if necessary.
/// The first assignment that consumes all unsaturations and whose formal charges sum
/// up to the given total charge is returned.
/// If no such assignment exists (or `max_iterations` is exceeded), the assignment with
/// the maximum total bond order found so far is returned as fallback.
///
/// Parameters
/// ----------
/// elements
///     The element of each atom in upper case.
/// bonds
///     The connectivity between the atoms.
///     The bond types are ignored.
/// total_charge
///     The total formal charge of the atoms, used as constraint.
/// max_iterations
///     The maximum number of bond order assignments to try.
///     If `None`, the search is only stopped, when all valence states are tried.
///
/// Returns
/// -------
/// bonds
///     The bonds with assigned bond orders.
/// charges
///     The formal charge of each atom, as implied by the assigned bond orders.
/// converged
///     Whether an assignment satisfying the unsaturation and charge conditions was
///     found.
#[pyfunction]
#[pyo3(signature = (elements, bonds, total_charge, max_iterations=None))]
pub fn infer_bond_types(
    py: Python<'_>,
    elements: Vec<String>,
    bonds: &BondList,
    total_charge: i32,
    max_iterations: Option<u64>,
) -> PyResult<(BondList, Py<PyAny>, bool)> {
    // No limit on the number of iterations means an exhaustive search,
    // which is already bounded by the finite number of valence states
    let max_iterations = max_iterations.unwrap_or(u64::MAX);
    let n_atoms = elements.len();
    if bonds.get_atom_count() != n_atoms {
        return Err(exceptions::PyValueError::new_err(format!(
            "Bond list represents {} atoms, but {} elements are given",
            bonds.get_atom_count(),
            n_atoms
        )));
    }

    let edges: Vec<(usize, usize)> = bonds
        .get_bonds_ref()
        .iter()
        .map(|bond| (bond.atom1, bond.atom2))
        .collect();
    // Adjacent atoms of each atom, given as `(partner_atom, edge_index)`
    let mut adjacency: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n_atoms];
    for (edge_i, &(atom1, atom2)) in edges.iter().enumerate() {
        adjacency[atom1].push((atom2, edge_i));
        adjacency[atom2].push((atom1, edge_i));
    }
    let degrees: Vec<u32> = adjacency.iter().map(|adj| adj.len() as u32).collect();

    // Every bond starts as a single bond
    let base_orders: Vec<u8> = vec![1; edges.len()];

    // Allowed valences of each atom, sorted in decreasing order
    // Atoms with unknown elements, no bonds at all or more bonds than any allowed
    // valence are treated as saturated
    let options: Vec<Vec<u32>> = (0..n_atoms)
        .map(|i| {
            if degrees[i] == 0 {
                return vec![degrees[i]];
            }
            match allowed_valences(&elements[i], degrees[i]) {
                Some(allowed) => {
                    let valid: Vec<u32> = allowed
                        .iter()
                        .copied()
                        .filter(|&valence| valence >= degrees[i])
                        .collect();
                    if valid.is_empty() {
                        vec![degrees[i]]
                    } else {
                        valid
                    }
                }
                None => vec![degrees[i]],
            }
        })
        .collect();
    // Only atoms with more than one allowed valence span the valence state space
    let multivalent: Vec<usize> = (0..n_atoms).filter(|&i| options[i].len() > 1).collect();

    let make_state = |option_indices: Vec<u8>| -> ValenceState {
        let mut valence_sum = 0;
        let mut charge_sum = 0;
        for (pos, &atom) in multivalent.iter().enumerate() {
            let valence = options[atom][option_indices[pos] as usize];
            valence_sum += valence;
            // The formal charge the atom would have if its valence is fully satisfied,
            // used to prefer states with fewer formal charges
            charge_sum += formal_charge(&elements[atom], 0, valence, 0).abs();
        }
        ValenceState {
            valence_sum,
            neg_charge_sum: -charge_sum,
            option_indices,
        }
    };

    // Enumerate the valence states lazily in decreasing order of total valence
    // to avoid materializing the potentially exponential state space
    let mut heap: BinaryHeap<ValenceState> = BinaryHeap::new();
    let mut visited: HashSet<Vec<u8>> = HashSet::new();
    let initial_indices = vec![0u8; multivalent.len()];
    visited.insert(initial_indices.clone());
    heap.push(make_state(initial_indices));

    let mut iteration: u64 = 0;
    let mut converged = false;
    let mut best_orders: Vec<u8> = base_orders.clone();
    let mut best_total: u64 = 0;

    'search: while let Some(state) = heap.pop() {
        // Enqueue all successor states,
        // i.e. states where a single atom has the next lower allowed valence
        for pos in 0..multivalent.len() {
            if (state.option_indices[pos] as usize) < options[multivalent[pos]].len() - 1 {
                let mut successor = state.option_indices.clone();
                successor[pos] += 1;
                if visited.insert(successor.clone()) {
                    heap.push(make_state(successor));
                }
            }
        }

        // Determine the unsaturation of each atom in this valence state
        let mut option_index_of_atom = vec![0usize; n_atoms];
        for (pos, &atom) in multivalent.iter().enumerate() {
            option_index_of_atom[atom] = state.option_indices[pos] as usize;
        }
        let initial_u: Vec<i32> = (0..n_atoms)
            .map(|i| options[i][option_index_of_atom[i]] as i32 - degrees[i] as i32)
            .collect();

        let unsaturated: Vec<usize> = (0..n_atoms).filter(|&i| initial_u[i] > 0).collect();

        // Quick rejection: If any unsaturated atom has only saturated neighbors,
        // no assignment in this valence state can consume all unsaturations
        if unsaturated
            .iter()
            .any(|&i| adjacency[i].iter().all(|&(j, _)| initial_u[j] <= 0))
        {
            iteration += 1;
            check_signals_periodically(py, iteration as usize)?;
            if iteration >= max_iterations {
                break 'search;
            }
            continue;
        }

        // The result of the greedy assignment depends on the starting atom
        // -> retry with different starting atoms until the conditions are met
        let starts: Vec<Option<usize>> = if unsaturated.is_empty() {
            vec![None]
        } else {
            unsaturated.iter().map(|&i| Some(i)).collect()
        };
        for start in starts {
            iteration += 1;
            check_signals_periodically(py, iteration as usize)?;
            let (orders, remaining_u) =
                assign_bond_orders(&adjacency, &base_orders, &initial_u, start);
            let charge_sum: i32 = (0..n_atoms)
                .map(|i| {
                    let bo_sum: u32 = adjacency[i].iter().map(|&(_, e)| orders[e] as u32).sum();
                    formal_charge(&elements[i], degrees[i], bo_sum, total_charge)
                })
                .sum();
            if remaining_u == 0 && charge_sum == total_charge {
                best_orders = orders;
                converged = true;
                break 'search;
            }
            let total: u64 = orders.iter().map(|&order| order as u64).sum();
            if total > best_total {
                best_total = total;
                best_orders = orders;
            }
            if iteration >= max_iterations {
                break 'search;
            }
        }
    }

    let mut bond_list = BondList::empty(n_atoms);
    for (edge_i, &(atom1, atom2)) in edges.iter().enumerate() {
        let bond_type = match best_orders[edge_i] {
            1 => BondType::Single,
            2 => BondType::Double,
            _ => BondType::Triple,
        };
        unsafe {
            // The input bond list guarantees `atom1 < atom2` and unique bonds
            bond_list.add(Bond {
                atom1,
                atom2,
                bond_type,
            })?;
        }
    }

    let charges: Vec<i64> = (0..n_atoms)
        .map(|i| {
            let bo_sum: u32 = adjacency[i]
                .iter()
                .map(|&(_, e)| best_orders[e] as u32)
                .sum();
            formal_charge(&elements[i], degrees[i], bo_sum, total_charge) as i64
        })
        .collect();

    Ok((
        bond_list,
        Array1::from_vec(charges)
            .into_pyarray(py)
            .into_any()
            .unbind(),
        converged,
    ))
}

/// A combination of chosen valences for all multivalent atoms.
///
/// Ordered such that the maximum element has the lowest sum of absolute formal
/// charges and among ties the highest total valence.
#[derive(PartialEq, Eq)]
struct ValenceState {
    valence_sum: u32,
    neg_charge_sum: i32,
    option_indices: Vec<u8>,
}

impl Ord for ValenceState {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.neg_charge_sum, self.valence_sum, &self.option_indices).cmp(&(
            other.neg_charge_sum,
            other.valence_sum,
            &other.option_indices,
        ))
    }
}

impl PartialOrd for ValenceState {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Allowed valences of the supported elements
/// (Table 1 in Kim & Kim, 2015, extended by further organic elements),
/// sorted in decreasing order.
/// `None` for unsupported elements.
///
/// In deviation from the reference, the allowed valences of sulfur and selenium are
/// restricted based on the number of bond partners:
/// The formal charge is zero for every sulfur valence, so valence states are ranked by
/// their total valence alone.
/// Hence, without this restriction the hypervalent state would be chosen for e.g.
/// aromatic sulfur, whose ring neighbors can accept the additional bond orders.
/// Furthermore, this allows thiolates (valence 1) and converging sulfoxides/sulfones
/// (formal charge 0 despite exceeding the octet).
fn allowed_valences(element: &str, degree: u32) -> Option<&'static [u32]> {
    match element {
        "H" | "F" | "CL" | "BR" | "I" => Some(&[1]),
        "B" => Some(&[3]),
        "C" | "SI" => Some(&[4]),
        "N" => Some(&[4, 3]),
        "O" => Some(&[2, 1]),
        "P" => Some(&[5, 4, 3]),
        "S" | "SE" => match degree {
            1 => Some(&[2, 1]),
            2 => Some(&[2]),
            3 => Some(&[4, 3]),
            _ => Some(&[6, 4]),
        },
        _ => None,
    }
}

/// Number of valence electrons of the supported elements.
fn valence_electrons(element: &str) -> Option<i32> {
    match element {
        "H" => Some(1),
        "B" => Some(3),
        "C" | "SI" => Some(4),
        "N" | "P" => Some(5),
        "O" | "S" | "SE" => Some(6),
        "F" | "CL" | "BR" | "I" => Some(7),
        _ => None,
    }
}

/// The formal charge of an atom based on the sum of its bond orders
/// (Table 2 in Kim & Kim, 2015).
/// The formal charge of unsupported elements is set to 0.
fn formal_charge(element: &str, degree: u32, bond_order_sum: u32, total_charge: i32) -> i32 {
    match element {
        "H" => 0,
        "B" => 3 - bond_order_sum as i32,
        // Carbon with three single bonds: +1/-1 depending on the total charge
        "C" | "SI" if degree == 3 && bond_order_sum == 3 => total_charge.signum(),
        // Carbene-like carbon with two single bonds
        "C" | "SI" if degree == 2 && bond_order_sum == 2 => 0,
        // Hypervalent phosphorus and sulfur have an expanded octet
        "P" if bond_order_sum == 5 => 0,
        "S" | "SE" if bond_order_sum == 4 || bond_order_sum == 6 => 0,
        _ => match valence_electrons(element) {
            // Octet rule
            Some(n_electrons) => n_electrons - 8 + bond_order_sum as i32,
            None => 0,
        },
    }
}

/// Greedily assign multiple bonds between adjacent unsaturated atoms
/// (Figure 1 in Kim & Kim, 2015):
/// In each step the bond order between the unsaturated atom with the fewest pairable
/// partners and its pairable partner with the fewest pairable partners is increased,
/// until no two adjacent unsaturated atoms are left.
/// `start` optionally forces the atom initiating the first pairing.
/// Returns the assigned bond orders and the total remaining unsaturation.
fn assign_bond_orders(
    adjacency: &[Vec<(usize, usize)>],
    base_orders: &[u8],
    initial_u: &[i32],
    start: Option<usize>,
) -> (Vec<u8>, i32) {
    let mut orders = base_orders.to_vec();
    let mut u = initial_u.to_vec();
    let n_atoms = u.len();
    let mut forced_start = start;

    loop {
        // Number of partners each unsaturated atom can still pair with
        let mut n_partners = vec![0u32; n_atoms];
        for i in 0..n_atoms {
            if u[i] <= 0 {
                continue;
            }
            n_partners[i] = adjacency[i]
                .iter()
                .filter(|&&(j, e)| u[j] > 0 && orders[e] < MAX_BOND_ORDER)
                .count() as u32;
        }

        let current = match forced_start
            .take()
            .filter(|&s| u[s] > 0 && n_partners[s] > 0)
        {
            Some(s) => s,
            None => {
                match (0..n_atoms)
                    .filter(|&i| u[i] > 0 && n_partners[i] > 0)
                    .min_by_key(|&i| n_partners[i])
                {
                    Some(i) => i,
                    // No pairable atoms left
                    None => break,
                }
            }
        };
        let &(partner, edge) = adjacency[current]
            .iter()
            .filter(|&&(j, e)| u[j] > 0 && orders[e] < MAX_BOND_ORDER)
            .min_by_key(|&&(j, _)| n_partners[j])
            // Guaranteed by `n_partners[current] > 0`
            .unwrap();
        orders[edge] += 1;
        u[current] -= 1;
        u[partner] -= 1;
    }

    (orders, u.iter().filter(|&&du| du > 0).sum())
}

/// Get indices to all atoms that are directly or indirectly connected
/// to the root atom indicated by the given index.
///
/// An atom is *connected* to the `root` atom, if that atom is reachable
/// by traversing an arbitrary number of bonds, starting from the
/// `root`.
/// Effectively, this means that all atoms are *connected* to `root`,
/// that are in the same molecule as `root`.
/// Per definition `root` is also *connected* to itself.
///
/// Parameters
/// ----------
/// bond_list : BondList
///     The reference bond list.
/// root : int
///     The index of the root atom.
/// as_mask : bool, optional
///     If true, the connected atom indices are returned as boolean
///     mask.
///     By default, the connected atom indices are returned as integer
///     array.
///
/// Returns
/// -------
/// connected : ndarray, dtype=int or ndarray, dtype=bool
///     Either a boolean mask or an integer array, representing the
///     connected atoms.
///     In case of a boolean mask: ``connected[i] == True``, if the atom
///     with index ``i`` is connected.
///
/// Examples
/// --------
/// Consider a system with 4 atoms, where only the last atom is not
/// bonded with the other ones (``0-1-2 3``):
///
/// >>> bonds = BondList(4)
/// >>> bonds.add_bond(0, 1)
/// >>> bonds.add_bond(1, 2)
/// >>> print(find_connected(bonds, 0))
/// [0 1 2]
/// >>> print(find_connected(bonds, 1))
/// [0 1 2]
/// >>> print(find_connected(bonds, 2))
/// [0 1 2]
/// >>> print(find_connected(bonds, 3))
/// [3]
///
/// The result can be output as a boolean mask:
///
/// >>> print(find_connected(bonds, 0, as_mask=True))
/// [ True  True  True False]
#[pyfunction]
#[pyo3(signature = (bond_list, root, as_mask=false))]
pub fn find_connected(
    py: Python<'_>,
    bond_list: &BondList,
    root: usize,
    as_mask: bool,
) -> PyResult<Py<PyAny>> {
    let atom_count = bond_list.get_atom_count();
    let bonds = bond_list.get_bonds_ref();

    if root >= atom_count {
        return Err(exceptions::PyValueError::new_err(format!(
            "Root atom index {} is out of bounds for bond list representing {} atoms",
            root, atom_count
        )));
    }

    // Build adjacency list
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); atom_count];
    for bond in bonds {
        adj[bond.atom1].push(bond.atom2);
        adj[bond.atom2].push(bond.atom1);
    }

    // Iterative breadth-first search
    let mut visited = vec![false; atom_count];
    let mut queue = VecDeque::new();
    queue.push_back(root);
    visited[root] = true;

    while let Some(current) = queue.pop_front() {
        for &neighbor in &adj[current] {
            if !visited[neighbor] {
                visited[neighbor] = true;
                queue.push_back(neighbor);
            }
        }
    }

    if as_mask {
        let mask = Array1::from_vec(visited);
        Ok(mask.into_pyarray(py).into_any().unbind())
    } else {
        let indices: Vec<i64> = visited
            .iter()
            .enumerate()
            .filter(|(_, &v)| v)
            .map(|(i, _)| i as i64)
            .collect();
        let arr = Array1::from_vec(indices);
        Ok(arr.into_pyarray(py).into_any().unbind())
    }
}

/// Squared Euclidean distance between two atoms in an Nx3 coordinate array.
#[inline(always)]
fn squared_distance(coord: &ArrayView2<f32>, a1: usize, a2: usize) -> f32 {
    let dx = coord[[a1, 0]] - coord[[a2, 0]];
    let dy = coord[[a1, 1]] - coord[[a2, 1]];
    let dz = coord[[a1, 2]] - coord[[a2, 2]];
    dx * dx + dy * dy + dz * dz
}

fn is_peptide_link(link_type: &str) -> bool {
    matches!(
        link_type,
        "PEPTIDE LINKING" | "L-PEPTIDE LINKING" | "D-PEPTIDE LINKING"
    )
}

fn is_nucleic_link(link_type: &str) -> bool {
    matches!(link_type, "RNA LINKING" | "DNA LINKING")
}

/// Find all indices where a value occurs in an array.
fn find<T: PartialEq>(array: &[T], value: &T) -> Vec<usize> {
    array
        .iter()
        .enumerate()
        .filter(|(_, x)| *x == value)
        .map(|(i, _)| i)
        .collect()
}
