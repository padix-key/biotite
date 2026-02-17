use numpy::ndarray::{Array1, Array2};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray2};
use pyo3::exceptions;
use pyo3::prelude::*;
use pyo3::types::{PyModule, PySet, PyTuple};
use std::collections::{HashMap, HashSet};
use std::convert::TryFrom;

/// This enum type represents the type of a chemical bond.
///
/// - ``ANY`` - Used if the actual type is unknown
/// - ``SINGLE`` - Single bond
/// - ``DOUBLE`` - Double bond
/// - ``TRIPLE`` - Triple bond
/// - ``QUADRUPLE`` - A quadruple bond
/// - ``AROMATIC_SINGLE`` - Aromatic bond with a single formal bond
/// - ``AROMATIC_DOUBLE`` - Aromatic bond with a double formal bond
/// - ``AROMATIC_TRIPLE`` - Aromatic bond with a triple formal bond
/// - ``COORDINATION`` - Coordination complex involving a metal atom
/// - ``AROMATIC`` - Aromatic bond without specification of the formal bond
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
#[repr(u8)]
pub enum BondType {
    Any = 0,
    Single = 1,
    Double = 2,
    Triple = 3,
    Quadruple = 4,
    AromaticSingle = 5,
    AromaticDouble = 6,
    AromaticTriple = 7,
    Coordination = 8,
    Aromatic = 9,
}

impl BondType {
    fn without_aromaticity(&self) -> Self {
        match self {
            BondType::AromaticSingle => BondType::Single,
            BondType::AromaticDouble => BondType::Double,
            BondType::AromaticTriple => BondType::Triple,
            BondType::Aromatic => BondType::Any,
            _ => *self,
        }
    }
}

#[pyfunction]
pub fn bond_type_members() -> HashMap<String, u8> {
    let mut map = HashMap::new();
    map.insert("ANY".to_string(), BondType::Any as u8);
    map.insert("SINGLE".to_string(), BondType::Single as u8);
    map.insert("DOUBLE".to_string(), BondType::Double as u8);
    map.insert("TRIPLE".to_string(), BondType::Triple as u8);
    map.insert("QUADRUPLE".to_string(), BondType::Quadruple as u8);
    map.insert(
        "AROMATIC_SINGLE".to_string(),
        BondType::AromaticSingle as u8,
    );
    map.insert(
        "AROMATIC_DOUBLE".to_string(),
        BondType::AromaticDouble as u8,
    );
    map.insert(
        "AROMATIC_TRIPLE".to_string(),
        BondType::AromaticTriple as u8,
    );
    map.insert("COORDINATION".to_string(), BondType::Coordination as u8);
    map.insert("AROMATIC".to_string(), BondType::Aromatic as u8);
    map
}

impl TryFrom<u32> for BondType {
    type Error = PyErr;

    fn try_from(value: u32) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(BondType::Any),
            1 => Ok(BondType::Single),
            2 => Ok(BondType::Double),
            3 => Ok(BondType::Triple),
            4 => Ok(BondType::Quadruple),
            5 => Ok(BondType::AromaticSingle),
            6 => Ok(BondType::AromaticDouble),
            7 => Ok(BondType::AromaticTriple),
            8 => Ok(BondType::Coordination),
            9 => Ok(BondType::Aromatic),
            _ => Err(exceptions::PyValueError::new_err(format!(
                "BondType {} is invalid",
                value
            ))),
        }
    }
}

/// Convert from Python's BondType IntEnum to Rust's BondType.
/// This allows accepting biotite.structure.BondType enum members directly.
impl<'a, 'py> FromPyObject<'a, 'py> for BondType {
    type Error = PyErr;

    fn extract(ob: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        // Extract the integer value from the Python enum (IntEnum supports __int__)
        let value: u32 = ob.extract()?;
        BondType::try_from(value)
    }
}

/// Convert from Rust's BondType to Python's biotite.structure.BondType IntEnum.
/// This ensures that BondType values returned to Python are proper enum members.
impl<'py> IntoPyObject<'py> for BondType {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> Result<Self::Output, Self::Error> {
        // Import biotite.structure.BondType and instantiate with the int value
        let bond_type_class = py.import("biotite.structure")?.getattr("BondType")?;
        let value = self as u8;
        bond_type_class.call1((value,))
    }
}

/// Internal representation of a single bond.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct Bond {
    /// First atom index (always <= atom2 after normalization)
    pub atom1: u32,
    /// Second atom index (always >= atom1 after normalization)
    pub atom2: u32,
    /// The type of the bond
    pub bond_type: BondType,
}

impl Bond {
    /// Create a new bond with normalized atom indices (smaller index first).
    fn new(atom1: u32, atom2: u32, bond_type: BondType) -> Self {
        let (a1, a2) = if atom1 <= atom2 {
            (atom1, atom2)
        } else {
            (atom2, atom1)
        };
        Bond {
            atom1: a1,
            atom2: a2,
            bond_type,
        }
    }

    /// Convert into an array for Python representation.
    fn as_array(&self) -> [u32; 3] {
        [self.atom1, self.atom2, self.bond_type as u32]
    }
}

/// __init__(atom_count, bonds=None)
///
/// A bond list stores indices of atoms
/// (usually of an :class:`AtomArray` or :class:`AtomArrayStack`)
/// that form chemical bonds together with the type (or order) of the
/// bond.
///
/// Internally the bonds are stored as *n x 3* :class:`ndarray`.
/// For each row, the first column specifies the index of the first
/// atom, the second column the index of the second atom involved in the
/// bond.
/// The third column stores an integer that is interpreted as member
/// of the the :class:`BondType` enum, that specifies the order of the
/// bond.
///
/// When indexing a :class:`BondList`, the index is not forwarded to the
/// internal :class:`ndarray`. Instead the indexing behavior is
/// consistent with indexing an :class:`AtomArray` or
/// :class:`AtomArrayStack`:
/// Bonds with at least one atom index that is not covered by the index
/// are removed, atom indices that occur after an uncovered atom index
/// move up.
/// Effectively, this means that after indexing an :class:`AtomArray`
/// and a :class:`BondList` with the same index, the atom indices in the
/// :class:`BondList` will still point to the same atoms in the
/// :class:`AtomArray`.
/// Indexing a :class:`BondList` with a single integer is equivalent
/// to calling :func:`get_bonds()`.
///
/// The same consistency applies to adding :class:`BondList` instances
/// via the '+' operator:
/// The atom indices of the second :class:`BondList` are increased by
/// the atom count of the first :class:`BondList` and then both
/// :class:`BondList` objects are merged.
///
/// Parameters
/// ----------
/// atom_count : int
///     A positive integer, that specifies the number of atoms the
///     :class:`BondList` refers to
///     (usually the length of an atom array (stack)).
///     Effectively, this value is the exclusive maximum for the indices
///     stored in the :class:`BondList`.
/// bonds : ndarray, shape=(n,2) or shape=(n,3), dtype=int, optional
///     This array contains the indices of atoms which are bonded:
///     For each row, the first column specifies the first atom,
///     the second row the second atom involved in a chemical bond.
///     If an *n x 3* array is provided, the additional column
///     specifies a :class:`BondType` instead of :attr:`BondType.ANY`.
///     By default, the created :class:`BondList` is empty.
///
/// Notes
/// -----
/// When initially providing the bonds as :class:`ndarray`, the input is
/// sanitized: Redundant bonds are removed, and each bond entry is
/// sorted so that the lower one of the two atom indices is in the first
/// column.
/// If a bond appears multiple times with different bond types, the
/// first bond takes precedence.
#[pyclass(module = "biotite.structure", subclass)]
#[derive(Clone)]
pub struct BondList {
    atom_count: u32,
    bonds: Vec<Bond>,
    max_bonds_per_atom: u32,
}

#[pymethods]
impl BondList {
    #[new]
    #[pyo3(signature = (atom_count, bonds=None))]
    fn new<'py>(
        py: Python<'py>,
        atom_count: u32,
        bonds: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Self> {
        match bonds {
            Some(bonds_obj) => {
                // Import numpy and convert to array with int64 dtype
                let np = PyModule::import(py, "numpy")?;
                let array_obj = np.call_method1("asarray", (bonds_obj,))?;

                // Get the shape
                let shape_attr = array_obj.getattr("shape")?;
                let shape: Vec<usize> = shape_attr.extract()?;

                // Validate the shape
                if shape.len() != 2 {
                    return Err(exceptions::PyValueError::new_err(
                        "Expected a 2D-ndarray for input bonds",
                    ));
                }

                let n_bonds = shape[0];
                let n_cols = shape[1];

                if n_cols != 2 && n_cols != 3 {
                    return Err(exceptions::PyValueError::new_err(
                        "Input array containing bonds must be either of shape (n,2) or (n,3)",
                    ));
                }

                // Now convert to int64 for processing
                let array_int64 = np.call_method1("asarray", (array_obj, np.getattr("int64")?))?;
                let bonds_array: PyReadonlyArray2<i64> = array_int64.extract()?;
                let array = bonds_array.as_array();

                let mut bond_vec: Vec<Bond> = Vec::with_capacity(n_bonds);
                for i in 0..n_bonds {
                    let idx1 = to_positive_index(array[[i, 0]], atom_count)?;
                    let idx2 = to_positive_index(array[[i, 1]], atom_count)?;
                    let bond_type = if n_cols == 3 {
                        let bt_val = array[[i, 2]];
                        if bt_val < 0 {
                            return Err(exceptions::PyValueError::new_err(format!(
                                "BondType {} is invalid",
                                bt_val
                            )));
                        }
                        BondType::try_from(bt_val as u32)?
                    } else {
                        BondType::Any
                    };
                    bond_vec.push(Bond::new(idx1, idx2, bond_type));
                }

                let mut bond_list = BondList {
                    atom_count,
                    bonds: bond_vec,
                    max_bonds_per_atom: 0,
                };
                bond_list.remove_redundant_bonds();
                bond_list.max_bonds_per_atom = bond_list.calculate_max_bonds_per_atom();
                Ok(bond_list)
            }
            None => Ok(BondList {
                atom_count,
                bonds: Vec::new(),
                max_bonds_per_atom: 0,
            }),
        }
    }

    /// Concatenate multiple :class:`BondList` objects into a single
    /// :class:`BondList`, respectively.
    ///
    /// Parameters
    /// ----------
    /// bonds_lists : iterable object of BondList
    ///     The bond lists to be concatenated.
    ///
    /// Returns
    /// -------
    /// concatenated_bonds : BondList
    ///     The concatenated bond lists.
    #[staticmethod]
    fn concatenate(bond_lists: Vec<BondList>) -> PyResult<BondList> {
        if bond_lists.is_empty() {
            return Ok(BondList {
                atom_count: 0,
                bonds: Vec::new(),
                max_bonds_per_atom: 0,
            });
        }

        let total_bonds: usize = bond_lists.iter().map(|bl| bl.bonds.len()).sum();
        let mut merged_bonds: Vec<Bond> = Vec::with_capacity(total_bonds);
        let mut cum_atom_count: u32 = 0;
        let mut max_bonds: u32 = 0;

        for bond_list in &bond_lists {
            for bond in &bond_list.bonds {
                merged_bonds.push(Bond {
                    atom1: bond.atom1 + cum_atom_count,
                    atom2: bond.atom2 + cum_atom_count,
                    bond_type: bond.bond_type,
                });
            }
            cum_atom_count += bond_list.atom_count;
            if bond_list.max_bonds_per_atom > max_bonds {
                max_bonds = bond_list.max_bonds_per_atom;
            }
        }

        Ok(BondList {
            atom_count: cum_atom_count,
            bonds: merged_bonds,
            max_bonds_per_atom: max_bonds,
        })
    }

    /// offset_indices(offset)
    ///
    /// Increase all atom indices in the :class:`BondList` by the given
    /// offset.
    ///
    /// Implicitly this increases the atom count.
    ///
    /// Parameters
    /// ----------
    /// offset : int
    ///     The atom indices are increased by this value.
    ///     Must be positive.
    fn offset_indices(&mut self, offset: i32) -> PyResult<()> {
        if offset < 0 {
            return Err(exceptions::PyValueError::new_err("Offset must be positive"));
        }
        let offset = offset as u32;
        for bond in &mut self.bonds {
            bond.atom1 += offset;
            bond.atom2 += offset;
        }
        self.atom_count += offset;
        Ok(())
    }

    /// as_array()
    ///
    /// Obtain a copy of the internal :class:`ndarray`.
    ///
    /// Returns
    /// -------
    /// array : ndarray, shape=(n,3), dtype=np.uint32
    ///     Copy of the internal :class:`ndarray`.
    ///     For each row, the first column specifies the index of the
    ///     first atom, the second column the index of the second atom
    ///     involved in the bond.
    ///     The third column stores the :class:`BondType`.
    fn as_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<u32>> {
        let n_bonds = self.bonds.len();
        let mut data: Vec<u32> = Vec::with_capacity(n_bonds * 3);
        for bond in &self.bonds {
            data.push(bond.atom1);
            data.push(bond.atom2);
            data.push(bond.bond_type as u32);
        }
        Array2::from_shape_vec((n_bonds, 3), data)
            .expect("Shape mismatch")
            .into_pyarray(py)
    }

    /// as_set()
    ///
    /// Obtain a set representation of the :class:`BondList`.
    ///
    /// Returns
    /// -------
    /// bond_set : set of tuple(int, int, int)
    ///     A set of tuples.
    ///     Each tuple represents one bond:
    ///     The first integer represents the first atom,
    ///     the second integer represents the second atom,
    ///     the third integer represents the :class:`BondType`.
    fn as_set<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PySet>> {
        let set = PySet::empty(py)?;
        for bond in &self.bonds {
            let tuple = PyTuple::new(py, [bond.atom1, bond.atom2, bond.bond_type as u32])?;
            set.add(tuple)?;
        }
        Ok(set)
    }

    /// as_graph()
    ///
    /// Obtain a graph representation of the :class:`BondList`.
    ///
    /// Returns
    /// -------
    /// bond_set : Graph
    ///     A *NetworkX* :class:`Graph`.
    ///     The atom indices are nodes, the bonds are edges.
    ///     Each edge has a ``"bond_type"`` attribute containing the
    ///     :class:`BondType`.
    ///
    /// Examples
    /// --------
    ///
    /// >>> bond_list = BondList(5, np.array([(1,0,2), (1,3,1), (1,4,1)]))
    /// >>> graph = bond_list.as_graph()
    /// >>> print(graph.nodes)
    /// [0, 1, 3, 4]
    /// >>> print(graph.edges)
    /// [(0, 1), (1, 3), (1, 4)]
    /// >>> for i, j in graph.edges:
    /// ...     print(i, j, graph.get_edge_data(i, j))
    /// 0 1 {'bond_type': <BondType.DOUBLE: 2>}
    /// 1 3 {'bond_type': <BondType.SINGLE: 1>}
    /// 1 4 {'bond_type': <BondType.SINGLE: 1>}
    fn as_graph<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let nx = PyModule::import(py, "networkx")?;
        let graph = nx.call_method0("Graph")?;

        // Build edges list with bond_type attribute
        let edges = pyo3::types::PyList::empty(py);
        for bond in &self.bonds {
            let edge_data = pyo3::types::PyDict::new(py);
            edge_data.set_item("bond_type", bond.bond_type)?;
            let edge_tuple = PyTuple::new(
                py,
                [
                    bond.atom1.into_pyobject(py)?.into_any(),
                    bond.atom2.into_pyobject(py)?.into_any(),
                    edge_data.into_any(),
                ],
            )?;
            edges.append(edge_tuple)?;
        }

        graph.call_method1("add_edges_from", (edges,))?;
        Ok(graph)
    }

    /// Remove aromaticity from the bond types.
    ///
    /// :attr:`BondType.AROMATIC_{ORDER}` is converted into
    /// :attr:`BondType.{ORDER}`.
    fn remove_aromaticity(&mut self) {
        for bond in &mut self.bonds {
            bond.bond_type = bond.bond_type.without_aromaticity();
        }
    }

    /// Remove the bond order information from aromatic bonds, i.e. convert all
    /// aromatic bonds to :attr:`BondType.ANY`.
    fn remove_kekulization(&mut self) {
        for bond in &mut self.bonds {
            if matches!(
                bond.bond_type,
                BondType::AromaticSingle | BondType::AromaticDouble | BondType::AromaticTriple
            ) {
                bond.bond_type = BondType::Aromatic;
            }
        }
    }

    /// Convert all bonds to :attr:`BondType.ANY`.
    fn remove_bond_order(&mut self) {
        for bond in &mut self.bonds {
            bond.bond_type = BondType::Any;
        }
    }

    /// convert_bond_type(original_bond_type, new_bond_type)
    ///
    /// Convert all occurences of a given bond type into another bond type.
    ///
    /// Parameters
    /// ----------
    /// original_bond_type : BondType
    ///     The bond type to convert.
    /// new_bond_type : BondType
    ///     The new bond type.
    fn convert_bond_type(&mut self, original_bond_type: BondType, new_bond_type: BondType) {
        for bond in &mut self.bonds {
            if bond.bond_type == original_bond_type {
                bond.bond_type = new_bond_type;
            }
        }
    }

    /// get_atom_count()
    ///
    /// Get the atom count.
    ///
    /// Returns
    /// -------
    /// atom_count : int
    ///     The atom count.
    fn get_atom_count(&self) -> u32 {
        self.atom_count
    }

    /// Private attribute for compatibility with tests.
    #[getter]
    fn _atom_count(&self) -> u32 {
        self.atom_count
    }

    /// Private attribute for compatibility with tests.
    #[getter]
    fn _max_bonds_per_atom(&self) -> u32 {
        self.max_bonds_per_atom
    }

    /// get_bond_count()
    ///
    /// Get the amount of bonds.
    ///
    /// Returns
    /// -------
    /// bond_count : int
    ///     The amount of bonds. This is equal to the length of the
    ///     internal :class:`ndarray` containing the bonds.
    fn get_bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// get_bonds(atom_index)
    ///
    /// Obtain the indices of the atoms bonded to the atom with the
    /// given index as well as the corresponding bond types.
    ///
    /// Parameters
    /// ----------
    /// atom_index : int
    ///     The index of the atom to get the bonds for.
    ///
    /// Returns
    /// -------
    /// bonds : np.ndarray, dtype=np.uint32, shape=(k,)
    ///     The indices of connected atoms.
    /// bond_types : np.ndarray, dtype=np.uint8, shape=(k,)
    ///     Array of integers, interpreted as :class:`BondType`
    ///     instances.
    ///     This array specifies the type (or order) of the bonds to
    ///     the connected atoms.
    #[allow(clippy::type_complexity)]
    fn get_bonds<'py>(
        &self,
        py: Python<'py>,
        atom_index: i64,
    ) -> PyResult<(Bound<'py, PyArray1<u32>>, Bound<'py, PyArray1<u8>>)> {
        let index = to_positive_index(atom_index, self.atom_count)?;

        let mut bonded_atoms: Vec<u32> = Vec::with_capacity(self.max_bonds_per_atom as usize);
        let mut bond_types: Vec<u8> = Vec::with_capacity(self.max_bonds_per_atom as usize);

        for bond in &self.bonds {
            if bond.atom1 == index {
                bonded_atoms.push(bond.atom2);
                bond_types.push(bond.bond_type as u8);
            } else if bond.atom2 == index {
                bonded_atoms.push(bond.atom1);
                bond_types.push(bond.bond_type as u8);
            }
        }

        Ok((
            Array1::from_vec(bonded_atoms).into_pyarray(py),
            Array1::from_vec(bond_types).into_pyarray(py),
        ))
    }

    /// get_all_bonds()
    ///
    /// For each atom index, give the indices of the atoms bonded to
    /// this atom as well as the corresponding bond types.
    ///
    /// Returns
    /// -------
    /// bonds : np.ndarray, dtype=np.int32, shape=(n,k)
    ///     The indices of connected atoms.
    ///     The first dimension represents the atoms,
    ///     the second dimension represents the indices of atoms bonded
    ///     to the respective atom.
    ///     Atoms can have have different numbers of atoms bonded to
    ///     them.
    ///     Therefore, the length of the second dimension *k* is equal
    ///     to the maximum number of bonds for an atom in this
    ///     :class:`BondList`.
    ///     For atoms with less bonds, the corresponding entry in the
    ///     array is padded with ``-1`` values.
    /// bond_types : np.ndarray, dtype=np.int8, shape=(n,k)
    ///     Array of integers, interpreted as :class:`BondType`
    ///     instances.
    ///     This array specifies the bond type (or order) corresponding
    ///     to the returned `bonds`.
    ///     It uses the same ``-1``-padding.
    fn get_all_bonds<'py>(
        &self,
        py: Python<'py>,
    ) -> (Bound<'py, PyArray2<i32>>, Bound<'py, PyArray2<i8>>) {
        let n_atoms = self.atom_count as usize;
        let k = self.max_bonds_per_atom as usize;

        let mut bonds_matrix = Array2::from_elem((n_atoms, k), -1i32);
        let mut types_matrix = Array2::from_elem((n_atoms, k), -1i8);
        let mut lengths: Vec<usize> = vec![0; n_atoms];

        for bond in &self.bonds {
            let i1 = bond.atom1 as usize;
            let i2 = bond.atom2 as usize;
            let bt = bond.bond_type as i8;

            bonds_matrix[[i1, lengths[i1]]] = bond.atom2 as i32;
            types_matrix[[i1, lengths[i1]]] = bt;
            lengths[i1] += 1;

            bonds_matrix[[i2, lengths[i2]]] = bond.atom1 as i32;
            types_matrix[[i2, lengths[i2]]] = bt;
            lengths[i2] += 1;
        }

        (bonds_matrix.into_pyarray(py), types_matrix.into_pyarray(py))
    }

    /// adjacency_matrix()
    ///
    /// Represent this :class:`BondList` as adjacency matrix.
    ///
    /// Returns
    /// -------
    /// matrix : ndarray, dtype=bool, shape=(n,n)
    ///     The created adjacency matrix.
    fn adjacency_matrix<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<bool>> {
        let n = self.atom_count as usize;
        let mut matrix = Array2::from_elem((n, n), false);
        for bond in &self.bonds {
            let i = bond.atom1 as usize;
            let j = bond.atom2 as usize;
            matrix[[i, j]] = true;
            matrix[[j, i]] = true;
        }
        matrix.into_pyarray(py)
    }

    /// bond_type_matrix()
    ///
    /// Represent this :class:`BondList` as a matrix depicting the bond
    /// type.
    ///
    /// Returns
    /// -------
    /// matrix : ndarray, dtype=np.int8, shape=(n,n)
    ///     The created bond type matrix.
    fn bond_type_matrix<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<i8>> {
        let n = self.atom_count as usize;
        let mut matrix = Array2::from_elem((n, n), -1i8);
        for bond in &self.bonds {
            let i = bond.atom1 as usize;
            let j = bond.atom2 as usize;
            let bt = bond.bond_type as i8;
            matrix[[i, j]] = bt;
            matrix[[j, i]] = bt;
        }
        matrix.into_pyarray(py)
    }

    /// add_bond(atom_index1, atom_index2, bond_type=BondType.ANY)
    ///
    /// Add a bond to the :class:`BondList`.
    ///
    /// If the bond is already existent, only the bond type is updated.
    ///
    /// Parameters
    /// ----------
    /// atom_index1, atom_index2 : int
    ///     The indices of the atoms to create a bond for.
    /// bond_type : BondType
    ///     The type of the bond. Default is :attr:`BondType.ANY`.
    #[pyo3(signature = (atom_index1, atom_index2, bond_type=BondType::Any))]
    fn add_bond(
        &mut self,
        atom_index1: i64,
        atom_index2: i64,
        bond_type: BondType,
    ) -> PyResult<()> {
        let idx1 = to_positive_index(atom_index1, self.atom_count)?;
        let idx2 = to_positive_index(atom_index2, self.atom_count)?;

        let new_bond = Bond::new(idx1, idx2, bond_type);

        // Check if bond already exists
        for bond in &mut self.bonds {
            if bond.atom1 == new_bond.atom1 && bond.atom2 == new_bond.atom2 {
                // Update bond type
                bond.bond_type = bond_type;
                return Ok(());
            }
        }

        // Bond doesn't exist, add it
        self.bonds.push(new_bond);
        self.max_bonds_per_atom = self.calculate_max_bonds_per_atom();
        Ok(())
    }

    /// remove_bond(atom_index1, atom_index2)
    ///
    /// Remove a bond from the :class:`BondList`.
    ///
    /// If the bond is not existent in the :class:`BondList`, nothing happens.
    ///
    /// Parameters
    /// ----------
    /// atom_index1, atom_index2 : int
    ///     The indices of the atoms whose bond should be removed.
    fn remove_bond(&mut self, atom_index1: i64, atom_index2: i64) -> PyResult<()> {
        let idx1 = to_positive_index(atom_index1, self.atom_count)?;
        let idx2 = to_positive_index(atom_index2, self.atom_count)?;
        let (min_idx, max_idx) = if idx1 <= idx2 {
            (idx1, idx2)
        } else {
            (idx2, idx1)
        };

        self.bonds
            .retain(|bond| !(bond.atom1 == min_idx && bond.atom2 == max_idx));
        // max_bonds_per_atom not recalculated (see Cython implementation comment)
        Ok(())
    }

    /// remove_bonds_to(self, atom_index)
    ///
    /// Remove all bonds from the :class:`BondList` where the given atom
    /// is involved.
    ///
    /// Parameters
    /// ----------
    /// atom_index : int
    ///     The index of the atom whose bonds should be removed.
    fn remove_bonds_to(&mut self, atom_index: i64) -> PyResult<()> {
        let index = to_positive_index(atom_index, self.atom_count)?;
        self.bonds
            .retain(|bond| bond.atom1 != index && bond.atom2 != index);
        Ok(())
    }

    /// remove_bonds(bond_list)
    ///
    /// Remove multiple bonds from the :class:`BondList`.
    ///
    /// All bonds present in `bond_list` are removed from this instance.
    /// If a bond is not existent in this instance, nothing happens.
    /// Only the bond indices, not the bond types, are relevant for
    /// this.
    ///
    /// Parameters
    /// ----------
    /// bond_list : BondList
    ///     The bonds in `bond_list` are removed from this instance.
    fn remove_bonds(&mut self, bond_list: &BondList) {
        let bonds_to_remove: HashSet<(u32, u32)> =
            bond_list.bonds.iter().map(|b| (b.atom1, b.atom2)).collect();
        self.bonds
            .retain(|bond| !bonds_to_remove.contains(&(bond.atom1, bond.atom2)));
    }

    /// merge(bond_list)
    ///
    /// Merge another :class:`BondList` with this instance into a new
    /// object.
    /// If a bond appears in both :class:`BondList`'s, the
    /// :class:`BondType` from the given `bond_list` takes precedence.
    ///
    /// The internal :class:`ndarray` instances containing the bonds are
    /// simply concatenated and the new atom count is the maximum of
    /// both bond lists.
    ///
    /// Parameters
    /// ----------
    /// bond_list : BondList
    ///     This bond list is merged with this instance.
    ///
    /// Returns
    /// -------
    /// bond_list : BondList
    ///     The merged :class:`BondList`.
    ///
    /// Notes
    /// -----
    /// This is not equal to using the `+` operator.
    fn merge(&self, bond_list: &BondList) -> BondList {
        let new_atom_count = self.atom_count.max(bond_list.atom_count);

        // Concatenate bonds: other bond_list takes precedence (put first)
        let total_bonds = self.bonds.len() + bond_list.bonds.len();
        let mut merged_bonds: Vec<Bond> = Vec::with_capacity(total_bonds);
        merged_bonds.extend(bond_list.bonds.iter().cloned());
        merged_bonds.extend(self.bonds.iter().cloned());

        let mut result = BondList {
            atom_count: new_atom_count,
            bonds: merged_bonds,
            max_bonds_per_atom: 0,
        };
        result.remove_redundant_bonds();
        result.max_bonds_per_atom = result.calculate_max_bonds_per_atom();
        result
    }

    fn __add__(&self, other: &BondList) -> BondList {
        BondList::concatenate(vec![self.clone(), other.clone()]).unwrap()
    }

    #[allow(deprecated)]
    fn __getitem__<'py>(&self, py: Python<'py>, index: &Bound<'py, PyAny>) -> PyResult<Py<PyAny>> {
        let np = PyModule::import(py, "numpy")?;

        // Handle single integer index
        if let Ok(int_index) = index.extract::<i64>() {
            let (bonds, types) = self.get_bonds(py, int_index)?;
            let tuple = PyTuple::new(py, [bonds.into_any(), types.into_any()])?;
            return Ok(tuple.into());
        }

        // Handle slices
        if index.is_instance_of::<pyo3::types::PySlice>() {
            // Use numpy's arange + slice to get the indices
            let arange = np.call_method1("arange", (self.atom_count,))?;
            let sliced_indices = arange.get_item(index)?;
            let index_array: Vec<i64> = sliced_indices.extract()?;
            return self.index_with_int_array(py, &index_array);
        }

        // Handle boolean mask
        if let Ok(bool_array) = index.extract::<Vec<bool>>() {
            if bool_array.len() != self.atom_count as usize {
                return Err(exceptions::PyIndexError::new_err(format!(
                    "Boolean index has length {}, expected {}",
                    bool_array.len(),
                    self.atom_count
                )));
            }
            let result = self.index_with_bool_mask(&bool_array)?;
            return Ok(result.into_pyobject(py)?.into_any().unbind());
        }

        // Try to convert to boolean array via numpy
        let result = np.call_method1("asarray", (index,))?;
        if let Ok(dtype) = result.getattr("dtype") {
            if let Ok(kind) = dtype.getattr("kind") {
                if let Ok(kind_str) = kind.extract::<String>() {
                    if kind_str == "b" {
                        // Boolean array
                        let bool_vec: Vec<bool> = result.extract()?;
                        if bool_vec.len() != self.atom_count as usize {
                            return Err(exceptions::PyIndexError::new_err(format!(
                                "Boolean index has length {}, expected {}",
                                bool_vec.len(),
                                self.atom_count
                            )));
                        }
                        let result = self.index_with_bool_mask(&bool_vec)?;
                        return Ok(result.into_pyobject(py)?.into_any().unbind());
                    }
                }
            }
        }

        // Handle integer array index
        let index_array: Vec<i64> = result.extract()?;
        self.index_with_int_array(py, &index_array)
    }

    fn __contains__(&self, item: &Bound<'_, PyTuple>) -> PyResult<bool> {
        if item.len() != 2 {
            return Err(exceptions::PyTypeError::new_err(
                "Expected a tuple of atom indices",
            ));
        }
        let idx1: u32 = item.get_item(0)?.extract()?;
        let idx2: u32 = item.get_item(1)?.extract()?;
        let (min_idx, max_idx) = if idx1 <= idx2 {
            (idx1, idx2)
        } else {
            (idx2, idx1)
        };

        for bond in &self.bonds {
            if bond.atom1 == min_idx && bond.atom2 == max_idx {
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn __eq__(&self, other: &BondList) -> bool {
        if self.atom_count != other.atom_count {
            return false;
        }
        let self_set: HashSet<_> = self.bonds.iter().map(|b| b.as_array()).collect();
        let other_set: HashSet<_> = other.bonds.iter().map(|b| b.as_array()).collect();
        self_set == other_set
    }

    fn __str__(&self, py: Python<'_>) -> PyResult<String> {
        let array = self.as_array(py);
        array.call_method0("__str__")?.extract()
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.__str__(py)
    }

    fn __iter__(&self) -> PyResult<()> {
        Err(exceptions::PyTypeError::new_err(
            "'BondList' object is not iterable",
        ))
    }

    fn __len__(&self) -> usize {
        self.bonds.len()
    }

    // Pickling support
    fn __reduce__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let struc = PyModule::import(py, "biotite.structure")?;
        let cls = struc.getattr("BondList")?;
        let from_state = cls.getattr("_from_state")?;

        // Serialize bonds as flat vector of [atom1, atom2, type] arrays
        let bonds_data: Vec<[u32; 3]> = self.bonds.iter().map(|b| b.as_array()).collect();

        let args = PyTuple::new(
            py,
            [
                self.atom_count.into_pyobject(py)?.into_any().unbind(),
                bonds_data.into_pyobject(py)?.into_any().unbind(),
                self.max_bonds_per_atom
                    .into_pyobject(py)?
                    .into_any()
                    .unbind(),
            ],
        )?;

        PyTuple::new(py, [from_state.unbind(), args.into_any().unbind()])
    }

    #[staticmethod]
    #[pyo3(name = "_from_state")]
    fn from_state(
        atom_count: u32,
        bonds_data: Vec<[u32; 3]>,
        max_bonds_per_atom: u32,
    ) -> PyResult<Self> {
        let bonds: Vec<Bond> = bonds_data
            .into_iter()
            .map(|[a1, a2, bt]| {
                Ok(Bond {
                    atom1: a1,
                    atom2: a2,
                    bond_type: BondType::try_from(bt)?,
                })
            })
            .collect::<PyResult<Vec<Bond>>>()?;

        Ok(BondList {
            atom_count,
            bonds,
            max_bonds_per_atom,
        })
    }

    // Copy support for biotite.Copyable interface
    fn copy(&self) -> BondList {
        self.clone()
    }

    fn __copy__(&self) -> BondList {
        self.clone()
    }

    fn __deepcopy__(&self, _memo: &Bound<'_, PyAny>) -> BondList {
        self.clone()
    }
}

impl BondList {
    /// Calculate the maximum number of bonds any single atom has.
    fn calculate_max_bonds_per_atom(&self) -> u32 {
        if self.atom_count == 0 {
            return 0;
        }

        let mut counts: Vec<u32> = vec![0; self.atom_count as usize];
        for bond in &self.bonds {
            counts[bond.atom1 as usize] += 1;
            counts[bond.atom2 as usize] += 1;
        }
        *counts.iter().max().unwrap_or(&0)
    }

    /// Remove redundant bonds, keeping the first occurrence.
    fn remove_redundant_bonds(&mut self) {
        let mut seen: HashSet<(u32, u32)> = HashSet::new();
        self.bonds.retain(|bond| {
            let key = (bond.atom1, bond.atom2);
            if seen.contains(&key) {
                false
            } else {
                seen.insert(key);
                true
            }
        });
    }

    /// Index the bond list with a boolean mask.
    fn index_with_bool_mask(&self, mask: &[bool]) -> PyResult<BondList> {
        // Calculate offsets: how much each index should be decreased
        let mut offsets: Vec<u32> = Vec::with_capacity(mask.len());
        let mut cumulative_offset: u32 = 0;
        for &is_selected in mask {
            if !is_selected {
                cumulative_offset += 1;
            }
            offsets.push(cumulative_offset);
        }

        let new_atom_count = mask.iter().filter(|&&x| x).count() as u32;

        let mut new_bonds: Vec<Bond> = Vec::new();
        for bond in &self.bonds {
            let i1 = bond.atom1 as usize;
            let i2 = bond.atom2 as usize;
            if mask[i1] && mask[i2] {
                let new_atom1 = bond.atom1 - offsets[i1];
                let new_atom2 = bond.atom2 - offsets[i2];
                new_bonds.push(Bond {
                    atom1: new_atom1,
                    atom2: new_atom2,
                    bond_type: bond.bond_type,
                });
            }
        }

        let mut result = BondList {
            atom_count: new_atom_count,
            bonds: new_bonds,
            max_bonds_per_atom: 0,
        };
        result.max_bonds_per_atom = result.calculate_max_bonds_per_atom();
        Ok(result)
    }

    /// Index the bond list with an integer array.
    fn index_with_int_array(&self, py: Python<'_>, indices: &[i64]) -> PyResult<Py<PyAny>> {
        // Build inverse index: for each original atom, what is its new index?
        // -1 means the atom is not in the new array
        let mut inverse_index: Vec<i32> = vec![-1; self.atom_count as usize];
        for (new_idx, &orig_idx) in indices.iter().enumerate() {
            let pos_idx = to_positive_index(orig_idx, self.atom_count)? as usize;
            if inverse_index[pos_idx] != -1 {
                return Err(exceptions::PyNotImplementedError::new_err(format!(
                    "Duplicate indices are not supported, but index {} appeared multiple times",
                    pos_idx
                )));
            }
            inverse_index[pos_idx] = new_idx as i32;
        }

        let new_atom_count = indices.len() as u32;
        let mut new_bonds: Vec<Bond> = Vec::new();

        for bond in &self.bonds {
            let new_idx1 = inverse_index[bond.atom1 as usize];
            let new_idx2 = inverse_index[bond.atom2 as usize];
            if new_idx1 != -1 && new_idx2 != -1 {
                new_bonds.push(Bond::new(new_idx1 as u32, new_idx2 as u32, bond.bond_type));
            }
        }

        let mut result = BondList {
            atom_count: new_atom_count,
            bonds: new_bonds,
            max_bonds_per_atom: 0,
        };
        result.max_bonds_per_atom = result.calculate_max_bonds_per_atom();
        Ok(result.into_pyobject(py)?.into_any().unbind())
    }
}

/// Convert a potentially negative index to a positive index.
fn to_positive_index(index: i64, array_length: u32) -> PyResult<u32> {
    let length = array_length as i64;
    if index < 0 {
        let pos_index = length + index;
        if pos_index < 0 {
            return Err(exceptions::PyIndexError::new_err(format!(
                "Index {} is out of range for an atom count of {}",
                index, array_length
            )));
        }
        Ok(pos_index as u32)
    } else {
        if index >= length {
            return Err(exceptions::PyIndexError::new_err(format!(
                "Index {} is out of range for an atom count of {}",
                index, array_length
            )));
        }
        Ok(index as u32)
    }
}
