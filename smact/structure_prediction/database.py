"""Tools for database interfacing for high throughput IO."""

import itertools
from multiprocessing import Pool
from operator import itemgetter

try:
    from pathos.pools import ParallelPool

    pathos_available = True
except ImportError:
    pathos_available = False

import sqlite3
from typing import Dict, Generator, List, Optional, Sequence, Tuple, Union

import pymatgen
from pymatgen.ext.matproj import MPRester

from . import logger
from .structure import SmactStructure
from .utilities import convert_next_gen_mprest_data, get_sign


class StructureDB:
    """SQLite Structure Database interface.

    Acts as a context manager for database interfacing
    and wraps several useful SQLite commands within
    methods.

    Attributes:
        db: The database name.
        conn: The database connection. Only open when
            used as a context manager.
        cur: The database connection cursor. Only usable
            when class implemented as context manager.

    Examples:
        Connecting to a database in memory:

        >>> DB = StructureDB(':memory:')
        >>> with DB as c:
        ...     _ = c.execute("CREATE TABLE test (id, val)")
        ...     c.execute("SELECT * FROM test").fetchall()
        []
        >>> DB.cur.execute("SELECT * FROM test").fetchall()
        Traceback (most recent call last):
            ...
        sqlite3.ProgrammingError: Cannot operate on a closed database.

    """

    def __init__(self, db: str):
        """Set database name.

        Args:
            db (str): The name of the database. Can also be ':memory:'
                to connect to a database in RAM.

        """
        self.db = db

    def __enter__(self) -> sqlite3.Cursor:
        """Initialize database connection.

        Returns:
            An SQLite cursor for interfacing with the database.

        """
        self.conn = sqlite3.connect(self.db)
        self.cur = self.conn.cursor()

        return self.cur

    def __exit__(self, exc_type, *args):
        """Close database connection.

        Commits all changes before closing.
        Alternatively, rolls back any changes if an exception
        was raised, causing the context to be exited.

        """
        if exc_type is not None:
            self.conn.rollback()
        else:
            self.conn.commit()

        self.conn.close()

    def add_mp_icsd(
        self,
        table: str,
        mp_data: Optional[
            List[Dict[str, Union[pymatgen.core.Structure, str]]]
        ] = None,
        mp_api_key: Optional[str] = None,
    ) -> int:
        """Add a table populated with Materials Project-hosted ICSD structures.

        Note:
            This is very computationally expensive for large datasets
            and will not likely run on a laptop.
            If possible, download a pre-constructed database.

        Args:
            table (str): The name of the table to add.
            mp_data: The Materials Project data to parse. If this is None, data
                will be downloaded. Downloading data needs `mp_api_key` to be set.
            mp_api_key (str): A Materials Project API key. Only needed if `mp_data`
                is None.

        Returns:
            The number of structs added.

        """
        if mp_data is None:  # pragma: no cover
            with MPRester(mp_api_key) as m:
                try:
                    data = m.query(
                        criteria={"icsd_ids.0": {"$exists": True}},
                        properties=["structure", "material_id"],
                    )
                except NotImplementedError:
                    docs = m.summary.search(
                        theoretical=False, fields=["structure", "material_id"]
                    )
                    data = [convert_next_gen_mprest_data(doc) for doc in docs]
        else:
            data = mp_data

        self.add_table(table)

        if pathos_available:
            pool = ParallelPool()
            parse_iter = pool.uimap(parse_mprest, data)
        else:  # pragma: no cover
            parse_iter = map(parse_mprest, data)

        return self.add_structs(parse_iter, table, commit_after_each=True)

    def add_table(self, table: str):
        """Add a table to the database.

        Args:
            table: The name of the table to add

        """
        with self as c:
            c.execute(
                f"""CREATE TABLE {table}
                (composition TEXT NOT NULL, structure TEXT NOT NULL)""",
            )

    def add_struct(self, struct: SmactStructure, table: str):
        """Add a SmactStructure to a table.

        Args:
            struct: The :class:`~.SmactStructure` to add.
            table: The name of the table to add the structure to.

        """
        entry = (struct.composition(), struct.as_poscar())

        with self as c:
            c.execute(f"INSERT into {table} VALUES (?, ?)", entry)

    def add_structs(
        self,
        structs: Sequence[SmactStructure],
        table: str,
        commit_after_each: Optional[bool] = False,
    ) -> int:
        """Add several SmactStructures to a table.

        Args:
            structs: Iterable of :class:`~.SmactStructure` s to add to table.
            table: The name of the table to add the structs to.
            commit_after_each (bool, optional): Whether to commit the addition
                after each structure is added.
                This is useful when adding a large number of structures over
                a long timeframe, as it ensures some structures are added,
                even if the program terminates before completion.
                Defaults to False.

        Returns:
            The number of structures added.

        """
        with self as c:
            num = 0
            for struct in structs:
                if struct is None:
                    # Handling for poorly decorated structures
                    continue
                entry = (struct.composition(), struct.as_poscar())
                c.execute(f"INSERT into {table} VALUES (?, ?)", entry)
                num += 1

                if commit_after_each:
                    self.conn.commit()

        return num

    def get_structs(
        self, composition: str, table: str
    ) -> List[SmactStructure]:
        """Get SmactStructures for a given composition.

        Args:
            composition: The composition to search for.
                See :meth:`SmactStructure.composition`.
            table: The name of the table in which to search.

        Returns:
            A list of :class:`~.SmactStructure` s.

        """
        with self as c:
            c.execute(
                f"SELECT structure FROM {table} WHERE composition = ?",
                (composition,),
            )
            structs = c.fetchall()
        return [SmactStructure.from_poscar(pos[0]) for pos in structs]

    def get_with_species(
        self,
        species: List[Tuple[str, int]],
        table: str,
    ) -> List[SmactStructure]:
        """Get SmactStructures containing given species.

        Args:
            species: A list of species as tuples, in (element, charge) format.
            table: The name of the table from which to get the species.

        Returns:
            A list of :class:`SmactStructure` s in the table that contain the species.

        """
        glob = "*".join("{}_*_{}{}" for _ in range(len(species)))
        glob = f"*{glob}*"

        species.sort(key=itemgetter(1), reverse=True)
        species.sort(key=itemgetter(0))

        # Generate a list of [element1, charge1, sign1, element2, ...]
        vals = list(
            itertools.chain.from_iterable(
                [x[0], abs(x[1]), get_sign(x[1])] for x in species
            )
        )

        glob_form = glob.format(*vals)

        with self as c:
            c.execute(
                f"SELECT structure FROM {table} WHERE composition GLOB ?",
                (glob_form,),
            )
            structs = c.fetchall()

        return [SmactStructure.from_poscar(pos[0]) for pos in structs]


def parse_mprest(
    data: Dict[str, Union[pymatgen.core.Structure, str]],
    determine_oxi: str = "BV",
) -> SmactStructure:
    """Parse MPRester query data to generate structures.

    Args:
        data: A dictionary containing the keys 'structure' and
            'material_id', with the associated values.
        determine_oxi (str): The method to determine the assignments oxidation states in the structure.
                Options are 'BV', 'comp_ICSD','both' for determining the oxidation states by bond valence,
                ICSD statistics or trial both sequentially, respectively.

    Returns:
        An oxidation-state-decorated :class:`SmactStructure`.

    """
    # Convert next gen query data to a dic
    # TODO check if the data is the same type as MPDataDoc
    if not isinstance(data, dict):
        data = convert_next_gen_mprest_data(data)

    try:
        return SmactStructure.from_py_struct(
            data["structure"], determine_oxi="BV"
        )
    except:
        # Couldn't decorate with oxidation states
        logger.warn(
            f"Couldn't decorate {data['material_id']} with oxidation states."
        )
