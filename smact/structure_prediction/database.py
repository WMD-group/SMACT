"""Tools for database interfacing for high throughput IO."""

import itertools
from multiprocessing import Pool
from operator import itemgetter
from typing import Generator, Sequence, List, Tuple, Dict, Union, Optional

import pymatgen
from pymatgen.ext.matproj import MPRester

from . import logger
from .structure import SmactStructure
from .utilities import get_sign
import sqlite3


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
        """Set database name."""
        self.db = db

    def __enter__(self) -> sqlite3.Cursor:
        """Initialize database connection."""
        self.conn = sqlite3.connect(self.db)
        self.cur = self.conn.cursor()

        return self.cur

    def __exit__(self, *args):
        """Close database connection."""
        self.conn.commit()
        self.conn.close()

    def add_mp_icsd(
      self,
      table: str,
      mp_data: Optional[List[Dict[str, Union[pymatgen.Structure, str]]]] = None,
      mp_api_key: Optional[str] = None,
      threads: Optional[int] = None
    ) -> int:
        """Add a table populated with Materials Project-hosted ICSD structures.

        Args:
            table (str): The name of the table to add.
            mp_data: The Materials Project data to parse. If this is None, data
                will be downloaded. Downloading data needs `mp_api_key` to be set.
            mp_api_key (str): A Materials Project API key. Only needed if `mp_data`
                is None.
            threads (int): The number of processing threads to use for `mp.Pool`.
        
        Returns:
            The number of structs added.

        """
        if mp_data is None:
            with MPRester(mp_api_key) as m:
                data = m.query(
                  criteria={
                    'icsd_ids': {
                      '$exists': True
                    }, },
                  properties=['structure', 'material_id'],
                )
        else:
            data = mp_data

        def parse_mprest(
          data: Dict[str, Union[pymatgen.Structure, str]],
        ) -> Generator[SmactStructure, None, None]:
            """Parse MPRester query data to generate structures."""
            try:
                yield SmactStructure.from_py_struct(data["structure"])
            except:
                # Couldn't decorate with oxidation states
                logger.warn(f"Couldn't decorate {data['material_id']} with oxidation states.")

        with Pool(threads) as p:
            return self.add_structs(p.imap_unordered(parse_mprest, data, 500), table)

    def add_table(self, table: str) -> bool:
        """Add a table to the database.

        Args:
            table: The name of the table to add

        Returns:
            bool: Whether the operation was successful.

        """
        with self as c:
            c.execute(
              f"""CREATE TABLE {table}
                (composition TEXT NOT NULL, structure TEXT NOT NULL)""",
            )
        return True

    def add_struct(self, struct: SmactStructure, table: str) -> bool:
        """Add a SmactStructure to a table.

        Args:
            struct: The :class:`~.SmactStructure` to add.
            table: The name of the table to add the structure to.

        Returns:
            bool: Whether the operation was successful.

        """
        entry = (struct.composition(), struct.as_poscar())

        with self as c:
            c.execute(f"INSERT into {table} VALUES (?, ?)", entry)

        return True

    def add_structs(self, structs: Sequence[SmactStructure], table: str) -> int:
        """Add several SmactStructures to a table.

        Args:
            structs: Iterable of :class:`~.SmactStructure`s to add to table.
            table: The name of the table to add the structs to.

        Returns:
            The number of structures added.

        """
        with self as c:
            for idx, struct in enumerate(structs):
                entry = (struct.composition(), struct.as_poscar())
                c.execute(f"INSERT into {table} VALUES (?, ?)", entry)

        return idx + 1

    def get_structs(self, composition: str, table: str) -> List[SmactStructure]:
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
              (composition, ), )
            structs = c.fetchall()
        return [SmactStructure.from_poscar(pos[0]) for pos in structs]

    def get_with_species(
      self,
      species: List[Tuple[str, int]],
      table: str, ) -> List[SmactStructure]:
        """Get SmactStructures containing given species."""
        glob = "*".join("{}_*_{}{}" for _ in range(len(species)))
        glob = f"*{glob}*"

        species.sort(key=itemgetter(1), reverse=True)
        species.sort(key=itemgetter(0))

        # Generate a list of [element1, charge1, sign1, element2, ...]
        vals = list(
          itertools.chain.from_iterable([x[0], abs(x[1]), get_sign(x[1])] for x in species)
        )

        glob_form = glob.format(*vals)

        with self as c:
            c.execute(
              f"SELECT structure FROM {table} WHERE composition GLOB ?",
              (glob_form, ), )
            structs = c.fetchall()

        return [SmactStructure.from_poscar(pos[0]) for pos in structs]
